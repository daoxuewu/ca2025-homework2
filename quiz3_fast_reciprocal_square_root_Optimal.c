#include <stdint.h>

/* ===================== Bit utilities ===================== */
/* Count leading zeros (portable, no builtins) */
static inline uint32_t clz32(uint32_t x) {
    if (!x) return 32u;
    uint32_t n = 0;
    if ((x >> 16) == 0) { n += 16; x <<= 16; }
    if ((x >> 24) == 0) { n += 8;  x <<= 8;  }
    if ((x >> 28) == 0) { n += 4;  x <<= 4;  }
    if ((x >> 30) == 0) { n += 2;  x <<= 2;  }
    if ((x >> 31) == 0) { n += 1; }
    return n;
}

/* ===================== 16x16 -> 32 multiply (no bit-by-bit loop) ===================== */
/* Nibble grouping: for each 4-bit chunk of b, accumulate a * (0..15) via shifts+adds. */
static inline uint32_t mul16x16_32(uint32_t a16, uint32_t b16) {
    uint32_t a   = a16 & 0xFFFFu;
    uint32_t acc = 0;
    for (unsigned shift = 0; shift < 16; shift += 4) {
        uint32_t nib = (b16 >> shift) & 0xFu;
        uint32_t p = 0;
        if (nib & 1u) p += a;
        if (nib & 2u) p += (a << 1);
        if (nib & 4u) p += (a << 2);
        if (nib & 8u) p += (a << 3);
        acc += (p << shift);
    }
    return acc; /* exact 16x16 product fits in 32 bits */
}

/* ===================== (x32 * y16) >> 16 without 64-bit multiply ===================== */
/* x = x_hi*2^16 + x_lo  =>  (x*y)>>16 = x_hi*y + ((x_lo*y)>>16) */
static inline uint32_t mul32x16_shr16(uint32_t x32, uint32_t y16) {
    uint32_t x_lo = x32 & 0xFFFFu;
    uint32_t x_hi = x32 >> 16;
    uint32_t p_hi = mul16x16_32(x_hi, y16);
    uint32_t p_lo = mul16x16_32(x_lo, y16) >> 16;
    return p_hi + p_lo;
}

/* ===================== Q16 LUT: 2^16 / sqrt(2^i), i in [0,31] ===================== */
static const uint16_t rsqrt_table[32] = {
    65535, 46341, 32768, 23170, 16384,
    11585,  8192,  5793,  4096,  2896,
     2048,  1448,  1024,   724,   512,
      362,   256,   181,   128,    90,
       64,    45,    32,    23,    16,
       11,     8,     6,     4,     3,
        2,     1
};

/* ===================== One Q16 Newton step (no 64-bit multiply) ===================== */
/*
   y <- y * (3 - x*y*y/2^16) / 2
   Compute xy2_q16 = (x * (y*y)) >> 16 using 16x16 pieces only.
*/
static inline uint32_t q16_newton_step(uint32_t y, uint32_t x) {
    y &= 0xFFFFu;                         /* Q16 */

    /* y^2 -> Q32 (in 32-bit container) */
    uint32_t y2    = mul16x16_32(y, y);
    uint32_t y2_lo = y2 & 0xFFFFu;
    uint32_t y2_hi = y2 >> 16;

    /* x = x_hi<<16 + x_lo */
    uint32_t x_lo = x & 0xFFFFu;
    uint32_t x_hi = x >> 16;

    /* (x*y^2) >> 16 = (x_hi*y2_hi)<<16 + (x_lo*y2_hi) + (x_hi*y2_lo) + ((x_lo*y2_lo)>>16) */
    uint64_t acc  = 0;
    acc += ((uint64_t)mul16x16_32(x_hi, y2_hi)) << 16;
    acc += (uint64_t)mul16x16_32(x_lo, y2_hi);
    acc += (uint64_t)mul16x16_32(x_hi, y2_lo);
    acc += ((uint64_t)mul16x16_32(x_lo, y2_lo)) >> 16;

    uint32_t xy2_q16 = (uint32_t)acc;     /* Q16 */
    uint32_t term    = (3u << 16) - xy2_q16;  /* Q16 */

    /* Split term to avoid general 32x32 multiply:
       (y*term)>>17 = ((y*(term>>16))>>1) + ((y*(term&0xFFFF))>>17) */
    uint32_t term_hi = term >> 16;        /* in {0,1,2,3} */
    uint32_t term_lo = term & 0xFFFFu;

    uint32_t y_mul_hi = (term_hi==3u) ? ((y<<1) + y)
                       : (term_hi==2u) ? (y<<1)
                       : (term_hi==1u) ? y : 0u;
    uint32_t t0 = y_mul_hi >> 1;
    uint32_t t1 = mul16x16_32(y, term_lo) >> 17;

    return t0 + t1;                       /* Q16 */
}

/* ===================== fast_rsqrt: LUT -> interpolate -> Newton (Q16) ===================== */
/* Returns floor(2^16 / sqrt(x)) in Q16; returns 0 for x == 0. */
uint32_t fast_rsqrt(uint32_t x) {
    if (x == 0u) return 0u;

    /* 1) Range bucket: find exp such that x in [2^exp, 2^(exp+1)) */
    uint32_t e = 31u - clz32(x);

    /* 2) Base and next estimates from LUT */
    uint32_t y0 = rsqrt_table[e];
    uint32_t y1 = (e < 31u) ? rsqrt_table[e + 1u] : 1u;

    /* 3) Linear interpolation:
          frac = ((x - 2^e) << 16) >> e  in [0, 2^16) without 64-bit ops */
    uint32_t base = 1u << e;              /* valid even for e=31 (0x80000000) */
    uint32_t diff = x - base;
    uint32_t frac = (e >= 16u) ? (diff >> (e - 16u)) : (diff << (16u - e));

    /* y = y0 - ((y0 - y1) * frac >> 16) */
    uint32_t dy = y0 - y1;
    uint32_t y  = y0 - (mul16x16_32(dy, frac) >> 16);

    /* 4) One Newton refinement (keeps behavior identical to your version) */
    y = q16_newton_step(y, x);

    return y; /* Q16 */
}
