#include <stdint.h>

/* -------------------- Utilities: CLZ and 32×32->64 shift-add multiply -------------------- */
/* Portable software CLZ (count leading zeros) */
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

/* Shift-add multiplication, RV32I-friendly (no hardware MUL). Returns a 64-bit product. */
static inline uint64_t mul32_shift_add(uint32_t a, uint32_t b) {
    uint64_t acc = 0;
    uint64_t aa  = (uint64_t)a;
    while (b) {
        if (b & 1u) acc += aa;
        aa <<= 1;
        b  >>= 1;
    }
    return acc;
}

/* -------------------- 32-entry Q16 lookup table: 2^16 / sqrt(2^i) -------------------- */
static const uint16_t rsqrt_table[32] = {
    UINT16_MAX, 46341, 32768, 23170, 16384,  /* 2^0..2^4  */
    11585,  8192,  5793,  4096,  2896,       /* 2^5..2^9  */
     2048,  1448,  1024,   724,   512,       /* 2^10..14  */
      362,   256,   181,   128,    90,       /* 2^15..19  */
       64,    45,    32,    23,    16,       /* 2^20..24  */
       11,     8,     6,     4,     3,       /* 2^25..29  */
        2,     1                            /* 2^30,31   */
};

/* One Q16 Newton step:
   y_{k+1} = y_k * (3 - x*y_k^2 / 2^16) / 2
   y and (3<<16) are Q16; all intermediate products use 64-bit accumulations to avoid overflow. */
static inline uint32_t q16_newton_step(uint32_t y, uint32_t x) {
    uint64_t y2 = mul32_shift_add(y, y);              /* y^2 (Q32 in 64-bit) */

    uint32_t y2_lo = (uint32_t)y2;
    uint32_t y2_hi = (uint32_t)(y2 >> 32);

    uint64_t prod_lo = mul32_shift_add(x, y2_lo);     /* x * y2_lo */
    uint64_t prod_hi = mul32_shift_add(x, y2_hi);     /* x * y2_hi */

    /* (x*y^2) >> 16 with rounding-to-nearest on the low part */
    uint64_t xy2_q16_64 = (prod_hi << 16) + ((prod_lo + (1ull << 15)) >> 16);
    uint32_t xy2_q16 = (uint32_t)xy2_q16_64;          /* Q16 */

    uint32_t term = (3u << 16) - xy2_q16;             /* Q16 */
    uint64_t prod = mul32_shift_add(y, term);         /* y * term */

    /* Final >>17 with rounding-to-nearest (add half LSB before shift) */
    return (uint32_t)((prod + (1ull << 16)) >> 17);   /* Q16 */
}

/* -------------------- fast_rsqrt core: LUT -> interpolate -> Newton -------------------- */
/* Returns y ≈ floor(2^16 / sqrt(x)); returns 0 for x == 0. */
uint32_t fast_rsqrt(uint32_t x) {
    if (x == 0u) return 0u;

    /* 1) Find MSB position to locate range [2^exp, 2^(exp+1)) */
    uint32_t exp = 31u - clz32(x);

    /* 2) Pick base and next LUT values (for exp=31, next is clamped to 1) */
    uint32_t y_base = rsqrt_table[exp];
    uint32_t y_next = (exp < 31u) ? rsqrt_table[exp + 1u] : 1u;

    /* 3) Linear interpolation:
          frac = ((x - 2^exp) << 16) >> exp  ∈ [0, 2^16) */
    uint32_t one_exp = (exp < 31u) ? (1u << exp) : 0u;  /* keep behavior unchanged */
    uint32_t frac = (uint32_t)((((uint64_t)x - (uint64_t)one_exp) << 16) >> exp);

    uint32_t delta = y_base - y_next;
    uint32_t y = y_base - (uint32_t)(mul32_shift_add(delta, frac) >> 16);

    /* 4) Two Newton iterations (as in the original) */
    y = q16_newton_step(y, x);
    y = q16_newton_step(y, x);

    return y;  /* Q16 fixed-point result representing 2^16 / sqrt(x) */
}
