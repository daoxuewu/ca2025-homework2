#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>   // for size_t

#define printstr(ptr, length)                   \
    do {                                        \
        asm volatile(                           \
            "add a7, x0, 0x40;"                 \
            "add a0, x0, 0x1;"  /* stdout */    \
            "add a1, x0, %0;"                   \
            "mv  a2, %1;"      /* byte length */\
            "ecall;"                            \
            :                                   \
            : "r"(ptr), "r"(length)             \
            : "a0", "a1", "a2", "a7");          \
    } while (0)

#define TEST_OUTPUT(msg, length) printstr(msg, length)

#define TEST_LOGGER(msg)                     \
    do {                                     \
        char _msg[] = msg;                   \
        TEST_OUTPUT(_msg, sizeof(_msg) - 1); \
    } while (0)

extern uint64_t get_cycles(void);
extern uint64_t get_instret(void);

/* ---------------- Bare-metal memcpy (byte-wise) ---------------- */
void *memcpy(void *dest, const void *src, size_t n)
{
    uint8_t *d = (uint8_t *)dest;
    const uint8_t *s = (const uint8_t *)src;
    while (n--) *d++ = *s++;
    return dest;
}

/* ---------------- Software div/mod/mul for RV32I (no M extension) ---------------- */
static unsigned long udiv(unsigned long dividend, unsigned long divisor)
{
    if (divisor == 0) return 0;

    unsigned long q = 0, r = 0;
    for (int i = 31; i >= 0; i--) {
        r <<= 1;
        r |= (dividend >> i) & 1UL;
        if (r >= divisor) {
            r -= divisor;
            q |= (1UL << i);
        }
    }
    return q;
}

static unsigned long umod(unsigned long dividend, unsigned long divisor)
{
    if (divisor == 0) return 0;

    unsigned long r = 0;
    for (int i = 31; i >= 0; i--) {
        r <<= 1;
        r |= (dividend >> i) & 1UL;
        if (r >= divisor) r -= divisor;
    }
    return r;
}

static uint32_t umul(uint32_t a, uint32_t b)
{
    uint32_t res = 0;
    while (b) {
        if (b & 1U) res += a;
        a <<= 1;
        b >>= 1;
    }
    return res;
}

/* GCC helper for soft-mul */
uint32_t __mulsi3(uint32_t a, uint32_t b) { return umul(a, b); }

/* ---------------- Printing helpers (no libc printf) ---------------- */
static void print_hex(unsigned long val)
{
    char buf[20];
    char *p = buf + sizeof(buf) - 1;
    *p-- = '\n';

    if (val == 0) {
        *p-- = '0';
    } else {
        while (val) {
            int d = (int)(val & 0xFUL);
            *p-- = (char)(d < 10 ? '0' + d : 'a' + (d - 10));
            val >>= 4;
        }
    }
    p++;
    printstr(p, (unsigned long)(buf + sizeof(buf) - p));
}

static void print_dec(unsigned long val)
{
    char buf[20];
    char *p = buf + sizeof(buf) - 1;
    *p-- = '\n';

    if (val == 0) {
        *p-- = '0';
    } else {
        while (val) {
            *p-- = (char)('0' + umod(val, 10));
            val  = udiv(val, 10);
        }
    }
    p++;
    printstr(p, (unsigned long)(buf + sizeof(buf) - p));
}

/* decimal printing without trailing newline */
static void print_dec_inline(unsigned long val)
{
    char buf[20];
    char *p = buf + sizeof(buf);
    do {
        *--p = (char)('0' + umod(val, 10));
        val   = udiv(val, 10);
    } while (val);
    printstr(p, (unsigned long)(buf + sizeof(buf) - p));
}

static inline void print_ch(char c) { printstr(&c, 1); }

/* print a C-string (null-terminated) */
static inline void print_str(const char *s) {
    const char *p = s;
    while (*p) p++;
    printstr(s, (unsigned long)(p - s));
}

/* Print a Q16 fixed-point value as "int.frac" with given digits (integer-only) */
static void print_q16_u(uint32_t y_q16, int frac_digits)
{
    print_dec_inline(y_q16 >> 16);
    if (frac_digits <= 0) { print_ch('\n'); return; }

    print_ch('.');
    uint32_t frac = y_q16 & 0xFFFFu;
    for (int i = 0; i < frac_digits; i++) {
        uint32_t frac10 = (frac << 3) + (frac << 1);  /* *10 in Q16 */
        uint32_t digit  = frac10 >> 16;               /* integer part */
        frac            = frac10 & 0xFFFFu;           /* remainder */
        print_ch((char)('0' + digit));
    }
    print_ch('\n');
}

/* ---------------- External test targets ---------------- */
extern void chacha20(uint8_t *out,
                     const uint8_t *in,
                     size_t inlen,
                     const uint8_t *key,
                     const uint8_t *nonce,
                     uint32_t ctr);

typedef uint8_t uf8;
extern uint32_t uf8_decode(uf8 fl);
extern uf8      uf8_encode(uint32_t value);

extern uint32_t fast_rsqrt(uint32_t x);

/* ---------------- Tests ---------------- */
static void test_UF8(void)
{
    int32_t previous_value = -1;

    for (int i = 0; i < 8; i++) {
        TEST_LOGGER("  Data: ");
        print_dec((unsigned long)i);

        uint8_t  fl    = (uint8_t)i;
        int32_t  value = (int32_t)uf8_decode(fl);
        uint8_t  fl2   = uf8_encode((uint32_t)value);

        if (fl != fl2) {
            TEST_LOGGER("  Mismatch!\n");
            return;
        }
        if (value <= previous_value) {
            TEST_LOGGER("  Mismatch!\n");
            return;
        }
        previous_value = value;
    }
    TEST_LOGGER("  PASSED\n");
}

extern void test_Hanoi(void);

void test_Fast_rsqrt(void)
{
    static const uint32_t tests[] = {
        1, 2, 4, 5, 10, 16, 20, 100, 1000, 0xFFFFFFFFu
    };
    for (unsigned i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        uint32_t x  = tests[i];
        uint32_t yq = fast_rsqrt(x);  /* y ≈ 2^16 / sqrt(x) */

        print_str("x=");
        print_dec_inline(x);
        print_str("  fast_rsqrt≈");
        print_q16_u(yq, 4);           /* show 4 fractional digits */
    }
}

/* RFC 7539 §2.4.2 test (kept for reference; not invoked by main) */
static void test_chacha20(void)
{
    const uint8_t key[32] = {
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
        16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    };
    const uint8_t nonce[12] = {0, 0, 0, 0, 0, 0, 0, 74, 0, 0, 0, 0};
    const uint32_t ctr = 1;

    uint8_t in[114] =
        "Ladies and Gentlemen of the class of '99: If I could offer you only "
        "one tip for the future, sunscreen would be it.";
    uint8_t out[114];

    static const uint8_t exp[] = {
        0x6e, 0x2e, 0x35, 0x9a, 0x25, 0x68, 0xf9, 0x80, 0x41, 0xba, 0x07, 0x28,
        0xdd, 0x0d, 0x69, 0x81, 0xe9, 0x7e, 0x7a, 0xec, 0x1d, 0x43, 0x60, 0xc2,
        0x0a, 0x27, 0xaf, 0xcc, 0xfd, 0x9f, 0xae, 0x0b, 0xf9, 0x1b, 0x65, 0xc5,
        0x52, 0x47, 0x33, 0xab, 0x8f, 0x59, 0x3d, 0xab, 0xcd, 0x62, 0xb3, 0x57,
        0x16, 0x39, 0xd6, 0x24, 0xe6, 0x51, 0x52, 0xab, 0x8f, 0x53, 0x0c, 0x35,
        0x9f, 0x08, 0x61, 0xd8, 0x07, 0xca, 0x0d, 0xbf, 0x50, 0x0d, 0x6a, 0x61,
        0x56, 0xa3, 0x8e, 0x08, 0x8a, 0x22, 0xb6, 0x5e, 0x52, 0xbc, 0x51, 0x4d,
        0x16, 0xcc, 0xf8, 0x06, 0x81, 0x8c, 0xe9, 0x1a, 0xb7, 0x79, 0x37, 0x36,
        0x5a, 0xf9, 0x0b, 0xbf, 0x74, 0xa3, 0x5b, 0xe6, 0xb4, 0x0b, 0x8e, 0xed,
        0xf2, 0x78, 0x5e, 0x42, 0x87, 0x4d
    };

    TEST_LOGGER("Test: ChaCha20\n");
    chacha20(out, in, sizeof(in), key, nonce, ctr);

    bool passed = true;
    for (size_t i = 0; i < sizeof(exp); i++) {
        if (out[i] != exp[i]) { passed = false; break; }
    }
    if (passed) {
        TEST_LOGGER("  ChaCha20 RFC 7539: PASSED\n");
    } else {
        TEST_LOGGER("  ChaCha20 RFC 7539: FAILED\n");
    }
}

/* ---------------- Main ---------------- */
int main(void)
{
    uint64_t start_cycles, end_cycles, cycles_elapsed;
    uint64_t start_instret, end_instret, instret_elapsed;

    /* Test 0: UF8 */
    TEST_LOGGER("\n=== Uf8 tests ===\n");
    start_cycles   = get_cycles();
    start_instret  = get_instret();
    test_UF8();
    end_cycles     = get_cycles();
    end_instret    = get_instret();
    cycles_elapsed   = end_cycles   - start_cycles;
    instret_elapsed  = end_instret  - start_instret;
    TEST_LOGGER("  Cycles: ");       print_dec((unsigned long)cycles_elapsed);
    TEST_LOGGER("  Instructions: "); print_dec((unsigned long)instret_elapsed);
    TEST_LOGGER("\n");

    /* Test 1: Hanoi */
    TEST_LOGGER("\n=== Hanoi tower tests ===\n\n");
    start_cycles   = get_cycles();
    start_instret  = get_instret();
    test_Hanoi();
    end_cycles     = get_cycles();
    end_instret    = get_instret();
    cycles_elapsed   = end_cycles   - start_cycles;
    instret_elapsed  = end_instret  - start_instret;
    TEST_LOGGER("  Cycles: ");       print_dec((unsigned long)cycles_elapsed);
    TEST_LOGGER("  Instructions: "); print_dec((unsigned long)instret_elapsed);
    TEST_LOGGER("\n");

    /* Test 2: Fast reciprocal square root */
    TEST_LOGGER("\n=== Fast reciprocal square root tests ===\n\n");
    start_cycles   = get_cycles();
    start_instret  = get_instret();
    test_Fast_rsqrt();
    end_cycles     = get_cycles();
    end_instret    = get_instret();
    cycles_elapsed   = end_cycles   - start_cycles;
    instret_elapsed  = end_instret  - start_instret;
    TEST_LOGGER("  Cycles: ");       print_dec((unsigned long)cycles_elapsed);
    TEST_LOGGER("  Instructions: "); print_dec((unsigned long)instret_elapsed);
    TEST_LOGGER("\n");

    TEST_LOGGER("\n=== All Tests Completed ===\n");
    return 0;
}
