    .text

# ------------------------------------------------------------
# uint32_t clz_branchless(uint32_t x)
# Branchless binary-search CLZ: uses sltu to form shift amounts
# ------------------------------------------------------------
    .globl clz_branchless
clz_branchless:
    li      t0, 32          # n = 32
    mv      t1, a0          # t1 = x
    srli    t2, t1, 16      # y = x >> 16
    sltu    t3, zero, t2    # b = (y != 0)
    slli    t4, t3, 4       # s = b * 16
    srl     t1, t1, t4      # x >>= s
    sub     t0, t0, t4      # n -= s

    srli    t2, t1, 8
    sltu    t3, zero, t2
    slli    t4, t3, 3       # s = b * 8
    srl     t1, t1, t4
    sub     t0, t0, t4

    srli    t2, t1, 4
    sltu    t3, zero, t2
    slli    t4, t3, 2       # s = b * 4
    srl     t1, t1, t4
    sub     t0, t0, t4

    srli    t2, t1, 2
    sltu    t3, zero, t2
    slli    t4, t3, 1       # s = b * 2
    srl     t1, t1, t4
    sub     t0, t0, t4

    srli    t2, t1, 1
    sltu    t3, zero, t2    # b = (x>>1) != 0
    mv      t4, t3          # s = b * 1
    srl     t1, t1, t4
    sub     t0, t0, t4

    sub     a0, t0, t1      # return n - x (x becomes 0 or 1)
    ret

# ------------------------------------------------------------
# uint32_t uf8_decode(uint8_t fl)
# D(b) = (m << e) + ((1<<e) - 1) << 4
# ------------------------------------------------------------
    .globl uf8_decode
uf8_decode:
    andi    t0, a0, 0x0F        # m = fl & 0x0f
    srli    t1, a0, 4           # e = fl >> 4
    li      t2, 1
    sll     t2, t2, t1
    addi    t2, t2, -1
    slli    t2, t2, 4           # offset = ((1<<e)-1) << 4
    sll     t0, t0, t1
    add     a0, t0, t2          # value = (m << e) + offset
    ret

# ------------------------------------------------------------
# uint8_t uf8_encode(uint32_t value)
# Uses branchless CLZ to locate msb -> initial exponent guess
# ------------------------------------------------------------
    .globl uf8_encode
uf8_encode:
    addi    sp, sp, -16
    sw      ra, 12(sp)          # store ra (because we call clz_branchless)
    sw      s0, 8(sp)           # s0 = exponent
    sw      s1, 4(sp)           # s1 = overflow (offset)

    sltiu   t0, a0, 16          # if (value < 16) return value;
    beqz    t0, enCode_normalVal
    andi    a0, a0, 0xFF        # keep in 8-bit just in case
    j       enCode_ret

enCode_normalVal:
    mv      t5, a0              # v = value
    jal     ra, clz_branchless  # a0 = clz(value)
    li      t0, 31
    sub     t1, t0, a0          # msb = 31 - lz
    mv      a0, t5              # restore value

    li      s0, 0               # exponent = 0
    li      s1, 0               # overflow = 0

    sltiu   t2, t1, 5           # if (msb >= 5) => !(t1<5)
    bnez    t2, enCode_find_up  # if msb < 5, skip initial guess
    addi    s0, t1, -4          # exponent = msb - 4
    sltiu   t3, s0, 16
    bnez    t3, enCode_initGuess_overflow_ok
    li      s0, 15              # clamp exponent to 15

# overflow = ((1<<e)-1) << 4
enCode_initGuess_overflow_ok:
    li      s1, 0               # overflow = 0
    li      t4, 0               # e = 0
enCode_overflow_loop:
    bge     t4, s0, enCode_adjust_down
    slli    s1, s1, 1
    addi    s1, s1, 16
    addi    t4, t4, 1
    j       enCode_overflow_loop

# if value < overflow, step exponent down until it fits
enCode_adjust_down:
    beqz    s0, enCode_find_up
    sltu    t4, a0, s1
    beqz    t4, enCode_find_up
    addi    s1, s1, -16
    srli    s1, s1, 1           # overflow = (overflow - 16) >> 1
    addi    s0, s0, -1
    j       enCode_adjust_down

# then go upward to the exact exponent
enCode_find_up:
    li      t4, 15
enCode_up_loop:
    bge     s0, t4, enCode_up_done
    slli    t1, s1, 1
    addi    t1, t1, 16          # next_overflow = (overflow << 1) + 16
    sltu    t2, a0, t1
    bnez    t2, enCode_up_done
    mv      s1, t1
    addi    s0, s0, 1
    j       enCode_up_loop

enCode_up_done:
    sub     t0, a0, s1          # num = value - overflow
    srl     t0, t0, s0          # mantissa = num >> exponent
    slli    t1, s0, 4
    andi    t0, t0, 0x0F        # mantissa &= 0x0F
    or      a0, t1, t0          # a0 = (exponent<<4) | mantissa
    andi    a0, a0, 0xFF

enCode_ret:
    lw      ra, 12(sp)
    lw      s0, 8(sp)
    lw      s1, 4(sp)
    addi    sp, sp, 16
    ret

# ------------------------------------------------------------
# int quiz1_uf8_test(void)
# Prints each case, checks decode/encode roundtrip and monotonicity
# return: a0 = 1 (pass) / 0 (fail)
# ------------------------------------------------------------
    .globl quiz1_uf8_test
quiz1_uf8_test:
    addi    sp, sp, -40
    sw      ra, 36(sp)
    sw      s0, 32(sp)          # previous_value
    sw      s1, 28(sp)          # passed
    sw      s2, 24(sp)          # i
    sw      s3, 20(sp)          # value
    sw      t0, 16(sp)
    sw      t1, 12(sp)
    sw      t2, 8(sp)
    sw      t3, 4(sp)

    li      s0, -1              # previous_value = -1
    li      s1, 1               # passed = true
    li      s2, 0               # i = 0

t_loop:
    li      t3, 256
    bge     s2, t3, t_done

    la      a0, msg_test        # "test data: "
    li      a7, 4
    ecall

    mv      a0, s2              # print i
    li      a7, 1
    ecall

    la      a0, msg_nl          # newline
    li      a7, 4
    ecall

    mv      a0, s2              # fl = i
    jal     ra, uf8_decode      # value = decode(fl)
    mv      s3, a0

    mv      a0, s3
    jal     ra, uf8_encode      # fl2 = encode(value)
    mv      t1, a0

    bne     s2, t1, t_fail_flag # if (fl != fl2) report

    slt     t2, s0, s3          # if (previous_value < value) OK
    bnez    t2, t_set_prev

# non-increasing
t_bad_inc:
    la      a0, msg_noninc_a
    li      a7, 4
    ecall
    mv      a0, s2              # fl
    li      a7, 1
    ecall
    la      a0, msg_noninc_b
    li      a7, 4
    ecall
    mv      a0, s3              # value
    li      a7, 1
    ecall
    la      a0, msg_noninc_c
    li      a7, 4
    ecall
    mv      a0, s0              # previous_value
    li      a7, 1
    ecall
    la      a0, msg_nl
    li      a7, 4
    ecall
    li      s1, 0               # passed = false
    j       t_set_prev

# mismatch fl vs fl2
t_fail_flag:
    la      a0, msg_mismatch_a
    li      a7, 4
    ecall
    mv      a0, s2              # fl
    li      a7, 1
    ecall
    la      a0, msg_mismatch_b
    li      a7, 4
    ecall
    mv      a0, s3              # value
    li      a7, 1
    ecall
    la      a0, msg_mismatch_c
    li      a7, 4
    ecall
    mv      a0, t1              # fl2
    li      a7, 1
    ecall
    la      a0, msg_nl
    li      a7, 4
    ecall
    li      s1, 0               # passed = false

t_set_prev:
    mv      s0, s3              # previous_value = value
    addi    s2, s2, 1           # i++
    j       t_loop

t_done:
    mv      a0, s1              # return passed
    lw      ra, 36(sp)
    lw      s0, 32(sp)
    lw      s1, 28(sp)
    lw      s2, 24(sp)
    lw      s3, 20(sp)
    lw      t0, 16(sp)
    lw      t1, 12(sp)
    lw      t2, 8(sp)
    lw      t3, 4(sp)
    addi    sp, sp, 40
    ret

# ------------------------------------------------------------
# 可選入口（示範用）：不與 main 衝突
# int quiz1_uf8_entry(void) { prints; return quiz1_uf8_test(); }
# ------------------------------------------------------------
    .globl quiz1_uf8_entry
quiz1_uf8_entry:
    addi    sp, sp, -16
    sw      ra, 12(sp)

    jal     ra, quiz1_uf8_test
    beqz    a0, q_entry_fail

    la      a0, msg_ok          # "All tests passed.\n"
    li      a7, 4
    ecall
    j       q_entry_exit

q_entry_fail:
    la      a0, msg_fail
    li      a7, 4
    ecall

q_entry_exit:
    lw      ra, 12(sp)
    addi    sp, sp, 16
    ret

# ------------------------------------------------------------
# Data section
# ------------------------------------------------------------
    .data
    .align 4
msg_ok:         .asciz "All tests passed.\n"
msg_fail:       .asciz "Failed.\n"
msg_test:       .asciz "test data: "
msg_nl:         .asciz "\n"
msg_mismatch_a: .asciz "mismatch: fl= "
msg_mismatch_b: .asciz " value= "
msg_mismatch_c: .asciz " fl2= "
msg_noninc_a:   .asciz "non-increasing: fl= "
msg_noninc_b:   .asciz " value= "
msg_noninc_c:   .asciz " prev= "
