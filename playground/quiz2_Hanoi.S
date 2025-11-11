    .attribute arch, "rv32i2p1_zicsr2p0"   # ISA: RV32I v2.1 + Zicsr v2.0
    .text
    .globl  test_Hanoi

# -----------------------------------------------------------------------------
# Iterative Tower of Hanoi (3 disks) using Gray code.
# - Stack frame (32 bytes): save x8,x9,x18,x19,x20 and 3 disk positions at [sp+20,24,28].
# - Disk positions are in {0,1,2}. Smallest disk (0) moves every step.
# - Printing via rv32emu SYS_write (a7=64).
# -----------------------------------------------------------------------------
test_Hanoi:
    addi    x2, x2, -32
    sw      x8, 0(x2)
    sw      x9, 4(x2)
    sw      x18, 8(x2)
    sw      x19, 12(x2)
    sw      x20, 16(x2)

    # Initialize disk positions to 0 (peg 'A')
    sw      x0, 20(x2)          # disk 0
    sw      x0, 24(x2)          # disk 1
    sw      x0, 28(x2)          # disk 2

    addi    x8, x0, 1           # step counter n = 1

game_loop:
    # Stop after 2^3 (=8) moves
    addi    x5, x0, 8
    beq     x8, x5, finish_game

    # Gray(n) = n ^ (n >> 1)
    srli    x5, x8, 1
    xor     x6, x8, x5          # x6 = gray(n)

    # Previous Gray: Gray(n-1) = (n-1) ^ ((n-1)>>1)
    addi    x7, x8, -1
    srli    x28, x7, 1
    xor     x7, x7, x28         # x7 = gray(n-1)

    # Changed bit between consecutive codes
    xor     x5, x6, x7          # bit mask of moved disk

    # Decode disk index x9 ∈ {0,1,2}
    addi    x9, x0, 0           # assume disk 0
    andi    x6, x5, 1           # test bit 0
    bne     x6, x0, disk_found
    addi    x9, x0, 1           # assume disk 1
    andi    x6, x5, 2           # test bit 1
    bne     x6, x0, disk_found
    addi    x9, x0, 2           # else disk 2

disk_found:
    # (Optional sanity check for multi-bit change) — no-op by design
    andi    x30, x5, 5
    addi    x31, x0, 5
    beq     x30, x31, 1f
1:

    # Load current position of disk x9
    slli    x5, x9, 2           # word index * 4
    addi    x5, x5, 20          # base offset of disk array
    add     x5, x2, x5
    lw      x18, 0(x5)          # x18 = pos(disk x9)

    # Compute next position x19
    bne     x9, x0, handle_large

    # Smallest disk moves +2 mod 3
    addi    x19, x18, 2
    addi    x6, x0, 3
    blt     x19, x6, display_move
    sub     x19, x19, x6
    jal     x0, display_move

handle_large:
    # Larger disk: target = (0+1+2) - pos(this) - pos(smallest)
    lw      x6, 20(x2)          # pos(smallest disk)
    addi    x19, x0, 3
    sub     x19, x19, x18
    sub     x19, x19, x6

display_move:
    # Map peg index {0,1,2} -> ASCII {'A','B','C'} via obfuscated table
    la      x20, obdata

    add     x5, x20, x18        # source peg
    lbu     x11, 0(x5)
    li      x6, 0x6F
    xor     x11, x11, x6
    addi    x11, x11, -0x12     # x11 = 'A'|'B'|'C'

    add     x7, x20, x19        # destination peg
    lbu     x12, 0(x7)
    xor     x12, x12, x6
    addi    x12, x12, -0x12     # x12 = 'A'|'B'|'C'

    # -------------------- Print: "Move Disk <n> from <X> to <Y>\n" --------------------
    mv      t2, x11             # src char
    mv      t3, x12             # dst char

    # "Move Disk "
    li      a7, 64              # SYS_write
    li      a0, 1               # STDOUT
    la      a1, str1
    la      a2, str1_end
    la      t0, str1
    sub     a2, a2, t0
    addi    a2, a2, -1          # strip NUL
    ecall

    # disk number '1'..'3'
    addi    t0, x9, 1
    addi    t0, t0, 48          # '0' + (1..3)
    la      t1, tmp_ch
    sb      t0, 0(t1)
    li      a7, 64
    li      a0, 1
    la      a1, tmp_ch
    addi    a2, x0, 1
    ecall

    # " from "
    li      a7, 64
    li      a0, 1
    la      a1, str2
    la      a2, str2_end
    la      t0, str2
    sub     a2, a2, t0
    addi    a2, a2, -1
    ecall

    # source peg
    la      t1, tmp_ch
    sb      t2, 0(t1)
    li      a7, 64
    li      a0, 1
    la      a1, tmp_ch
    addi    a2, x0, 1
    ecall

    # " to "
    li      a7, 64
    li      a0, 1
    la      a1, str3
    la      a2, str3_end
    la      t0, str3
    sub     a2, a2, t0
    addi    a2, a2, -1
    ecall

    # destination peg
    la      t1, tmp_ch
    sb      t3, 0(t1)
    li      a7, 64
    li      a0, 1
    la      a1, tmp_ch
    addi    a2, x0, 1
    ecall

    # newline
    li      t0, 10              # '\n'
    la      t1, tmp_ch
    sb      t0, 0(t1)
    li      a7, 64
    li      a0, 1
    la      a1, tmp_ch
    addi    a2, x0, 1
    ecall

    # Store updated position and continue
    slli    x5, x9, 2
    addi    x5, x5, 20
    add     x5, x2, x5
    sw      x19, 0(x5)

    addi    x8, x8, 1           # n++
    jal     x0, game_loop

finish_game:
    lw      x8, 0(x2)
    lw      x9, 4(x2)
    lw      x18, 8(x2)
    lw      x19, 12(x2)
    lw      x20, 16(x2)
    addi    x2, x2, 32
    ret

    .data
obdata:     .byte   0x3c, 0x3b, 0x3a
str1:       .asciz  "Move Disk "
str1_end:
str2:       .asciz  " from "
str2_end:
str3:       .asciz  " to "
str3_end:
tmp_ch:     .byte   0
