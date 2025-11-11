    .attribute arch, "rv32i2p1_zicsr2p0"     # ISA: RV32I v2.1 + Zicsr v2.0
    .text
    .globl  test_Hanoi

# -----------------------------------------------------------------------------
# Iterative Tower of Hanoi (3 disks) using Gray code
# Registers:
#   x8  = step counter (n)
#   x9  = moved disk index (0..2)
#   x18 = current position of the selected disk (0..2)
#   x19 = next position of the selected disk (0..2)
#   x20 = &outbuf (persisted across loop)
# Stack layout (32 bytes):
#   [sp+00]=x8, [sp+04]=x9, [sp+08]=x18, [sp+12]=x19, [sp+16]=x20
#   [sp+20]=pos(disk0), [sp+24]=pos(disk1), [sp+28]=pos(disk2)
# Output per move: "Move Disk <n> from <X> to <Y>\n" via rv32emu SYS_write (a7=64)
# -----------------------------------------------------------------------------
test_Hanoi:
    addi    x2, x2, -32
    sw      x8,  0(x2)
    sw      x9,  4(x2)
    sw      x18, 8(x2)
    sw      x19, 12(x2)
    sw      x20, 16(x2)

    # Keep &outbuf in x20 across the loop
    la      x20, outbuf

    # Initialize all disk positions to peg 0 ('A')
    sw      x0, 20(x2)                     # disk 0
    sw      x0, 24(x2)                     # disk 1
    sw      x0, 28(x2)                     # disk 2

    addi    x8, x0, 1                      # n = 1

game_loop:
    # Stop after 2^3 (=8) states (n == 8)
    addi    x5, x0, 8
    beq     x8, x5, finish_game

    # Gray(n)   = n ^ (n >> 1)
    # Gray(n-1) = (n-1) ^ ((n-1) >> 1)
    srli    x5, x8, 1
    xor     x6, x8, x5                     # x6 = Gray(n)

    addi    x7, x8, -1
    srli    x28, x7, 1
    xor     x7, x7, x28                    # x7 = Gray(n-1)

    # Which bit changed between consecutive Gray codes?
    xor     x5, x6, x7                     # bitmask of moved disk

    # Decode moved disk index into x9 ∈ {0,1,2}
    addi    x9, x0, 0
    andi    x6, x5, 1
    bne     x6, x0, disk_found
    addi    x9, x0, 1
    andi    x6, x5, 2
    bne     x6, x0, disk_found
    addi    x9, x0, 2

disk_found:
    # Load current position pos[disk]
    slli    x5, x9, 2                      # index * 4
    addi    x5, x5, 20                     # base offset of pos[]
    add     x5, x2, x5
    lw      x18, 0(x5)                     # x18 = pos[disk]

    # Compute next position x19
    bne     x9, x0, handle_large           # disk 0 (smallest) moves every step
    # Smallest disk moves +2 (mod 3)
    addi    x19, x18, 2
    addi    x6,  x0, 3
    blt     x19, x6, display_move
    sub     x19, x19, x6
    jal     x0, display_move

handle_large:
    # Larger disk: target = (0+1+2) - pos(this) - pos(smallest)
    lw      x6, 20(x2)                     # pos(smallest)
    addi    x19, x0, 3
    sub     x19, x19, x18
    sub     x19, x19, x6

display_move:
    # Map peg index -> ASCII via table
    la      x6, peg_chars
    add     x5,  x6, x18
    lbu     x11, 0(x5)                     # src char 'A'|'B'|'C'
    add     x7,  x6, x19
    lbu     x12, 0(x7)                     # dst char 'A'|'B'|'C'

    # Patch message template in outbuf:
    #   "Move Disk 0 from X to Y\n"
    #               ^    ^    ^
    #             [10] [17]  [22]
    addi    x28, x9, 1
    addi    x28, x28, 48                   # '1'..'3'
    sb      x28, 10(x20)                   # digit
    sb      x11, 17(x20)                   # source peg
    sb      x12, 22(x20)                   # destination peg

    # write(fd=1, buf=&outbuf, len=24)
    li      a0, 1
    li      a7, 64
    li      a2, 24
    addi    a1, x20, 0
    ecall

    # Store updated position and iterate
    slli    x5, x9, 2
    addi    x5, x5, 20
    add     x5, x2, x5
    sw      x19, 0(x5)

    addi    x8, x8, 1                      # n++
    jal     x0, game_loop

finish_game:
    lw      x8,  0(x2)
    lw      x9,  4(x2)
    lw      x18, 8(x2)
    lw      x19, 12(x2)
    lw      x20, 16(x2)
    addi    x2,  x2, 32
    ret

    .data
    .balign 4
outbuf:     .ascii  "Move Disk 0 from X to Y\n"
peg_chars:  .byte   'A','B','C'
