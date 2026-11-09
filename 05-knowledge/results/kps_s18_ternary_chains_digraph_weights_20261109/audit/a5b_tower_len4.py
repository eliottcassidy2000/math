#!/usr/bin/env python3
"""AUDIT extension of THM-5070 (5): products of length exactly 4 of the 57 base-3 order-<=4 blocks; list the distinct
products N with N (3,1,2) parallel to (3,1,2) and test whether each is a power of M_T. (Not claimed by the packet;
checks how far the 'no digit system' gloss is supported.)  Usage: a5b_tower_len4.py MATS_X3"""
import sys
mats = []
for line in open(sys.argv[1]):
    kmin, mat, rep = line.split("|")
    if int(kmin) <= 4:
        mats.append(tuple(map(int, mat.split())))
def mul(A, B):
    return tuple(sum(A[3*i+k]*B[3*k+j] for k in range(3)) for i in range(3) for j in range(3))
MT = (1, 3, 0, 0, 0, 1, 1, 1, 0)
powers, P = {}, (1, 0, 0, 0, 1, 0, 0, 0, 1)
for k in range(1, 30):
    P = mul(MT, P); powers[P] = k
L3 = set(mats)
for _ in range(2):
    L3 = {mul(M, S) for S in L3 for M in mats}
v = (3, 1, 2)
found = {}
for S in L3:
    w = (S[0]*3+S[1]+S[2]*2, S[3]*3+S[4]+S[5]*2, S[6]*3+S[7]+S[8]*2)
    for M in mats:
        x = (M[0]*w[0]+M[1]*w[1]+M[2]*w[2], M[3]*w[0]+M[4]*w[1]+M[5]*w[2], M[6]*w[0]+M[7]*w[1]+M[8]*w[2])
        if x[0]*1 == x[1]*3 and x[2]*1 == x[1]*2:
            N = mul(M, S)
            found[N] = powers.get(N, -1)
print(f"length-3 layer: {len(L3)} distinct products; length-4 products preserving (3,1,2): {len(found)} distinct;"
      f" tower powers {sorted(p for p in found.values() if p > 0)}; non-tower {sum(1 for p in found.values() if p < 0)}")
