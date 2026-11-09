#!/usr/bin/env python3
"""AUDIT: THM-5070 (5). Products of length <= L of my 57 base-3 order-<= 4 block matrices (a1_blocks.py, mats_x3.txt,
kmin <= 4); keep the distinct products N with N v parallel to v = (3, 1, 2) (exact cross-product test) and say
which are powers of M_T = [[1,3,0],[0,0,1],[1,1,0]]. Also lists the length-1 blocks equal to M_T^j.
Usage: a5_tower_direction.py MATS_X3 L
"""
import sys


def mul(A, B):
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)) for i in range(3))


def main():
    mats = []
    for line in open(sys.argv[1]):
        kmin, mat, rep = line.split("|")
        if int(kmin) <= 4:
            t = list(map(int, mat.split()))
            mats.append((t[0:3], t[3:6], t[6:9]))
    mats = [tuple(tuple(r) for r in M) for M in mats]
    L = int(sys.argv[2])
    MT = ((1, 3, 0), (0, 0, 1), (1, 1, 0))
    powers, P = {}, ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    for k in range(1, 4 * L + 4):
        P = mul(MT, P)
        powers[P] = k
    v = (3, 1, 2)

    def parallel(N):
        w = [sum(N[r][c] * v[c] for c in range(3)) for r in range(3)]
        return all(w[i] * v[j] == w[j] * v[i] for i in range(3) for j in range(3))

    print(f"{len(mats)} matrices; single blocks equal to tower powers:",
          sorted(powers[M] for M in mats if M in powers))
    layer = set(mats)
    allp = set(mats)
    for length in range(1, L + 1):
        if length > 1:
            layer = {mul(M, S) for S in layer for M in mats}
            allp |= layer
        found = sorted((powers.get(N, -1), N) for N in allp if parallel(N))
        tags = [p for p, N in found]
        print(f"length <= {length}: {len(allp)} distinct products, {len(found)} preserve the direction (3,1,2);"
              f" tower powers among them: {tags}; non-tower: {[N for p, N in found if p < 0]}")


if __name__ == "__main__":
    main()
