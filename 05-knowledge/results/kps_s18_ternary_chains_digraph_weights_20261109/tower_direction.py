#!/usr/bin/env python3
"""THM-5070 (5): which port-chain blocks (and short block products) preserve the tower's dominant direction?

In (h, E, O) coordinates the tower block (one new vertex, w = u') is M_T = [[1, x, 0], [0, 0, 1], [1, 1, 0]].
At x = 3 its dominant eigenvector is v = (3, 1, 2) (eigenvalue 2); the other eigenvalues are the primitive cube
roots of unity.  A finite-state "digit system" along the tower would need block products N with N v = lambda v
other than powers of M_T.  This script lists, for the blocks in MATFILE and their products of length <= L, every
distinct matrix N with N v parallel to v, and reports whether N is a power of M_T.

Usage: tower_direction.py MATFILE L
"""
import sys
from fractions import Fraction as Fr


def mul(A, B):
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)) for i in range(3))


def main():
    mats = []
    for line in open(sys.argv[1]):
        t = list(map(int, line.split("|")[0].split()))
        mats.append(((t[0], t[1], t[2]), (t[3], t[4], t[5]), (t[6], t[7], t[8])))
    L = int(sys.argv[2])
    MT = ((1, 3, 0), (0, 0, 1), (1, 1, 0))
    powers = {}
    P = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    for k in range(1, 4 * L + 2):
        P = mul(MT, P)
        powers[P] = k
    v = (3, 1, 2)
    found = {}
    cur = {(): ((1, 0, 0), (0, 1, 0), (0, 0, 1))}
    for length in range(1, L + 1):
        nxt = {}
        for word, S in cur.items():
            for i, M in enumerate(mats):
                S2 = mul(M, S)
                nxt[word + (i,)] = S2
                Sv = [sum(S2[r][c] * v[c] for c in range(3)) for r in range(3)]
                lam = Fr(Sv[0], v[0])
                if all(Sv[r] == lam * v[r] for r in range(3)):
                    found.setdefault(S2, word + (i,))
        cur = nxt
        print(f"products of length <= {length}: {len(found)} distinct direction-preserving matrices")
    for S, word in sorted(found.items(), key=lambda t: powers.get(t[0], 10 ** 9)):
        tag = f"M_T^{powers[S]}" if S in powers else "NOT a tower power"
        print(word, S, tag)


if __name__ == "__main__":
    main()
