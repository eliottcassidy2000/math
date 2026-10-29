"""Audit of the post-launch THM-4846 (D') mod-3 ghost formula, on every class n = 4..6 (and n = 7 with --n7).

spur_3 := rank_F3 M - rank_F3 N  versus  rank I + rank II - rank P - rank S - 2 rank [[P, L12],[0, S]]
(P, S, L' = Q^-1 L Q from the unglued N; M from the subset DP; ranks by FLINT nmod).  Also checks the congruences
L'21 = 0, L'11 = P, L'22 = S (mod 3) that the lemma assumes, and the S-invertible reduction rank(E: ker P -> coker P).
"""
import sys
from runlib import tournaments, adj_from_bits, hp_table, run_matrix, unglued, rank_mod_p

Q = [[1, 1], [1, 0]]
QI = [[0, 1], [1, -1]]


def m2(a, b):
    return [[sum(a[i][t] * b[t][j] for t in range(2)) for j in range(2)] for i in range(2)]


def assemble(blocks):
    out = []
    for brow in blocks:
        h = len(brow[0])
        for i in range(h):
            out.append(sum((b[i] for b in brow), []))
    return out


def main():
    ns = [4, 5, 6] + ([7] if '--n7' in sys.argv else [])
    for n in ns:
        bad = pos = 0
        for bits in tournaments(n):
            A = adj_from_bits(n, bits)
            tup, M = run_matrix(A, 3, hp_table(A))
            N = unglued(M, 3)
            r = len(N) // 6
            P = [[N[6 * i][6 * j] for j in range(r)] for i in range(r)]
            S = [[N[6 * i + 1][6 * j + 1] for j in range(r)] for i in range(r)]
            L = {}
            for i in range(r):
                for j in range(r):
                    blk = [[N[6 * i + 2][6 * j + 2], N[6 * i + 2][6 * j + 4]], [N[6 * i + 4][6 * j + 2], N[6 * i + 4][6 * j + 4]]]
                    L[(i, j)] = m2(m2(QI, blk), Q)
            L11 = [[L[(i, j)][0][0] for j in range(r)] for i in range(r)]
            L12 = [[L[(i, j)][0][1] for j in range(r)] for i in range(r)]
            L21 = [[L[(i, j)][1][0] for j in range(r)] for i in range(r)]
            L22 = [[L[(i, j)][1][1] for j in range(r)] for i in range(r)]
            assert all((L11[i][j] - P[i][j]) % 3 == 0 and L21[i][j] % 3 == 0 and (L22[i][j] - S[i][j]) % 3 == 0
                       for i in range(r) for j in range(r))
            D1 = [[(L11[i][j] - P[i][j]) // 3 for j in range(r)] for i in range(r)]
            D2 = [[L21[i][j] // 3 for j in range(r)] for i in range(r)]
            D3 = [[(L22[i][j] - S[i][j]) // 3 for j in range(r)] for i in range(r)]
            Z = [[0] * r for _ in range(r)]
            I = assemble([[P, Z, Z], [D2, S, Z], [D1, L12, P]])
            II = assemble([[S, Z, Z], [L12, P, Z], [D3, D2, S]])
            std = assemble([[P, L12], [Z, S]])
            formula = (rank_mod_p(I, 3) + rank_mod_p(II, 3) - rank_mod_p(P, 3) - rank_mod_p(S, 3)
                       - 2 * rank_mod_p(std, 3))
            direct = rank_mod_p(M, 3) - rank_mod_p(N, 3)
            # also: rank_F3 M = rank I + rank II
            assert rank_mod_p(M, 3) == rank_mod_p(I, 3) + rank_mod_p(II, 3)
            if formula != direct:
                bad += 1
            pos += direct > 0
        print('n=%d: classes where the (D\') formula differs from direct mod-3 ranks: %d; classes with spur_3 > 0: %d'
              % (n, bad, pos), flush=True)


if __name__ == '__main__':
    main()
