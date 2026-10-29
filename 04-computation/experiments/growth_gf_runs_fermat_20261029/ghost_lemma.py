"""p = 2 ghost lemma: spurious 2-factors = rank of the trivial-sign connecting map on ker/coker of the plus sector.

For an S_k-equivariant integer M whose 2-local glue is the single parity condition u = v (mod 2) between the
trivial (plus, P) and sign (minus, S) sectors (k = 2: the whole glue; k = 3: the 2-part of the index-54 glue,
the standard sector splitting off 2-locally), the glued principal block in the basis u = x, v = x + 2y is
[[P, 0], [B, S]] with B = (S - P)/2, so over F_2
    spur_2 := rank_F2(M) - rank_F2(N) = rank_F2 [[P,0],[B,P]] - 2 rank_F2 P = rank(B: ker P -> coker P).
Checked here against direct mod-2 ranks of M and of the unglued N for every class, n <= 6, k = 2, 3.
"""
import sys, itertools
from runs import tournaments, tour_from_bits, run_matrix, unglued_k3, matmul

def rank2(M):
    rows = [int(''.join(str(x & 1) for x in row), 2) for row in M]
    r = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol - 1, -1, -1):
        bit = 1 << c
        piv = next((i for i in range(r, len(rows)) if rows[i] & bit), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] & bit:
                rows[i] ^= rows[r]
        r += 1
    return r

def sectors_k2(M):
    m = len(M) // 2
    P = [[M[2*a][2*b] + M[2*a][2*b+1] for b in range(m)] for a in range(m)]
    S = [[M[2*a][2*b] - M[2*a][2*b+1] for b in range(m)] for a in range(m)]
    return P, S

tot = 0
for n in (4, 5, 6):
    for k in (2, 3):
        bad = 0
        nonzero = 0
        for bits in tournaments(n):
            A = tour_from_bits(n, bits)
            tuples, M = run_matrix(A, k)
            if k == 3:
                N, P, S, L0, L1 = unglued_k3(M, len(tuples) // 6)
                rN = rank2(N)
            else:
                P, S = sectors_k2(M)
                rN = rank2(P) + rank2(S)
                # sanity: the k=2 sector relation (equivariance) M[(a,b)(c,d)] symmetric under joint swap
            r = len(P)
            B = [[(S[i][j] - P[i][j]) // 2 for j in range(r)] for i in range(r)]
            assert all((S[i][j] - P[i][j]) % 2 == 0 for i in range(r) for j in range(r))
            block = [P[i] + [0] * r for i in range(r)] + [B[i] + P[i] for i in range(r)]
            conn = rank2(block) - 2 * rank2(P)
            spur = rank2(M) - rN
            tot += 1
            if conn != spur:
                bad += 1
            if spur:
                nonzero += 1
        print('n=%d k=%d: classes with spur_2 != connecting rank: %d; classes with spur_2 > 0: %d' % (n, k, bad, nonzero))
print('checked', tot, 'class-observer pairs')
