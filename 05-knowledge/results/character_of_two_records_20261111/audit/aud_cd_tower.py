"""AUDIT: is the Cayley-Dickson sign tournament CD_k the doubling tower of THM-448 (as the note says in section 2)?
Compare CD_k (audit's own construction, rule (a,b)(c,d) = (ac - conj(d) b, d a + b conj(c))) with the double
D(CD_(k-1)) (tournament form of [[M, M], [-M^T, M^T]]) and with its converse, by isomorphism invariants:
the number of transitive 5-subsets and the sorted vertex profile of transitive 4-subsets through a vertex.
Also: anti-automorphisms of CD_5 (backtracking with degree-profile pruning is too weak here, so only the
invariant comparison T vs converse is reported)."""
import itertools, sys
from aud_common import cd_tournament, double_tournament, is_doubly_regular


def tt_counts(A):
    n = len(A)
    out = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
    # transitive k-subsets: count linear orders: number of directed "transitive chains" v1 > v2 > ... (v_i beats all later)
    # TT4 through each vertex and total TT5 via chains: choose v1, then v2 in out(v1), v3 in out(v1)&out(v2), ...
    tt5 = 0
    tt4_vertex = [0] * n
    for a in range(n):
        A1 = out[a]
        x = A1
        while x:
            b = (x & -x).bit_length() - 1
            x &= x - 1
            A2 = A1 & out[b]
            y = A2
            while y:
                c = (y & -y).bit_length() - 1
                y &= y - 1
                A3 = A2 & out[c]
                z = A3
                while z:
                    d = (z & -z).bit_length() - 1
                    z &= z - 1
                    for v in (a, b, c, d):
                        tt4_vertex[v] += 1
                    tt5 += bin(A3 & out[d]).count("1")
    return tt5, tuple(sorted(tt4_vertex))


for k in (4, 5, 6):
    C = cd_tournament(k)
    D = double_tournament(cd_tournament(k - 1))
    Dc = [[D[j][i] for j in range(len(D))] for i in range(len(D))]
    Cc = [[C[j][i] for j in range(len(C))] for i in range(len(C))]
    tc, td, tdc, tcc = tt_counts(C), tt_counts(D), tt_counts(Dc), tt_counts(Cc)
    print("k = %d (order %d): CD_k doubly regular: %s; D(CD_(k-1)) doubly regular: %s" % (k, len(C), is_doubly_regular(C), is_doubly_regular(D)))
    print("   transitive 5-subsets: CD_k %d, converse of CD_k %d, D(CD_(k-1)) %d, its converse %d" % (tc[0], tcc[0], td[0], tdc[0]))
    print("   TT4-through-vertex profiles equal: CD_k vs D(CD_(k-1)): %s; CD_k vs converse of D: %s; CD_k vs its own converse: %s" % (
        tc[1] == td[1], tc[1] == tdc[1], tc[1] == tcc[1]))
    print("   distinct values of TT4 through a vertex: CD_k %s ; D(CD_(k-1)) %s" % (sorted(set(tc[1])), sorted(set(td[1]))))
