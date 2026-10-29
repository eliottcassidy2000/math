"""Audit of THM-4846 (C), (D) at T5 and of the boundary-killer uniqueness for further k."""
import itertools
import flint
from runlib import (circulant, hp_table, run_matrix, plus_sector, snf_diag, snf_summary, fmt, rank_mod_p, unglued,
                    sector, tournaments, adj_from_bits)


def compound(A):
    n = len(A)
    pairs = list(itertools.combinations(range(n), 2))
    return pairs, [[A[b[0]][a[0]] * A[b[1]][a[1]] + A[b[0]][a[1]] * A[b[1]][a[0]] for b in pairs] for a in pairs]


def main():
    out = []
    n = 5
    A = circulant(5, {1, 2})
    pairs, Pc = compound(A)
    pidx = {p: i for i, p in enumerate(pairs)}
    stars = []
    for w in range(5):
        Q = [y for y in range(5) if y != w]
        c = {y: (-1) ** pos for pos, y in enumerate(Q)}
        v = [0] * len(pairs)
        for y in Q:
            v[pidx[tuple(sorted((w, y)))]] = c[y]
        assert all(sum(Pc[i][j] * v[j] for j in range(len(v))) == 0 for i in range(len(Pc)))
        stars.append(v)
        # T5 - w is the strong 4-tournament with Hamiltonian cycle w+1..w+4 and diagonals w+1->w+3, w+2->w+4
        S = [(w + i) % 5 for i in range(1, 5)]
        sub = [[A[a][b] for b in S] for a in S]
        assert sorted(sum(r) for r in sub) == [1, 1, 2, 2]
        hcs = [p for p in itertools.permutations(S[1:]) if all(A[x][y] for x, y in zip((S[0],) + p, p + (S[0],)))]
        assert hcs == [tuple(S[1:])], hcs
        assert A[S[0]][S[2]] and A[S[1]][S[3]]
        assert c[S[0]] == c[S[2]] == -c[S[1]] == -c[S[3]]
    rk = flint.fmpz_mat(Pc).rank()
    span = flint.fmpz_mat(stars).rank()
    rel = [sum((-1) ** w * stars[w][j] for w in range(5)) for j in range(len(pairs))]
    out.append('T5: rank Pc = %d (kernel %d); the five stars are in the kernel, span %d; relation sum (-1)^w v_w = 0: %s; '
               'each T5-w is S4 with unique HC w+1..w+4, diagonals w+1->w+3, w+2->w+4, and c_w = +-1 on the diagonal ends'
               % (rk, 10 - rk, span, all(x == 0 for x in rel)))
    # the ordered three-run observer and the lifted tetrahedron boundaries (THREE_RUN T4)
    hp = hp_table(A)
    tup, M = run_matrix(A, 3, hp)
    tidx = {t: i for i, t in enumerate(tup)}
    killed = True
    for Q in itertools.combinations(range(5), 4):
        v = [0] * len(tup)
        for pos, y in enumerate(Q):
            tri = tuple(x for x in Q if x != y)
            for perm in itertools.permutations(tri):
                v[tidx[perm]] += -1 if pos % 2 else 1
        if any(sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))):
            killed = False
    out.append('T5: M_3 kills every lifted tetrahedron boundary: %s; rank M_3 = %d; Smith M_3 = %s'
               % (killed, flint.fmpz_mat(M).rank(), fmt(snf_summary(snf_diag(M)))))
    # ghost lemma pieces mod 2
    N = unglued(M, 3)
    P = sector(N, 3, [0])
    S = sector(N, 3, [1])
    sets, Pd = plus_sector(A, 3, hp)
    assert P == Pd
    # star indicators on end SETS (triangles): star at w = triangles V-{w,y}
    sidx = {s: i for i, s in enumerate(sets)}
    ind = []
    for w in range(5):
        v = [0] * 10
        for y in range(5):
            if y != w:
                v[sidx[tuple(sorted(set(range(5)) - {w, y}))]] = 1
        ind.append(v)
    kerP2 = 10 - rank_mod_p(P, 2)
    in_ker = all(sum(P[i][j] * v[j] for j in range(10)) % 2 == 0 for v in ind for i in range(10))
    B = [[(S[i][j] - P[i][j]) // 2 for j in range(10)] for i in range(10)]
    # image of the star indicators under B, modulo im P (mod 2): rank of [P | B*ind] - rank P
    Bv = [[sum(B[i][j] * v[j] for j in range(10)) % 2 for i in range(10)] for v in ind]
    Pcols = [[P[i][j] % 2 for i in range(10)] for j in range(10)]
    r_im = rank_mod_p(Pcols + Bv, 2) - rank_mod_p(Pcols, 2)
    out.append('T5 mod 2: rank P = %d, dim ker P = %d, star indicators in ker: %s, their rank %d; '
               'B(ker) has rank %d modulo im P (dim coker = %d)'
               % (rank_mod_p(P, 2), kerP2, in_ker, rank_mod_p(ind, 2), r_im, 10 - rank_mod_p(P, 2)))
    sm = snf_summary(snf_diag(M))
    sn = snf_summary(snf_diag(N))
    out.append('T5: Smith N = %s; spur_2 = %d, spur_3 = %d (factor counts)'
               % (fmt(sn), sum(m for d, m in sn if d == 0 or d % 2 == 0) - sum(m for d, m in sm if d == 0 or d % 2 == 0),
                  sum(m for d, m in sn if d == 0 or d % 3 == 0) - sum(m for d, m in sm if d == 0 or d % 3 == 0)))
    # boundary killers for further (n, k)
    for n2, k2 in ((6, 4), (7, 4), (7, 5), (8, 2)):
        cnt = 0
        for bits in tournaments(n2):
            A2 = adj_from_bits(n2, bits)
            sets2, P2 = plus_sector(A2, k2)
            idx2 = {s: i for i, s in enumerate(sets2)}
            ok = True
            for Q in itertools.combinations(range(n2), k2 + 1):
                v = [0] * len(sets2)
                for pos, y in enumerate(Q):
                    v[idx2[tuple(x for x in Q if x != y)]] += -1 if pos % 2 else 1
                if any(sum(P2[i][j] * v[j] for j in range(len(v))) for i in range(len(P2))):
                    ok = False
                    break
            cnt += ok
        out.append('boundary killers at n=%d, k=%d: %d classes' % (n2, k2, cnt))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
