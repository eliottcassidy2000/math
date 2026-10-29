"""Audit of THM-4846 (F): plus/sign sectors as [x_V] per/det of the path matrix, and the Jacobi complementary-minor
form of the sign sector, for k = 2 AND k = 3 (every class n <= 6) and k = 4 (n = 6).

Independent of the producer's jacobi_sign_sector.py: the path matrix P is built from the Hamiltonian-path table of
induced subtournaments (P_st = sum_S hp[S][s,t] x^S), not by inverting B; determinants by Leibniz expansion over S_k;
M_k from the subset DP of runlib; sectors S_k[alpha,beta] = sum_g sign(g) M_k[(alpha sorted),(beta permuted by g)].
"""
import itertools
from runlib import tournaments, adj_from_bits, hp_table, run_matrix


def mul(f, g):
    out = {}
    for a, c in f.items():
        for b, d in g.items():
            if a & b:
                continue
            out[a | b] = out.get(a | b, 0) + c * d
    return {k: v for k, v in out.items() if v}


def add(f, g, s=1):
    out = dict(f)
    for k, v in g.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


def sgn(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


def leibniz(Mat, perm_sign=True):
    k = len(Mat)
    tot = {0: 1} if k == 0 else {}
    for p in itertools.permutations(range(k)):
        term = {0: 1}
        for i in range(k):
            term = mul(term, Mat[i][p[i]])
            if not term:
                break
        if term:
            tot = add(tot, term, sgn(p) if perm_sign else 1)
    return tot


def inv_unit(f, n):
    N = {k: -v for k, v in f.items() if k}
    assert f.get(0) == 1
    out, term = {0: 1}, {0: 1}
    for _ in range(n + 1):
        term = mul(term, N)
        if not term:
            break
        out = add(out, term)
    return out


def main():
    total = {2: 0, 3: 0, 4: 0}
    for n in (4, 5, 6):
        full = (1 << n) - 1
        for bits in tournaments(n):
            A = adj_from_bits(n, bits)
            hp = hp_table(A)
            Pm = [[{} for _ in range(n)] for _ in range(n)]
            for S, d in hp.items():
                for (s, t), c in d.items():
                    Pm[s][t] = add(Pm[s][t], {S: c})
            Bm = [[({0: 1} if i == j else {}) for j in range(n)] for i in range(n)]
            for i in range(n):
                for j in range(n):
                    if A[i][j]:
                        Bm[i][j] = add(Bm[i][j], {1 << j: -1})
            detB = leibniz(Bm)
            inv = inv_unit(detB, n)
            for k in (2, 3, 4):
                if k > n - 1 or (k == 4 and n < 6):
                    continue
                tup, M = run_matrix(A, k, hp)
                idx = {t: i for i, t in enumerate(tup)}
                for al in itertools.combinations(range(n), k):
                    for be in itertools.combinations(range(n), k):
                        Ssec = sum(sgn(g) * M[idx[al]][idx[tuple(be[g[i]] for i in range(k))]]
                                   for g in itertools.permutations(range(k)))
                        Psec = sum(M[idx[al]][idx[tuple(be[g[i]] for i in range(k))]]
                                   for g in itertools.permutations(range(k)))
                        sub = [[Pm[a][b] for b in be] for a in al]
                        assert leibniz(sub).get(full, 0) == Ssec
                        assert leibniz(sub, perm_sign=False).get(full, 0) == Psec
                        rows = [v for v in range(n) if v not in be]
                        cols = [v for v in range(n) if v not in al]
                        comp = leibniz([[Bm[r][c] for c in cols] for r in rows])
                        xa = 0
                        for a in al:
                            xa |= 1 << a
                        val = mul(mul({xa: 1}, comp), inv).get(full, 0)
                        assert (-1) ** (sum(al) + sum(be)) * val == Ssec, (n, bits, k, al, be)
                        total[k] += 1
        print('n=%d: plus = [x_V] per P, sign = [x_V] det P and the Jacobi form hold for every class and every (alpha,beta)'
              % n, flush=True)
    print('entries checked by k:', total)


if __name__ == '__main__':
    main()
