"""Audit item 2 (independent, exact dense linear algebra with python-flint; no class machinery).

For p in {3, 5, 7, 1, -1} and small K:
  C = 2^(K-1) P_K (integer counting matrix, built from scratch), J = all-ones, so 2^(K-1) N = C - J (N = P - Pi).
  (i)   exact minimal polynomial of C (flint) vs (x - 2^(K-1)) x^(K-1)  [<=> P = Pi + N, N^(K-1) = 0 != N^(K-2)]
  (ii)  exact ranks of (C - J)^n, n = 0..K-1 (modular lower bound + exact distinct-row upper bound), Jordan type of N
        on the whole space (= note's type on W plus one 1x1 block from the constants)
  (iii) rows of C^n: supports are classes; precision histogram -> w_m(n); hitting times h(a); worst/avg TV; Dobrushin
  (iv)  Gram identity (P^n)^T P^n = sum_m w_m(n) E_m exactly (E_m = conditional expectation on classes mod 2^m);
        hence the exact singular spectrum of P^n and ||(P-Pi)^n||^2 = sum_{m>=2} w_m(n) = P(Bin(K-2,1/2) >= n)
"""
import math
import sys
import time
from fractions import Fraction

import flint

Q1 = 2305843009213693951      # 2^61 - 1
Q2 = 4611686018427387847      # 2^62 - 57


def v2(n):
    return (n & -n).bit_length() - 1


def build_C(K, p):
    mod = 1 << K
    odd = list(range(1, mod, 2))
    idx = {a: i for i, a in enumerate(odd)}
    n = len(odd)
    rows = [[0] * n for _ in range(n)]
    for a in odd:
        t = (p * a + 1) % mod
        if t == 0:
            for b in odd:
                rows[idx[a]][idx[b]] = 1
            continue
        v = v2(t)
        r = (t >> v) % (1 << (K - v))
        for b in range(r, mod, 1 << (K - v)):
            rows[idx[a]][idx[b]] = 1 << (K - 1 - v)
    return odd, rows


def nmod_rank(M, q):
    n, m = M.nrows(), M.ncols()
    return flint.nmod_mat(n, m, [int(x) % q for x in M.entries()], q).rank()


def distinct_rows(M):
    return len(set(tuple(int(x) for x in r) for r in M.tolist()))


def analyse(K, p, do_minpoly=True, do_gram=True):
    t0 = time.time()
    odd, rows = build_C(K, p)
    N = len(odd)
    idx_of = {a: i for i, a in enumerate(odd)}
    C = flint.fmpz_mat(rows)
    J = flint.fmpz_mat(N, N, [1] * (N * N))
    # double stochastic in integers
    assert all(sum(r) == N for r in rows)
    assert all(sum(rows[i][j] for i in range(N)) == N for j in range(N))
    res = {}
    if do_minpoly:
        mp = C.minpoly()
        x = flint.fmpz_poly([0, 1])
        target = (x - (1 << (K - 1))) * x ** (K - 1)
        res['minpoly_ok'] = (mp == target)
        assert res['minpoly_ok'], (K, p, mp)
    # powers
    Cn = [flint.fmpz_mat(N, N, [1 if i == j else 0 for i in range(N) for j in range(N)])]
    for s in range(1, K):
        Cn.append(Cn[-1] * C)
    # C^(K-1) = 2^((K-1)(K-2)) J ; C^(K-2) != 2^((K-1)(K-3)) J
    assert Cn[K - 1] == J * (1 << ((K - 1) * (K - 2)))
    # P^(K-2) = Pi  <=>  2^(K-1) C^(K-2) = 2^((K-1)(K-2)) J   (valid also at K = 2)
    assert Cn[K - 2] * (1 << (K - 1)) != J * (1 << ((K - 1) * (K - 2)))
    # (C - J)^n directly (no use of P Pi = Pi): exact
    D = C - J
    Dn = [Cn[0]]
    for s in range(1, K):
        Dn.append(Dn[-1] * D)
    assert Dn[K - 1].is_zero() and not Dn[K - 2].is_zero()
    # ranks of N^n: lower bound mod two primes, upper bound rank(P^n) - 1 <= distinct rows of C^n - 1 (rows of N^n
    # lie in rowspace(P^n) intersected with the sum-zero hyperplane, and rowspace(P^n) contains pi, which is not sum-zero)
    rN, rP = [], []
    for s in range(K):
        lowN = max(nmod_rank(Dn[s], Q1), nmod_rank(Dn[s], Q2))
        lowP = max(nmod_rank(Cn[s], Q1), nmod_rank(Cn[s], Q2))
        # upper bound: every row of C^s is constant on the classes mod 2^(K-s) (exact check), so the row space lies in
        # the span of the 2^(K-1-s) indicators of those classes  (distinct rows would overcount: nested classes)
        Ls = Cn[s].tolist()
        step = 1 << (K - s)
        for r in Ls:
            for j, b in enumerate(odd):
                j2 = idx_of[(b + step) % (1 << K)]
                assert r[j] == r[j2]
        upP = (1 << (K - 1 - s))
        upN = upP - 1 if s >= 1 else N
        assert lowP == upP == (1 << (K - 1 - s)), (K, p, s, lowP, upP)
        if s >= 1:
            assert lowN == upN, (K, p, s, lowN, upN)
        else:
            assert lowN == N
        rN.append(lowN)
        rP.append(lowP)
    rN.append(0)
    ge = {s: rN[s - 1] - rN[s] for s in range(1, K)}
    exact_blocks = {s: ge[s] - ge.get(s + 1, 0) for s in range(1, K)}
    # note's claim on W: one block K-1 and 2^(K-2-s) blocks of size s ; whole space: + one extra 1x1 block
    claim = {s: (2 ** (K - 2 - s) if s <= K - 2 else 0) for s in range(1, K)}
    claim[K - 1] += 1
    claim[1] += 1
    assert exact_blocks == claim, (K, p, exact_blocks, claim)
    res['jordan_full'] = exact_blocks
    # (iii) rows of C^n: supports are classes (residues congruent mod 2^m), values constant
    prec_hist = []
    hit = [None] * N
    worst_tv, avg_tv, dob = [], [], []
    for s in range(K):
        L = Cn[s].tolist()
        tot = 1 << ((K - 1) * s)
        hist = {}
        tvs = []
        for i, r in enumerate(L):
            sup = [odd[j] for j, x in enumerate(r) if x != 0]
            vals = set(int(r[j]) for j, x in enumerate(r) if x != 0)
            assert len(vals) == 1
            size = len(sup)
            assert size & (size - 1) == 0
            m = K - (size.bit_length() - 1)               # class of precision m has 2^(K-m) odd residues
            c0 = sup[0]
            assert all((b - c0) % (1 << m) == 0 for b in sup)   # support is the full class c0 mod 2^m
            assert int(list(vals)[0]) * size == tot
            hist[m] = hist.get(m, 0) + 1
            if m == 1 and hit[i] is None:
                hit[i] = s
            tvs.append(1 - Fraction(size, N))
        prec_hist.append(hist)
        worst_tv.append(max(tvs))
        avg_tv.append(sum(tvs) / N)
        # Dobrushin exactly: max over pairs of TV between rows; rows are uniform laws on classes, so TV(row_i,row_j)
        # = 1 - |Ci cap Cj|/max(|Ci|,|Cj|) ; compute by brute force over distinct rows
        distinct = {}
        for r in L:
            key = tuple(int(x) for x in r)
            distinct[key] = r
        dk = list(distinct.keys())
        best = Fraction(0)
        if K <= 8:          # brute force over all pairs of distinct rows
            for i in range(len(dk)):
                for j in range(i + 1, len(dk)):
                    a, b = dk[i], dk[j]
                    d = sum(abs(x - y) for x, y in zip(a, b))
                    best = max(best, Fraction(d, 2 * tot))
        else:               # rows are uniform on classes (checked above): TV = 1 iff some two supports are disjoint
            sups = [frozenset(j for j, x in enumerate(key) if x) for key in dk]
            for i in range(len(sups)):
                for j in range(i + 1, len(sups)):
                    if not (sups[i] & sups[j]):
                        best = Fraction(1)
                        break
                    small, big = (sups[i], sups[j]) if len(sups[i]) <= len(sups[j]) else (sups[j], sups[i])
                    best = max(best, 1 - Fraction(len(small & big), len(big)))
                if best == 1:
                    break
        dob.append(best)
    assert all(h is not None for h in hit)
    from collections import Counter
    hc = Counter(hit)
    assert all(hc[h] == 2 * math.comb(K - 2, h - 1) for h in range(1, K)) and hc[0] == 0
    assert Fraction(sum(hit), N) == Fraction(K, 2)
    assert worst_tv == [1 - Fraction(1, 2 ** (K - 1 - s)) for s in range(K)]
    assert dob == [Fraction(1)] * (K - 1) + [Fraction(0)], dob
    # level-only multiplicities: each class of precision m >= 2 appears the same number of times among the rows
    for s in range(1, K):
        L = Cn[s].tolist()
        mult = {}
        for r in L:
            key = tuple(int(x) for x in r)
            mult[key] = mult.get(key, 0) + 1
        bylevel = {}
        for key, c in mult.items():
            size = sum(1 for x in key if x)
            m = K - (size.bit_length() - 1)
            bylevel.setdefault(m, set()).add(c)
        assert all(len(v) == 1 for m, v in bylevel.items() if m >= 2), (K, p, s, bylevel)
        assert sorted(bylevel) == list(range(1, K - s + 1))
        # number of distinct classes at each level = 2^(m-1)
        cnt = Counter()
        for key in mult:
            size = sum(1 for x in key if x)
            cnt[K - (size.bit_length() - 1)] += 1
        assert all(cnt[m] == 2 ** (m - 1) for m in range(1, K - s + 1))
    res['hit'] = [hc[h] for h in range(1, K)]
    res['worst_tv'] = worst_tv
    res['avg_tv'] = avg_tv
    # (iv) Gram identity and l^2 norms
    if do_gram:
        l2 = []
        for s in range(K):
            w = {m: Fraction(c, N) for m, c in prec_hist[s].items()}
            # build sum_m w_m E_m scaled by S = 2^(2(K-1)s) * 2^(K-1) to integers:
            # (C^s)^T C^s = 2^(2(K-1)s) (P^s)^T P^s ;  E_m(a,b) = [a = b mod 2^m] / 2^(K-m)
            G = Cn[s].transpose() * Cn[s]
            scale = 1 << (2 * (K - 1) * s)
            ok = True
            Gl = G.tolist()
            for i, a in enumerate(odd):
                for j, b in enumerate(odd):
                    val = Fraction(0)
                    for m, wm in w.items():
                        if (a - b) % (1 << m) == 0:
                            val += wm * Fraction(1, 1 << (K - m))
                    if Fraction(int(Gl[i][j]), scale) != val:
                        ok = False
                        break
                if not ok:
                    break
            assert ok, (K, p, s)
            # eigenvalues of sum_m w_m E_m on level-m Haar wavelets: tail sums sum_{m' >= m} w_m'
            norm2_W = sum(wm for m, wm in w.items() if m >= 2)
            binom_tail = sum(Fraction(math.comb(K - 2, j), 2 ** (K - 2)) for j in range(s, K - 1))
            assert norm2_W == binom_tail, (K, p, s, norm2_W, binom_tail)
            l2.append(norm2_W)
        res['l2'] = l2
    res['secs'] = time.time() - t0
    return res


def main():
    print('A2: exact dense checks (python-flint), no class machinery')
    for p in (3, 5, 7, 1, -1):
        Kmax = 11 if p == 3 else 9
        for K in range(2, Kmax + 1):
            gram = K <= (8 if p == 3 else 7)
            mp = K <= 10
            r = analyse(K, p, do_minpoly=mp, do_gram=gram)
            line = (f'p = {p:2d} K = {K:2d}: ' + ('minpoly(C) = (x - 2^(K-1)) x^(K-1) exactly; ' if mp else '') +
                    f'C^(K-1) = 2^((K-1)(K-2)) J, C^(K-2) != ..., (C - J)^(K-1) = 0 != (C - J)^(K-2); '
                    f'rank P^n = 2^(K-1-n) and rank N^n = 2^(K-1-n) - 1 exact; Jordan type of N on the whole space '
                    f'{r["jordan_full"]} (= note + one 1x1 block); hitting counts {r["hit"]}; worst TV, Dobrushin as claimed')
            if gram:
                line += f'; Gram identity exact; ||(P-Pi)^n||^2 = {[str(x) for x in r["l2"]]} = binomial tails'
            line += f'  [{r["secs"]:.1f}s]'
            print(line)
            sys.stdout.flush()
    print('A2 ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
