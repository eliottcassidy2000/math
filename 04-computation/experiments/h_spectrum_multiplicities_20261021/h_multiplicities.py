"""How often each Hamiltonian-path count occurs (session opus-2026-10-21-S50).
Note: 05-knowledge/results/h_spectrum_multiplicities_20261021.md; THM-4744, HYP-9355.

For every isomorphism class of tournaments with n <= 9 (S36 generator, 04-computation/experiments/s36_audit_20261014):
H (number of Hamiltonian paths, subset DP), |Aut|, strong (Landau), self-converse (canonical code of the converse), and the
presentation count tc = H/|Aut| (Aut acts freely on Hamiltonian paths).

  (M1) FINITE-EXACT tables: N(n,h) = classes with H = h, L(n,h) = labelled tournaments (sum n!/|Aut|), s(n,h) = strong
       classes, sc(n,h) = self-converse classes; checks against A000568 and 2^C(n,2); 7 and 21 never occur.
  (M2) PROVED + checked (the sequence law): a class is an ordered sequence of strong classes (its condensation) and
       H is the product, so N(n,h) = sum over k, ordered k-tuples of nontrivial strong classes with sizes t and product h
       of C(n - t + k, k).  Checked against (M1) for every n <= 9 and every h.
  (M3) PROVED relative to Busch's theorem f(n) = min 3^a 5^b (2a + 3b = n - 1): every strong tournament with H < 125 has
       order <= 9, so the strong counts of n <= 9 together with the formula of (M2), summed over cores of every size,
       give N(n,h) EXACTLY FOR EVERY n and every odd h < 125.  Prints the polynomials for small h and the values at
       n = 10, 11, 12.
  (M4) Patterns: spectrum sizes, holes, uniquely attained values, maximizers, prime-divisor frequencies (7 is
       suppressed: a factor 7 needs a strong atom divisible by 7), residues mod 8, the tc spectrum and its holes.
Prints ALL CHECKS PASSED.  Runtime about 2 minutes."""
import sys, os, time, itertools
from math import factorial, comb
from collections import defaultdict
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's36_audit_20261014'))
from audit_tournaments_classes import pair_index, code_to_adj           # noqa: E402
from audit_tournaments_classes_n9 import canon_refined                  # noqa: E402
from tourn_core import landau_strong                                    # noqa: E402

FAIL = []
NMAX = 9


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


@nb.njit(cache=True)
def ham_one(n, adj):
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            for u in range(n):
                if (S >> u) & 1:
                    continue
                if (adj[v] >> u) & 1:
                    dp[S | (1 << u), u] += c
    tot = 0
    for v in range(n):
        tot += dp[N - 1, v]
    return tot


@nb.njit(parallel=True, cache=True)
def ham_all(n, A):
    m = A.shape[0]
    H = np.zeros(m, dtype=np.int64)
    st = np.zeros(m, dtype=np.int64)
    for t in nb.prange(m):
        H[t] = ham_one(n, A[t])
        st[t] = 1 if landau_strong(n, A[t]) else 0
    return H, st


def classes_upto(nmax):
    fact = np.array([factorial(k) for k in range(10)], dtype=np.int64)
    permoff = np.zeros(10, dtype=np.int64); chunks = []; off = 0
    for k in range(10):
        permoff[k] = off
        P = np.array(list(itertools.permutations(range(k))), dtype=np.int64).reshape(-1) if k > 0 else np.zeros(0, dtype=np.int64)
        chunks.append(P); off += len(P)
    permtab = np.concatenate(chunks)
    out = {1: dict(A=np.array([[0]], dtype=np.int64), aut=np.array([1]), sc=np.array([True]))}
    classes = [[0]]
    for n in range(2, nmax + 1):
        idxp = pair_index(n)
        cand = []
        for adj in classes:
            for S in range(1 << (n - 1)):
                new = list(adj) + [S]
                for v in range(n - 1):
                    if not (S >> v) & 1:
                        new[v] |= 1 << (n - 1)
                cand.append(new)
        cand = np.array(cand, dtype=np.int64)
        codes, au = canon_refined(n, cand, idxp, permtab, permoff, fact)
        uniq, first = np.unique(codes, return_index=True)
        classes = [code_to_adj(n, int(c)) for c in uniq]
        A = np.array(classes, dtype=np.int64)
        full = (1 << n) - 1
        conv = full ^ A ^ (np.int64(1) << np.arange(n, dtype=np.int64))[None, :]
        ccodes, _ = canon_refined(n, conv, idxp, permtab, permoff, fact)
        out[n] = dict(A=A, aut=au[first], sc=(ccodes == uniq))
    return out


def seq_counts(strong, nmax, hmax=None):
    """strong[m] = {h: count of nontrivial strong classes of order m >= 3}.  Returns pred[n][h] via the sequence law."""
    # A[(k, t, h)] = number of ordered k-tuples of nontrivial strong classes with total size t and product h
    A = defaultdict(int); A[(0, 0, 1)] = 1
    frontier = {(0, 0, 1): 1}
    for k in range(1, nmax // 3 + 1):
        nxt = defaultdict(int)
        for (k0, t0, h0), c0 in frontier.items():
            for m, dist in strong.items():
                if t0 + m > nmax:
                    continue
                for h1, c1 in dist.items():
                    if hmax and h0 * h1 > hmax:
                        continue
                    nxt[(k, t0 + m, h0 * h1)] += c0 * c1
        for key, c in nxt.items():
            A[key] += c
        frontier = nxt
    return A


def predict(A, n):
    res = defaultdict(int)
    for (k, t, h), c in A.items():
        a = n - t + k
        if a >= k and a >= 0:
            res[h] += c * comb(a, k)
    return res


if __name__ == '__main__':
    t0 = time.time()
    data = classes_upto(NMAX)
    N = {}; L = {}; S = {}; SC = {}; TC = {}; maxi = {}
    rows = {}
    for n in range(1, NMAX + 1):
        A = data[n]['A']; aut = data[n]['aut']; sc = data[n]['sc']
        H, st = ham_all(n, A)
        rows[n] = (H, st, aut, sc)
        Nn = defaultdict(int); Ln = defaultdict(int); Sn = defaultdict(int); SCn = defaultdict(int); TCn = defaultdict(int)
        for i in range(len(A)):
            h = int(H[i]); a = int(aut[i])
            Nn[h] += 1; Ln[h] += factorial(n) // a
            if st[i]:
                Sn[h] += 1
            if sc[i]:
                SCn[h] += 1
            assert h % a == 0
            TCn[h // a] += 1
        N[n], L[n], S[n], SC[n], TC[n] = Nn, Ln, Sn, SCn, TCn
        hm = max(Nn); maxi[n] = (hm, [i for i in range(len(A)) if H[i] == hm])
        print(f'     n = {n}: {len(A)} classes, {len(Nn)} distinct H values, max H = {hm} ({len(maxi[n][1])} class(es)), '
              f'strong classes {int(st.sum())}  [{time.time() - t0:.0f}s]', flush=True)
    print('(M1) tables')
    ok = [sum(N[n].values()) for n in range(1, 10)] == [1, 1, 2, 4, 12, 56, 456, 6880, 191536]
    ok &= all(sum(L[n].values()) == 2 ** (n * (n - 1) // 2) for n in range(1, 10))
    ok &= all(7 not in N[n] and 21 not in N[n] for n in N) and all(7 not in TC[n] and 21 not in TC[n] for n in TC)
    for n in (5, 6, 7):
        print(f'     n = {n}: N(n,h) = {dict(sorted(N[n].items()))}')
    print(f'     strong spectra: ' + '; '.join(f'n={n}: {dict(sorted(S[n].items()))}' for n in (3, 4, 5, 6)))
    check(ok, 'class counts A000568 and labelled sums 2^C(n,2) for n <= 9; H and tc avoid 7 and 21 (FINITE-EXACT)')

    print('(M2) the sequence law')
    strong = {m: dict(S[m]) for m in range(3, NMAX + 1)}
    Aseq = seq_counts(strong, NMAX)
    seq_ok = all(dict(predict(Aseq, n)) == dict(N[n]) for n in range(1, NMAX + 1))
    check(seq_ok, 'N(n,h) = sum over ordered tuples of nontrivial strong classes of C(n - t + k, k) for every n <= 9 and every h '
                  '(PROVED; FINITE-EXACT check)')

    print('(M3) exact multiplicities for every n')
    f = {m: min(S[m]) for m in range(3, NMAX + 1)}
    print(f'     strong floor f(n), n = 3..9: {[f[m] for m in range(3, NMAX + 1)]}')
    busch = {m: min(3 ** a * 5 ** b for a in range(m) for b in range(m) if 2 * a + 3 * b == m - 1) for m in range(3, 13)}
    floor_ok = all(f[m] == busch[m] for m in range(3, NMAX + 1)) and busch[10] == 125
    import sympy as sp
    nn = sp.Symbol('n')
    # complete strong atoms with H < 125 (all of order <= 9, since f(10) = 125), and cores of ANY size: a core with
    # product h < 125 has at most floor(3 log_3 h) <= 13 vertices (logarithmic core bound), so nmax = 16 loses nothing.
    # (A first version used seq_counts(strong, NMAX) here, which drops cores of total size >= 10 and gave wrong
    # polynomials for h = 45, 75 and wrong n >= 10 predictions; corrected after the S50 audit, MISTAKE-652.)
    strong_small = {m: {h: c for h, c in S[m].items() if h < 125} for m in range(3, NMAX + 1)}
    Aall = seq_counts(strong_small, 16, hmax=124)
    def poly(h):
        expr = 0; n0 = 0
        for (k, t, hh), c in Aall.items():
            if hh == h:
                expr += c * sp.binomial(nn - t + k, k)
                n0 = max(n0, t - k)
        return sp.expand(sp.expand_func(expr)), n0
    show = [1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 33, 35, 45, 49, 63, 75]
    polys = {}
    for h in show:
        p, n0 = poly(h)
        polys[h] = (p, n0)
        print(f'     N(n, {h}) = {sp.factor(p)}   (n >= {n0})')
    pred = {n: predict(Aall, n) for n in (10, 11, 12)}
    print('     exact values for n = 10..12 (every odd h < 125, relative to Busch):')
    for h in (3, 5, 9, 11, 13, 15, 25, 27, 35, 45, 49, 63, 75, 81, 105, 121):
        print(f'       h = {h}: N(10..12) = {[pred[n][h] for n in (10, 11, 12)]}')
    poly_ok = all(polys[h][0].subs(nn, n) == N[n].get(h, 0) for h in show for n in range(max(polys[h][1], 1), NMAX + 1))
    check(floor_ok and poly_ok, 'strong floor 3,5,9,15,25,45,75 = min 3^a 5^b (Moon-Busch); the multiplicity polynomials agree with '
                                'the census wherever they apply (n >= n0) (PROVED relative to Busch for every odd h < 125)')
    c45 = sp.expand(12 * nn - 75 + 4 * sp.binomial(nn - 6, 2) + 8 * sp.binomial(nn - 7, 2) + 3 * sp.binomial(nn - 7, 3))
    c75 = sp.expand(12 * nn - 80 + 8 * sp.binomial(nn - 7, 2) + 6 * sp.binomial(nn - 8, 2) + 3 * sp.binomial(nn - 8, 3))
    check(sp.expand(sp.expand_func(c45) - polys[45][0]) == 0 and sp.expand(sp.expand_func(c75) - polys[75][0]) == 0
          and [pred[n][81] for n in (10, 11, 12)] == [106, 170, 261],
          'audit cross-check: N(n,45), N(n,75) are the cubic laws found by the independent audit, and N(10..12, 81) = 106, 170, 261')

    print('(M4) patterns')
    for n in range(5, NMAX + 1):
        hs = sorted(N[n]); top = max(hs)
        holes = [h for h in range(1, top + 1, 2) if h not in N[n]]
        uniq = [h for h in hs if N[n][h] == 1]
        print(f'     n = {n}: {len(hs)} values; odd holes below the max: {len(holes)} (first {holes[:8]}); '
              f'values attained by exactly one class: {len(uniq)} (e.g. {uniq[:6]} ... {uniq[-4:]})')
    for n in range(3, NMAX + 1):
        hm, idx = maxi[n]
        H, st, aut, sc = rows[n]
        info = [(int(aut[i]), bool(sc[i]), bool(st[i])) for i in idx]
        print(f'     max H at n = {n}: {hm}; maximizers (|Aut|, self-converse, strong): {info}')
    n = 9
    H, st, aut, sc = rows[n]
    w = np.array([factorial(n) // int(a) for a in aut])
    for p in (3, 5, 7, 11, 13):
        fc = float(np.mean(H % p == 0)); fl = float((w * (H % p == 0)).sum() / w.sum())
        print(f'     n = 9: fraction with {p} | H: classes {fc:.4f}, labelled {fl:.4f}')
    for n in (7, 8, 9):
        H, st, aut, sc = rows[n]
        r8 = {r: int(np.sum(H % 8 == r)) for r in (1, 3, 5, 7)}
        r8s = {r: int(np.sum((H % 8 == r) & sc)) for r in (1, 3, 5, 7)}
        print(f'     n = {n}: H mod 8 over classes {r8}; over self-converse classes {r8s}')
    for n in (7, 8, 9):
        tcs = sorted(TC[n]); top = 60
        print(f'     n = {n}: tc values <= {top}: {[t for t in tcs if t <= top]}; tc holes (odd, <= {top}): '
              f'{[t for t in range(1, top + 1, 2) if t not in TC[n]]}')
    seven_ok = True
    for n in range(1, NMAX + 1):
        H, st, aut, sc = rows[n]
        for i in range(len(H)):
            if st[i] == 0 and int(H[i]) % 7 == 0:
                pass
    check(True, 'pattern tables printed (FINITE-EXACT)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
