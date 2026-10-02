#!/usr/bin/env python3
"""procgen_edim2_20261001_run.py -- runner for the edim2 lane (2026-10-01).

Re-verifies every computational claim of
05-knowledge/results/procgen_edim2_20261001_growth_and_uniform_bounds.md
(Allikvere, arXiv:2608.09983, Open Problems 2-4; THM-4525 settled Open Problem 1).

Sections
  S1  cells, pair types and orbit counts (brute force)                                      [VERIFIED]
  S2  forest lemma for general q: brute force over all 2^16 subsets of Q_4                    [VERIFIED]
  S3  Lemma A (Fourier atom bound), Lemma A2, Lemma B (binomial bounds), Lemma S (star lemma),
      conditional mean/variance identities, path/zigzag forests and their cell bounds         [checks of PROVED lemmas]
  S4  Lemma FH (Fourier-Hoelder multi-forest lemma): exact check on random small level graphs   [check of PROVED lemma]
  S5  Theorem B: exact entropy (L5) values and the constant c* = (3 ln2/(2 sqrt2))^(2/3)        [illustration]
  S6  Theorem D: Lemma PI, crs/par cell bounds of the closed-form forests, endpoint structure,
      interval evaluation at d = 17, 18, 19, crude bounds E1..E5 for n >= 19                    [PROVED; numeric inputs]
  S7  certificates: explicit resolving sets for Q_6 .. Q_16, exact verification               [VERIFIED]
  S8  certified union bounds at q = 1/2 with Lemma FH: U_10 < 1, U_11                          [FINITE-EXACT (intervals)]
  S9  certified sparse union bounds edim_m(Q_d) <= M_d (Lemma FH, intervals)                   [FINITE-EXACT (intervals)]
  S10 Monte Carlo first moments for density 1/2 (fixed seed)                                   [EMPIRICAL]
  S11 Q_7: exhaustive search (complete normal form, C), validated against Python counts and
      against THM-4525 on Q_6 and Burnside counts; no resolving set of size <= 10            [FINITE-EXACT]
Prints to stdout only; ends with ALL CHECKS PASSED.
Usage: python3 -u procgen_edim2_20261001_run.py [--quick] [--q7k11]
       (--quick: smaller S4/S9/S10, Q_7 up to k = 9 only; --q7k11: also the (unfinished in this lane) Q_7 k = 11 search)
"""
import os, sys, time, random, itertools
from fractions import Fraction as Fr
from math import comb, log, sqrt, pi, exp, lgamma, ceil
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_edim2_20261001_lib as L

NCHECK = 0
def ok(cond, msg):
    global NCHECK
    if not cond:
        print('[FAIL]', msg, flush=True)
        sys.exit(1)
    NCHECK += 1
    print('[OK]', msg, flush=True)
def section(t):
    print('\n==== %s ====' % t, flush=True)
QUICK = '--quick' in sys.argv
Q7K = ((9, 1), (10, 2), (11, 2)) if '--q7k11' in sys.argv else ((9, 1), (10, 2))
T0 = time.time()
rng = random.Random(20261001)
iv = L.iv_ctx(160)

# ======================================================================================= S1
section('S1 cells, pair types, orbit counts')
good = all(L.brute_cells(d, *L.representative(t, d, h)) == L.cells(t, d, h)
           for d in range(3, 8) for (t, h) in L.pair_types(d))
ok(good, 'cell formulas (paper Prop. 14) equal brute force for every pair type, 3 <= d <= 7')
def classify_pairs(d):
    E = [(u, u | (1 << i)) for u in range(1 << d) for i in range(d) if not (u >> i) & 1]
    cnt = {}
    for x in range(len(E)):
        for y in range(x + 1, len(E)):
            (u, v), (w, z) = E[x], E[y]
            par = (u ^ v) == (w ^ z)
            h = min(bin(a ^ b).count('1') for a in (u, v) for b in (w, z))
            key = ('par' if par else 'crs', h)
            cnt[key] = cnt.get(key, 0) + 1
    return cnt
for d in (4, 5):
    c = classify_pairs(d)
    ok(c == {(t, h): L.type_count(t, d, h) for (t, h) in L.pair_types(d)},
       'Q_%d: every unordered edge pair classified by (parallel?, d(e,f)); counts = d 2^(d-2) C(n,h), C(d,2) 2^d C(n-1,h)' % d)
ok(all(sum(L.type_count(t, d, h) for (t, h) in L.pair_types(d)) == (d * 2 ** (d - 1)) * (d * 2 ** (d - 1) - 1) // 2
       for d in range(2, 60)), 'sum of type counts = C(d 2^(d-1), 2) for 2 <= d < 60')
# probabilistic description: row sums 2 C(n,a) (A ~ Bin(n,1/2)), column sums 2 C(n,b)
ok(all(all(sum(M[a]) == 2 * comb(d - 1, a) and sum(M[r][a] for r in range(d)) == 2 * comb(d - 1, a) for a in range(d))
       for d in range(3, 30) for (t, h) in L.pair_types(d) for M in [L.cells(t, d, h)]),
   'A = d(e,w) and B = d(f,w) are Bin(n,1/2) for uniform w (row/column sums 2C(n,a)), all types, d < 30')

# ======================================================================================= S2
section('S2 forest lemma for general q (brute force, Q_4, all 2^16 subsets)')
import numpy as np
d = 4; V = 1 << d
masks = np.arange(1 << V, dtype=np.int64)
pcm = np.zeros(1 << V, dtype=np.int64)
for w in range(V): pcm += (masks >> w) & 1
def dist_e(edge, w):
    u, v = edge
    return min(bin(u ^ w).count('1'), bin(v ^ w).count('1'))
allok = True; ncase = 0
for q in (Fr(1, 2), Fr(1, 4), Fr(1, 10)):
    wts = np.array([float(q ** k * (1 - q) ** (V - k)) for k in range(V + 1)])[pcm]
    for (t, h) in L.pair_types(d):
        e, f = L.representative(t, d, h)
        He = np.zeros((1 << V, d), dtype=np.int64); Hf = np.zeros((1 << V, d), dtype=np.int64)
        for w in range(V):
            bit = (masks >> w) & 1
            He[:, dist_e(e, w)] += bit; Hf[:, dist_e(f, w)] += bit
        Pex = wts[np.all(He == Hf, axis=1)].sum()
        M = L.cells(t, d, h); Ne = L.Nmat(M)
        F = L.kruskal(Ne, d, lambda x: Ne[x])
        bound = 1.0
        for (a, b) in F: bound *= float(L.maxatom_exact(M[a][b], M[b][a], q))
        allok &= Pex <= bound * (1 + 1e-12); ncase += 1
ok(allok, 'P(H_e = H_f) <= prod over a max-weight forest of exact max atoms, all %d (q, type) cases' % ncase)

# ======================================================================================= S3
section('S3 atom lemmas, binomial bounds, star lemma, forests')
bad = 0; nt = 0
for q in (Fr(1, 2), Fr(1, 3), Fr(1, 7), Fr(1, 20), Fr(2, 5)):
    for N in list(range(1, 50)) + [80, 120]:
        for n1 in sorted({0, N // 3, N // 2, N}):
            A = float(L.maxatom_exact(n1, N - n1, q)); x = float(N * q * (1 - q))
            if A > ((1 + 1 / (4 * x)) / sqrt(2 * pi * x) + exp(-x) / 2) * (1 + 1e-12): bad += 1
            if A < (12 * x + 1) ** -0.5 * (1 - 1e-12): bad += 1
            # e^{-x} I_0(x) bound (interval implementation)
            if A > float(L.expI0_iv(iv, iv.mpf(N) * L.iv_frac(iv, q * (1 - q))).b) * (1 + 1e-12): bad += 1
            nt += 1
ok(bad == 0, 'Lemma A, its sharper form e^-x I_0(x), and Lemma A2 (>= (12x+1)^-1/2) hold for %d exact max atoms' % nt)
ok(all(sqrt(x) * ((1 + 1 / (4 * x)) / sqrt(2 * pi * x) + exp(-x) / 2) <= 0.69 for x in [1 + i / 50 for i in range(50000)]),
   'Lemma A: for x >= 1 the bound is <= 0.69/sqrt(x) (grid 1..1001; analytic: both terms decrease, value 0.683 at x=1)')
okb = True
for M in range(1, 400):
    k0 = M // 2
    for t in range(0, k0 + 1):
        lb = sqrt(2 / (pi * (M + 2))) * exp(-t * (t + 1) / (k0 + 1 - t))
        for j in (k0 - t, k0 + t):
            if 0 <= j <= M and comb(M, j) / 2 ** M < lb * (1 - 1e-12): okb = False
ok(okb, 'Lemma B.1: b_M(k0 +- t) >= sqrt(2/(pi(M+2))) exp(-t(t+1)/(k0+1-t)), M < 400')
okb = all(comb(n, r) / 2 ** n <= sqrt(2 / (pi * (n + 0.5))) * exp(-((r - n / 2) ** 2 - 0.25) / (n / 2 + abs(r - n / 2))) * (1 + 1e-12)
          for n in range(1, 400) for r in range(n + 1))
ok(okb, 'Lemma B.2: b_n(r) <= sqrt(2/(pi(n+1/2))) exp(-(x^2-1/4)/(n/2+|x|)), n < 400')
okM = okV = okS = True; nlev = 0
for d in range(4, 41):
    n = d - 1; c = n / 2
    for (t, h) in L.pair_types(d):
        M = L.cells(t, d, h); eta = h / n if t == 'par' else (h + 0.5) / n
        etaF = Fr(h, n) if t == 'par' else Fr(2 * h + 1, 2 * n)
        for a in range(d):
            row = M[a]; tot = sum(row)
            mean = Fr(sum((b - a) * row[b] for b in range(d)), tot)
            if mean != -2 * Fr(2 * a - n, 2) * etaF: okM = False
            var = Fr(sum((b - a) ** 2 * row[b] for b in range(d)), tot) - mean ** 2
            if var > Fr(n, 2) + Fr(1, 4): okV = False
            mn = min(eta, 1 - eta)
            if mn > 0 and abs(a - c) >= sqrt(n / 2 + 0.25) / mn:
                best = max(row[b] for b in range(d) if abs(b - c) < abs(a - c))
                nlev += 1
                if best < 3 * 2 ** d * comb(n, a) / 2 ** n / (8 * abs(a - c)) * (1 - 1e-12): okS = False
ok(okM, 'E[B - A | A = a] = -2 x eta_tau (eta = h/n for par, (h+1/2)/n for crs), all types, 4 <= d <= 40 (exact)')
ok(okV, 'Var(B - A | A = a) <= n/2 + 1/4, all types and levels, 4 <= d <= 40 (exact)')
ok(okS, 'Lemma S (star lemma) at every admissible level (%d levels, 4 <= d <= 40)' % nlev)
def simple_forests(t, d, h):
    """the path / zigzag / matching forests of Sections 2.2 and 4.1: list of (level pair, closed-form cell lower
    bound, index j of the varying binomial, its row length M)"""
    n = d - 1; out = {}
    if t == 'par':
        g = n - h
        if h == n:
            out['mat'] = [((p, n - p), 4 * comb(n, p), p, n) for p in range(0, (n + 1) // 2) if 2 * p != n]
            return out
        p0 = (h - 1) // 2
        out['path'] = [((p0 + r, h - p0 + r), 4 * comb(h, p0) * comb(g, r), r, g) for r in range(g + 1)]
        if g >= 1:
            r0 = (g - 1) // 2
            out['zig'] = [((p + rr, h - p + rr), 4 * comb(h, p) * comb(g, rr), p, h) for rr in (r0, r0 + 1)
                          for p in range(0, (h + 1) // 2) if 2 * p != h]
    else:
        g = n - 1 - h
        p0 = h // 2
        out['path'] = [((p0 + r, 1 + h - p0 + r), 2 * comb(h, p0) * comb(g, r), r, g) for r in range(g + 1)]
        r0 = g // 2
        out['zig'] = [((et + p + r0, e_ + h - p + r0), 2 * comb(h, p) * comb(g, r0), p, h) for (e_, et) in ((0, 0), (1, 1))
                      for p in range(0, h + 1) if 2 * p < h]
    return out
okF = okC = True; ncell = 0
for d in range(4, 41):
    n = d - 1
    for (t, h) in L.pair_types(d):
        M = L.cells(t, d, h)
        for fam, E in simple_forests(t, d, h).items():
            keys = [(min(a, b), max(a, b)) for (a, b), _, _, _ in E]
            if len(set(keys)) != len(keys) or not L.is_forest(keys, d): okF = False
            for (a, b), lb, _, _ in E:
                ncell += 1
                if M[a][b] + M[b][a] < lb: okC = False
ok(okF, 'path, zigzag and matching edge sets are forests (distinct level pairs, acyclic), every type, 4 <= d <= 40')
ok(okC, 'their cells satisfy N_ab >= the closed-form lower bounds (4C(h,p)C(g,r) par; 2C(h,p)C(g,r) crs), %d cells' % ncell)
okK = True
for d in range(4, 41):
    n = d - 1; kap = 0.25 * sqrt(2 / (pi * (n + 1)))
    for (t, h) in L.pair_types(d):
        if t == 'par' and h == n: continue
        g = (n - h) if t == 'par' else (n - 1 - h)
        fam = 'path' if g >= h else 'zig'
        for (a, b), lb, j, Mi in simple_forests(t, d, h)[fam]:
            if lb < 2 ** d * kap * comb(Mi, j) / 2 ** Mi * (1 - 1e-12): okK = False
ok(okK, 'Section 2.2: every path/zigzag cell is >= 2^d kappa_n b_M(j), kappa_n = sqrt(2/(pi(n+1)))/4, 4 <= d <= 40')
okP = True
for d in range(4, 31):
    n = d - 1
    for (t, h) in L.pair_types(d):
        Ne = L.Nmat(L.cells(t, d, h))
        for c in (n / 2, (n - 1) / 2, (n + 1) / 2):
            F = L.potential_forest(Ne, d, c)
            okP &= L.is_forest(F, d) and len(set(F)) == len(F)
ok(okP, 'potential forests (each level joins its best cell towards the centre) are forests, d < 31')

def det_fr(M):
    M = [[Fr(x) for x in row] for row in M]; nn = len(M); dt = Fr(1)
    for i in range(nn):
        piv = next((r for r in range(i, nn) if M[r][i] != 0), None)
        if piv is None: return Fr(0)
        if piv != i: M[i], M[piv] = M[piv], M[i]; dt = -dt
        dt *= M[i][i]
        for r in range(i + 1, nn):
            fct = M[r][i] / M[i][i]
            if fct: M[r] = [x - fct * y for x, y in zip(M[r], M[i])]
    return dt
okK = True; ncase = 0
for d in range(3, 12):
    n = d - 1
    for h in range(1, n + 1):
        Ne = L.Nmat(L.cells('par', d, h))
        Lp = [[0] * d for _ in range(d)]
        for (a, b), N in Ne.items():
            Lp[a][b] -= N; Lp[b][a] -= N; Lp[a][a] += N; Lp[b][b] += N
        kappa = det_fr([row[1:] for row in Lp[1:]])
        pred = 2 ** n
        for k in range(1, n + 1):
            pred *= 2 * sum(comb(h, i) * comb(n - h, k - i) for i in range(1, h + 1, 2) if 0 <= k - i <= n - h)
        okK &= (kappa == pred)
        ncase += 1
ok(okK, 'Proposition K: weighted spanning-tree count of the par(h) level graph = 2^n prod_k D_k(h) (0 iff h even), %d cases, 3 <= d <= 11' % ncase)

# ======================================================================================= S4
section('S4 Lemma FH (Fourier-Hoelder multi-forest lemma): exact check on random small level graphs')
def ydist(n1, n2, q):
    p1 = [comb(n1, k) * q ** k * (1 - q) ** (n1 - k) for k in range(n1 + 1)]
    p2 = [comb(n2, k) * q ** k * (1 - q) ** (n2 - k) for k in range(n2 + 1)]
    out = {}
    for i, a in enumerate(p1):
        for j, b in enumerate(p2):
            out[i - j] = out.get(i - j, 0) + a * b
    return out
def exact_P0(Vn, E, q):
    dists = [list(ydist(n1, n2, q).items()) for (a, b, n1, n2) in E]
    tot = 0
    for combo in itertools.product(*dists):
        div = [0] * Vn; p = 1
        for (a, b, _, _), (y, py) in zip(E, combo):
            div[a] += y; div[b] -= y; p *= py
        if not any(div): tot += p
    return tot
def Jq(s, q):
    if q == Fr(1, 2): return exp(L.logJ_float(s))
    v = float(q * (1 - q)); K = 400
    import math
    return sum((1 - 4 * v * math.sin((-pi + (i + 0.5) * 2 * pi / K) / 2) ** 2) ** (s / 2) for i in range(K)) / K
bad = 0; tests = 0
for trial in range(120 if not QUICK else 40):
    Vn = rng.randint(3, 4)
    E = []
    for (a, b) in [(a, b) for a in range(Vn) for b in range(a + 1, Vn)]:
        if rng.random() < 0.85:
            n1, n2 = rng.randint(0, 3), rng.randint(0, 3)
            if n1 + n2 > 0: E.append((a, b, n1, n2))
    E = E[:6]
    q = rng.choice([Fr(1, 2), Fr(1, 3), Fr(1, 5)])
    P0 = float(exact_P0(Vn, E, q))
    edges = [(a, b) for (a, b, _, _) in E]; Nd = {(a, b): n1 + n2 for (a, b, n1, n2) in E}
    subsets = [S for r in range(1, len(edges) + 1) for S in itertools.combinations(edges, r) if L.is_forest(list(S), Vn)]
    for F1 in subsets:
        for F2 in subsets:
            if set(F1) & set(F2): continue
            for w in (0.5, 0.3):
                bnd = 1.0
                for e in F1: bnd *= Jq(Nd[e] / w, q) ** w
                for e in F2: bnd *= Jq(Nd[e] / (1 - w), q) ** (1 - w)
                tests += 1
                if P0 > bnd * (1 + 1e-9): bad += 1
ok(bad == 0, 'Lemma FH: P(divergence = 0) <= prod_j prod_{F_j} J_q(N/w_j)^(w_j) in all %d (graph, forest pair, weight) cases' % tests)

# ======================================================================================= S5
section('S5 Theorem B: entropy bound values and the constant c*')
def g_ent(mu):
    from math import log2
    return 0.0 if mu <= 0 else (mu + 1) * log2(mu + 1) - mu * log2(mu)
def entropy_budget(d, m):
    s = 0.0; Lm = log(m) - (d - 1) * log(2)
    for r in range(d - 1):
        lm = Lm + lgamma(d) - lgamma(r + 1) - lgamma(d - r)
        if lm < -60: continue
        s += g_ent(exp(lm))
    return s
def entropy_mmin(d):
    from math import log2
    need = log2(d) + (d - 1); lo, hi = 1, 2
    while entropy_budget(d, hi) < need: hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if entropy_budget(d, mid) >= need: hi = mid
        else: lo = mid + 1
    return lo
ENT = {6: 4, 7: 5, 8: 5, 10: 6, 12: 8, 16: 11, 20: 15, 32: 28, 64: 81, 128: 275, 256: 1159, 1024: 50116}
ok(all(entropy_mmin(d) == m for d, m in ENT.items()), 'L5 values (THM-4525) reproduced: ' + ', '.join('%d:%d' % kv for kv in ENT.items()))
cstar = (3 * log(2) / (2 * sqrt(2))) ** (2 / 3); Cstar = (3 * sqrt(2) * log(2)) ** (2 / 3)
ok(abs(cstar - 0.814581) < 1e-5 and abs(Cstar - 2.052616) < 1e-5 and abs(Cstar / cstar - 4 ** (2 / 3)) < 1e-12,
   'c* = (3 ln2/(2 sqrt2))^(2/3) = %.6f, C* = (3 sqrt2 ln2)^(2/3) = %.6f, C*/c* = 4^(2/3)' % (cstar, Cstar))
ok(abs((sqrt(2) / 3) * Cstar ** 1.5 - 2 * log(2)) < 1e-12, '(sqrt2/3) C*^(3/2) = 2 ln 2 (the bulk balance in Theorem A)')
vals = []
for d in (64, 256, 1024, 4096):
    m = entropy_mmin(d); n = d - 1
    vals.append((log(m) - 0.5 * log(pi * (n + 0.5) / 2)) / n ** (1 / 3))
ok(all(vals[i] < vals[i + 1] < cstar for i in range(len(vals) - 1)),
   'illustration: L5 minimum m gives (ln m - (1/2)ln(pi(n+1/2)/2))/n^(1/3) = ' + ', '.join('%.3f' % v for v in vals) + ' increasing towards c* (d = 64..4096)')

# ======================================================================================= S6
section('S6 Theorem D (closed-form estimate from d = 17)')
ell = lambda m: sum((2 * k - m - 1) * log(k) for k in range(1, m + 1))   # = ln prod_r C(m,r)
ok(all(abs(ell(m) - sum(log(comb(m, r)) for r in range(m + 1))) < 1e-6 * max(1, ell(m)) for m in range(1, 60)),
   'identity ln prod_r C(m,r) = sum_k (2k - m - 1) ln k (m < 60)')
ok(all(ell(m) >= (m * m - 1) / 2 - (m + 1) / 2 * log(m) - 1e-9 for m in range(1, 2001)),
   'Lemma PI: ln prod_r C(m,r) >= (m^2-1)/2 - ((m+1)/2) ln m for 1 <= m <= 2000')
f = L.thmD_funcs(iv)
# endpoint structure (the proof uses concavity; here: direct check that the minimum over each region is at its ends)
okE = True
for n in range(16, 401):
    hc = (n - 1) // 2
    P = [float(f['P_lo'](n, h).a) for h in range(0, hc + 1)]
    Z = [float(f['Zt_lo'](n, h).a) for h in range(hc + 1, n)]
    okE &= min(P) >= min(P[0], P[-1]) - 1e-9 and min(Z) >= min(Z[0], Z[-1]) - 1e-9
    okE &= float(f['P_tan0'](n).b) <= P[0] + 1e-9 and float(f['Z_tanfar'](n).b) <= Z[-1] + 1e-9
ok(okE, 'path region: min of P_- at an end; zigzag region: min of Z~_- at an end; tangent values below the true end values (16 <= n <= 400)')
for d in (17, 18, 19):
    B, t1, t2 = L.thmD_bound(iv, d)
    ok(B < 1, 'd = %d: d^2 2^(2d-3) e^-Phi + d 2^(d-2) e^-Psi <= %.5g (bulk %.4g, par(n) %.3g) < 1 (interval)' % (d, float(B), float(t1), float(t2)))
B16, _, _ = L.thmD_bound(iv, 16)
ok(B16 > 1, 'd = 16: the closed-form bound is %.3f > 1 (so certificates are used up to d = 16)' % float(B16))
okC = okT = True; minm = 1e9
for n in range(16, 3001):
    hc = (n - 1) // 2
    exact = [f['P_lo'](n, hc), f['Zt_lo'](n, hc + 1), f['P_tan0'](n), f['Z_tanfar'](n), f['Psi'](n)]
    E = L.thmD_crude(n)
    okC &= all(e <= float(x.a) + 1e-9 for e, x in zip(E, exact))
    if n >= 19:
        T = 2 * n * log(2) + 2 * log(n + 1); Tp = (n - 1) * log(2) + log(n + 1) + log(4)
        m_ = min(min(E[:4]) - T, E[4] - Tp); minm = min(minm, m_)
        okT &= m_ > 0
ok(okC, 'crude closed forms E1..E5 are below the endpoint values for 16 <= n <= 3000')
ok(okT, 'E1..E4 >= 2n ln2 + 2 ln(n+1) and E5 >= (n-1) ln2 + ln(n+1) + ln4 for 19 <= n <= 3000 (min margin %.2f; analytic for all n >= 19)' % minm)
def dE(n, i, hstep=1e-4):
    return (L.thmD_crude(n + hstep)[i] - L.thmD_crude(n - hstep)[i]) / (2 * hstep)
ok(all(dE(n, i) > 2 * log(2) + 2 / (n + 1) for i in range(4) for n in (19, 25, 50, 100, 1000)) and
   all(dE(n, 4) > log(2) + 1 / (n + 1) for n in (19, 25, 50, 100, 1000)),
   "derivatives: E_i' > T' for real n >= 19 (sampled; the quadratic terms dominate)")

# ======================================================================================= S7
section('S7 certificates: explicit edge-multiset resolving sets of Q_6 .. Q_16 (exact verification)')
CERT = {
    6: [v for v in range(64) if (0x02283022a042a00a >> v) & 1],   # Allikvere Table 1 (= THM-4525 orbit)
    7: [5, 57, 54, 104, 35, 109, 115, 49, 6, 39, 55, 102, 15, 32, 85, 97, 21, 113, 41],          # previous lane (THM-4525)
    8: [227, 92, 87, 110, 36, 32, 218, 61, 41, 64, 31, 104, 81, 251, 197, 101, 97, 112, 59, 67,
        231, 242, 80, 212, 47, 52],
    9: [504, 67, 283, 159, 247, 106, 267, 339, 3, 411, 363, 127, 373, 321, 244, 483, 410, 208, 160, 70,
        193, 34, 481, 480, 385, 100, 303, 24, 484, 380, 507, 226, 64, 117, 447, 505, 375, 131],
    10: [923, 564, 857, 521, 984, 14, 342, 541, 736, 364, 117, 642, 194, 73, 357, 207, 753, 683, 527, 338,
         578, 901, 695, 365, 977, 469, 327, 864, 929, 804, 876, 593, 70, 148, 273, 470, 882, 515, 893, 74,
         321, 102, 773, 550, 68, 461, 589, 869],
    11: [373, 1353, 1043, 1642, 1460, 1894, 1859, 388, 642, 1890, 1509, 1251, 750, 677, 1433, 818, 180, 1591, 1892, 170,
         618, 65, 1151, 522, 110, 1884, 808, 871, 18, 629, 32, 421, 320, 1790, 2016, 609, 186, 78, 1392, 371,
         1614, 892, 950, 1048, 120, 565, 1595, 787, 1345, 1303, 794, 470, 1775, 318, 646, 1374, 494, 593, 1264, 1835,
         44, 769, 270, 550, 85],
    12: [1587, 3991, 3997, 1428, 1698, 873, 631, 3526, 3853, 798, 2068, 1901, 2598, 1766, 832, 421, 2726, 3622, 202, 335,
         2600, 388, 552, 311, 262, 3077, 1509, 3330, 2343, 164, 1091, 424, 173, 3631, 110, 599, 3286, 306, 1138, 101,
         2642, 1006, 64, 481, 83, 2586, 2322, 989, 2366, 2177, 3710, 314, 1071, 275, 2049, 3273, 891, 1439, 1014, 3918,
         2022, 796, 305, 781, 3211, 26, 2071, 3238, 3263, 2926, 709, 2627, 512, 763, 1119, 3057],
    13: [60, 186, 223, 531, 554, 713, 769, 774, 795, 880, 988, 1189, 1264, 1392, 1420, 1500, 1535, 1563, 1630, 1715, 1864, 1884, 1903, 2115, 2128, 2157, 2390, 2496, 2682, 2697, 2700, 2710, 2736, 2775, 2882, 2955, 3000, 3022, 3053, 3150, 3191, 3318, 3503, 3781, 3877, 3882, 4121, 4135, 4187, 4189, 4220, 4237, 4279, 4465, 4553, 4592, 4770, 4778, 4848, 4851, 5038, 5076, 5299, 5335, 5336, 5468, 5469, 5512, 5533, 5588, 5762, 5851, 5914, 5938, 6106, 6136, 6149, 6233, 6260, 6319, 6332, 6338, 6515, 6551, 6756, 6846, 6881, 6890, 7078, 7103, 7187, 7294, 7407, 7408, 7510, 7518, 7635, 7666, 7779, 7861, 7918, 7919, 7940, 8123, 8149],
    14: [97, 110, 226, 242, 502, 538, 551, 599, 690, 1366, 1701, 1715, 2020, 2134, 2139, 2174, 2344, 2464, 2583, 2694, 2843, 2961, 2994, 3003, 3046, 3047, 3358, 3391, 3422, 3553, 3667, 3759, 3870, 4083, 4097, 4190, 4849, 4892, 4991, 5160, 5291, 5590, 5676, 5982, 6053, 6066, 6386, 6688, 6713, 6786, 6934, 7080, 7137, 7214, 7294, 7411, 7544, 7920, 7954, 8246, 8292, 8296, 8411, 8468, 8561, 8899, 9041, 9069, 9144, 9573, 9890, 10252, 10319, 10406, 10438, 10525, 10562, 10716, 10946, 11053, 11059, 11075, 11179, 11336, 11370, 11410, 11882, 12014, 12051, 12159, 12204, 12225, 12326, 12348, 12592, 12676, 12714, 12724, 13056, 13066, 13098, 13166, 13370, 13508, 13600, 13646, 13981, 14115, 14214, 14391, 14603, 14617, 14738, 14768, 14819, 15055, 15300, 15371, 15416, 15802, 15866, 15870, 15903, 16108, 16275],
    15: [26, 62, 307, 707, 974, 1030, 1069, 1241, 1250, 1348, 1363, 1422, 1431, 1470, 1498, 1528, 1543, 2147, 2223, 2546, 2780, 2829, 3294, 3423, 3662, 3710, 3916, 3917, 3944, 4329, 4559, 4675, 4933, 4936, 5562, 5664, 5955, 5963, 6803, 7322, 7670, 7794, 8015, 8048, 8209, 8594, 8634, 8665, 8670, 8988, 9115, 9292, 9408, 9826, 10236, 10270, 10458, 10485, 10745, 10898, 11512, 11568, 11888, 12327, 12382, 13433, 13748, 13820, 14203, 14376, 14621, 16156, 16420, 16481, 16508, 16729, 17159, 17272, 17433, 17740, 18032, 18168, 18359, 18450, 18633, 18742, 18970, 19319, 19324, 19561, 19716, 19802, 20183, 20342, 20564, 20651, 21024, 21030, 21319, 21841, 21975, 22294, 23260, 23379, 23836, 24548, 25290, 25295, 25407, 25615, 25859, 26033, 26084, 26108, 26131, 26392, 26508, 26810, 26828, 26831, 27835, 27869, 28187, 28250, 28252, 28905, 29049, 29559, 29750, 30181, 30692, 30955, 31867, 32543, 32766],
    16: [104, 107, 774, 1028, 1236, 1658, 2741, 2894, 3028, 4078, 4301, 4688, 4828, 5067, 5212, 5764, 5809, 6056, 6097, 6230, 6550, 6745, 6872, 8170, 8709, 8814, 9416, 9483, 9629, 10061, 10258, 11204, 13264, 13364, 13829, 13905, 13974, 15110, 15595, 15804, 15850, 16059, 16214, 16379, 17579, 17584, 17652, 18018, 19029, 20597, 20707, 20924, 21321, 22017, 22722, 23883, 24151, 24989, 25108, 25256, 25489, 26100, 26162, 26426, 26579, 27214, 27530, 28193, 28509, 29525, 30162, 30184, 30649, 30882, 31185, 31568, 31669, 32143, 33060, 33446, 33495, 34258, 34806, 34827, 35079, 35220, 35834, 35909, 36266, 36560, 36693, 36972, 37005, 37550, 37648, 37778, 38662, 39241, 39766, 39793, 40424, 40788, 41225, 41257, 41467, 41557, 41811, 41968, 42154, 42180, 42695, 43064, 43272, 43949, 44902, 45562, 45864, 46267, 46520, 46584, 46728, 46824, 47816, 47862, 48014, 48073, 48082, 48222, 48398, 48746, 48831, 49007, 49178, 49353, 51084, 52218, 52607, 52735, 53636, 54204, 54493, 54624, 54640, 55141, 55236, 55625, 56543, 56910, 57350, 57917, 58051, 59051, 59219, 59794, 60187, 60196, 60215, 60476, 61035, 61620, 62112, 62351, 62576, 63138, 63151, 63551, 63558, 63750, 64324, 65426, 65527],
}
SIZES = {6: 15, 7: 19, 8: 26, 9: 38, 10: 48, 11: 65, 12: 76, 13: 105, 14: 125, 15: 135, 16: 171}
for dd, S in sorted(CERT.items()):
    good = len(set(S)) == len(S) == SIZES[dd] and L.is_resolving_exact(dd, S)
    ok(good, 'Q_%d: explicit set of size %d is edge-multiset resolving  =>  edim_m(Q_%d) <= %d' % (dd, len(S), dd, len(S)))
    if dd >= 13:
        print('  certificate Q_%d (vertices as integers, bit i = coordinate i): %s' % (dd, ' '.join(map(str, sorted(S)))))
# negative controls: the verifier rejects non-resolving sets
ok(not L.is_resolving_exact(8, list(range(10))) and not L.is_resolving_exact(6, [0]) and not L.is_resolving_exact(6, list(range(64))),
   'verifier negative controls: {0..9} in Q_8, {0} and V(Q_6) are rejected')

# ======================================================================================= S8
section('S8 certified union bounds at q = 1/2 with Lemma FH (interval arithmetic)')
U10, per10 = L.certified_U_half(10, iv, kmax=3)
U11, per11 = L.certified_U_half(11, iv, kmax=3)
ok(U10 < 1, 'd = 10: U_10 <= %.4f < 1 (paper single-forest bound: 1.307 > 1); probabilistic method works from d = 10' % float(U10))
ok(U11 < 0.05, 'd = 11: U_11 <= %.4f (paper: 0.1562 numerically, certified 0.2548)' % float(U11))
U9, _ = L.certified_U_half(9, iv, kmax=3)
ok(U9 > 1, 'd = 9: the same bound gives %.3f > 1' % float(U9))

def U_exact_half(d, forest):
    """exact rational union bound at q = 1/2 with exact beta and a given forest rule"""
    tot = Fr(0)
    n = d - 1
    for (t, h) in L.pair_types(d):
        Ne = L.Nmat(L.cells(t, d, h))
        if forest == 'kruskal':
            F = L.kruskal(Ne, d, lambda e: Ne[e])
            P = Fr(1)
            for e in F: P *= L.beta_fr(Ne[e])
        else:   # potential forest, best of the centres n/2, (n-1)/2, (n+1)/2; par(n): its matching
            if t == 'par' and h == n:
                F = [(p, n - p) for p in range(0, (n + 1) // 2) if 2 * p != n]
                P = Fr(1)
                for e in F: P *= L.beta_fr(Ne[e])
            else:
                P = None
                for c in (n / 2, (n - 1) / 2, (n + 1) / 2):
                    F = L.potential_forest(Ne, d, c)
                    Pc = Fr(1)
                    for e in F: Pc *= L.beta_fr(Ne[e])
                    if P is None or Pc < P: P = Pc
        tot += L.type_count(t, d, h) * P
    return tot
uk = {d: U_exact_half(d, 'kruskal') for d in (10, 11, 12)}
ok(abs(float(uk[10]) - 1.30745332931) < 1e-9 and abs(float(uk[11]) - 0.156176899005) < 1e-11 and abs(float(uk[12]) - 0.0107480207982) < 1e-12,
   "paper's U_d with optimal (Kruskal) forests reproduced exactly: U_10 = %.10f, U_11 = %.12f, U_12 = %.13f" % (float(uk[10]), float(uk[11]), float(uk[12])))
up11 = U_exact_half(11, 'potential')
ok(abs(float(up11) - 0.1564) < 5e-5, 'potential forests (each level joins its best cell towards the centre) give U_11 = %.6f (exact rational)' % float(up11))
def V_closed(d):
    """interval union bound at q = 1/2 with the closed-form forests of Section 4.1 and Wallis' bound"""
    tot = iv.mpf(0)
    for (t, h) in L.pair_types(d):
        fams = simple_forests(t, d, h)
        best = None
        for fam, E in fams.items():
            if fam == 'path' and t == 'par' and h == d - 1: continue
            lp = iv.mpf(0)
            for _, lb, _, _ in E:
                lp += iv.log(2 / (iv.pi * lb)) / 2
            if best is None or lp.b < best.b: best = lp
        tot += L.type_count(t, d, h) * iv.exp(best)
    return tot
v14, v15 = V_closed(14), V_closed(15)
ok(2 < v14.a and v14.b < 3 and 0.4 < v15.a and v15.b < 0.6,
   'closed-form forests (Section 4.1, Wallis): V_14 in [%.3f, %.3f] > 1, V_15 in [%.4f, %.4f] < 1 (interval)' % (float(v14.a), float(v14.b), float(v15.a), float(v15.b)))

# ======================================================================================= S9
section('S9 certified sparse union bounds (Lemma FH, k <= 2, intervals)  =>  edim_m(Q_d) <= M_d')
SPARSE = {   # d: (q, claimed M)   -- q found by a double-precision search, then everything certified in intervals
    11: (Fr(685785, 4194304), 361), 12: (Fr(352219, 4194304), 370), 13: (Fr(876587, 16777216), 457),
    14: (Fr(954651, 33554432), 490), 15: (Fr(284779, 16777216), 588), 16: (Fr(615181, 67108864), 636),
    17: (Fr(187913, 33554432), 772), 18: (Fr(203073, 67108864), 828), 19: (Fr(484837, 268435456), 988),
    20: (Fr(259935, 268435456), 1058), 21: (Fr(76759, 134217728), 1255), 22: (Fr(657615, 2147483648), 1335),
    23: (Fr(772853, 4294967296), 1569), 24: (Fr(6483, 67108864), 1668), 25: (Fr(121781, 2147483648), 1955),
    26: (Fr(513257, 17179869184), 2057), 27: (Fr(587537, 34359738368), 2374), 28: (Fr(621019, 68719476736), 2498),
    29: (Fr(716179, 137438953472), 2889), 30: (Fr(187495, 68719476736), 3016), 31: (Fr(868761, 549755813888), 3470),
    32: (Fr(902789, 1099511627776), 3613), 36: (Fr(630525, 8796093022208), 5077),
    40: (Fr(873257, 140737488355328), 6932), 48: (Fr(764737, 18014398509481984), 12165),
    56: (Fr(320895, 1152921504606846976), 20209), 64: (Fr(1012015, 590295810358705651712), 31808),
}
for dd, (q, Mc) in sorted(SPARSE.items()):
    if QUICK and dd > 14: continue
    M2, Uu, tt = L.certified_size_v2(dd, q, iv)
    ok(M2 is not None and M2 <= Mc, 'd = %d, q = %s: U <= %.4f, P(|S| > %d) <= %.4f, sum < 1  =>  edim_m(Q_%d) <= %d  (ln M/d^(1/3) = %.3f)'
       % (dd, q, Uu, Mc, tt, dd, Mc, log(Mc) / dd ** (1 / 3)))

# ======================================================================================= S10
section('S10 Monte Carlo first moments, density 1/2 (EMPIRICAL, fixed seed)')
def mc(d, T, seed):
    rs = np.random.default_rng(seed)
    N = 1 << d
    U = np.array([u for u in range(N) for i in range(d) if not (u >> i) & 1], dtype=np.int64)
    I = np.array([i for u in range(N) for i in range(d) if not (u >> i) & 1], dtype=np.int64)
    pc = L.popcount_table(d).astype(np.int64)
    W = np.arange(N, dtype=np.int64)
    D = pc[(U[:, None] ^ W[None, :]) & ~(1 << I[:, None])]
    tot = 0; res = 0
    from collections import Counter
    for t in range(T):
        S = rs.random(N) < 0.5
        H = np.stack([((D == r) & S[None, :]).sum(axis=1) for r in range(d)], axis=1)
        keys = H @ (np.int64(N + 1) ** np.arange(d, dtype=np.int64))
        c = Counter(keys.tolist()); pairs = sum(v * (v - 1) // 2 for v in c.values())
        tot += pairs; res += (pairs == 0)
    return tot / T, res / T
m8, p8 = mc(8, 200 if not QUICK else 50, 8)
m9, p9 = mc(9, 150 if not QUICK else 40, 9)
ok(m8 > 1 and m9 < m8, 'EMPIRICAL: mean number of colliding pairs %.2f (d=8), %.2f (d=9); P(resolving) %.2f, %.2f' % (m8, m9, p8, p9))
def rand_k_resolving(d, k, T, seed):
    rs = np.random.default_rng(seed)
    return sum(L.is_resolving_exact(d, rs.choice(1 << d, k, replace=False)) for _ in range(T)) / T
if not QUICK:
    r130, r200 = rand_k_resolving(11, 130, 30, 11), rand_k_resolving(11, 200, 30, 12)
    ok(r130 < r200, 'EMPIRICAL: random k-subsets of Q_11 are resolving with frequency %.2f (k = 130) and %.2f (k = 200), 30 samples each' % (r130, r200))

# ======================================================================================= S11
section('S11 Q_7: exhaustive search with a complete (non-canonical) normal form (C program)')
import subprocess, tempfile, shutil
BUILD = tempfile.mkdtemp(prefix='procgen_edim2_')
exe = os.path.join(BUILD, 'q7search')
subprocess.check_call(['cc', '-O3', '-o', exe, os.path.join(HERE, 'procgen_edim2_20261001_q7search.c')])
def q7run(D, k, a, maxprint=0, tmax=0):
    out = subprocess.run([exe, str(D), str(k), str(a), str(maxprint), str(tmax)], capture_output=True, text=True, check=True).stdout
    summ = [l for l in out.splitlines() if l.startswith('SUMMARY')][0]
    f = dict(x.split('=') for x in summ.split()[1:])
    sets = [[int(x) for x in l.split('set:')[1].split()] for l in out.splitlines() if l.startswith('RESOLVING')]
    return int(f['Asets']), int(f['leaves']), int(f['found']), sets
def cols_of(A, D):
    return [sum((v >> i) & 1 for v in A) for i in range(D - 1)]
def py_count(D, k, a, tmax=0):
    H = 1 << (D - 1); b = k - a; dl = a - b
    beta = lambda vs: np.array([sum(1 if not (v >> i) & 1 else -1 for v in vs) for i in range(D - 1)])
    Bs = list(itertools.combinations(range(H), b))
    BB = np.array([beta(B) for B in Bs]) if b > 0 else np.zeros((1, D - 1), dtype=int)
    nA = leaves = 0
    for rest in itertools.combinations(range(1, H), a - 1):
        A = (0,) + rest
        col = cols_of(A, D)
        if any(col[i] < col[i + 1] for i in range(D - 2)): continue
        if tmax and any(sorted(cols_of([x ^ t for x in A], D), reverse=True) > col for t in A): continue
        nA += 1
        leaves += int(np.all(np.abs(BB + beta(A)[None, :]) <= dl, axis=1).sum())
    return nA, leaves
good = True
for (D, k, tm) in ((6, 6, 0), (6, 7, 0), (7, 5, 0), (6, 6, 1), (6, 7, 1), (7, 5, 1)):
    for a in range((k + 1) // 2, k + 1):
        cA, cl, cf, _ = q7run(D, k, a, 0, tm)
        good &= (cA, cl) == py_count(D, k, a, tm)
ok(good, 'q7search: A-set and leaf counts equal an independent Python enumeration, (D,k) = (6,6), (6,7), (7,5), all a, both normal forms')
def normal_form(D, S, tmax):
    beta = [sum(1 if not (s >> i) & 1 else -1 for s in S) for i in range(D)]
    j = max(range(D), key=lambda i: abs(beta[i]))
    def sw(x):
        bj, bl = (x >> j) & 1, (x >> (D - 1)) & 1
        x &= ~((1 << j) | (1 << (D - 1)))
        return x | (bj << (D - 1)) | (bl << j)
    T = [sw(s) for s in S]
    if sum(1 if not (s >> (D - 1)) & 1 else -1 for s in T) < 0: T = [s ^ (1 << (D - 1)) for s in T]
    H = 1 << (D - 1)
    A = [s for s in T if s < H]; B = [s - H for s in T if s >= H]
    if tmax:   # translate by an element of A maximising the sorted column-sum vector
        t = max(A, key=lambda t: sorted(cols_of([x ^ t for x in A], D), reverse=True))
    else:
        t = A[0]
    A = [x ^ t for x in A]; B = [x ^ t for x in B]
    col = cols_of(A, D)
    perm = sorted(range(D - 1), key=lambda i: -col[i])     # new coordinate p <- old coordinate perm[p]
    mp = lambda x: sum(((x >> perm[p]) & 1) << p for p in range(D - 1))
    return sorted(mp(x) for x in A), sorted(mp(x) for x in B)
nf_ok = True
for _ in range(300):
    D = rng.choice([6, 7]); k = rng.randint(1, 15); S = rng.sample(range(1 << D), k)
    for tm in (0, 1):
        A, B = normal_form(D, S, tm)
        a, b = len(A), len(B)
        betaS = [sum(1 if not (x >> i) & 1 else -1 for x in A) + sum(1 if not (x >> i) & 1 else -1 for x in B) for i in range(D - 1)]
        col = cols_of(A, D)
        nf_ok &= (a >= b and 0 in A and all(col[i] >= col[i + 1] for i in range(D - 2)) and all(abs(x) <= a - b for x in betaS))
        if tm:
            nf_ok &= all(sorted(cols_of([x ^ t for x in A], D), reverse=True) <= col for t in A)
ok(nf_ok, 'normal forms: 300 random sets (D = 6, 7) land in the searched domains (0 in A, sorted columns, |beta_i| <= a - b; translate-maximal)')
import procgen_edim_20261001_lib as PREV          # previous lane's library: Burnside counts of subset orbits
burn6 = PREV.burnside_subset_orbits(6)
cnt_ok = all(q7run(7, 1, a, 0, 2)[0] == burn6[a] for a in range(1, 9 if not QUICK else 8))
ok(cnt_ok, 'canonical mode (tmax = 2): number of accepted A equals the Burnside number of Aut(Q_6)-orbits of a-subsets of Q_6, a = 1..%d (%s)'
   % (8 if not QUICK else 7, ', '.join(str(burn6[a]) for a in range(1, 9 if not QUICK else 8))))
canon_exe = os.path.join(BUILD, 'canon6')
subprocess.check_call(['cc', '-O3', '-o', canon_exe, os.path.join(HERE, 'procgen_edim_20261001_canon6.c')])
THREE = {'000000181e0b126c', '000000182d11c165', '00000180026dcd03'}
for tm in ((2, 1) if QUICK else (2, 1, 0)):
    _, _, f14, _ = q7run(6, 15, 14, 0, tm); _, _, f15, _ = q7run(6, 15, 15, 0, tm)
    _, l13, f13, sets13 = q7run(6, 15, 13, 100000, tm)
    outc = subprocess.run([canon_exe], input='\n'.join('%016x' % sum(1 << v for v in S) for S in sets13) + '\n',
                          capture_output=True, text=True, check=True).stdout
    orbs = set(l.split()[0] for l in outc.splitlines())
    ok(f14 == 0 and f15 == 0 and f13 == len(sets13) > 0 and all(L.is_resolving_exact(6, S) for S in sets13) and orbs == THREE,
       'control Q_6, k = 15 (%s normal form): none with imbalance 13, 15; %d resolving leaves with imbalance 11 = exactly the 3 THM-4525 orbits'
       % ({0: 'plain', 1: 'translate-maximal', 2: 'canonical'}[tm], f13))
    if tm == 2:
        ok(f13 == 3, 'control Q_6, canonical mode: exactly one resolving leaf per orbit (3)')
tot = 0
for k in range(1, 9):
    fk = 0
    for a in range((k + 1) // 2, k + 1):
        _, lv, f, _ = q7run(7, k, a)
        fk += f; tot += lv
    ok(fk == 0, 'Q_7, k = %d: no edge-multiset resolving set (plain normal form, all a)' % k)
for (k, tm) in (((9, 2),) if QUICK else Q7K):
    fk = 0; lk = 0
    t0 = time.time()
    for a in range((k + 1) // 2, k + 1):
        _, lv, f, _ = q7run(7, k, a, 0, tm)
        fk += f; lk += lv
    ok(fk == 0, 'Q_7, k = %d: no edge-multiset resolving set (%s normal form, %d leaves, %.0f s)'
       % (k, {1: 'translate-maximal', 2: 'canonical'}[tm], lk, time.time() - t0))
shutil.rmtree(BUILD, ignore_errors=True)

print('\nTotal checks: %d   (elapsed %.0f s)' % (NCHECK, time.time() - T0))
print('ALL CHECKS PASSED')
