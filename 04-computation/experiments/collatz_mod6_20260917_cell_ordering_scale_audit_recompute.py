#!/usr/bin/env python3
"""
collatz_mod6_20260917_cell_ordering_scale_audit_recompute.py

Independent recomputation (adversarial audit, lens = recompute) of the numbers in
04-computation/experiments/collatz_mod6_20260917_cell_ordering_scale.py.  Nothing
is imported from the explorer's script.  Different algorithms where possible:

  * Omega via a smallest-prime-factor sieve + chunked repeated division
    (explorer: prime-power strided increments).
  * Selberg-Delange coefficients via explicit power-series arithmetic with the
    Abramowitz-Stegun 6.1.34 coefficients of 1/Gamma(z) and prime-zeta values
    from the prime list with tail correction (explorer: mpmath taylor/rgamma).
  * Mertens identity checked numerically: sum_p [log(1-1/p)+1/p] over the
    sieved primes vs M - gamma.
  * Ranking-change census via sign patterns of all pairwise differences in
    chunks; stability / crossover / mirror statistics via global cumsums.
  * Tournament score sequences by my own enumeration.

Universe: centers 6k, 1 <= k <= 10^7, no filter.  Rows = Omega(6k-1), columns =
Omega(6k+1), tail >= 4 collapsed.  Explicit raise on every load-bearing check.
"""
import math
import sys
import time
from fractions import Fraction

import numpy as np

T0 = time.time()
KMAX = 10**7
NMAX = 6 * KMAX + 1
DECADES = [10**3, 10**4, 10**5, 10**6, 10**7]
FIT_SCALES = [10**3, 2 * 10**3, 5 * 10**3, 10**4, 2 * 10**4, 5 * 10**4, 10**5,
              2 * 10**5, 5 * 10**5, 10**6, 2 * 10**6, 5 * 10**6, 10**7]
LAB = {0: "P", 1: "S", 2: "C", 3: "4+"}


def req(c, m):
    if not c:
        raise RuntimeError("AUDIT CHECK FAILED: " + m)


def t():
    return "[t=%.1fs]" % (time.time() - T0)


print("=== 0. smallest-prime-factor sieve on [0, %d] and chunked Omega" % NMAX)
spf = np.zeros(NMAX + 1, dtype=np.int32)
r = math.isqrt(NMAX)
for p in range(2, r + 1):
    if spf[p] == 0:
        seg = spf[p * p::p]
        seg[seg == 0] = p
idx = np.flatnonzero(spf == 0)
spf[idx] = idx            # primes (and 0,1 map to themselves; handled below)
spf[0] = 1
spf[1] = 1
om = np.zeros(NMAX + 1, dtype=np.uint8)
CH = 5 * 10**6
for s in range(0, NMAX + 1, CH):
    e = min(s + CH, NMAX + 1)
    cur = np.arange(s, e, dtype=np.int64)
    cur[cur < 2] = 1
    cnt = np.zeros(e - s, dtype=np.uint8)
    while True:
        alive = np.flatnonzero(cur > 1)
        if len(alive) == 0:
            break
        cnt[alive] += 1
        cur[alive] //= spf[cur[alive]]
    om[s:e] = cnt
del spf
primes = np.flatnonzero(om == 1)
primes = primes[np.array([True] * len(primes))]
# om==1 <=> prime; check a few
req(om[1] == 0 and om[2] == 1 and om[4] == 2 and om[12] == 3 and om[2**25] == 25 and om[3**15] == 15, "edge values")
req(len(primes) == 3562115, "pi(6e7+1) = 3562115 (explorer)")
req(int((primes <= 10**6).sum()) == 78498 and int((primes <= 10**7).sum()) == 664579, "pi(10^6), pi(10^7)")


def omega_trial(n):
    c = 0
    d = 2
    while d * d <= n:
        while n % d == 0:
            n //= d
            c += 1
        d += 1
    return c + (1 if n > 1 else 0)


rng = np.random.default_rng(1)
for n in rng.integers(2, NMAX + 1, size=2000).tolist() + list(range(2, 5001)):
    req(int(om[n]) == omega_trial(n), "trial division mismatch at %d" % n)
print(t(), "sieve OK; 2000 random + 2..5000 trial-division audit PASS")

L = np.minimum(om[5:NMAX + 1:6], 4).astype(np.int64) - 1     # Omega(6k-1) capped, 0..3
R = np.minimum(om[7:NMAX + 1:6], 4).astype(np.int64) - 1     # Omega(6k+1)
req(len(L) == KMAX and len(R) == KMAX, "lengths")
req(int(om[5:NMAX + 1:6].min()) >= 1 and int(om[7:NMAX + 1:6].min()) >= 1, "Omega>=1")
cellid = (L * 4 + R).astype(np.int8)


def mat(K):
    return np.bincount(cellid[:K], minlength=16).reshape(4, 4)


EXPLORER = {
    10**3: [[142, 182, 60, 13], [167, 192, 71, 7], [65, 71, 7, 2], [10, 9, 2, 0]],
    10**4: [[810, 1338, 687, 206], [1323, 1942, 915, 231], [685, 935, 367, 62], [196, 233, 61, 9]],
    10**5: [[5330, 10091, 6537, 2614], [10050, 17794, 10719, 3851], [6561, 10751, 5600, 1664], [2583, 3844, 1666, 345]],
    10**6: [[37915, 78689, 59706, 30192], [78277, 157420, 112416, 52182], [59992, 112305, 72363, 28755], [30161, 52125, 28736, 8766]],
    10**7: [[280557, 632905, 539140, 328636], [632766, 1380659, 1118019, 632028], [539061, 1118744, 837851, 416682], [328491, 631688, 416902, 165871]],
}
print("=== 1. 4x4 matrices vs explorer's frozen values")
for K in DECADES:
    M = mat(K)
    ok = M.tolist() == EXPLORER[K]
    print("  K=%-9d match=%s" % (K, ok))
    req(ok, "matrix at K=%d differs from explorer" % K)
M7 = mat(KMAX)
print("  rows(10^7) =", M7.sum(axis=1).tolist(), " cols =", M7.sum(axis=0).tolist())
req(M7.sum(axis=1).tolist() == [1781238, 3763472, 2912338, 1542952], "row marginals 10^7")
req(M7.sum(axis=0).tolist() == [1780875, 3763996, 2911912, 1543217], "col marginals 10^7")
M272 = mat(272)
print("  K=272: N33=%d N22=%d ; K=270: N33=%d N22=%d" % (M272[2][2], M272[1][1], mat(270)[2][2], mat(270)[1][1]))
print("  e^{e^2}/6 = %.2f  (note says 'K ~ 272')" % (math.exp(math.exp(2)) / 6))

# ---------------------------------------------------------------- 2. SD coefficients
print("=== 2. Selberg-Delange coefficients by explicit power series (independent of mpmath.taylor)")
# 1/Gamma(z) = sum c_k z^k, Abramowitz-Stegun 6.1.34
CG = [0.0, 1.0, 0.5772156649015329, -0.6558780715202538, -0.0420026350340952,
      0.1665386113822915, -0.0421977345555443, -0.0096219715278770]
EG = 0.5772156649015329
M_MERTENS = 0.2614972128476428   # Meissel-Mertens constant (reference value)
pf = primes.astype(np.float64)
# Mertens identity numeric check: sum_p [log(1-1/p)+1/p] (+ tail) vs M - gamma
s1 = float(np.sum(np.log1p(-1.0 / pf) + 1.0 / pf))
# tail: sum_{p>N} [log(1-1/p)+1/p] ~ -sum_{p>N} 1/(2p^2) ~ -1/(2 N log N)
tail1 = -1.0 / (2 * NMAX * math.log(NMAX))
print("  sum_{p<=6e7}[log(1-1/p)+1/p] = %.12f, +tail = %.12f, M - gamma = %.12f, diff = %.2e"
      % (s1, s1 + tail1, M_MERTENS - EG, s1 + tail1 - (M_MERTENS - EG)))
req(abs(s1 + tail1 - (M_MERTENS - EG)) < 1e-8, "Mertens identity numeric (f_2 = M route)")
# prime zeta P(m), m>=2, with tail ~ 1/((m-1) N^{m-1} log N)
NT = 8
P = [0.0, 0.0]
for m in range(2, NT + 1):
    val = float(np.sum(pf ** (-m)))
    val += 1.0 / ((m - 1) * NMAX ** (m - 1) * math.log(NMAX))
    P.append(val)
# log H(z) = h_1 z + sum_{m>=2} P(m)/m z^m  with h_1 = M - gamma (use the numerically verified value)
h = [0.0, M_MERTENS - EG] + [P[m] / m for m in range(2, NT + 1)]


def series_exp(a, n):
    """exp of power series a (a[0]=0), n terms."""
    b = [0.0] * n
    b[0] = 1.0
    for k in range(1, n):
        b[k] = sum(j * a[j] * b[k - j] for j in range(1, k + 1)) / k
    return b


def series_mul(a, b, n):
    return [sum(a[i] * b[k - i] for i in range(k + 1)) for k in range(n)]


NS = 7
Hs = series_exp(h[:NS], NS)
f = series_mul(Hs, CG[:NS], NS)
print("  f_m (mine):    ", ["%.9f" % v for v in f])
FEXP = [0.0, 1.0, 0.261497212847643, -0.562152927240038, 0.305977961616984, 0.0262972632346772, -0.0644501098005552]
for m in range(NS):
    print("   f_%d mine=%.10f explorer=%.10f diff=%.1e" % (m, f[m], FEXP[m], f[m] - FEXP[m]))
req(all(abs(f[m] - FEXP[m]) < 2e-7 for m in range(NS)), "f_m agree with explorer to 2e-7")
req(abs(f[2] - M_MERTENS) < 1e-12, "f_2 = M to 1e-12 (algebraically forced in this construction)")
g = series_mul(f, [1.0, -5.0 / 6.0, 1.0 / 6.0, 0, 0, 0, 0], NS)
GEXP = [0.0, 1.0, -0.571836120485691, -0.613400604613073, 0.818021603124956, -0.322376525986149, -0.0353681688932889]
for m in range(NS):
    print("   g_%d mine=%.10f explorer=%.10f diff=%.1e" % (m, g[m], GEXP[m], g[m] - GEXP[m]))
req(all(abs(g[m] - GEXP[m]) < 2e-7 for m in range(NS)), "g_m agree")
req(abs(g[2] - (M_MERTENS - 5 / 6)) < 1e-12, "g_2 = M - 5/6")
print("  M - 5/6 = %.9f ; g_3 = f_3 - 5/6 f_2 + 1/6 f_1 = %.9f" % (M_MERTENS - 5 / 6, FEXP[3] - 5 / 6 * FEXP[2] + FEXP[1] / 6))


def sdpoly(coef, k, Lv):
    return sum(coef[m] * Lv ** (k - m) / math.factorial(k - m) for m in range(0, k + 1))


def landau(k, x):
    Lv = math.log(math.log(x))
    return x / math.log(x) * Lv ** (k - 1) / math.factorial(k - 1)


print("=== 3. Global pi_k ratios (actual/Landau, actual/SD-k-term)")
for x in [10**6, 10**7, 6 * 10**7]:
    bc = np.bincount(om[1:x + 1], minlength=8)
    row = []
    for k in [1, 2, 3, 4]:
        a = int(bc[k])
        row.append("k=%d act=%d La=%.4f SD=%.4f" % (k, a, a / landau(k, x), a / (x / math.log(x) * sdpoly(f, k, math.log(math.log(x))))))
    print("  x=%d  " % x + " | ".join(row))
req(int(np.bincount(om[1:10**6 + 1], minlength=8)[2]) == 210035, "pi_2(10^6)=210035")

print("=== 4. Class marginals N_i. vs Landau and SD-F6 (i=1,2,3)")
for K in DECADES:
    M = mat(K)
    x = 6 * K
    Lv = math.log(math.log(x))
    out = []
    for i in [1, 2, 3]:
        la = 3 * K / math.log(x) * Lv ** (i - 1) / math.factorial(i - 1)
        sd = 3 * K / math.log(x) * sdpoly(g, i, Lv)
        out.append("i=%d N=%d La=%.4f SD=%.4f" % (i, M.sum(axis=1)[i - 1], M.sum(axis=1)[i - 1] / la, M.sum(axis=1)[i - 1] / sd))
    print("  K=%-8d " % K + " | ".join(out))

# ---------------------------------------------------------------- 5. rankings and Landau merged order
print("=== 5. Rankings (nine-cell, six merged) and the Landau six-class order")


def rank_str(items):
    items = sorted(items, key=lambda t: -t[1])
    s = items[0][0]
    for p, q in zip(items, items[1:]):
        s += (" = " if p[1] == q[1] else " > ") + q[0]
    return s


def six(M):
    return [("PP", int(M[0][0])), ("PS+SP", int(M[0][1] + M[1][0])), ("SS", int(M[1][1])),
            ("PC+CP", int(M[0][2] + M[2][0])), ("SC+CS", int(M[1][2] + M[2][1])), ("CC", int(M[2][2]))]


def nine(M):
    return [(LAB[i] + LAB[j], int(M[i][j])) for i in range(3) for j in range(3)]


def merged_w(w):
    return [("PP", w[1] * w[1]), ("PS+SP", 2 * w[1] * w[2]), ("SS", w[2] * w[2]),
            ("PC+CP", 2 * w[1] * w[3]), ("SC+CS", 2 * w[2] * w[3]), ("CC", w[3] * w[3])]


for K in DECADES:
    M = mat(K)
    Lv = math.log(math.log(6 * K))
    wl = {i: Lv ** (i - 1) / math.factorial(i - 1) for i in [1, 2, 3]}
    ws = {i: sdpoly(g, i, Lv) for i in [1, 2, 3]}
    nine_s = rank_str(nine(M))
    six_s = rank_str(six(M))
    la9 = rank_str([(LAB[i] + LAB[j], wl[i + 1] * wl[j + 1]) for i in range(3) for j in range(3)])
    la6 = rank_str(merged_w(wl))
    sd6 = rank_str(merged_w(ws))
    row = M.sum(axis=1)
    col = M.sum(axis=0)
    mg6 = rank_str([("PP", row[0] * col[0]), ("PS+SP", row[0] * col[1] + row[1] * col[0]), ("SS", row[1] * col[1]),
                    ("PC+CP", row[0] * col[2] + row[2] * col[0]), ("SC+CS", row[1] * col[2] + row[2] * col[1]),
                    ("CC", row[2] * col[2])])
    pos_cc_9 = [n for n, _ in sorted(nine(M), key=lambda t: -t[1])].index("CC") + 1
    pos_cc_6 = [n for n, _ in sorted(six(M), key=lambda t: -t[1])].index("CC") + 1
    pos_ss_6 = [n for n, _ in sorted(six(M), key=lambda t: -t[1])].index("SS") + 1
    print("  K=%d L=%.3f" % (K, Lv))
    print("    nine actual : %s   (CC rank %d of 9)" % (nine_s, pos_cc_9))
    print("    six  actual : %s   (CC rank %d of 6, SS rank %d of 6)" % (six_s, pos_cc_6, pos_ss_6))
    print("    Landau nine : %s" % la9)
    print("    Landau six  : %s   %s" % (la6, "MATCH" if la6 == six_s else "differs"))
    print("    SD-F6 six   : %s   %s" % (sd6, "MATCH" if sd6 == six_s else "differs"))
    print("    marg six    : %s   %s" % (mg6, "MATCH" if mg6 == six_s else "differs"))
print("  Landau merged: SC+CS (L^3) > CC (L^4/4) iff L < 4, i.e. CC is FIRST in the six-class Landau order only for L > 4"
      " (x = e^{e^4} ~ 10^%.1f); for 2<L<4 Landau puts SC+CS first, CC second, and SS = PC+CP tied." % (math.exp(4) / math.log(10)))
for K in [100, 270, 500] + DECADES:
    M = mat(K)
    print("  2x2 K=%-8d %s" % (K, rank_str([("PP", int(M[0][0])), ("PS", int(M[0][1])), ("SP", int(M[1][0])), ("SS", int(M[1][1]))])))

# ---------------------------------------------------------------- 6. running statistics via global cumsum
print("=== 6. Running statistics over all K <= 10^7 (global cumsums)")
onehot_cum = {}
for c in range(16):
    onehot_cum[c] = np.cumsum(cellid == c, dtype=np.int32)
req(all(int(onehot_cum[c][-1]) == int(M7.reshape(16)[c]) for c in range(16)), "cumsum totals")
Ks = np.arange(1, KMAX + 1)
cPP, cPS, cSP, cSS = onehot_cum[0], onehot_cum[1], onehot_cum[4], onehot_cum[5]
cPC, cCP, cSC, cCS, cCC = onehot_cum[2], onehot_cum[8], onehot_cum[6], onehot_cum[9], onehot_cum[10]
v1 = ~(cSS > np.maximum(cPS, cSP))
v2 = ~(np.minimum(cPS, cSP) > cPP)
print("  #K with not(SS>max(PS,SP)) = %d, last = %d, all K<=last? %s" % (v1.sum(), Ks[v1][-1], bool(v1[:Ks[v1][-1]].all())))
print("  #K with not(min(PS,SP)>PP) = %d, last = %d ; the K<=160 where it HOLDS: %s"
      % (v2.sum(), Ks[v2][-1], Ks[:160][~v2[:160]].tolist()))
req(int(v1.sum()) == 761 and int(Ks[v1][-1]) == 761, "SS>max(PS,SP) violations 761/761")
req(int(v2.sum()) == 151 and int(Ks[v2][-1]) == 160, "min(PS,SP)>PP violations 151/160")
for name, (a, b) in {"PS-SP": (cPS, cSP), "PC-CP": (cPC, cCP), "SC-CS": (cSC, cCS)}.items():
    D = a.astype(np.int64) - b
    nz = D[D != 0]
    flips = int(np.sum(np.sign(nz[1:]) != np.sign(nz[:-1])))
    zeros = int((D == 0).sum())
    am = int(np.argmax(np.abs(D)))
    # returns to zero = number of maximal runs of D==0 (entered from nonzero)
    z = (D == 0)
    runs = int(np.sum(z[1:] & ~z[:-1])) + int(z[0])
    print("  %s final=%+d zeros(#K with D=0)=%d zero-runs(returns to 0)=%d flips=%d max|D|=%d at K=%d sqrt(tot)=%.0f"
          % (name, D[-1], zeros, runs, flips, abs(D[am]), Ks[am], math.sqrt(a[-1] + b[-1])))
D = cPS.astype(np.int64) - cSP
req(int(D[-1]) == 139 and int((D == 0).sum()) == 14947 and int(np.max(np.abs(D))) == 650 and int(Ks[np.argmax(np.abs(D))]) == 5512997, "PS-SP stats")

# crossovers 'for good': last K with A<=B
def last_le(A, B):
    le = A <= B
    return (int(Ks[le][-1]) if le.any() else 0, int(le.sum()), bool((~le).any()))


cross = {("SS", "PS"): (cSS, cPS), ("SS", "SP"): (cSS, cSP), ("SS", "PS+SP"): (cSS, cPS + cSP),
         ("SC+CS", "SS"): (cSC + cCS, cSS), ("SC+CS", "PS+SP"): (cSC + cCS, cPS + cSP),
         ("CC", "PS"): (cCC, cPS), ("CC", "SP"): (cCC, cSP), ("CC", "PS+SP"): (cCC, cPS + cSP),
         ("CC", "PC+CP"): (cCC, cPC + cCP), ("PC+CP", "PP"): (cPC + cCP, cPP), ("CC", "PP"): (cCC, cPP),
         ("CC", "SS"): (cCC, cSS), ("SC", "SS"): (cSC, cSS), ("CS", "SS"): (cCS, cSS)}
EXP_CROSS = {("SS", "PS"): 761, ("SS", "SP"): 672, ("SS", "PS+SP"): 920807, ("SC+CS", "SS"): 16815,
             ("SC+CS", "PS+SP"): 66959, ("CC", "PS"): 1591506, ("CC", "SP"): 1559904, ("PC+CP", "PP"): 1445, ("CC", "PP"): 83683}
for key, (A, B) in cross.items():
    lk, n, ever = last_le(A, B)
    print("  %-6s vs %-6s: last K with A<=B = %d (#=%d); A>B ever: %s" % (key[0], key[1], lk, n, ever))
    if key in EXP_CROSS:
        req(lk == EXP_CROSS[key], "crossover %s" % (key,))
    else:
        req(not ever, "%s never exceeds %s" % key)

# ranking-change census via sign patterns (chunked)
def census(groups):
    G = [sum(onehot_cum[c].astype(np.int64) for c in gl) for gl in groups]
    ng = len(G)
    pairs = [(a, b) for a in range(ng) for b in range(a + 1, ng)]
    S = np.stack([np.sign(G[a] - G[b]).astype(np.int8) for (a, b) in pairs], axis=1)
    ch = np.any(S[1:] != S[:-1], axis=1)
    idx = np.flatnonzero(ch) + 2      # K index (1-based) of the changed state
    return len(idx), (int(idx[-1]) if len(idx) else 0), idx[:6].tolist()


SYS = {"nine": [[0], [1], [2], [4], [5], [6], [8], [9], [10]],
       "six": [[0], [1, 4], [5], [2, 8], [6, 9], [10]],
       "2x2": [[0], [1], [4], [5]],
       "2x2m": [[0], [1, 4], [5]]}
EXP_CENSUS = {"nine": (8918, 8714402), "six": (259, 920808), "2x2": (3885, 8714402), "2x2m": (33, 920808)}
for name, groups in SYS.items():
    n, last, first = census(groups)
    print("  census %-5s changes=%d last=%d first=%s" % (name, n, last, first))
    req((n, last) == EXP_CENSUS[name], "census %s" % name)
del onehot_cum
print(t(), "running statistics done")

# ---------------------------------------------------------------- 7. (3,3) vs (2,2), R matrix, models
print("=== 7. N33/N22 table, independence ratios, crossover models")
ROWS = []
for K in FIT_SCALES:
    M = mat(K).astype(np.float64)
    Lv = math.log(math.log(6 * K))
    row = M.sum(axis=1)
    col = M.sum(axis=0)
    N33, N22 = M[2][2], M[1][1]
    c33 = N33 * K / (row[2] * col[2])
    c22 = N22 * K / (row[1] * col[1])
    ROWS.append((K, Lv, N33, N22, N33 / N22, row[2] / row[1], col[2] / col[1], c33 / c22))
    print("  K=%-8d L=%.4f N33=%d N22=%d ratio=%.5f r-=%.5f r+=%.5f c33/c22=%.5f" % ROWS[-1])
M = mat(KMAX).astype(np.float64)
Rm = M * KMAX / np.outer(M.sum(axis=1), M.sum(axis=0))
print("  R(10^7):")
for i in range(4):
    print("    " + " ".join("%.4f" % Rm[i][j] for j in range(4)))
C2 = 1.0
for p in primes[1:].tolist():
    if p > 10**7:
        break
    C2 *= 1 - 1 / (p - 1) ** 2
print("  C_2 (p<=10^7) = %.8f  4C_2/3 = %.5f ; R_PP(10^7) = %.4f" % (C2, 4 * C2 / 3, Rm[0][0]))
req(abs(C2 - 0.6601618158) < 1e-7, "C2")
for K in DECADES:
    Mk = mat(K).astype(np.float64)
    Rk = Mk * K / np.outer(Mk.sum(axis=1), Mk.sum(axis=0))
    print("  K=%-8d R_PP=%.4f R_SS=%.4f R_CC=%.4f R_44=%.4f" % (K, Rk[0][0], Rk[1][1], Rk[2][2], Rk[3][3]))
# shuffle with explorer's seed (same RNG stream) and with a different seed
for seed in [20260917, 7]:
    rg = np.random.default_rng(seed)
    Rs = R[rg.permutation(KMAX)]
    Msh = np.bincount((L * 4 + Rs), minlength=16).reshape(4, 4).astype(np.float64)
    Rsh = Msh * KMAX / np.outer(Msh.sum(axis=1), Msh.sum(axis=0))
    z = np.max(np.abs(Rsh - 1) * np.sqrt(Msh))
    print("  shuffle seed=%d: max|R-1|sqrt(N) = %.2f ; R_PP = %.4f" % (seed, z, Rsh[0][0]))
    req(z < 5, "shuffle control")


def bis(fn, lo, hi):
    flo = fn(lo)
    for _ in range(200):
        mid = (lo + hi) / 2
        fm = fn(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def rth(Lv):
    return sdpoly(g, 3, Lv) / sdpoly(g, 2, Lv)


def Kof(Lv):
    return math.exp(math.exp(Lv)) / 6


c7 = ROWS[-1][7]
c6 = ROWS[-4][7]
LB = bis(lambda Lv: c7 * rth(Lv) ** 2 - 1, 2, 8)
sel = [rw for rw in ROWS if rw[0] >= 10**5]
Ls = np.array([rw[1] for rw in sel])
ac, bc_ = np.polyfit(Ls, np.array([math.sqrt(rw[5] * rw[6]) for rw in sel]), 1)
LC = bis(lambda Lv: c7 * (ac * Lv + bc_) ** 2 - 1, 1, 12)
ad, bd = np.polyfit(Ls, np.array([math.sqrt(rw[4]) for rw in sel]), 1)
LD = (1 - bd) / ad
L6, L7 = ROWS[-4][1], ROWS[-1][1]
slope = (c7 - c6) / (1 / L7 - 1 / L6)
cinf = c7 + slope * (0 - 1 / L7)
LE = bis(lambda Lv: (c7 + slope * (1 / Lv - 1 / L7)) * rth(Lv) ** 2 - 1, 2, 8)
for tag, Lv in [("B", LB), ("C", LC), ("D", LD), ("E", LE)]:
    print("  Model %s: L*=%.4f K*=%.3e log10K*=%.2f" % (tag, Lv, Kof(Lv), math.log10(Kof(Lv))))
print("  fit C: sqrt(r-r+) = %.3f L %+.3f ; fit D: sqrt(N33/N22) = %.3f L %+.3f ; E c_inf = %.4f" % (ac, bc_, ad, bd, cinf))
print("  r_theory at decades:", ["%.4f" % rth(math.log(math.log(6 * K))) for K in DECADES])
req(abs(LB - 3.1589) < 2e-4 and abs(LC - 3.2197) < 2e-4 and abs(LD - 3.1867) < 2e-4 and abs(LE - 3.1268) < 2e-4, "model L*")
w = lambda i, Lv: sdpoly(g, i, Lv)
for tag, fn in [("CC>SS", lambda Lv: w(3, Lv) - w(2, Lv)), ("CC>SC+CS", lambda Lv: w(3, Lv) - 2 * w(2, Lv)), ("SS>PS+SP", lambda Lv: w(2, Lv) - 2 * w(1, Lv))]:
    Lv = bis(fn, 2, 9)
    print("  SD-F6 c=1 %-9s at L=%.4f, log10 x = %.2f" % (tag, Lv, math.exp(Lv) / math.log(10)))
# Landau (pure) thresholds for comparison
print("  pure Landau merged thresholds: CC>SC+CS at L=4 (log10 x=%.1f); CC>SS at L=2; SS>PS+SP at L=2" % (math.exp(4) / math.log(10)))
print("  Landau claim 'x~10^60 for CC>SC+CS' is the SD-F6 number; the pure-Landau value is 10^%.1f" % (math.exp(4) / math.log(10)))

# ---------------------------------------------------------------- 8. tournaments
print("=== 8. Tournament score sequences n=4 (own enumeration)")
seqs = set()
edges = [(a, b) for a in range(4) for b in range(a + 1, 4)]
for mask in range(1 << 6):
    sc = [0] * 4
    for e, (a, b) in enumerate(edges):
        sc[a if (mask >> e) & 1 else b] += 1
    seqs.add(tuple(sorted(sc)))
seqs = sorted(seqs)
print("  ", seqs)
req(seqs == [(0, 1, 2, 3), (0, 2, 2, 2), (1, 1, 1, 3), (1, 1, 2, 2)], "score sequences")
req(not any(len(set(s)) == 3 and s[1] == s[2] for s in seqs), "no max>mid=mid>min")
# Landau's theorem check: nondecreasing s with partial sums >= C(k,2) and total 6
land = sorted({s for s in __import__("itertools").combinations_with_replacement(range(4), 4)
               if sum(s) == 6 and all(sum(s[:k]) >= k * (k - 1) // 2 for k in range(1, 5))})
req(land == seqs, "Landau 1953 characterization agrees with enumeration")
print("  Landau-condition sequences:", land, "-> agree")
# weak order a > {b,c} > d: arcs of the strict comparability digraph
arcs = [(x, y) for x in "abcd" for y in "abcd" if {"a": 3, "b": 2, "c": 2, "d": 1}[x] > {"a": 3, "b": 2, "c": 2, "d": 1}[y]]
print("  strict arcs of weak order a>{b,c}>d:", len(arcs))
req(len(arcs) == 5, "5 arcs")

# ---------------------------------------------------------------- 9. 3x3 block
print("=== 9. 3x3 block fractions")
for K in DECADES:
    M = mat(K)
    blk = M[:3, :3]
    tot = int(blk.sum())
    fr = {"PP": Fraction(int(blk[0][0]), tot), "PS+SP": Fraction(int(blk[0][1] + blk[1][0]), tot), "SS": Fraction(int(blk[1][1]), tot),
          "PC+CP": Fraction(int(blk[0][2] + blk[2][0]), tot), "SC+CS": Fraction(int(blk[1][2] + blk[2][1]), tot), "CC": Fraction(int(blk[2][2]), tot)}
    req(sum(fr.values()) == 1, "fractions sum")
    print("  K=%-8d block=%d (%.4f)  PP=%.5f PS=%.5f PC=%.5f SP=%.5f SS=%.5f SC=%.5f CP=%.5f CS=%.5f CC=%.5f"
          % ((K, tot, tot / K) + tuple(blk[i][j] / tot for i in range(3) for j in range(3))))
req(mat(10**7)[:3, :3].sum() == 7079702 and mat(10**6)[:3, :3].sum() == 769083, "block totals")

print(t(), "ALL AUDIT CHECKS PASSED")
