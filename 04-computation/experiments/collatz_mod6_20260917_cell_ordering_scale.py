#!/usr/bin/env python3
"""
collatz_mod6_20260917_cell_ordering_scale.py -- sandwich cell ordering across scales.

Lane: cell_ordering_scale (session collatz-mod6-20260917).

Object.  Centers W = 6k (k = 1..K), endpoints 6k-1 and 6k+1.  Omega counts prime
factors with multiplicity.  Cell (i,j) = #{k <= K : Omega(6k-1) = i, Omega(6k+1) = j}
with the tail i or j >= 4 collapsed.  P = Omega 1 (prime), S = Omega 2 (semiprime),
C = Omega 3 (3-almost-prime).  Rows describe W-1, columns W+1 (same convention as
the inherited SW3 table).

Inheritance (read, not re-derived; cited by path in the note):
  05-knowledge/results/arithmetic_braids_20260917_divisors.md  (SW1 CRT identity,
     SW2 k -> -k mirror and class-parity law, SW3 K = 10^6 matrix, Chen audit)
  04-computation/experiments/collatz_mod6_20260917_sandwich_bias.py  (sibling lane:
     owns the antisymmetric part N_ij - N_ji; here only sign-change counts)

What this script does (sections printed):
  0  Omega sieve on [1, 6*10^7+1] (uint8, ~60 MB) + independent trial-division
     and sympy.factorint audits (positive controls); K = 10^6 matrix must equal
     the inherited SW3 table (raise otherwise).
  1  Landau / Selberg-Delange marginals: pi_k(x) actual vs Landau leading term vs
     the k-term Selberg-Delange polynomial (coefficients computed with mpmath from
     the prime zeta function; f_2 must equal the Mertens constant).  Predicted
     asymptotic order of the nine cells under the product heuristic.
  2  FINITE-EXACT 4x4 matrices at K = 10^3..10^7, nine-cell rankings, ranking
     change census (four ranking systems), 2x2 stability (last violating K).
  3  (3,3) vs (2,2): ratio across scales, independence-ratio matrix
     N_ij K / (N_i. N_.j), shuffle hostile control, crossover extrapolation under
     four models with the spread reported as the uncertainty.
  4  Tournament audit: all 64 labelled 4-tournaments, their score sequences
     (Landau's theorem, FINITE-EXACT), the mirror involution on cells.
  5  3x3 (P,S,C) block: conditional probabilities and symmetry ties.

All load-bearing checks use explicit `raise` (active under python -O).
Run:  python3 04-computation/experiments/collatz_mod6_20260917_cell_ordering_scale.py
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
LABEL = {1: "P", 2: "S", 3: "C", 4: "4+"}

INHERITED_K1E6 = [
    [37915, 78689, 59706, 30192],
    [78277, 157420, 112416, 52182],
    [59992, 112305, 72363, 28755],
    [30161, 52125, 28736, 8766],
]


def require(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def elapsed():
    return "[t=%.1fs]" % (time.time() - T0)


def hdr(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


# ----------------------------------------------------------------------------
# 0. Omega sieve and audits
# ----------------------------------------------------------------------------
def omega_sieve(n_max):
    is_p = np.ones(n_max + 1, dtype=bool)
    is_p[:2] = False
    r = math.isqrt(n_max)
    for p in range(2, r + 1):
        if is_p[p]:
            is_p[p * p::p] = False
    primes = np.flatnonzero(is_p).astype(np.int64)
    del is_p
    om = np.zeros(n_max + 1, dtype=np.uint8)
    small = primes[primes <= r]
    big = primes[primes > r]
    for p in small.tolist():
        q = p
        while q <= n_max:
            om[q::q] += 1
            q *= p
    # primes p > sqrt(n_max) divide n <= n_max at most once: n = m*p, m <= sqrt(n_max)
    for m in range(1, r + 1):
        cnt = int(np.searchsorted(big, n_max // m, side="right"))
        if cnt == 0:
            break
        om[m * big[:cnt]] += 1
    return om, primes


def omega_trial(n):
    c = 0
    d = 2
    while d * d <= n:
        while n % d == 0:
            n //= d
            c += 1
        d += 1
    if n > 1:
        c += 1
    return c


hdr("0. Omega sieve on [1, %d] and audits" % NMAX)
om, primes = omega_sieve(NMAX)
print("sieve done", elapsed(), "primes up to", NMAX, "=", len(primes))
require(om[1] == 0 and om[2] == 1 and om[4] == 2 and om[2**25] == 25, "edge Omega values")
require(om[6 * 10**7 + 1] == omega_trial(6 * 10**7 + 1), "top endpoint")
for n in range(2, 30001):
    if om[n] != omega_trial(n):
        raise RuntimeError("Omega mismatch at n=%d" % n)
print("trial-division audit 2..30000: PASS")
try:
    from sympy import factorint
    rng0 = np.random.default_rng(20260917)
    samp = rng0.integers(2, NMAX + 1, size=300)
    for n in samp.tolist():
        if om[n] != sum(factorint(n).values()):
            raise RuntimeError("sympy Omega mismatch at n=%d" % n)
    print("sympy.factorint audit on 300 random n <= NMAX: PASS")
except ImportError:
    print("sympy not available; skipped factorint audit (trial division audit stands)")

Lraw = om[5:NMAX + 1:6]          # Omega(6k-1), k = 1..KMAX
Rraw = om[7:NMAX + 1:6]          # Omega(6k+1), k = 1..KMAX
require(len(Lraw) == KMAX and len(Rraw) == KMAX, "endpoint array lengths")
require(int(Lraw.min()) >= 1 and int(Rraw.min()) >= 1, "Omega >= 1 on endpoints")
Lc = np.minimum(Lraw, 4).astype(np.uint8)
Rc = np.minimum(Rraw, 4).astype(np.uint8)
cell = ((Lc - 1) * 4 + (Rc - 1)).astype(np.uint8)   # 0..15, row-major (i-1, j-1)


def matrix_at(K, cells=cell):
    return np.bincount(cells[:K], minlength=16).reshape(4, 4).astype(np.int64)


M6 = matrix_at(10**6)
require(M6.tolist() == INHERITED_K1E6, "K=10^6 matrix must equal inherited SW3 table")
print("K=10^6 matrix equals inherited SW3 table: PASS (positive control)")

# independent K<=1000 matrix by trial division
ind = [[0] * 4 for _ in range(4)]
for k in range(1, 1001):
    a = min(omega_trial(6 * k - 1), 4)
    b = min(omega_trial(6 * k + 1), 4)
    ind[a - 1][b - 1] += 1
require(ind == matrix_at(1000).tolist(), "independent K<=1000 matrix")
print("independent trial-division matrix at K=1000: PASS")


def print_matrix(M, title):
    print(title)
    print("  L\\R    " + "".join("%10s" % LABEL[j] for j in range(1, 5)) + "%10s" % "row")
    for i in range(4):
        print("  %-5s  " % LABEL[i + 1] + "".join("%10d" % M[i][j] for j in range(4))
              + "%10d" % sum(M[i]))
    print("  %-5s  " % "col" + "".join("%10d" % sum(M[i][j] for i in range(4)) for j in range(4))
          + "%10d" % int(np.sum(M)))


# ----------------------------------------------------------------------------
# 1. Landau / Selberg-Delange
# ----------------------------------------------------------------------------
hdr("1. Landau / Selberg-Delange marginals and the predicted asymptotic order")
import mpmath
from mpmath import mp

mp.dps = 30
MERTENS = mp.mertens
# F(z) = H(1,z)/Gamma(z), H(1,z) = prod_p (1-1/p)^z (1-z/p)^{-1};
# log H(1,z) = sum_{m>=2} (z^m - z) P(m)/m, P = prime zeta.
PZ = [None, None] + [mp.primezeta(m) for m in range(2, 90)]


def logH(z):
    return mp.fsum((z**m - z) * PZ[m] / m for m in range(2, 90))


def F_all(z):
    return mp.exp(logH(z)) * mp.rgamma(z)


f = mp.taylor(F_all, 0, 6)     # f[0..6]
require(abs(f[0]) < mp.mpf(10) ** -20, "f_0 = 0 (1/Gamma(0)=0)")
require(abs(f[1] - 1) < mp.mpf(10) ** -20, "f_1 = 1")
require(abs(f[2] - MERTENS) < mp.mpf(10) ** -18, "f_2 must equal the Mertens constant")
print("Selberg-Delange coefficients of F(z)=H(1,z)/Gamma(z):")
for m in range(0, 7):
    print("  f_%d = %s" % (m, mp.nstr(f[m], 15)))
print("  (f_2 == Mertens constant M = %s : PASS)" % mp.nstr(MERTENS, 15))
# coprime-to-6 version: F_6(z) = F(z)(1 - z/2)(1 - z/3)
g = mp.taylor(lambda z: F_all(z) * (1 - z / 2) * (1 - z / 3), 0, 6)
print("Coefficients of F_6(z) = F(z)(1-z/2)(1-z/3) (integers coprime to 6):")
for m in range(0, 7):
    print("  g_%d = %s" % (m, mp.nstr(g[m], 15)))
require(abs(g[2] - (MERTENS - mp.mpf(5) / 6)) < mp.mpf(10) ** -18, "g_2 = M - 5/6")


def sd_poly(coef, k, Lval):
    """[z^k] exp(z L) * sum coef_m z^m  =  sum_{m<=k} coef_m L^{k-m}/(k-m)!"""
    return sum(coef[m] * mp.mpf(Lval) ** (k - m) / mp.factorial(k - m) for m in range(0, k + 1))


def landau(k, x):
    Lv = math.log(math.log(x))
    return x / math.log(x) * Lv ** (k - 1) / math.factorial(k - 1)


def sd_full(coef, k, x, factor=1.0):
    Lv = math.log(math.log(x))
    return float(factor * x / math.log(x) * sd_poly(coef, k, Lv))


print()
print("Global pi_k(x) = #{n<=x: Omega(n)=k}: actual vs Landau leading term vs")
print("k-term Selberg-Delange main term (x/log x) * sum_{m<=k} f_m L^{k-m}/(k-m)!, L=loglog x")
print("  %10s %2s %10s %12s %8s %12s %8s" % ("x", "k", "actual", "Landau", "ratio", "SD-k-term", "ratio"))
GLOBAL_PI = {}
for x in [10**6, 10**7, 6 * 10**7]:
    bc = np.bincount(om[1:x + 1], minlength=8)
    for k in [1, 2, 3, 4]:
        act = int(bc[k])
        GLOBAL_PI[(x, k)] = act
        la = landau(k, x)
        sd = sd_full(f, k, x)
        print("  %10d %2d %10d %12.0f %8.4f %12.0f %8.4f" % (x, k, act, la, act / la, sd, act / sd))
require(GLOBAL_PI[(10**6, 1)] == 78498, "pi(10^6) = 78498")
require(GLOBAL_PI[(10**7, 1)] == 664579, "pi(10^7) = 664579")
require(GLOBAL_PI[(10**6, 2)] == 210035, "pi_2(10^6) = 210035 (OEIS A066265)")
print("  pi(10^6), pi(10^7), pi_2(10^6) match known values: PASS")

print()
print("HEURISTIC product model: N_ij(K) ~ c_ij * N_i.(K) N_.j(K) / K with Landau marginals")
print("  N_i.(K) ~ (3K/log 6K) L^{i-1}/(i-1)!  =>  N_ij ~ c * (9K/log^2 6K) L^{i+j-2}/((i-1)!(j-1)!)")
print("  Landau weights w_ij = L^{i+j-2}/((i-1)!(j-1)!) for the nine cells:")
print("   (1,1):1  (1,2)=(2,1):L  (2,2):L^2  (1,3)=(3,1):L^2/2  (2,3)=(3,2):L^3/2  (3,3):L^4/4")
print("  Pairs (i,j),(j,i) tie exactly (symmetric weight).  (2,2) vs (1,3): factor 1 vs 1/2.")
print("  Threshold crossings under Landau weights: (3,3)=(2,3)=(2,2) all at L=2 (K=e^{e^2}/6~272);")
print("  (2,2)>(1,2) for L>1; (3,3)>(1,1) for L>sqrt(2); CC>PS for L^3>4 i.e. L>1.587.")
print("  Landau predicted order for L>2:  CC > CS=SC > SS > CP=PC > PS=SP > PP")
print("  Landau predicted order for 1<L<2: SS > CS=SC > CC ...  (see actual data below)")
print("  STATUS: HEURISTIC (Hardy-Littlewood-type product with a common singular series);")
print("          the marginal Landau term itself is very inaccurate at these scales (table above).")


# ----------------------------------------------------------------------------
# 2. Exact matrices, rankings, ranking-change census
# ----------------------------------------------------------------------------
hdr("2. FINITE-EXACT 4x4 matrices, nine-cell rankings, ranking-change census")
MATS = {}
for K in DECADES:
    MATS[K] = matrix_at(K)
    print_matrix(MATS[K], "K = %d  (centers 6..%d)" % (K, 6 * K))
    print()

NINE = [(i, j) for i in range(1, 4) for j in range(1, 4)]


def cname(i, j):
    return LABEL[i] + LABEL[j]


def ranking_string(items):
    """items: list of (name, count). Return 'a > b = c > d' string (desc)."""
    items = sorted(items, key=lambda t: -t[1])
    out = items[0][0]
    for p, q in zip(items, items[1:]):
        out += (" = " if p[1] == q[1] else " > ") + q[0]
    return out


print("Nine-cell (i,j<=3) ranking at each decade (rows W-1 / columns W+1):")
for K in DECADES:
    M = MATS[K]
    items = [(cname(i, j), int(M[i - 1][j - 1])) for (i, j) in NINE]
    print("  K=%-9d %s" % (K, ranking_string(items)))
print("Mirror-merged six-class ranking (XY+YX pooled):")
for K in DECADES:
    M = MATS[K]
    items = [("PP", int(M[0][0])), ("PS+SP", int(M[0][1] + M[1][0])), ("SS", int(M[1][1])),
             ("PC+CP", int(M[0][2] + M[2][0])), ("SC+CS", int(M[1][2] + M[2][1])), ("CC", int(M[2][2]))]
    print("  K=%-9d %s" % (K, ranking_string(items)))
print()
print("Class marginals N_i.(K) = #{k<=K: Omega(6k-1)=i} vs Landau (3K/log 6K) L^{i-1}/(i-1)! vs")
print("Selberg-Delange with F_6 coefficients g_m: (3K/log 6K) sum_{m<=i} g_m L^{i-m}/(i-m)!")
print("  %9s %2s %9s %9s %11s %8s %11s %8s" % ("K", "i", "N_i.", "N_.i", "Landau", "ratio", "SD-F6", "ratio"))
for K in DECADES:
    M = MATS[K]
    x = 6 * K
    Lv = math.log(math.log(x))
    for i in [1, 2, 3]:
        la = 3 * K / math.log(x) * Lv ** (i - 1) / math.factorial(i - 1)
        sd = float(3 * K / math.log(x) * sd_poly(g, i, Lv))
        print("  %9d %2d %9d %9d %11.0f %8.4f %11.0f %8.4f" % (K, i, int(M.sum(axis=1)[i - 1]),
              int(M.sum(axis=0)[i - 1]), la, M.sum(axis=1)[i - 1] / la, sd, M.sum(axis=1)[i - 1] / sd))
print("  (Landau is off by up to a factor ~2 for i=3 on the coprime-to-6 classes: g_2 = M - 5/6 < 0.)")
print()
print("Predicted six-class rankings from three symmetric product models vs actual (FINITE-EXACT):")
print("  [Landau]  w_i = L^{i-1}/(i-1)!;  [SD-F6]  w_i = sum g_m L^{i-m}/(i-m)!;")
print("  [marg]    w_ij = N_i. * N_.j measured (c_ij := 1).  Merged classes double the mixed weights.")
n_sd_match = 0
n_marg_match = 0
for K in DECADES:
    M = MATS[K]
    x = 6 * K
    Lv = math.log(math.log(x))
    actual = ranking_string([("PP", int(M[0][0])), ("PS+SP", int(M[0][1] + M[1][0])), ("SS", int(M[1][1])),
                             ("PC+CP", int(M[0][2] + M[2][0])), ("SC+CS", int(M[1][2] + M[2][1])),
                             ("CC", int(M[2][2]))])

    def merged(w):
        return [("PP", w[1] * w[1]), ("PS+SP", 2 * w[1] * w[2]), ("SS", w[2] * w[2]),
                ("PC+CP", 2 * w[1] * w[3]), ("SC+CS", 2 * w[2] * w[3]), ("CC", w[3] * w[3])]

    wl = {i: Lv ** (i - 1) / math.factorial(i - 1) for i in [1, 2, 3]}
    ws = {i: float(sd_poly(g, i, Lv)) for i in [1, 2, 3]}
    row = M.sum(axis=1)
    col = M.sum(axis=0)
    wm = [("PP", row[0] * col[0]), ("PS+SP", row[0] * col[1] + row[1] * col[0]), ("SS", row[1] * col[1]),
          ("PC+CP", row[0] * col[2] + row[2] * col[0]), ("SC+CS", row[1] * col[2] + row[2] * col[1]),
          ("CC", row[2] * col[2])]
    r_l, r_s, r_m = ranking_string(merged(wl)), ranking_string(merged(ws)), ranking_string(wm)
    n_sd_match += int(r_s == actual)
    n_marg_match += int(r_m == actual)
    print("  K=%d" % K)
    print("    actual : %s" % actual)
    print("    Landau : %s   %s" % (r_l, "MATCH" if r_l == actual else "differs"))
    print("    SD-F6  : %s   %s" % (r_s, "MATCH" if r_s == actual else "differs"))
    print("    marg   : %s   %s" % (r_m, "MATCH" if r_m == actual else "differs"))
print("  SD-F6 matches the actual six-class ranking at %d of 5 decades; measured-marginal product at %d of 5."
      % (n_sd_match, n_marg_match))
print("  STATUS: the observed ordering is explained by the marginals (HEURISTIC product model with c~1);")
print("  the naive Landau ordering is REFUTED at every decade (it puts CC first for L>2).")
print()
print("2x2 block ranking:")
for K in [100, 270, 500] + DECADES:
    M = matrix_at(K)
    items = [("PP", int(M[0][0])), ("PS", int(M[0][1])), ("SP", int(M[1][0])), ("SS", int(M[1][1]))]
    print("  K=%-9d %s" % (K, ranking_string(items)))

# Ranking systems: list of (name, groups) where a group is a list of cell ids.
cid = lambda i, j: (i - 1) * 4 + (j - 1)
SYSTEMS = {
    "nine cells": [(cname(i, j), [cid(i, j)]) for (i, j) in NINE],
    "six mirror-merged": [("PP", [cid(1, 1)]), ("PS+SP", [cid(1, 2), cid(2, 1)]), ("SS", [cid(2, 2)]),
                          ("PC+CP", [cid(1, 3), cid(3, 1)]), ("SC+CS", [cid(2, 3), cid(3, 2)]),
                          ("CC", [cid(3, 3)])],
    "2x2 block": [("PP", [cid(1, 1)]), ("PS", [cid(1, 2)]), ("SP", [cid(2, 1)]), ("SS", [cid(2, 2)])],
    "2x2 mirror-merged": [("PP", [cid(1, 1)]), ("PS+SP", [cid(1, 2), cid(2, 1)]), ("SS", [cid(2, 2)])],
}
CHUNK = 2 * 10**5
state = {}
for name, groups in SYSTEMS.items():
    ng = len(groups)
    pairs = [(a, b) for a in range(ng) for b in range(a + 1, ng)]
    state[name] = {"pairs": pairs, "prev": None, "changes": [], "n_changes": 0}
viol_ss = []      # K with not (SS > max(PS,SP))
viol_mid = []     # K with not (min(PS,SP) > PP)
n_viol_ss = 0
n_viol_mid = 0
# mirror differences: D12 = N12 - N21, D13 = N13 - N31, D23 = N23 - N32
mirror = {"PS-SP": (cid(1, 2), cid(2, 1)), "PC-CP": (cid(1, 3), cid(3, 1)), "SC-CS": (cid(2, 3), cid(3, 2))}
mstat = {m: {"zeros": 0, "maxabs": 0, "argmax": 0, "flips": 0, "prev_sign": 0} for m in mirror}
# CC vs PS+SP and CC vs each: track last K where CC <= PS (per cell) for the crossover report
PSSP = [cid(1, 2), cid(2, 1)]
PCCP = [cid(1, 3), cid(3, 1)]
SCCS = [cid(2, 3), cid(3, 2)]
cross_track = {("SS", "PS"): (cid(2, 2), cid(1, 2)), ("SS", "SP"): (cid(2, 2), cid(2, 1)),
               ("SS", "PS+SP"): (cid(2, 2), PSSP), ("SC+CS", "SS"): (SCCS, cid(2, 2)),
               ("SC+CS", "PS+SP"): (SCCS, PSSP), ("CC", "PS"): (cid(3, 3), cid(1, 2)),
               ("CC", "SP"): (cid(3, 3), cid(2, 1)), ("CC", "PS+SP"): (cid(3, 3), PSSP),
               ("CC", "PC+CP"): (cid(3, 3), PCCP), ("PC+CP", "PP"): (PCCP, cid(1, 1)),
               ("CC", "PP"): (cid(3, 3), cid(1, 1)), ("CC", "SS"): (cid(3, 3), cid(2, 2)),
               ("SC", "SS"): (cid(2, 3), cid(2, 2)), ("CS", "SS"): (cid(3, 2), cid(2, 2))}
# per pair: (last K with A<=B, #K with A<=B, last K with A>B, #K with A>B)
cross_last = {key: [0, 0, 0, 0] for key in cross_track}
offset = np.zeros(16, dtype=np.int64)
for start in range(0, KMAX, CHUNK):
    stop = min(start + CHUNK, KMAX)
    onehot = np.zeros((stop - start, 16), dtype=np.int32)
    onehot[np.arange(stop - start), cell[start:stop]] = 1
    C = np.cumsum(onehot, axis=0) + offset          # (chunk, 16) cumulative counts at K=start+1..stop
    offset = C[-1].copy()
    Ks = np.arange(start + 1, stop + 1)
    for name, groups in SYSTEMS.items():
        G = np.stack([C[:, gl].sum(axis=1) for (_, gl) in groups], axis=1)
        st = state[name]
        S = np.stack([np.sign(G[:, a] - G[:, b]).astype(np.int8) for (a, b) in st["pairs"]], axis=1)
        prev = st["prev"]
        if prev is None:
            prevS = np.vstack([np.zeros((1, S.shape[1]), dtype=np.int8), S[:-1]])
            prevS[0] = S[0]           # K=1 is not a change
        else:
            prevS = np.vstack([prev, S[:-1]])
        ch = np.any(S != prevS, axis=1)
        idx = np.flatnonzero(ch)
        st["n_changes"] += int(len(idx))
        if len(idx):
            st["changes"].extend(Ks[idx].tolist())
        st["prev"] = S[-1:].copy()
    PP, PS, SP, SS = C[:, cid(1, 1)], C[:, cid(1, 2)], C[:, cid(2, 1)], C[:, cid(2, 2)]
    v1 = ~(SS > np.maximum(PS, SP))
    v2 = ~(np.minimum(PS, SP) > PP)
    n_viol_ss += int(v1.sum())
    n_viol_mid += int(v2.sum())
    if v1.any():
        viol_ss.append(int(Ks[np.flatnonzero(v1)[-1]]))
    if v2.any():
        viol_mid.append(int(Ks[np.flatnonzero(v2)[-1]]))
    for m, (a, b) in mirror.items():
        D = C[:, a] - C[:, b]
        ms = mstat[m]
        ms["zeros"] += int((D == 0).sum())
        am = int(np.argmax(np.abs(D)))
        if abs(int(D[am])) > ms["maxabs"]:
            ms["maxabs"] = abs(int(D[am]))
            ms["argmax"] = int(Ks[am])
        sg = np.sign(D)
        nz = sg[sg != 0]
        seq = np.concatenate([[ms["prev_sign"]], nz]) if ms["prev_sign"] != 0 else nz
        if len(seq) > 1:
            ms["flips"] += int(np.sum(seq[1:] != seq[:-1]))
        if len(nz):
            ms["prev_sign"] = int(nz[-1])
    for key, (a, b) in cross_track.items():
        A = C[:, a] if isinstance(a, int) else C[:, a].sum(axis=1)
        B = C[:, b] if isinstance(b, int) else C[:, b].sum(axis=1)
        le = A <= B
        rec = cross_last[key]
        if le.any():
            rec[0] = int(Ks[np.flatnonzero(le)[-1]])
        rec[1] += int(le.sum())
        gt = ~le
        if gt.any():
            rec[2] = int(Ks[np.flatnonzero(gt)[-1]])
        rec[3] += int(gt.sum())
require(offset.tolist() == matrix_at(KMAX).reshape(16).tolist(), "chunked cumulative totals")
print(elapsed(), "ranking-change census done")
print()
print("Ranking-change census over all K = 1..%d (a change = the weak order incl. ties differs from K-1):" % KMAX)
for name, groups in SYSTEMS.items():
    st = state[name]
    ch = st["changes"]
    print("  %-20s changes=%-7d last change at K=%d" % (name, st["n_changes"], ch[-1] if ch else 0))
    if ch:
        print("     first changes: %s" % ch[:12])
        print("     last changes : %s" % ch[-12:])
print()
print("2x2 stability (running check over all K <= %d):" % KMAX)
print("  #K violating SS > max(PS,SP): %d ; last violating K = %d" % (n_viol_ss, max(viol_ss) if viol_ss else 0))
print("  #K violating min(PS,SP) > PP: %d ; last violating K = %d" % (n_viol_mid, max(viol_mid) if viol_mid else 0))
K_STABLE = max(max(viol_ss) if viol_ss else 0, max(viol_mid) if viol_mid else 0) + 1
print("  => SS > {PS,SP} > PP holds for every K in [%d, %d]  (FINITE-EXACT)" % (K_STABLE, KMAX))
print("  (PS vs SP itself keeps flipping; see mirror statistics.)")
print()
print("Mirror differences D(K)=N_ij-N_ji (sibling lane owns the mechanism; here only census):")
for m, ms in mstat.items():
    a, b = mirror[m]
    tot = int(offset[a] + offset[b])
    print("  %-6s final D=%+d  sqrt(N_ij+N_ji)=%.0f  #K with D=0: %d  strict sign flips: %d  max|D|=%d at K=%d"
          % (m, int(offset[a] - offset[b]), math.sqrt(tot), ms["zeros"], ms["flips"], ms["maxabs"], ms["argmax"]))
print()
print("Crossovers A vs B (FINITE-EXACT; 'A>B for good' means for every K in (K_last, 10^7]):")
for key, (lastle, nle, lastgt, ngt) in cross_last.items():
    A, B = key
    if lastle < KMAX:
        print("  %-7s > %-7s for good after K=%-9d (#K with A<=B: %d)" % (A, B, lastle, nle))
    elif lastgt == 0:
        print("  %-7s <= %-6s at every K<=10^7 (A never exceeds B)" % (A, B))
    else:
        print("  %-7s <= %-6s at K=10^7; A>B held for %d values of K, last at K=%d" % (A, B, ngt, lastgt))


# ----------------------------------------------------------------------------
# 3. (3,3) vs (2,2), independence ratios, shuffle control, crossover
# ----------------------------------------------------------------------------
hdr("3. (3,3) vs (2,2) across scales; independence ratios; crossover extrapolation")
print("  %9s %8s %9s %9s %9s %9s %9s %9s %9s" % ("K", "L", "N33", "N22", "N33/N22", "r-", "r+", "c33/c22", "Landau L^2/4"))
ROWS = []
for K in FIT_SCALES:
    M = matrix_at(K)
    x = 6 * K
    Lv = math.log(math.log(x))
    row = M.sum(axis=1)
    col = M.sum(axis=0)
    N33, N22 = int(M[2][2]), int(M[1][1])
    rminus = row[2] / row[1]      # N_3./N_2.  (6k-1 side)
    rplus = col[2] / col[1]       # N_.3/N_.2  (6k+1 side)
    c33 = N33 * K / (row[2] * col[2])
    c22 = N22 * K / (row[1] * col[1])
    ROWS.append((K, Lv, N33, N22, N33 / N22, rminus, rplus, c33 / c22, c33, c22))
    print("  %9d %8.4f %9d %9d %9.5f %9.5f %9.5f %9.5f %9.4f" % (K, Lv, N33, N22, N33 / N22, rminus, rplus, c33 / c22, Lv**2 / 4))
print("  STATUS: counts FINITE-EXACT; L=loglog(6K).  Naive Landau ratio L^2/4 predicts N33>N22 already")
print("  at K~272: REFUTED (N33(272)=%d, N22(272)=%d)." % (int(matrix_at(272)[2][2]), int(matrix_at(272)[1][1])))

print()
print("Independence-ratio matrix  R_ij = N_ij * K / (N_i. * N_.j)  (1 = product of marginals):")
C2 = 1.0
for p in primes[1:].tolist():
    if p > 10**6:
        break
    C2 *= 1 - 1 / (p - 1) ** 2
print("  twin-prime constant C_2 (primes<=10^6) = %.7f ; Hardy-Littlewood prediction for R_PP -> 4C_2/3 = %.5f"
      % (C2, 4 * C2 / 3))
require(abs(C2 - 0.6601618158) < 2e-6, "C_2 value")
RATIOS = {}
for K in DECADES:
    M = MATS[K].astype(np.float64)
    row = M.sum(axis=1)
    col = M.sum(axis=0)
    Rm = M * K / np.outer(row, col)
    RATIOS[K] = Rm
    print("  K=%d" % K)
    print("    L\\R  " + "".join("%9s" % LABEL[j] for j in range(1, 5)))
    for i in range(4):
        print("    %-4s " % LABEL[i + 1] + "".join("%9.4f" % Rm[i][j] for j in range(4)))
print("  Row/column marginals at K=10^7: rows N_i. =", MATS[KMAX].sum(axis=1).tolist(),
      " cols N_.j =", MATS[KMAX].sum(axis=0).tolist())
print("  Hostile control: randomly re-pair right endpoints (seed 20260917) at K=10^7 and recompute R_ij:")
rng = np.random.default_rng(20260917)
Rsh = Rc[rng.permutation(KMAX)]
cell_sh = ((Lc - 1) * 4 + (Rsh - 1)).astype(np.uint8)
Msh = matrix_at(KMAX, cell_sh).astype(np.float64)
Rm_sh = Msh * KMAX / np.outer(Msh.sum(axis=1), Msh.sum(axis=0))
print("    L\\R  " + "".join("%9s" % LABEL[j] for j in range(1, 5)))
for i in range(4):
    print("    %-4s " % LABEL[i + 1] + "".join("%9.4f" % Rm_sh[i][j] for j in range(4)))
zmax = 0.0
for i in range(4):
    for j in range(4):
        n = Msh[i][j]
        zmax = max(zmax, abs(Rm_sh[i][j] - 1) * math.sqrt(n))
print("    max |R_ij - 1| * sqrt(N_ij) over shuffled cells = %.2f (should be O(1))" % zmax)
require(zmax < 5.0, "shuffled control must be independent within 5 sigma")
require(RATIOS[KMAX][0][0] < 0.95 and abs(RATIOS[KMAX][0][0] - 4 * C2 / 3) < 0.02,
        "true PP ratio near 4C_2/3 and clearly below 1")
print("    true data at K=10^7: R_PP = %.4f (HL prediction 4C_2/3 = %.4f); shuffled R_PP = %.4f"
      % (RATIOS[KMAX][0][0], 4 * C2 / 3, Rm_sh[0][0]))

print()
print("Crossover extrapolation for N33/N22 = 1 (all HEURISTIC; exact data end at K=10^7):")
# Model A: naive Landau: L^2/4 = 1 -> L = 2
KA = math.exp(math.exp(2.0)) / 6
print("  Model A (naive Landau L^2/4):            L*=2.000  K* = %.0f  -> REFUTED by exact data above" % KA)


def r_theory(coef, Lv):
    return float(sd_poly(coef, 3, Lv) / sd_poly(coef, 2, Lv))


def bisect(fn, lo, hi, tol=1e-9):
    flo = fn(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = fn(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def report_model(tag, Lstar):
    Kstar = math.exp(math.exp(Lstar)) / 6
    print("  %-42s L*=%.4f  K* = %.3e  (log10 K* = %.2f)" % (tag, Lstar, Kstar, math.log10(Kstar)))
    return Kstar


# measured c-ratio at the largest scales
c_last = ROWS[-1][7]
c_prev = ROWS[-4][7]   # at 10^6
# Model B: Selberg-Delange coprime-to-6 marginal ratio, c fixed at K=10^7 value
LB = bisect(lambda Lv: c_last * r_theory(g, Lv) ** 2 - 1, 2.0, 8.0)
KB = report_model("Model B (SD F_6 marginals, c=c(10^7)=%.4f)" % c_last, LB)
print("      SD ratio check: r_theory(L) vs measured r-, r+ at the decades:")
for (K, Lv, N33, N22, rat, rm, rp, cr, c33, c22) in ROWS:
    if K in DECADES:
        print("        K=%-8d L=%.3f  r_theory=%.4f  r-=%.4f  r+=%.4f  (Landau L/2=%.4f)" % (K, Lv, r_theory(g, Lv), rm, rp, Lv / 2))
# Model C: empirical linear fit of sqrt(r- r+) vs L on K >= 10^5, c fixed
sel = [r for r in ROWS if r[0] >= 10**5]
Ls = np.array([r[1] for r in sel])
ys = np.array([math.sqrt(r[5] * r[6]) for r in sel])
a_c, b_c = np.polyfit(Ls, ys, 1)
LC = bisect(lambda Lv: c_last * (a_c * Lv + b_c) ** 2 - 1, 1.0, 12.0)
KC = report_model("Model C (linear fit sqrt(r- r+)=%.3fL%+.3f, c fixed)" % (a_c, b_c), LC)
# Model D: direct linear fit of sqrt(N33/N22) vs L on K >= 10^5
yd = np.array([math.sqrt(r[4]) for r in sel])
a_d, b_d = np.polyfit(Ls, yd, 1)
LD = bisect(lambda Lv: (a_d * Lv + b_d) - 1, 1.0, 12.0)
KD = report_model("Model D (direct fit sqrt(N33/N22)=%.3fL%+.3f)" % (a_d, b_d), LD)
# Model E: Model B but with c drifting: linear extrapolation of c in 1/L, using 10^6 and 10^7
L6, L7 = ROWS[-4][1], ROWS[-1][1]
slope_c = (c_last - c_prev) / (1 / L7 - 1 / L6)


def c_of_L(Lv):
    return c_last + slope_c * (1 / Lv - 1 / L7)


LE = bisect(lambda Lv: c_of_L(Lv) * r_theory(g, Lv) ** 2 - 1, 2.0, 8.0)
KE = report_model("Model E (SD marginals, c drifting in 1/L -> c_inf=%.4f)" % c_of_L(1e9), LE)
Ks_all = [KB, KC, KD, KE]
print("  Spread of the non-refuted models: log10 K* in [%.2f, %.2f]  (x = 6K: log10 x in [%.2f, %.2f])"
      % (min(map(math.log10, Ks_all)), max(map(math.log10, Ks_all)),
         min(map(math.log10, Ks_all)) + math.log10(6), max(map(math.log10, Ks_all)) + math.log10(6)))
print("  STATUS: HEURISTIC extrapolation by 2-3 orders of magnitude beyond the data; the models span")
print("  about one order of magnitude in K*, all assume c33/c22 converges, and none is a theorem.")
print()
print("When does the SD-F6 product model (c=1) reach the Landau asymptotic order CC > SC+CS > SS > ...?")
w = lambda i, Lv: float(sd_poly(g, i, Lv))
L_cc_ss = bisect(lambda Lv: w(3, Lv) - w(2, Lv), 2.0, 9.0)
L_cc_sccs = bisect(lambda Lv: w(3, Lv) - 2 * w(2, Lv), 2.0, 9.0)
L_ss_pssp = bisect(lambda Lv: w(2, Lv) - 2 * w(1, Lv), 2.0, 9.0)
for tag, Lv in [("CC > SS  (w3 > w2)", L_cc_ss), ("CC > SC+CS (w3 > 2 w2)", L_cc_sccs), ("SS > PS+SP (w2 > 2 w1)", L_ss_pssp)]:
    print("  %-26s at L = %.4f, i.e. x = 6K = exp(exp(L)) ~ 10^%.1f" % (tag, Lv, math.exp(Lv) / math.log(10)))
print("  (actual SS > PS+SP for good after K=%d, x~10^%.2f; the c=1 model puts it at x~10^%.2f: about one decade early.)"
      % (cross_last[("SS", "PS+SP")][0], math.log10(6 * cross_last[("SS", "PS+SP")][0]), math.exp(L_ss_pssp) / math.log(10)))
print("  STATUS: HEURISTIC; the CC-on-top order of the naive heuristic is a ~10^60 phenomenon under SD-F6.")


# ----------------------------------------------------------------------------
# 4. Tournament audit
# ----------------------------------------------------------------------------
hdr("4. Tournament audit (repo guardrail)")
scores = set()
for mask in range(64):
    s = [0, 0, 0, 0]
    e = 0
    for a in range(4):
        for b in range(a + 1, 4):
            if (mask >> e) & 1:
                s[a] += 1
            else:
                s[b] += 1
            e += 1
    scores.add(tuple(sorted(s)))
scores = sorted(scores)
print("Score sequences of all 64 labelled tournaments on 4 vertices:", scores)
require(scores == [(0, 1, 2, 3), (0, 2, 2, 2), (1, 1, 1, 3), (1, 1, 2, 2)], "Landau score sequences n=4")
has_max_two_equal_min = any(len(set(s)) == 3 and s[1] == s[2] for s in scores)
require(not has_max_two_equal_min, "no 4-tournament has score profile max > two equal middles > min")
print("PROVED/FINITE-EXACT: no 4-vertex tournament has score profile (max, mid, mid, min) with")
print("  max > mid > min.  The user's 'max, min, two equal middles' is therefore not the score")
print("  profile of any tournament; it is a weak order (total preorder) with one tie, i.e. a")
print("  transitive comparability digraph on 4 vertices with 5 arcs, not 6.")
print("Intrinsic relations between cells: only the involution iota: k -> -k, which maps cell (i,j)")
print("  to (j,i) (inherited SW2).  On the nine cells iota has 3 fixed points (PP,SS,CC) and 3")
print("  2-cycles; it is symmetric, so it orients nothing.  An orientation such as (i,j)->(j,i) for")
print("  i<j is a gauge (left endpoint < right endpoint), not data.")
print("Frequency ordering = a scale-dependent total preorder on the cells; at K=%d:" % KMAX)
M = MATS[KMAX]
print("  " + ranking_string([(cname(i, j), int(M[i - 1][j - 1])) for (i, j) in NINE]))
print("VERDICT: the '4-tournament' reading is a cosmetic ranking.  It carries no content beyond the")
print("  four counts, has no intrinsic pairwise observable, and its claimed tie is not a tournament")
print("  feature but a heuristic asymptotic symmetry (mirror cells) that exact counts violate at")
print("  every finite K where N_12 != N_21.")


# ----------------------------------------------------------------------------
# 5. 3x3 block
# ----------------------------------------------------------------------------
hdr("5. 3x3 (P,S,C) block: conditional probabilities given both endpoints have Omega<=3")
for K in DECADES:
    M = MATS[K]
    blk = M[:3, :3]
    tot = int(blk.sum())
    print("K=%d  block total=%d of %d centers (%.4f)" % (K, tot, K, tot / K))
    print("    L\\R  " + "".join("%9s" % LABEL[j] for j in range(1, 4)) + "%9s" % "row")
    for i in range(3):
        print("    %-4s " % LABEL[i + 1] + "".join("%9.5f" % (blk[i][j] / tot) for j in range(3))
              + "%9.5f" % (blk[i].sum() / tot))
    print("    %-4s " % "col" + "".join("%9.5f" % (blk[:, j].sum() / tot) for j in range(3)))
    # exact counts for the six mirror-merged classes (denominator = block total)
    cnts = {"PP": int(blk[0][0]), "PS+SP": int(blk[0][1] + blk[1][0]), "SS": int(blk[1][1]),
            "PC+CP": int(blk[0][2] + blk[2][0]), "SC+CS": int(blk[1][2] + blk[2][1]), "CC": int(blk[2][2])}
    require(sum(Fraction(v, tot) for v in cnts.values()) == 1, "conditional probabilities sum to 1 exactly")
    print("    mirror-merged exact: " + ", ".join("%s=%d/%d" % (k, v, tot) for k, v in cnts.items()))
print()
print("How C changes the picture (HEURISTIC asymptotic symmetry, FINITE-EXACT order):")
print("  asymptotic ties by mirror symmetry: PS~SP, PC~CP, SC~CS (3 tied pairs); PP,SS,CC unpaired.")
print("  Landau leading order for L>2 would put CC on top; at every computed K<=10^7 SS is on top,")
print("  SC/CS second, PS/SP third, CC fourth (CC overtakes PC/CP and PP early; see crossover table).")
print("  With C the 2x2 'max, min, two middles' becomes a 9-cell preorder with three symmetric pairs")
print("  and three singletons: 6 classes, of which only the pairs are asymptotically tied.")

print()
print(elapsed(), "ALL CHECKS PASSED")
