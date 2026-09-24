#!/usr/bin/env python3
"""collatz_procgen_20260924_tree_counts.py

Lane: inverse tree mod 192 (session collatz-procgen-20260922, 2026-09-24).

Part 3: what the type automaton CAN prove, and where it stops.

  K1  Predecessor-count functions N_d(a) = |T^-d(a)| on Z/3^d (d <= 14): exact Haar means
      E[N_{d,e}] = C(d,e) 3^-e, per-type means, spread (min over units, second moment).
  K2  The Haar-averaged size series 1/(1 - 2^-s - (3/2)^s/3): poles s = 1, 2; residue 1/log(2/sqrt3)
      (the reciprocal AM-GM gap); dyadic window means M_n / 2^n -> 1/log(2/sqrt3); qx+1 comparison.
  K3  Actual tree densities of integer roots on both sheets (X = 10^7) against the averaged model.
  K4  Incongruent classes mod 3*2^k: forward-only, owner-refined (one backward E at 2 mod 3 points),
      and minimal-chain refined (all class-determined 3-adic precision), with growth exponents.
  K5  Rational periodic points: sign = sign(2^p - 3^a), integrality = divisibility by 2^p - 3^a;
      the integer ones for p <= 18 and those lying in Bad.
  K6  The 16-point trap F of E_{6 mod 8} is the negative image of the 3n-1 cycles {1,2}, {5,7,10}.
  K7  The exact renewal identity for basins and its neutral main term (AM = 1).

Every check raises on failure.  Peak memory about 400 MB; runtime about 15 s.
"""
import sys
import time
from fractions import Fraction
from math import comb, log, log2, lgamma, exp, sqrt

import numpy as np

T0 = time.time()
LOG2_3 = log2(3.0)
R_AMGM = 1.0 / log(2.0 / sqrt(3.0))


def check(cond, msg):
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)


def T(n, b=1):
    return n // 2 if n % 2 == 0 else (3 * n + b) // 2


def mem(tag):
    import resource
    print(f"[mem] {tag}: max RSS so far {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB", file=sys.stderr)


def hdr(s):
    mem("before " + s[:3])
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


# ---------------------------------------------------------------------------
hdr("K1  Predecessor-count functions on Z/3^d (the 3-adic half of the owner's type)")
# ---------------------------------------------------------------------------
print("N_d(a) = number of depth-d T-predecessors of a (words in D, E legal at a).  Legality is 3-adic only;")
print("N_d is a function of a mod 3^d:  N_d(a) = N_(d-1)(2a) + [a = 2 mod 3] N_(d-1)((2a-1)/3).")
print()
print("  d | mean over Z/3^d | (4/3)^d  | mean(a=1)  mean(a=2) | min over units (arg a mod 3^d) | max    | E[N^2]/E[N]^2 (units)")
Nprev = np.ones(1, dtype=np.int32)  # d = 0 on Z/1 (int32 throughout: counts <= 2^14, indices < 3^14)
mins = []
for d in range(1, 15):
    Md = 3 ** d
    a = np.arange(Md, dtype=np.int32)
    prevM = 3 ** (d - 1)
    Nd = Nprev[(2 * a) % prevM].copy()
    br = (a % 3) == 2
    ab = a[br]
    Nd[br] += Nprev[((2 * ab - 1) // 3) % prevM]
    tot = int(Nd.sum(dtype=np.int64))
    check(tot * 3 ** 0 == sum(comb(d, e) * 3 ** (d - e) for e in range(d + 1)), f"Haar mean identity d={d}")
    check(tot == 4 ** d, "sum over Z/3^d of N_d = 4^d")
    units = (a % 3) != 0
    Nu = Nd[units]
    m1 = Nd[a % 3 == 1].mean()
    m2 = Nd[a % 3 == 2].mean()
    imin = int(np.argmin(np.where(units, Nd, np.iinfo(np.int32).max)))
    mn = int(Nd[imin])
    mx = int(Nd.max())
    sec = float((Nu.astype(np.float64) ** 2).mean() / Nu.mean() ** 2)
    mins.append(mn)
    print(f" {d:2d} | {tot / Md:15.6f} | {(4/3)**d:8.4f} | {m1:9.4f}  {m2:9.4f} | {mn:6d} ({imin:>8d})           | {mx:6d} | {sec:.4f}")
    Nprev = Nd
check(all(mins[i] <= mins[i + 1] for i in range(len(mins) - 1)), "min over units nondecreasing")
print("  Exact: sum_a N_d(a) = 4^d over Z/3^d, so the Haar mean is (4/3)^d.  Per-type means (m0 = 1 on 3Z):")
print("  m1_d = m2_(d-1) (the D-child of a = 1 is = 2 mod 3), m2_d = m1_(d-1) + (4/3)^(d-1) (the E-child is Haar-")
print("  uniform on Z_3); ratio m2/m1 -> 4/3, the Perron vector of A_6.")
m1, m2 = Fraction(1), Fraction(1)
for dd in range(1, 15):
    m1, m2 = m2, m1 + Fraction(4, 3) ** (dd - 1)
    check((1 + m1 + m2) / 3 == Fraction(4, 3) ** dd, "per-type recursion consistent with the Haar mean")
print("  (recursion checked exactly: (m0 + m1_d + m2_d)/3 = (4/3)^d for d <= 14)")
# bivariate identity for d <= 9
for d in range(1, 10):
    Md = 3 ** d
    # enumerate words by DP on (a mod 3^d) with e-count polynomials
    cnt = np.zeros(d + 1, dtype=np.int64)
    # recursive enumeration over residues: count legal words with e E-steps, summed over a
    # use DP over depth with arrays indexed by residue mod 3^(d-t) and e
    # layer t: residues mod 3^(d) at start
    cur = {}
    # simple exact enumeration: for each residue class mod 3^d, walk all words
    tot_e = [0] * (d + 1)
    for a0 in range(Md):
        stack = [(a0, 0, 0)]
        while stack:
            x, t, e = stack.pop()
            if t == d:
                tot_e[e] += 1
                continue
            stack.append((2 * x, t + 1, e))
            if x % 3 == 2:
                stack.append(((2 * x - 1) // 3, t + 1, e + 1))
    for e in range(d + 1):
        check(tot_e[e] == comb(d, e) * 3 ** (d - e), f"E[N_(d,e)] = C(d,e) 3^-e at d={d}, e={e}")
print("  Bivariate identity sum_a N_(d,e)(a) = C(d,e) 3^(d-e), i.e. E[N_(d,e)] = C(d,e) 3^-e, checked d <= 9;")
print("  generating function E[sum_(d,e) N_(d,e) z^d u^e] = 1/(1 - z(1 + u/3)).")
growth = [mins[i + 1] / mins[i] for i in range(len(mins) - 1)]
print(f"  The minimum over 3-adic units grows by factors {', '.join(f'{g:.3f}' for g in growth[-6:])} (last six):")
print(f"  min_d^(1/d) at d = 14: {mins[-1] ** (1/14):.4f} (vs 4/3 = 1.3333).  The spread is the object of Wirsching's")
print("  'heuristic principle' (individual counts comparable to the mean); the automaton computes it exactly")
print("  level by level but does not bound it uniformly.  (The minimising unit classes are those of the small")
print("  integers 5, 7, 25, 34: thin trees are a property of specific 3-adic points, visible level by level.)")
del a, Nd, Nprev, Nu, br, ab, units

# ---------------------------------------------------------------------------
hdr("K2  The Haar-averaged size series: poles at s = 1, 2 and the reciprocal AM-GM gap")
# ---------------------------------------------------------------------------
print("A depth-d predecessor with e E-steps has size ~ a 2^d/3^e (carry ignored).  Haar average over a in Z_3:")
print("  F(s) = sum_(d,e) C(d,e) 3^-e (2^d/3^e)^-s = 1/(1 - g(s)),   g(s) = 2^-s + (3/2)^s/3.")
g = lambda s: 2.0 ** (-s) + (1.5 ** s) / 3.0
check(abs(g(1.0) - 1) < 1e-15 and abs(g(2.0) - 1) < 1e-15, "g(1) = g(2) = 1")
gp1 = -log(2) / 2 + log(1.5) / 2
check(abs(-gp1 - log(2 / sqrt(3))) < 1e-15, "-g'(1) = log(2/sqrt3)")
print(f"  g(1) = 1/2 + 1/2 = 1 (arithmetic mean of the step factors 1/2, 3/2 is 1);  g(2) = 1/4 + 3/4 = 1;")
print(f"  g < 1 exactly on (1, 2).  -g'(1) = (1/2) log(4/3) = log(2/sqrt 3) = {log(2/sqrt(3)):.6f}, the Collatz drift.")
print(f"  Residue of F at s = 1:  R = 1/log(2/sqrt 3) = {R_AMGM:.6f}.  (No other pole on Re s = 1: |g(1+it)| = 1")
print("  would need 2^-it = (3/2)^it = 1, impossible for t != 0.)  By Wiener-Ikehara the averaged number of")
print("  predecessors of relative size <= X is ~ R X, and in each dyadic window [2^n, 2^(n+1)) it is ~ R 2^n.")
# per-type residues via the linear system
def F_types(s):
    x = 2.0 ** (-s)
    y = (1.5 ** s) / 3.0
    F0 = 1 / (1 - x)
    F2 = (1 + x + y * (F0 + 1)) / ((1 + x) * (1 - x - y))
    F1 = 1 + x * F2
    return F0, F1, F2
eps = 1e-7
F0, F1, F2 = F_types(1 + eps)
res1, res2 = F1 * eps, F2 * eps
check(abs(res1 - R_AMGM) < 1e-4 and abs(res2 - 2 * R_AMGM) < 1e-4, "per-type residues R and 2R")
print(f"  Per-type residues (root uniform on a class mod 3): a = 0: 0 (the tree is the D-chain);")
print(f"  a = 1 mod 3: {res1:.5f} = R;  a = 2 mod 3: {res2:.5f} = 2R.  Average over Z_3: (0 + R + 2R)/3 = R.")
print("  Reading: orbits visit 2 mod 3 twice as often as 1 mod 3 (stationary law 2/3, 1/3 after the first odd")
print("  step), and the visit rate at size y is 1/(y log(2/sqrt3)): renewal theory with the drift as mean step.")
# dyadic window means M_n
LN2, LN3 = log(2), log(3)
def lbin(d, e):
    return lgamma(d + 1) - lgamma(e + 1) - lgamma(d - e + 1)
Mn = []
for n in range(0, 41):
    tot = 0.0
    for dd in range(0, 6000):
        # e with n <= dd - e*log2(3) < n+1  ->  (dd-n-1)/log2(3) < e <= (dd-n)/log2(3)
        elo = int((dd - n - 1) / LOG2_3) + 1
        ehi = int((dd - n) / LOG2_3 + 1e-12)
        for e in range(max(elo, 0), min(ehi, dd) + 1):
            val = dd - e * LOG2_3
            if n <= val < n + 1:
                tot += exp(lbin(dd, e) - e * LN3)
    Mn.append(tot)
print("  Dyadic window means M_n = sum over n <= d - e log2 3 < n+1 of C(d,e) 3^-e (d < 6000):")
print("   n  : " + "  ".join(f"{n:>6d}" for n in (0, 5, 10, 15, 20, 25, 30, 35, 40)))
print("  M/2^n: " + "  ".join(f"{Mn[n] / 2 ** n:6.3f}" for n in (0, 5, 10, 15, 20, 25, 30, 35, 40)))
tail = [Mn[n] / 2 ** n for n in range(30, 41)]
avg = sum(tail) / len(tail)
print(f"  mean of M_n/2^n over n = 30..40: {avg:.4f}   (R = {R_AMGM:.4f})")


def A_cum(n, dmax=20000):
    """sum of C(d,e) 3^-e over d - e log2 3 < n (all d < dmax; negligible tail)."""
    tot = 0.0
    for dd in range(dmax):
        elo = max(0, int((dd - n) / LOG2_3) + 1)
        for e in range(elo, dd + 1):
            v = lbin(dd, e) - e * LN3
            if v < -60 and e > dd / 4 + 5:
                break
            tot += exp(v)
    return tot


cum = {n: A_cum(n) / 2 ** n for n in (20, 40, 60, 80)}
print("  cumulative A(2^n)/2^n (Wiener-Ikehara predicts -> R): " +
      ", ".join(f"n={n}: {v:.4f}" for n, v in cum.items()))
print("  The ratio oscillates slowly around R (the non-lattice renewal theorem gives convergence without a")
print("  rate; the near-lattice structure comes from the convergents 19/12, 84/53, ... of log2 3).")
check(all(abs(v - R_AMGM) / R_AMGM < 0.01 for v in cum.values()), "cumulative means within 1% of R for n = 20..80")
# qx+1 comparison
def roots_q(q):
    gq = lambda s: 2.0 ** (-s) * (1 + q ** (s - 1))
    rts = []
    grid = [(i + 0.5) / 1000 for i in range(0, 4000)]
    for u, v in zip(grid, grid[1:]):
        if (gq(u) - 1) * (gq(v) - 1) < 0:
            lo, hi = u, v
            for _ in range(60):
                mid = (lo + hi) / 2
                if (gq(lo) - 1) * (gq(mid) - 1) <= 0:
                    hi = mid
                else:
                    lo = mid
            rts.append(round((lo + hi) / 2, 6))
    return rts
print("  Same series for the shortcut qx+1 map (branch probability 1/q): g_q(s) = 2^-s (1 + q^(s-1)).")
for q in (3, 5, 7, 9):
    rts = roots_q(q)
    print(f"    q = {q}: roots of g_q = 1 at s = {rts}; tree-count exponent = {min(rts)}")
rq3 = roots_q(3)
check(len(rq3) == 2 and abs(rq3[0] - 1) < 1e-5 and abs(rq3[1] - 2) < 1e-5, "q=3 roots 1 and 2")
rq5 = roots_q(5)
check(len(rq5) == 2 and abs(rq5[1] - 1) < 1e-5 and 0.64 < rq5[0] < 0.66, "q=5 roots 0.651 and 1")
print("  s = 1 is always a root (arithmetic mean of 1/2 and q/2 is (q+1)/4, and the s = 1 identity reads")
print("  1/2 + (1/q)(q/2) = 1).  For q = 3 the other root is 2 > 1, so the averaged tree has POSITIVE DENSITY;")
print("  for q >= 5 the other root is < 1 (0.651 for q = 5): averaged trees have density zero.  The automaton")
print("  sees DRIFT exactly.  It is blind to the sheet: g does not contain the offset b.")

# ---------------------------------------------------------------------------
hdr("K3  Actual tree densities of integer roots on both sheets (X = 10^7) against the averaged model")
# ---------------------------------------------------------------------------
X = 10 ** 7
targets = [4, 8, 13, 16, 20, 21, 22, 26, 29, 32, 35, 40, 44, 49, 53, 58, 100, 101, 103, 104, 105, 106,
           107, 109, 110, 112, 113, 1000, 1001, 1002, 1003, 1004, 1006, 1007, 1009, 1010, 1012, 1013, 1015, 1016,
           2000, 2001, 2002, 2003, 2005, 2006, 2008, 2009, 2011, 2012, 2014, 2015, 2017, 2018, 2020]
targets = sorted(set(targets))
check(len(targets) <= 64, "fits in a 64-bit mask")
targetsB = [a for a in range(10001, 10200) if a % 3 != 0][:64]


BUF = {}


def tree_counts(X, b, targets, cycle_pts):
    maxA = max(targets)
    bit = np.zeros(maxA + 1, dtype=np.uint64)
    for i, a in enumerate(targets):
        bit[a] = np.uint64(1) << np.uint64(i)
    if X not in BUF:  # one reusable buffer pair (keeps the peak memory low)
        BUF[X] = (np.zeros(X + 1, dtype=np.uint64), np.zeros(X + 1, dtype=bool))
    H, res = BUF[X]
    H.fill(0)
    res.fill(False)
    for c in cycle_pts:
        res[c] = True  # cycle points: their orbit (the cycle) contains no target
    for c in cycle_pts:
        check(c > maxA or bit[c] == 0, "targets avoid cycle points")
    chunk = 500_000
    for lo in range(1, X + 1, chunk):
        hi = min(X, lo + chunk - 1)
        act = np.arange(lo, hi + 1, dtype=np.int64)
        act = act[~res[lo:hi + 1]]
        x = act.copy()
        m = np.zeros(act.size, dtype=np.uint64)
        small = x <= maxA
        m[small] |= bit[x[small]]
        rounds = 0
        while act.size:
            rounds += 1
            odd = (x & 1).astype(bool)
            xo = x[odd] * 3 + b
            x >>= 1
            x[odd] = xo >> 1
            small = x <= maxA
            if small.any():
                m[small] |= bit[x[small]]
            below = x < act
            if below.any():
                idx = np.nonzero(below)[0]
                ok = res[x[idx]]
                if ok.any():
                    good = idx[ok]
                    H[act[good]] = m[good] | H[x[good]]
                    res[act[good]] = True
                    keep = np.ones(act.size, dtype=bool)
                    keep[good] = False
                    act, x, m = act[keep], x[keep], m[keep]
            check(rounds < 5000, "tree loop terminates")
    counts = [0] * len(targets)
    for lo in range(1, X + 1, 10 ** 6):
        blk = H[lo:min(X, lo + 10 ** 6 - 1) + 1]
        for i in range(len(targets)):
            counts[i] += int(((blk >> np.uint64(i)) & np.uint64(1)).sum())
    return counts


cyc_plus = [1, 2]
cyc_minus = [1, 5, 7, 10, 17, 25, 34, 37, 41, 55, 61, 68, 82, 91, 136]
tp = [a for a in targets if a not in cyc_minus]
cp = tree_counts(X, 1, tp, cyc_plus)
cm = tree_counts(X, -1, tp, cyc_minus)
print(f"  Model (Haar average over the root's 3-adic class, carries ignored): density ~ c R / a, with c = 2 at")
print("  branch roots (a = 2 mod 3 plus, a = 1 mod 3 minus), c = 1 at the other units, and a tree of")
print("  ~log2(X/a) elements (density 0) at multiples of 3.  Ratio rho = a * density / (c R).")
print("   a     | plus: count      rho   | minus: count     rho")
rows = []
for a, x1, x2 in zip(tp, cp, cm):
    cpl = 0 if a % 3 == 0 else (2 if a % 3 == 2 else 1)
    cmi = 0 if a % 3 == 0 else (2 if a % 3 == 1 else 1)
    r1 = a * x1 / X / (cpl * R_AMGM) if cpl else float('nan')
    r2 = a * x2 / X / (cmi * R_AMGM) if cmi else float('nan')
    rows.append((a, x1, r1, x2, r2))
    if a % 3 == 0:
        dchain = sum(1 for j in range(0, 64) if a * 2 ** j <= X)
        check(x1 == dchain and x2 == dchain, "multiple of 3: the tree is its D-chain")
    if a <= 40 or a in (100, 101, 103, 105, 1000, 1001, 1002, 2000, 2001, 2002):
        print(f"   {a:5d} | {x1:10d}  {r1:7.4f}   | {x2:10d}  {r2:7.4f}")
def group(lo, hi, k):
    vals = [r[k] for r in rows if lo <= r[0] <= hi and r[0] % 3 != 0]
    return sum(vals) / len(vals), min(vals), max(vals), len(vals)
for lo, hi in ((4, 20), (100, 113), (1000, 1016), (2000, 2020)):
    gp_ = group(lo, hi, 2)
    gm_ = group(lo, hi, 4)
    print(f"   roots in [{lo},{hi}] (units): mean rho plus {gp_[0]:.3f} (range {gp_[1]:.3f}..{gp_[2]:.3f}), "
          f"minus {gm_[0]:.3f} (range {gm_[1]:.3f}..{gm_[2]:.3f}), {gp_[3]} roots")
cpB = tree_counts(X, 1, targetsB, cyc_plus)
cmB = tree_counts(X, -1, targetsB, cyc_minus)
rB_p, rB_m = [], []
for a, x1, x2 in zip(targetsB, cpB, cmB):
    cpl = 2 if a % 3 == 2 else 1
    cmi = 2 if a % 3 == 1 else 1
    rB_p.append(a * x1 / X / (cpl * R_AMGM))
    rB_m.append(a * x2 / X / (cmi * R_AMGM))
def msd(v):
    mu = sum(v) / len(v)
    sd = (sum((t - mu) ** 2 for t in v) / (len(v) - 1)) ** 0.5
    return mu, sd, sd / len(v) ** 0.5
mp_, sp_, ep_ = msd(rB_p)
mm_, sm_, em_ = msd(rB_m)
print(f"   64 consecutive units a in [10001, {targetsB[-1]}]: mean rho plus {mp_:.3f} +- {ep_:.3f} (sd {sp_:.3f}),"
      f" minus {mm_:.3f} +- {em_:.3f} (sd {sm_:.3f})")
print(f"   rho range plus {min(rB_p):.3f}..{max(rB_p):.3f}, minus {min(rB_m):.3f}..{max(rB_m):.3f}")
check(abs(mp_ - 1) < 4 * ep_ + 0.1 and abs(mm_ - 1) < 4 * em_ + 0.1, "large roots: model right on average, both sheets")
print("  Reading: for large roots the averaged model is right on average on BOTH sheets (it is sheet-blind);")
print("  individual roots scatter strongly with their 3-adic expansion (Wirsching's s_n(a) is unbounded on a")
print("  dense set of 3-adic a).  Small minus-sheet roots sit inside one of the three basins, so their trees")
print("  are capped by 0.327, 0.3245 or 0.348: SHEET appears only at the small end.  Near the root cycle the model")
print("  fails by an order-one factor (tree(4) = all n >= 3 on the plus sheet: rho = 4/R = 0.575): that is the")
print("  archimedean boundary, where carries are as large as the values.")

# ---------------------------------------------------------------------------
hdr("K4  Incongruent (undecided) classes mod 3*2^k")
# ---------------------------------------------------------------------------
print("A class (n mod 3, parity word of length k) is DECIDED if the class determines a certificate of")
print("multiplier < 1.  Three certificate sets:")
print("  (I)   forward only: some prefix j <= k with 3^(a_j) < 2^j (Terras).  The mod-3 digit plays no role.")
print("  (II)  owner-refined: (I), or n = 2 mod 3 (E(n) = (2n-1)/3 < n), or an even step j-1 -> j with")
print("        T^j(n) = 2 mod 3 and (2/3) 3^(a_j)/2^j < 1 (the 'A/2 branch' E(T^j n) = (T^(j-1) n - 1)/3).")
print("  (III) minimal-chain refined: (II), plus at every window point the non-retracing predecessor branch")
print("        followed along its minimal-child chain for as many E-steps as the class determines")
print("        (T^j(n) is known mod 3^(1+a_j): odd steps MANUFACTURE 3-adic precision).")


def count_I(k):
    # words of length k with 3^(a_j) > 2^j for all j <= k
    cur = {0: 1}
    for j in range(1, k + 1):
        nxt = {}
        for a, c in cur.items():
            for bitv in (0, 1):
                a2 = a + bitv
                if 3 ** a2 > 2 ** j:
                    nxt[a2] = nxt.get(a2, 0) + c
        cur = nxt
    return sum(cur.values())


def count_II(k):
    # state: (a, t, last) with t = T^j(n) mod 3, last = parity of the previous step (None at j = 0)
    cur = {}
    for t0 in (0, 1, 2):
        if t0 == 2:
            continue  # decided at j = 0 by E(n) < n
        cur[(0, t0, -1)] = cur.get((0, t0, -1), 0) + 1
    for j in range(1, k + 1):
        nxt = {}
        for (a, t, last), c in cur.items():
            for bitv in (0, 1):
                a2 = a + bitv
                t2 = 2 if bitv else (2 * t) % 3
                if 3 ** a2 < 2 ** j:
                    continue  # forward certificate
                if bitv == 0 and t2 == 2 and 2 * 3 ** a2 < 3 * 2 ** j:
                    continue  # sibling branch certificate
                key = (a2, t2, bitv)
                nxt[key] = nxt.get(key, 0) + c
        cur = nxt
    return sum(cur.values())


# brute-force cross-check of (I) and (II) for small k
def brute_I_II(k):
    M = 3 * 2 ** k
    undec_I = undec_II = 0
    for r in range(M):
        n = r if r > 0 else M
        x, a = n, 0
        decI = decII = False
        if n % 3 == 2:
            decII = True
        prev_even = None
        for j in range(1, k + 1):
            bitv = x % 2
            a += bitv
            x = T(x)
            if 3 ** a < 2 ** j:
                decI = decII = True
                break
            if bitv == 0 and x % 3 == 2 and 2 * 3 ** a < 3 * 2 ** j:
                decII = True
        undec_I += (not decI)
        undec_II += (not decII)
    return undec_I, undec_II


for k in range(1, 13):
    bi, bii = brute_I_II(k)
    check(bi == 3 * count_I(k) and bii == count_II(k), f"DP = brute force at k={k}")
print("  DP counts agree with brute force over all classes mod 3*2^k for k <= 12.")


def chain_decided(n, k):
    """(III): class-determined certificates using the minimal-child chain of each non-retracing branch."""
    # forward window
    xs = [n]
    bits = []
    for _ in range(k):
        bits.append(xs[-1] % 2)
        xs.append(T(xs[-1]))
    a = 0
    for j in range(0, k + 1):
        if j >= 1:
            a += bits[j - 1]
            if 3 ** a < 2 ** j:
                return True
        y = xs[j]
        M = Fraction(3 ** a, 2 ** j)
        branches = []
        if j == 0:
            branches.append(("D", 2 * y, M * 2))
            if y % 3 == 2:
                branches.append(("E", (2 * y - 1) // 3, M * Fraction(2, 3)))
        elif bits[j - 1] == 0:
            if y % 3 == 2:
                branches.append(("E", (2 * y - 1) // 3, M * Fraction(2, 3)))
        else:
            branches.append(("D", 2 * y, M * 2))
        prec = 1 + a  # y is known mod 3^(1+a) from the class
        for kind, c, mult in branches:
            used = 1 if kind == "E" else 0  # trits consumed so far
            if mult < 1:
                return True
            z = c
            while True:
                # next minimal block needs z mod 3, available iff used < prec
                if used >= prec:
                    break
                tz = z % 3
                if tz == 0:
                    break
                if tz == 2:
                    z = (2 * z - 1) // 3
                    mult *= Fraction(2, 3)
                else:
                    z = (4 * z - 1) // 3
                    mult *= Fraction(4, 3)
                used += 1
                if mult < 1:
                    return True
    return False


def count_III(k, verify_lifts=False):
    M = 3 * 2 ** k
    und = 0
    for r in range(M):
        n = r if r > 0 else M
        # quick exit: forward-decided classes
        x, a, fwd = n, 0, False
        for j in range(1, k + 1):
            a += x % 2
            x = T(x)
            if 3 ** a < 2 ** j:
                fwd = True
                break
        if fwd:
            continue
        dec = chain_decided(n, k)
        if verify_lifts:
            for t in (1, 2, 5):
                check(chain_decided(n + M * t, k) == dec, "(III) decision is a class invariant")
        und += (not dec)
    return und


print()
print("   k | 3*2^k      | (I) 3|Bad_k| | (II) owner  | (III) chains | (II)/(I) | (III)/(I) | log2(I)/k log2(II)/k log2(III)/k")
rowsK = []
for k in list(range(2, 17, 2)) + [18]:
    cI = 3 * count_I(k)
    cII = count_II(k)
    cIII = count_III(k, verify_lifts=(k <= 10)) if k <= 18 else None
    rowsK.append((k, cI, cII, cIII))
    print(f"  {k:2d} | {3 * 2 ** k:10d} | {cI:12d} | {cII:11d} | {cIII:12d} | {cII / cI:8.4f} | {cIII / cI:9.4f} | "
          f"{log2(cI) / k:9.4f} {log2(cII) / k:9.4f} {log2(max(cIII,1)) / k:9.4f}")
for k in (24, 32, 40, 48, 64, 96, 128, 192, 256, 384):
    cI = 3 * count_I(k)
    cII = count_II(k)
    print(f"  {k:3d} |            | {cI:.4e} | {cII:.4e} |              | {cII / cI:8.4f} |           | "
          f"{log2(cI) / k:9.4f} {log2(cII) / k:9.4f}")
print("  (III) was checked to be a class invariant on three lifts per class for k <= 10.")
# the owner's modulus 192 explicitly
undI6, undII6 = [], []
for r in range(192):
    n = r if r > 0 else 192
    x, a, dI, dII = n, 0, False, n % 3 == 2
    for j in range(1, 7):
        bitv = x % 2
        a += bitv
        x = T(x)
        if 3 ** a < 2 ** j:
            dI = dII = True
            break
        if bitv == 0 and x % 3 == 2 and 2 * 3 ** a < 3 * 2 ** j:
            dII = True
    if not dI:
        undI6.append(r)
    if not dII:
        undII6.append(r)
check(len(undI6) == 24 and len(undII6) == 14 and 191 in undI6 and 191 not in undII6, "mod-192 lists")
print(f"  mod 192, (I): {len(undI6)} classes = 3 lifts of {sorted(set(r % 64 for r in undI6))} mod 64")
print(f"  mod 192, (II): {len(undII6)} classes {undII6}")
print("   all begin with two odd steps and are 0 or 1 mod 3.  The class 191 = -1 is DECIDED in (II): for n >= 1")
print("   in it, E(n) = (2n-1)/3 < n.  The 3-adic/2-adic point -1 itself is E's fixed point (E(-1) = -1), where")
print("   the multiplier 2/3 never becomes descent: the class is decided, the hostile point is not.")
print("  Reading: the mod-3 refinement removes a fraction of the undecided classes (about 45-50% for (II), with")
print("  the ratio drifting slowly, and about 60% for (III) at k <= 18), never the exponent: log2(count)/k")
print("  creeps up towards h(log_3 2) = 0.9500 in every column.  PROVED for (II) in the note (shifted barrier);")
print("  (III) is FINITE-EXACT only.")

# ---------------------------------------------------------------------------
hdr("K5  Rational periodic points: sign, integrality, and which lie in Bad")
# ---------------------------------------------------------------------------
print("A primitive parity word w of length p with a odd steps has the unique 2-adic periodic point")
print("x_w = c_w/(2^p - 3^a) (plus sheet; c_w > 0).  Sign(x_w) = sign(2^p - 3^a) [sign law]; x_w is an integer")
print("iff (2^p - 3^a) | c_w, a congruence modulo a number PRIME TO 6, invisible to the 3*2^k tower.")


def word_point(w, b=1):
    A, C, den = 1, 0, 1
    for bitv in w:
        if bitv:
            A, C = 3 * A, 3 * C + b * den
            den *= 2
        else:
            den *= 2
    return Fraction(C, den - A)  # x = (A x + C)/den


def is_primitive(w):
    p = len(w)
    for q in range(1, p):
        if p % q == 0 and w == w[q:] + w[:q]:
            return False
    return True


int_points = {}
bad_rational = {}
nwords = {}
for p in range(1, 19):
    seen = set()
    for bits in range(2 ** p):
        w = tuple((bits >> t) & 1 for t in range(p))
        if not is_primitive(list(w)):
            continue
        a = sum(w)
        x = word_point(w)
        # check periodicity of the point under T in Z_2: parity sequence of x is w (for x rational, 2-adic)
        nwords[p] = nwords.get(p, 0) + 1
        if x.denominator == 1:
            int_points.setdefault(p, set()).add(int(x))
        if 3 ** a > 2 ** p:
            # in Bad iff every prefix of w (one period) is above the line
            ok = True
            aa = 0
            for j in range(1, p + 1):
                aa += w[j - 1]
                if 3 ** aa < 2 ** j:
                    ok = False
                    break
            if ok:
                bad_rational.setdefault(p, []).append(x)
allint = sorted({v for s in int_points.values() for v in s})
print(f"  integer periodic points for p <= 18 (all primitive words): {allint}")
exp_int = {0, 1, 2, -1, -5, -7, -10, -17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34}
check(set(allint) == exp_int, "integer periodic points = the five known cycles (p <= 18)")
print("  = {0}, {1,2}, {-1}, {-5,-7,-10}, {-17,...,-136}: positive ones only on contracting words (2^p > 3^a),")
print("  negative ones only on expanding words.")
bad_int = sorted({int(x) for L in bad_rational.values() for x in L if x.denominator == 1})
check(bad_int == [-17, -5, -1], "integer points of Bad among periodic points: -17, -5, -1")
print(f"  periodic points lying in Bad (every prefix of the period above the line): integer ones {bad_int};")
print("   count of rational ones per period p: " + ", ".join(f"p={p}: {len(bad_rational.get(p, []))}" for p in range(1, 19)))
print("  All of them are negative.  Under x -> -x they are the positive rational cycle points of 3x-1 that lie")
print("  in Bad_-, among them the minima 1, 5, 17 of its three integer cycles.")

# ---------------------------------------------------------------------------
hdr("K6  The 16-point trap of E_{6 mod 8} is the negative image of the 3n-1 cycles")
# ---------------------------------------------------------------------------
F = {-1, -2, -4, -5, -7, -8, -10, -14, -16, -20, -29, -32, -43, -64, -86, -128}


def ES_moves(x, sheet=1, S=((6, 8),)):
    out = []
    if x % 2:
        out.append(3 * x + sheet)
    else:
        out.append(x // 2)
        if any(x % mod == r % mod for (r, mod) in S):
            out.append(3 * x + sheet)
    return out


for x in F:
    for y in ES_moves(x, 1, ((6, 8),)):
        check(y in F, "F closed under E_{6 mod 8} (plus)")
nuF = {-x for x in F}
for x in nuF:
    for y in ES_moves(x, -1, ((2, 8),)):
        check(y in nuF, "-F closed under E^-_{2 mod 8} (minus)")
# orbit of 1 under minus E_S equals -F
orb = set()
st = [1]
while st:
    x = st.pop()
    if x in orb:
        continue
    orb.add(x)
    st.extend(ES_moves(x, -1, ((2, 8),)))
check(orb == nuF, "-F = E^-_{2 mod 8}-orbit of 1")
# C-form minus cycles inside -F
check({1, 2} <= nuF and {5, 14, 7, 20, 10} <= nuF, "-F contains the 3n-1 cycles {1,2} and {5,14,7,20,10} (C-form)")
print(f"  -F = {sorted(nuF)}")
print("  is closed under the minus-sheet relaxation (x odd -> 3x-1; x even -> x/2, and 3x-1 also at x = 2 mod 8),")
print("  equals the orbit of 1 there, and contains the 3n-1 cycles 1 -> 2 -> 1 and 5 -> 14 -> 7 -> 20 -> 10 -> 5")
print("  (C-form), glued by the excursions at 2 and 10.  Its rate 2/3 is the rate of the 3n-1 cycle {5,7,10}")
print("  (2 odd steps, 3 halvings; multiplier 9/8).  So the trap that refutes HYP-9121 is SHEET, seen from")
print("  inside the plus sheet: the hostile neighbourhood of -1 is the mirror image of 3n-1's small cycles.")
# E (full) orbit of -1 reaches the -17 cycle along an explicit path
path = [-1, -2, -5, -14, -41]
for u, v in zip(path, path[1:]):
    check(v in ES_moves(u, 1, ((0, 2),)), "legal E-move on the path")
cyc17 = [-17]
y = -17
while True:
    y = y // 2 if y % 2 == 0 else 3 * y + 1
    if y == -17:
        break
    cyc17.append(y)
check(-41 in cyc17 and len(cyc17) == 18, "-41 lies on the C-form -17 cycle (18 points)")
print("  In the full relaxation E (3x+1 allowed at every even) the orbit of -1 also meets the third cycle:")
print("  -1 -> -2 -> -5 -> -14 -> -41 (3x+1 at the even -14), and -41 lies on the C-form -17 cycle.")
print("  So all three 3n-1 cycles sit, mirrored, inside the plus sheet's hostile neighbourhood of -1.")

# ---------------------------------------------------------------------------
hdr("K7  The renewal identity for basins and its neutral main term")
# ---------------------------------------------------------------------------
print("For any set R with T^-1(R) = R (a basin, or a union of basins):  R = D(R) disjoint-union E(R cap G),")
print("G = guard class (2 mod 3 plus, 1 mod 3 minus).  Hence, exactly, for every X:")
print("   N_R(X) = N_R(X/2) + #{n in R cap G : E(n) <= X}.")
# verify on census labels at 10^6 for the minus basins and the plus basin
sys.path.insert(0, ".")
Nc = 10 ** 6


def labels(N, b, cycles):
    lab = np.full(N + 1, -1, dtype=np.int8)
    for ci, c in enumerate(cycles):
        for v in c:
            if v <= N:
                lab[v] = ci
    for n0 in range(1, N + 1):
        if lab[n0] >= 0:
            continue
        x = n0
        while not (x < n0 and lab[x] >= 0):
            x = T(x, b)
        lab[n0] = lab[x]
    return lab


labm = labels(Nc, -1, [[1], [5, 7, 10], [17, 25, 37, 55, 82, 41, 61, 91, 136, 68, 34]])
labp = labels(Nc, 1, [[1, 2]])
for (b, lab, nb) in ((1, labp, 1), (-1, labm, 3)):
    guard = 2 if b == 1 else 1
    for ci in range(nb):
        inR = lab == ci
        pref = np.cumsum(inR)
        for Xv in (1000, 12345, 200000, 600000):
            lhs = int(pref[Xv])
            # E(n) = (2n - b)/3 <= X  <=>  n <= (3X + b)/2
            nmax = (3 * Xv + b) // 2
            nn = np.arange(1, nmax + 1)
            rhs = int(pref[Xv // 2]) + int((inR[1:nmax + 1] & (nn % 3 == guard)).sum())
            check(lhs == rhs, f"renewal identity b={b} basin {ci} X={Xv}")
print("  verified exactly for the plus basin of {1,2} and each of the three minus basins, X in {10^3,")
print("  12345, 2*10^5, 6*10^5}.  Main terms: if R has density delta, equidistributed mod 3 (as the census")
print("  shows for all four basins), then delta = delta/2 + (3/2)(delta/3) = delta for EVERY delta.")
print("  The equation is neutral because the arithmetic mean of the step factors is 1: the automaton's")
print("  renewal structure cannot fix a basin's density.  Density 1 (plus) versus 0.327/0.3245/0.348 (minus)")
print("  is boundary data from the small numbers, i.e. non-residue information.")

print()
print(f"[tree_counts] ALL CHECKS PASSED   ({time.time() - T0:.1f} s)", file=sys.stderr)
print("[tree_counts] ALL CHECKS PASSED")
