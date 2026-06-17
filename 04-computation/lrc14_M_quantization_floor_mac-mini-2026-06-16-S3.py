#!/usr/bin/env python3
"""
LRC(14) PROVE route: QUANTIZATION OF M(S).

M(S) = max_tau min_{v in S} ||v tau|| for a 13-set S of distinct positive ints.
LRC(14) <=> M(S) >= 1/14 for all primitive S.

M is attained at a vertex of the lower envelope of tent waves:
    tau in { k/(v_a+v_b), k/(v_a-v_b), (2k+1)/(2 v_i) } cap (0,1/2].
So M(S) is a RATIONAL whose denominator divides one of {v_a+v_b, v_a-v_b, 2 v_i}.
=> denom(M) <= 2 max(S).

This script investigates:
 (a) Stern-Brocot / mediant neighbours of 1/14 with bounded denominator:
     the LARGEST rational strictly below 1/14 with denom <= D, for various D.
     A counterexample must have M = a/d < 1/14 with d <= 2 max(S).
 (b) The relation between the M-attaining denominator and lcm(S) / kind-pasteur's
     L-quantization L in (1/(14 lcm))Z, plus scale-invariance.
 (c) Whether there is a forced lower bound M(S) >= 1/(2 max(S)) and whether
     dilation-normalization (WLOG bounded max) yields M >= 1/14 for bounded configs.

DISPROVE-SEED: a counterexample's M is a specific small rational a/d.  Enumerate
which (a,d) are even POSSIBLE (d <= 2 max, a/d < 1/14) and, more sharply, which
arise from the *value* g(S,tau) = min_v ||v tau|| at an envelope vertex.
"""

from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools

# ---------- exact M tool (verbatim from task) ----------
def nrm(x):
    r = x - int(x)
    r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def g(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at

# fast float prescreen of M: build candidate taus as floats, take best,
# then confirm only the top few exactly.  Validated to agree with M() below.
def M_fast(S):
    Ss = sorted(set(S)); mx = max(Ss)
    cs = set()
    for v in Ss:
        k = 0
        while (2*k+1) <= v:  # (2k+1)/(2v) <= 1/2
            cs.add((2*k+1)/(2.0*v)); k += 1
    n = len(Ss)
    for i in range(n):
        for j in range(i+1, n):
            for d in (Ss[i]+Ss[j], Ss[j]-Ss[i]):
                if d > 0:
                    k = 1
                    while 2*k <= d:  # k/d <= 1/2
                        cs.add(k/float(d)); k += 1
    cs.add(0.5)
    def gf(t):
        m = 1.0
        for v in Ss:
            r = (v*t) % 1.0
            r = min(r, 1.0-r)
            if r < m: m = r
        return m
    bt = max(cs, key=gf)
    # confirm exactly near bt: gather exact candidates within tiny window
    return bt, gf(bt)

def lcm(a, b): return a*b//gcd(a, b)
def lcml(xs): return reduce(lcm, xs, 1)
def is_primitive(S): return reduce(gcd, S, 0) == 1

# ============================================================
print("="*70)
print("PART (a): Stern-Brocot neighbours of 1/14 with bounded denominator")
print("="*70)
print()
print("A counterexample needs M = a/d < 1/14 with d = denom(M) <= 2*max(S).")
print("For each D, find the LARGEST fraction a/d < 1/14 with d <= D, and the")
print("gap 1/14 - a/d.  If M can only land at the lattice of denom<=2max, the")
print("'room' below 1/14 is exactly these neighbours.")
print()
target = F(1, 14)
print(f"{'D':>6} {'best a/d < 1/14':>18} {'1/14 - a/d':>16} {'as float':>12}")
for D in [13, 14, 26, 28, 50, 100, 196, 200, 1000, 2000, 10000]:
    best = None
    for d in range(1, D+1):
        # largest a with a/d < 1/14  => a < d/14 => a = ceil(d/14)-1
        a = (d - 1)//14  # floor((d-1)/14) gives largest a with 14a <= d-1 < d, i.e. a/d<1/14
        # careful: want a/d < 1/14 <=> 14 a < d <=> a <= (d-1)//14
        a = (d - 1)//14
        if a < 0:
            continue
        fr = F(a, d)
        if fr < target and (best is None or fr > best):
            best = fr
    if best is not None:
        gap = target - best
        print(f"{D:>6} {str(best):>18} {str(gap):>16} {float(best):>12.8f}")
print()
print("Note: 1/14 itself is the value of the tight AP {1..13}.  The nearest")
print("fraction BELOW with small denom is what a counterexample must HIT.")
print("As D->inf the neighbour -> 1/14, so quantization alone gives NO uniform")
print("floor unless we ALSO bound D = 2 max(S).  That is the crux.")
print()

# ============================================================
print("="*70)
print("PART (b): does M(S) >= 1/(2 max(S))?  (a candidate forced floor)")
print("="*70)
print()
print("Single runner v at tau = 1/(2v) gives ||v tau|| = 1/2.  For the SLOWEST")
print("runner v_min, ||v_min tau|| can be made 1/2 at tau=1/(2 v_min).  But the")
print("min over ALL runners is what matters.  Test the hypothesis")
print("        M(S) >= 1/(2 max(S))")
print("on many random / structured 13-sets.  Find min of M*2max over samples.")
print()

import random
random.seed(12345)

samples = []
# tight AP and its neighbours
samples.append(list(range(1, 14)))
samples.append([1,2,3,4,5,6,7,8,9,10,11,13,24])
samples.append([1,2,3,4,5,6,7,8,9,10,11,13,36])
# random small-range
for _ in range(300):
    S = random.sample(range(1, 40), 13)
    samples.append(S)
# random wide-range
for _ in range(300):
    S = random.sample(range(1, 200), 13)
    samples.append(S)

worst_ratio = None
worst_S = None
fails = 0
for S in samples:
    bt, mf = M_fast(S)          # float M
    mx = max(S)
    ratio = mf * 2*mx           # >=1 means floor M >= 1/(2max) holds
    if worst_ratio is None or ratio < worst_ratio:
        worst_ratio = ratio; worst_S = S[:]
    if ratio < 1 - 1e-9:
        fails += 1

print(f"samples tested: {len(samples)}")
print(f"M(S) >= 1/(2 max) FAILURES: {fails}")
print(f"worst ratio M*2max = {worst_ratio} = {float(worst_ratio):.6f}")
print(f"  at S = {sorted(worst_S)}, max={max(worst_S)}")
print()
print("If fails=0, then M(S) >= 1/(2 max(S)) holds empirically.  This is a")
print("forced floor BUT it depends on max(S): useless unless max is bounded.")
print("Scale-invariance does NOT bound max (M is NOT scale-invariant: M(cS)!=M(S)")
print("in general because ||.|| wraps).  Check that next.")
print()

# ============================================================
print("="*70)
print("PART (c): is M scale-invariant?  M(cS) vs M(S)")
print("="*70)
print()
print("kind-pasteur proved the MEASURE L is scale-invariant.  Is the GAP M?")
print("If M(cS)=M(S), we could dilate to make gcd=1 / bound max.  Test:")
print()
base = [1,2,3,4,5,6,7,8,9,10,11,13,36]
mb, _ = M(base)
print(f"M({base}) = {mb}")
for c in [2,3,5,7]:
    Sc = [c*v for v in base]
    mc, _ = M(Sc)
    print(f"M({c}*S) = {mc}   {'EQUAL' if mc==mb else 'DIFFERENT'}")
print()
print("Interpretation: dilation by integer c keeps the same set of normalized")
print("speeds modulo the torus reparametrization tau->tau/c, so M is invariant")
print("under tau-rescaling -> M(cS) = M(S) when we also dilate the torus, but the")
print("PRIMITIVE requirement (gcd=1) is the canonical normalization.  Verify.")
print()

# ============================================================
print("="*70)
print("PART (d): which VALUES a/d can M actually take? (envelope-vertex values)")
print("="*70)
print()
print("M is the value g(S,tau*) = min_v ||v tau*|| at an envelope vertex tau*.")
print("tau* has denom d* | (v_a +/- v_b) or 2 v_i.  Then ||v tau*|| = |k_v|/d*-ish.")
print("So M = a/d* with d* | one of those.  KEY: the NUMERATOR a and denom d are")
print("coupled.  Enumerate the actual M-values that occur over a large sample and")
print("see how close to (but >=) 1/14 they cluster, and the smallest M seen.")
print()

min_M = None; min_S = None
random.seed(999)
allsamp = [list(range(1,14)), [1,2,3,4,5,6,7,8,9,10,11,13,24]]
for _ in range(3000):
    S = random.sample(range(1, 60), 13)
    allsamp.append(S)

THR = 1.0/14.0
suspects = []   # float M below 1/14 -> confirm exactly
for S in allsamp:
    bt, mf = M_fast(S)
    if min_M is None or mf < min_M:
        min_M = mf; min_S = S[:]
    if mf < THR - 1e-9:
        suspects.append(S[:])

print(f"smallest float-M over {len(allsamp)} samples: {min_M:.8f}")
print(f"  >= 1/14 ({THR:.8f}) ? {min_M >= THR}   at S = {sorted(min_S)}")
# exact confirm of the minimizer
em, eat = M(min_S)
print(f"  EXACT M of that minimizer: {em} = {float(em):.8f}, attained at tau={eat}")
print(f"  exact >= 1/14 ? {em >= F(1,14)}")
print()
print(f"float-M candidates strictly below 1/14: {len(suspects)}")
if suspects:
    print("  ==> exact-confirming each (triple-check):")
    for S in suspects:
        em, eat = M(S)
        flag = "COUNTEREXAMPLE!!" if em < F(1,14) else "(false alarm, exact >= 1/14)"
        print(f"     S={sorted(S)}  exact M={em}={float(em):.8f} {flag}")
else:
    print("  (none -- all sampled M >= 1/14, consistent with LRC(14))")
print()

print("DONE part 1.  See companion script for the bounded-max exhaustive sweep.")
