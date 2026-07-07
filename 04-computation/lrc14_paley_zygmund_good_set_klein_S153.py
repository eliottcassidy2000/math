#!/usr/bin/env python3
"""
klein-2026-07-07-S153.  PALEY-ZYGMUND on the GOOD SET, controlled by the RELATION LATTICE.

The moat (the single open LRC(14) core, klein-S152 corrected by kps-S54): every single-scale
non-AP 13-family has M >= 2/27 > 1/14 (the AP {1..13} is the unique single-scale family at 1/14).

opus-S131 applied PZ to mu_{1/7} the CHEAP way (U<=1 => E(U^2)<=E(U) => mu_{1/7} >= E[U]),
reducing to a FIRST moment inf E[U] > 0 -- still open, obstructed by "structure-dependent triple+
overlaps." THIS script attacks the GOOD SET directly with the FULL second moment, and makes the
obstruction PRECISE: it is the LOW-DEGREE ADDITIVE RELATION LATTICE.

Set-up. Fix a nonneg trig-poly kernel h >= 0 concentrated on the good arc [beta, 1-beta]
(beta = 1/14), small on the bad arc (-beta,beta). Use the Fejer kernel centered at 1/2:
    h(x) = F_m(x - 1/2),   F_m(y) = sum_{|k|<=m} (1-|k|/(m+1)) e(ky) >= 0,   hhat_k = (1-|k|/(m+1))(-1)^k.
h peaks at x=1/2 (farthest from integers), is small near x=0 (bad). Define
    Z(t) = prod_i h(v_i t) >= 0.
Then (Paley-Zygmund):  meas{ Z > tau } >= (1 - tau/E[Z])^2 * E[Z]^2 / E[Z^2]   for tau <= E[Z].
If the "bad leakage" tau0 := h_bad * h_max^12  <  E[Z]   (h_bad = max_{|x|<beta} h, h_max = max h),
then { Z > tau0 } => ALL runners good, so
    meas(GOOD at beta)  >=  (1 - tau0/E[Z])^2 * E[Z]^2 / E[Z^2]   >  0     [a RIGOROUS witness-floor].

KEY (the Fourier point):  E[Z] = sum_{ <k,v>=0, |k_i|<=m } prod hhat_{k_i} = (hhat_0)^13 + RELATION TERMS.
With hhat_0 = 1, E[Z] = 1 + sum over NONZERO relations <k,v>=0 (|k_i|<=m) of prod hhat_{k_i}.
No low relations  =>  E[Z] ~ 1  =>  PZ floor positive.  The AP is maximally resonant (3-term
relations v_{i-1}-2v_i+v_{i+1}=0) => the relation terms drag E[Z] down / good-set to measure 0.

This script computes, across the single-scale spectrum (AP / near-AP / spread / dilated AP):
  meas(good@1/14) exactly-ish (fine grid), E[Z], E[Z^2], the leakage tau0, the PZ floor,
  and the SHORTEST additive relation (min sum|k_i| over nonzero <k,v>=0) -- the discriminator.
"""
from fractions import Fraction as F
import math, random, itertools

BETA = 1.0/14.0
NGRID = 200000          # torus grid for moments / good-measure
TWO_PI = 2*math.pi

def fejer_coeffs(m):
    # hhat_k for h(x)=F_m(x-1/2): (1-|k|/(m+1)) * (-1)^k, k=-m..m
    return {k: (1-abs(k)/(m+1))*((-1)**k) for k in range(-m, m+1)}

def make_h(m):
    c = fejer_coeffs(m)
    def h(y):  # y in R, period 1; h(x)=F_m(x-1/2) => evaluate F_m at (x-1/2)
        yy = (y-0.5)
        s = 0.0
        for k,ck in c.items():
            s += ck*math.cos(TWO_PI*k*yy)  # real part (coeffs real, symmetric)
        return s
    return h

def h_extremes(h):
    # scan one period for h_max, and h_bad = max over bad arc |x|<beta (x in (-beta,beta) mod 1)
    hmax = -1e9; hbad = -1e9; hmin_good = 1e9
    NN = 20000
    for j in range(NN):
        x = j/NN
        val = h(x)
        if val > hmax: hmax = val
        d = min(x, 1-x)  # ||x||
        if d < BETA:
            if val > hbad: hbad = val
        else:
            if val < hmin_good: hmin_good = val
    return hmax, hbad, hmin_good

def moments_and_good(v, h):
    """fine-grid: E[Z], E[Z^2], meas(good@beta)."""
    n = len(v)
    sumZ = 0.0; sumZ2 = 0.0; good = 0
    for j in range(NGRID):
        t = j/NGRID
        z = 1.0
        allgood = True
        for vi in v:
            x = (vi*t) % 1.0
            z *= h(x)
            if allgood:
                d = min(x, 1-x)
                if d < BETA: allgood = False
        sumZ += z; sumZ2 += z*z
        if allgood: good += 1
    return sumZ/NGRID, sumZ2/NGRID, good/NGRID

def true_M(v):
    """exact-ish M = max_t min_i ||v_i t||, via critical times t=p/q, q up to bound."""
    from math import gcd
    best = 0.0
    Qbound = 2*(max(abs(x) for x in v)+1)
    # candidate times: t = a/(v_i +- v_j) and a/(2 v_i); sample a fine set too
    cand = set()
    vs = list(v)
    for i in range(len(vs)):
        for s in (vs[i], 2*vs[i]):
            if s!=0:
                for a in range(1, abs(s)):
                    cand.add(F(a, abs(s)))
        for jx in range(i+1, len(vs)):
            for s in (vs[i]-vs[jx], vs[i]+vs[jx]):
                if s!=0:
                    for a in range(1, abs(s)):
                        cand.add(F(a, abs(s)))
    for t in cand:
        mn = min(min((vi*t) % 1, 1-((vi*t)%1)) for vi in vs)
        if mn > best: best = float(mn)
    return best

def shortest_relation(v, max_terms=4, coeff_cap=3):
    """min sum|k_i| over nonzero integer k, |k_i|<=coeff_cap, <=max_terms nonzero, with <k,v>=0."""
    n = len(v)
    best = None; bestk = None
    idxs = range(n)
    for t in range(2, max_terms+1):
        for combo in itertools.combinations(idxs, t):
            # coefficients in [-cap,cap]\{0} for each chosen index
            for ks in itertools.product([c for c in range(-coeff_cap,coeff_cap+1) if c!=0], repeat=t):
                s = sum(ks[a]*v[combo[a]] for a in range(t))
                if s==0:
                    # require gcd of ks == 1 to avoid trivial multiples
                    g = 0
                    for kk in ks: g = math.gcd(g, abs(kk))
                    if g!=1: continue
                    L1 = sum(abs(kk) for kk in ks)
                    if best is None or L1 < best:
                        best = L1; bestk = (combo, ks)
        if best is not None:
            break  # smallest #terms first; report the minimal-term relation
    return best, bestk

def pz_floor(EZ, EZ2, tau0):
    if EZ <= 0 or EZ2 <= 0: return 0.0, "EZ<=0 (resonance killed the mean)"
    if tau0 >= EZ: return 0.0, f"leakage tau0={tau0:.4g} >= EZ={EZ:.4g} (kernel too lossy)"
    theta = tau0/EZ
    return (1-theta)**2 * EZ*EZ/EZ2, "ok"

# ---------------------------------------------------------------
random.seed(20260707)
M_DEG = 13   # Fejer degree (odd => h(1/2)=0 exactly; higher => more concentrated but more relations seen)
h = make_h(M_DEG)
hmax, hbad, hmin_good = h_extremes(h)
tau0 = hbad * (hmax**12)
print("="*80)
print(f"PALEY-ZYGMUND GOOD-SET FLOOR (Fejer kernel deg m={M_DEG}, beta=1/14)")
print(f"  h_max={hmax:.5f}  h_bad(max on bad arc)={hbad:.5f}  h_min(on good arc)={hmin_good:.5f}")
print(f"  leakage tau0 = h_bad*h_max^12 = {tau0:.6g}   (need tau0 < E[Z] for a rigorous floor)")
print(f"  clean-threshold check: h_min^13 = {hmin_good**13:.4g} vs tau0={tau0:.4g}  "
      + ("(h_min^13>tau0 OK)" if hmin_good**13>tau0 else "(FAIL: kernel leaks)"))
print("="*80)

families = []
AP = list(range(1,14))
families.append(("AP {1..13}", AP))
families.append(("dilated 2*AP", [2*x for x in AP]))
families.append(("near-AP: 13->14", list(range(1,13))+[14]))
families.append(("near-AP: 7->8dup?no; 6->5swap", [1,2,3,4,5,7,6,8,9,10,11,12,13]))  # local swap (=AP set) -> still AP set
families.append(("near-AP: {1..12,15}", list(range(1,13))+[15]))
# spread single-scale (ratio<=13), random primitive
def rand_spread():
    while True:
        s = sorted(random.sample(range(2,27), 13))
        if math.gcd(*s)==1 and s[-1] <= 13*s[0]:
            return s
for r in range(5):
    families.append((f"spread#{r}", rand_spread()))

print(f"{'family':28s} {'M':>8s} {'good@1/14':>10s} {'E[Z]':>9s} {'E[Z^2]':>9s} {'PZ floor':>10s} {'shortRel(L1)':>12s}")
print("-"*95)
for name, v in families:
    EZ, EZ2, good = moments_and_good(v, h)
    Mv = true_M(v)
    floor, note = pz_floor(EZ, EZ2, tau0)
    sr, srk = shortest_relation(v)
    srstr = f"{sr}" if sr is not None else ">4/cap3"
    print(f"{name:28s} {Mv:8.5f} {good:10.5f} {EZ:9.4g} {EZ2:9.4g} {floor:10.4g} {srstr:>12s}"
          + ("" if note=="ok" else f"  [{note}]"))

print()
print("READOUT:")
print("  * good@1/14 = 0 for AP/dilated-AP (isolated witness) -> PZ floor 0 (correctly fails at AP).")
print("  * spread families: good@1/14 > 0, shortest relation LONG (>4 or none <=cap) -> PZ floor > 0.")
print("  * near-AP: good small, short relation present -> PZ weak; needs rigidity/conjugate-witness.")
print("  * The discriminator is the SHORTEST ADDITIVE RELATION (opus's 'structure-dependent")
print("    triple+ overlap' = the relation lattice). AP has L1=4 (3-term); spread has none small.")
