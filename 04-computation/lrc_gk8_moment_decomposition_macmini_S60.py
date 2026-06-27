#!/usr/bin/env python3
"""gK8/gK9/gK11 concentration bound = which MOMENT drives the binding?

S59 redirect: the proof needs the BOUND max_E L_y(E) <= scale*cap (CRUX 1, the
covering-moment / bounded-core positivity), NOT the census. The Delsarte duals in
moment form (LRCFactorialAtom.lean):
   L_yK8  = 10 S0 - 10 S1 + 10 S2 -  9 S3 + 6 S4   (k=8)
   L_yK9  = 18 S0 - 13 S1 +  8 S2 -  3 S3           (k=9,10)
   L_yK11 =  6 S0 -  3 S1 +    S2                   (k=11,12,13)
where S_r = E[C(N,r)] = sum_t C(t,r) q_t is the r-th factorial (binomial) moment of
the miss-count N = #{inner sectors j=1..6 that are empty}.  S2 = sum over the 15=C(6,2)
sector-PAIRS of meas{both empty} -- the pairwise co-emptiness (covariance) layer.

CLAIM under test (mac-mini-S60): the *binding direction* of every L_y dual is the
+S2 (pairwise) term; the concentration GAP (consec binding minus wide) is driven by S2.
This places CRUX 1 on the Clebsch "support-six" pairwise carrier (HYP-2890/2892), which
is the cut-space Cayley graph of K5 (S40) -- the 2-adic / H4 face of the tower.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(7)

def sector_of(p): return int((p % 1) * 7)

def missdist(E):
    """q[t] = meas{tau : exactly t of the 6 inner sectors are empty}, t=0..6."""
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): b.add(F(m, 7 * e))
    b = sorted(b); q = [F(0)] * 7
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0: continue
        covered = set(sector_of(e * ((x0 + x1) / 2)) for e in E)
        t = 7 - len(covered)          # inner sectors empty (sector 0 always hit by e=0)
        if 0 <= t <= 6: q[t] += x1 - x0
    return q

def fac_moments(q):
    """S_r = sum_t C(t,r) q_t for r=0..6."""
    from math import comb
    return [sum(comb(t, r) * q[t] for t in range(7)) for r in range(7)]

# Delsarte duals in moment form (S_0..S_6 coefficients)
DUALS = {
    8:  [10, -10, 10, -9, 6, 0, 0],
    9:  [18, -13,  8, -3, 0, 0, 0],
    10: [18, -13,  8, -3, 0, 0, 0],
    11: [6,   -3,  1,  0, 0, 0, 0],
}
READOUT = {  # direct q-form, for cross-check
    8:  lambda q: 10*q[0] + q[3] + 10*q[6],
    9:  lambda q: 18*q[0] + 5*q[1] + 2*q[4] + 3*q[5],
    10: lambda q: 18*q[0] + 5*q[1] + 2*q[4] + 3*q[5],
    11: lambda q: 6*q[0] + 3*q[1] + q[2] + q[5] + 3*q[6],
}

def Ly(coeffs, S): return sum(c * s for c, s in zip(coeffs, S))

print("=" * 78)
print(" gK8/gK9/gK11: WHICH MOMENT DRIVES THE CONCENTRATION BINDING?")
print("=" * 78)

for k in (8, 9, 10, 11):
    coeffs = DUALS[k]
    consk = tuple(range(k))                       # binding bounded config (consec)
    qc = missdist(consk); Sc = fac_moments(qc)
    # cross-check moment form == readout
    assert Ly(coeffs, Sc) == READOUT[k](qc), f"moment-form mismatch k={k}"
    Lc = Ly(coeffs, Sc)

    # comparison configs: single-far, wide doublet, even-AP (all same size k)
    far = consk[:-1] + (21,)                       # single-far (r=1), the wide binding
    wide = consk[:-1] + (40,)                       # wider single-far
    evenAP = tuple(2 * i for i in range(k))         # dilated AP
    comps = {"single-far(21)": far, "single-far(40)": wide, "even-AP(2i)": evenAP}

    print(f"\n--- k={k}   dual S0..S6 = {coeffs} ---")
    print(f"  consec_{k}: L_y = {float(Lc):.4f}   S = [{', '.join(f'{float(s):.3f}' for s in Sc)}]")
    for name, cfg in comps.items():
        qx = missdist(cfg); Sx = fac_moments(qx); Lx = Ly(coeffs, Sx)
        # per-term contribution to the GAP  Lc - Lx = sum_r coeff_r (S^consec_r - S^x_r)
        gap_terms = [coeffs[r] * (Sc[r] - Sx[r]) for r in range(7)]
        gap = sum(gap_terms)
        # which term dominates the gap?
        driver = max(range(7), key=lambda r: abs(gap_terms[r]))
        frac_S2 = float(gap_terms[2] / gap) if gap != 0 else float('nan')
        print(f"  vs {name:16s}: L_y={float(Lx):.4f}  gap={float(gap):+.4f}  "
              f"per-term ΔS·c = [{', '.join(f'{float(g):+.3f}' for g in gap_terms)}]")
        print(f"      -> driver = S{driver} (coeff {coeffs[driver]:+d}); "
              f"S2-share of gap = {frac_S2:+.1%}")

print("\n" + "=" * 78)
print(" PAIRWISE STRUCTURE of S2 (the Clebsch 'support-six' layer)")
print("=" * 78)
# S2 = sum over the 15 = C(6,2) inner-sector pairs of meas{both empty}.
# Verify: does consec MAXIMIZE the pairwise co-emptiness S2 among size-k configs?
def pair_coemptiness(E):
    """For each inner-sector pair (i,j) in C(6,2), meas{tau: both i,j empty}. Returns dict+sum."""
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): b.add(F(m, 7 * e))
    b = sorted(b)
    pairs = list(itertools.combinations(range(1, 7), 2))   # 15 inner-sector pairs
    acc = {p: F(0) for p in pairs}
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0: continue
        covered = set(sector_of(e * ((x0 + x1) / 2)) for e in E)
        empty = set(range(1, 7)) - covered
        for p in pairs:
            if p[0] in empty and p[1] in empty: acc[p] += x1 - x0
    return acc, sum(acc.values())

for k in (8, 9, 10):
    consk = tuple(range(k))
    _, S2c = pair_coemptiness(consk)
    # random bounded + single-far competitors
    best = S2c; bestcfg = "consec"
    trials = [("single-far21", consk[:-1] + (21,)), ("even-AP", tuple(2*i for i in range(k)))]
    for _ in range(40):
        cfg = tuple(sorted([0] + random.sample(range(1, 15), k - 1)))
        trials.append(("rand", cfg))
    n_beat = 0
    for name, cfg in trials:
        _, s2 = pair_coemptiness(cfg)
        if s2 > best + F(1, 10**9): best = s2; bestcfg = f"{name}{cfg}"; n_beat += 1
    print(f"  k={k}: S2(consec)={float(S2c):.4f}   "
          f"configs beating consec's S2: {n_beat}/{len(trials)}   "
          f"{'consec MAXIMIZES pairwise co-emptiness' if n_beat==0 else 'max at '+bestcfg}")
print("\n(15 = C(6,2) inner-sector pairs = 2^4 - 1 nonzero vectors of the Clebsch cut-space (Z/2)^4 = K5 cut-space, S40.)")
