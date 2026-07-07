#!/usr/bin/env python3
"""
lrc_tail_diameter_floor_kps_S59.py   (kind-pasteur-2026-07-07-S59, HYP-4797)

THE TAIL DIAMETER FLOOR: pointwise monotonicity feeds the IRREDUCIBLE tail mu_{1/7}
(death-star-S1 capstone: mean detours are dead; the tail must be bounded directly).

ONE-LINE LEMMA (pointwise, rigorous):  E ⊆ F (integer sets)  ⟹  {frac(ex)} ⊆ {frac(fx)}
  ⟹  maxgap(E,x) ≥ maxgap(F,x) for EVERY x  ⟹  the superlevel sets nest:
      mu_theta(E) ≥ mu_theta(F)   for every theta.
With F = {0..D} (D = primitive diameter of E, after translate+gcd-reduce, both
mu-invariant) and {0..D} ≅ AP_{D+1} = {1..D+1} (translate):

      mu_{1/7}(E) ≥ mu_{1/7}(AP_{D+1})        [TAIL DIAMETER FLOOR]
      E[maxgap(E)] ≥ A(D+1) = E[maxgap(AP_{D+1})]   [mean version]

and the RHS are EXACT rationals via the opus-S134 Farey roof.  So the open lemma (A')
["mu_{1/7}(E) ≥ mu_{1/7}(AP_13) = 477/1078"] becomes a THEOREM-with-loss on bounded
diameter: the loss is mu(AP_13) -> mu(AP_{D+1}), and the floor stays useful as long as
mu_{1/7}(AP_{D+1}) ≥ the needed constant (m_P for k=13; thr_k + m_P for the k=8..12
union-bound legs).  This carves the bounded-diameter bulk off the open (A')/decorrelation
core with EXACT constants, k-uniformly (monotonicity never uses |E|).

This script:
  PART 1: exact mu_{1/7}(AP_n) for n = 13..130 via the Farey roof; crossings vs
          m_P = 14249/252252, vs 477/1078, and vs the six union-bound bars
          thr_k(+m_P) of THM-530/monad-HYP-4787.
  PART 2: verify tail monotonicity mu_theta(E) ≥ mu_theta(AP_{D+1}) numerically on the
          known-minimizer zoo + random families (sanity; the lemma is one line).
  PART 3: per-k bites: for each cluster size k = 8..12, the diameter range covered by
          the tail floor against its union-bound bar.
  PART 4: LEG-B groundwork (the residual): the exact pair-decorrelation constant
          |meas{frac(px)∈I, frac(qx)∈J} − |I||J|| ≤ 1/(3·pq/gcd²) — verified numerically;
          this is the explicit-constant tool for the diam > D1 residual (Chung-Erdos /
          Bonferroni windows, monad redirect target 4).
"""
from fractions import Fraction as F
import math, random

def farey(n):
    fr = set()
    for q in range(1, n + 1):
        for p in range(0, q + 1):
            fr.add(F(p, q))
    return sorted(fr)

def mu_theta_AP_roof(n, theta):
    """Exact mu_theta(AP_n) via the Farey roof (opus-S134): per-cell linear tail."""
    Fs = farey(n)
    tot = F(0)
    for a, b in zip(Fs[:-1], Fs[1:]):
        q, qp = a.denominator, b.denominator
        L = b - a
        vl, vr = F(1, q), F(1, qp)
        if vl <= theta and vr <= theta:
            continue
        if vl > theta and vr > theta:
            tot += L
            continue
        t_star = (theta - vl) * L / (vr - vl)
        tot += t_star if vl > theta else L - t_star
    return tot

def maxgap_float(E, x):
    ph = sorted((e * x) % 1.0 for e in E)
    g = max(ph[i + 1] - ph[i] for i in range(len(ph) - 1)) if len(ph) > 1 else 0.0
    return max(g, ph[0] + 1.0 - ph[-1])

def mu_numeric(E, theta, res=40000):
    return sum(1 for r in range(res) if maxgap_float(E, (r + .5) / res) > theta) / res

# ------------------------------------------------------------------ PART 1
print("=" * 92)
print("PART 1 -- exact mu_{1/7}(AP_n) (Farey roof), n = 13..130; crossings vs the bars")
print("=" * 92)
mP = F(14249, 252252)
bar_ap13 = F(477, 1078)
th = F(1, 7)
# six union-bound bars: mu needed = 1 - min meas(G_P) + m_P  (monad HYP-4787 quant form)
minGP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7), 13: F(1, 1)}
bars = {k: 1 - minGP[k] + mP for k in range(8, 14)}
print(f"  m_P = {mP} = {float(mP):.6f}")
for k in range(8, 14):
    print(f"    union-bound bar (k={k}): mu >= {float(bars[k]):.6f}")
mu_vals = {}
marks = {}
for n in list(range(13, 61)) + list(range(62, 131, 4)):
    mu_vals[n] = mu_theta_AP_roof(n, th)
crossings = {}
for label, bar in [("m_P", mP)] + [(f"thr_{k}+m_P", bars[k]) for k in range(8, 13)]:
    ns = sorted(mu_vals)
    last_ok = None
    for n in ns:
        if mu_vals[n] >= bar:
            last_ok = n
        else:
            break
    crossings[label] = last_ok
sel = [13,14,15,16,18,20,22,24,26,28,30,34,38,42,46,50,54,58,66,78,90,102,114,126,130]
for n in sel:
    if n in mu_vals:
        v = mu_vals[n]
        tags = []
        if v >= bar_ap13: tags.append(">=477/1078")
        if v >= mP: tags.append(">=m_P")
        print(f"  n={n:3d} (diam {n-1:3d}): mu_1/7(AP_n) = {float(v):.6f}   {' '.join(tags)}")
print()
for label, n_ok in crossings.items():
    print(f"  TAIL FLOOR covers diam <= {n_ok-1 if n_ok else '??'} against bar {label}"
          f"   (mu(AP_{n_ok})={float(mu_vals[n_ok]):.5f})" if n_ok else f"  bar {label}: not met even at n=13")

# ------------------------------------------------------------------ PART 2
print()
print("=" * 92)
print("PART 2 -- verify tail monotonicity mu_theta(E) >= mu_theta(AP_{D+1}) (numeric sanity)")
print("=" * 92)
zoo = {
    "AP {1..13}":                list(range(1, 14)),
    "GW {1..11,13,24}":          list(range(1, 12)) + [13, 24],
    "opus stretch":              [0] + list(range(2, 13)) + [17, 28],
    "death-star prim-sat":       [2,4,6,8,10,12,13,14,16,18,20,22,24],
    "monad record":              [2,4,6,8,10,11,12,13,14,16,18,20,22],
    "kps-S57 adversarial":       [2,6,8,10,11,12,14,16,18,20,22,26,42],
}
viol = 0
for nm, E in zoo.items():
    E0 = [e - min(E) for e in E]
    g = math.gcd(*[d for d in E0 if d != 0])
    E0 = [e // g for e in E0]
    D = max(E0)
    lhs = mu_numeric(E, 1/7)
    rhs = float(mu_theta_AP_roof(D + 1, th))
    ok = lhs >= rhs - 2e-3
    if not ok: viol += 1
    print(f"  {nm:26s} primdiam={D:3d}: mu(E)={lhs:.4f} >= mu(AP_{D+1})={rhs:.4f}  {'OK' if ok else '*** VIOL ***'}")
rng = random.Random(4797)
for _ in range(6):
    D = rng.randrange(14, 45)
    E = sorted(rng.sample(range(1, D), 11)) + [0, D]
    E = sorted(set(E))
    if len(E) != 13: continue
    lhs = mu_numeric(E, 1/7, 20000)
    rhs = float(mu_theta_AP_roof(D + 1, th))
    ok = lhs >= rhs - 3e-3
    if not ok: viol += 1
    print(f"  random diam {D:3d}:          mu(E)={lhs:.4f} >= mu(AP_{D+1})={rhs:.4f}  {'OK' if ok else '*** VIOL ***'}")
print(f"  violations: {viol}")

# ------------------------------------------------------------------ PART 3
print()
print("=" * 92)
print("PART 3 -- per-k bites: diameter range the tail floor covers vs each union-bound bar")
print("=" * 92)
print("  (monotonicity never uses |E|: a k-point cluster E of primitive diameter D has")
print("   mu_1/7(E) >= mu_1/7(AP_{D+1}); note D >= k-1 always, so the bite starts at D=k-1.)")
for k in range(8, 14):
    bar = bars[k]
    best = None
    for n in sorted(mu_vals):
        if mu_vals[n] >= bar:
            best = n
        else:
            break
    if best is None:
        print(f"  k={k}: bar {float(bar):.4f} not met at any n -- no bite")
    else:
        note = "COVERS ALL sizes (D>=k-1 start)" if best - 1 >= k - 1 else "no room"
        cov = f"diam {k-1}..{best-1}" if best - 1 >= k - 1 else "empty"
        print(f"  k={k}: bar {float(bar):.4f}  bite = {cov}   (mu(AP_n)>=bar up to n={best})")

# ------------------------------------------------------------------ PART 4
print()
print("=" * 92)
print("PART 4 -- LEG-B groundwork: exact pair-decorrelation constant (verified)")
print("=" * 92)
print("  LEMMA (Fourier, rigorous): for gcd(p,q)=1, intervals I,J:")
print("    |meas{x: frac(px) in I, frac(qx) in J} - |I||J|| <= sum_{t>=1} 2/(pi^2 t^2 pq) = 1/(3pq).")
print("  For speeds e_i, e_j with g = gcd: bound = g^2/(3 e_i e_j).")
rng = random.Random(59)
worst = 0.0
for _ in range(300):
    p = rng.randrange(2, 40); q = rng.randrange(2, 40)
    g = math.gcd(p, q); pp, qq = p // g, q // g
    if pp == qq: continue
    a, s1 = rng.random(), rng.random() * 0.3
    b, s2 = rng.random(), rng.random() * 0.3
    res = 30011
    cnt = sum(1 for r in range(res)
              if ((p * (r + .5) / res) % 1.0 - a) % 1.0 < s1 and ((q * (r + .5) / res) % 1.0 - b) % 1.0 < s2)
    joint = cnt / res
    err = abs(joint - s1 * s2)
    bound = 1 / (3 * pp * qq)
    ratio = err / bound if bound > 0 else 0
    worst = max(worst, ratio)
print(f"  300 random (p,q,I,J): max |err|/bound = {worst:.3f}  (<= 1 + numeric noise expected)")
print()
print("  COMPOSITION (the emerging 3-leg tail program for the k=13 floor):")
print("   LEG A (this session): primitive diam <= D1  =>  mu_1/7 >= mu_1/7(AP_{D+1}) >= m_P  [PROVED]")
print("   LEG B (residual):     primitive diam > D1   =>  decorrelation via the 1/(3pq) pair bound")
print("                          + staggered-window Chung-Erdos (monad target 4) -- most pairs are")
print("                          large when the diameter is large; small-pair clusters re-enter the")
print("                          multi-scale (P,E)-recursion / peel (canon GREEN).")
print("   LEG C (upstream):     the bounded-spread compact reduction (THM-527-D) bounds the honest")
print("                          domain; if its S0 <= D1 the k=13 floor is CLOSED by LEG A alone.")
print("DONE.")
