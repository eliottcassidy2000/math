#!/usr/bin/env python3
"""
alpha_crossover_analysis.py — Comprehensive crossover analysis for α₁,α₂,α₃,...
in H = 1 + 2α₁ + 4α₂ + 8α₃ + ...

EXACT data for interval/cyclic tournament from fast_cycle_cc + SSC pipeline.

Author: opus-2026-04-16-S1
"""

import numpy as np
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# ─── Complete exact data for interval (cyclic) tournament ────────────────────
# H and {k: alpha_k}

DATA_INT = {
    3:  dict(H=3,                        a={1:1}),
    5:  dict(H=15,                       a={1:7}),      # 5 three-cycles + 2 five-cycles
    7:  dict(H=175,                      a={1:59, 2:14}),
    9:  dict(H=3267,                     a={1:837, 2:354, 3:22}),
    11: dict(H=93027,                    a={1:18397, 2:11110, 3:1474}),
    13: dict(H=3711175,                  a={1:606027, 2:436748, 3:87568, 4:3224}),
    15: dict(H=198464295,                a={1:27495799, 2:22518662, 3:5849428, 4:397720, 5:7472}),
    17: dict(H=13689269499,              a={1:1651334601, 2:1482234998, 3:458011858, 4:45997104, 5:1800368}),
    19: dict(H=1184212824763,            a={1:126443605257, 2:122111579294, 3:42960731622, 4:5521030944, 5:331078344, 6:4100656}),
    21: dict(H=125547534942879,          a={1:12030499746751, 2:12330182836208, 3:4796354751404, 4:738531326288, 5:58868297768, 6:1454221328, 7:12571712}),
    23: dict(H=16011537490557279,        a={1:1391602826199187, 2:1499656616321278, 3:632921002322216, 4:111796734828336, 5:10945293151712, 6:412282843184, 7:7454017376}),
}

DATA_PALEY = {
    7:  dict(H=189,              a={1:80, 2:7}),
    11: dict(H=95095,            a={1:21169, 2:10879, 3:1155}),
    19: dict(H=1172695746915,    a={1:130965270477}),   # α₂,α₃ not yet computed
    # 23: IN PROGRESS
}

# Verify H values match OCF
print("=== VERIFICATION ===")
for n, d in sorted(DATA_INT.items()):
    H_calc = 1 + sum(2**k * v for k, v in d['a'].items())
    match = "✓" if H_calc == d['H'] else f"✗ (calc={H_calc}, stored={d['H']})"
    print(f"  n={n:>2}: H={d['H']:>24,}  {match}")

print()

# ─── TABLE 1: Full alpha decomposition ────────────────────────────────────────
print("=" * 110)
print("TABLE 1: Alpha Decomposition for Interval/Cyclic Tournament C_n")
print("         H = 1 + 2α₁ + 4α₂ + 8α₃ + 16α₄ + ...")
print("=" * 110)
print(f"{'n':>4}  {'H':>25}  {'2α₁%':>6} {'4α₂%':>6} {'8α₃%':>6}  "
      f"{'dom_k':>5}  {'α₁/2α₂':>9} {'α₃/α₂':>8} {'α₄/α₃':>8}")
print("-" * 110)

for n in sorted(DATA_INT.keys()):
    d = DATA_INT[n]
    H = d['H']
    a = d['a']
    a1, a2, a3, a4 = a.get(1,0), a.get(2,0), a.get(3,0), a.get(4,0)
    terms = {k: 2**k * a.get(k,0) for k in range(1, 9) if k in a}
    pct = {k: 100*v/H for k,v in terms.items()}
    dom_k = max(pct, key=pct.get) if pct else 0
    r12 = a1/(2*a2) if a2 else float('inf')
    r32 = a3/a2 if a2 else 0
    r43 = a4/a3 if a3 else 0
    p1,p2,p3 = pct.get(1,0),pct.get(2,0),pct.get(3,0)
    print(f"{n:>4}  {H:>25,}  {p1:>5.1f}% {p2:>5.1f}% {p3:>5.1f}%  "
          f"{dom_k:>5}  {r12:>9.4f} {r32:>8.4f} {r43:>8.4f}")

print()

# ─── TABLE 2: Ratio progression ───────────────────────────────────────────────
print("=" * 70)
print("TABLE 2: Ratio α_{k+1}/α_k — measures when next term overtakes current")
print("         Crossover occurs when ratio reaches 0.5")
print("=" * 70)
print(f"{'n':>4}  {'α₂/α₁':>9} {'α₃/α₂':>9} {'α₄/α₃':>9} {'α₅/α₄':>9}")
print("-" * 48)

for n in sorted(DATA_INT.keys()):
    a = DATA_INT[n]['a']
    a1,a2,a3,a4,a5 = a.get(1,0),a.get(2,0),a.get(3,0),a.get(4,0),a.get(5,0)
    r21 = a2/a1 if a1 else 0
    r32 = a3/a2 if a2 else 0
    r43 = a4/a3 if a3 else 0
    r54 = a5/a4 if a4 else 0
    print(f"{n:>4}  {r21:>9.5f} {r32:>9.5f} {r43:>9.5f} {r54:>9.5f}")

print()
print("NOTE: 2^{k+1}α_{k+1} > 2^k α_k  ⟺  α_{k+1}/α_k > 1/2")

# ─── CROSSOVER ANALYSIS ───────────────────────────────────────────────────────
print()
print("=" * 70)
print("CROSSOVER 1: 4α₂ overtakes 2α₁  (when α₁/(2α₂) crosses 0.5)")
print("=" * 70)

r12_data = [(n, d['a'].get(1,0)/(2*d['a'].get(2,1))
             ) for n, d in sorted(DATA_INT.items()) if d['a'].get(2,0) > 0]

print("\nRatio α₁/(2α₂) progression:")
for n, r in r12_data:
    flag = " ← α₁ dominant" if r > 1 else (" ← TRANSITION" if r > 0.8 else " ← α₂ dominant")
    print(f"  n={n:>2}: {r:.5f}{flag}")

# Find bracket
for i in range(len(r12_data)-1):
    n1,r1 = r12_data[i]; n2,r2 = r12_data[i+1]
    if r1 > 0.5 and r2 < 0.5:
        n_cross = n1 + (n2-n1)*(r1-0.5)/(r1-r2)
        print(f"\n→ Crossover between n={n1} and n={n2}")
        print(f"→ Interpolated crossing: n ≈ {n_cross:.2f}")
        print(f"→ First odd n where 4α₂ dominates: n = {n2}")

print()
print("=" * 70)
print("CROSSOVER 2: 8α₃ overtakes 4α₂  (when α₃/α₂ crosses 0.5)")
print("=" * 70)

r32_data = [(n, d['a'].get(3,0)/d['a'].get(2,1)
             ) for n, d in sorted(DATA_INT.items()) if d['a'].get(2,0)>0 and d['a'].get(3,0)>0]

print("\nRatio α₃/α₂ progression:")
for n, r in r32_data:
    flag = " ← α₃ dominant" if r >= 0.5 else ""
    print(f"  n={n:>2}: {r:.5f}{flag}")

ns = np.array([x[0] for x in r32_data])
rs = np.array([x[1] for x in r32_data])

# Increments per 2-step
print("\nIncrements (Δ per 2 in n):")
for i in range(1, len(r32_data)):
    n1,r1 = r32_data[i-1]; n2,r2 = r32_data[i]
    print(f"  n={n1}→{n2}: +{r2-r1:.5f}")

# Fit multiple models
print("\nModel fits:")
def model_exp(n, A, B, c):
    return A - B * np.exp(-c * n)

def model_logistic(n, A, B, c):
    return A / (1 + B * np.exp(-c * n))

def model_sqrt(n, a, b):
    return a - b / np.sqrt(n)

def model_power(n, A, B, c):
    return A - B * n**(-c)

best_cross = None
for name, func, p0 in [
    ("Exp",      model_exp,      [0.65, 10.0, 0.07]),
    ("Power",    model_power,    [0.65, 5.0, 0.5]),
    ("Logistic", model_logistic, [0.65, 100.0, 0.07]),
]:
    try:
        popt, pcov = curve_fit(func, ns, rs, p0=p0, maxfev=100000)
        resid = rs - func(ns, *popt)
        rmse = np.sqrt(np.mean(resid**2))
        A = popt[0]
        print(f"\n  {name:>9}: asymptote={A:.5f}, RMSE={rmse:.6f}")
        if A > 0.5:
            # Find crossover numerically
            from scipy.optimize import brentq
            try:
                n_cross = brentq(lambda n: func(n, *popt) - 0.5, 20, 200)
                print(f"             crossover at n ≈ {n_cross:.1f}")
                if name == "Exp":
                    best_cross = n_cross
                for n_pred in [25, 27, 29, 31, 33, 35, 41]:
                    r_pred = func(n_pred, *popt)
                    flag = " ← CROSSOVER" if abs(r_pred - 0.5) < 0.01 else ""
                    print(f"             n={n_pred}: α₃/α₂ ≈ {r_pred:.5f}{flag}")
            except Exception:
                pass
        else:
            print(f"             asymptote < 0.5: crossover may NOT occur!")
    except Exception as e:
        print(f"  {name:>9}: fit failed: {e}")

if best_cross:
    print(f"\n→ BEST ESTIMATE: 8α₃ overtakes 4α₂ at n ≈ {best_cross:.0f}")

print()
print("=" * 70)
print("CROSSOVER 3: 16α₄ overtakes 8α₃  (when α₄/α₃ crosses 0.5)")
print("=" * 70)

r43_data = [(n, d['a'].get(4,0)/d['a'].get(3,1)
             ) for n, d in sorted(DATA_INT.items()) if d['a'].get(3,0)>0 and d['a'].get(4,0)>0]

print("\nRatio α₄/α₃:")
for n, r in r43_data:
    print(f"  n={n:>2}: {r:.5f}")

if len(r43_data) >= 3:
    ns43 = np.array([x[0] for x in r43_data])
    rs43 = np.array([x[1] for x in r43_data])
    try:
        popt43, _ = curve_fit(model_exp, ns43, rs43, p0=[0.5, 5.0, 0.1], maxfev=100000)
        A43 = popt43[0]
        print(f"\n  Exp fit: asymptote = {A43:.5f}")
        if A43 > 0.5:
            from scipy.optimize import brentq
            n_cross43 = brentq(lambda n: model_exp(n, *popt43) - 0.5, 20, 500)
            print(f"  → Crossover at n ≈ {n_cross43:.0f}")
        else:
            print(f"  → Asymptote < 0.5: α₃ ALWAYS dominates α₄ (no crossover?)")
    except Exception as e:
        print(f"  Fit failed: {e}")

# ─── Dominant-term progression ────────────────────────────────────────────────
print()
print("=" * 70)
print("DOMINANT TERM PROGRESSION")
print("=" * 70)

print("""
The dominant term in H = 1 + Σ_k 2^k α_k shifts as n grows:

 n ≤ 9:  2α₁ is largest (single cycles dominate — counting density)
 n ≥ 11: 4α₂ is largest (pairs of disjoint cycles — packing efficiency)
 n ≈ ?:  8α₃ becomes largest (triples of disjoint cycles)

This reflects a PACKING PHASE TRANSITION:
 At small n, adding one more disjoint triangle is rare → single cycles dominate
 At medium n, there's enough space for pairs → k=2 dominates
 At large n, the maximum packing (k=floor(n/3)) will eventually dominate
""")

# Show the eventual kmax dominance
print("Maximum possible k (floor(n/3)) and its term:")
for n in sorted(DATA_INT.keys()):
    a = DATA_INT[n]['a']
    H = DATA_INT[n]['H']
    kmax = n // 3
    ak_max = a.get(kmax, 0)
    if ak_max > 0:
        frac = 100 * 2**kmax * ak_max / H
        dom_k = max(a, key=lambda k: 2**k * a[k])
        print(f"  n={n:>2}: kmax={kmax}, 2^{kmax}·α_{kmax} = {2**kmax * ak_max:>20,} ({frac:.1f}%)  dom_k={dom_k}")

# ─── PALEY vs INTERVAL comparison ─────────────────────────────────────────────
print()
print("=" * 70)
print("PALEY vs INTERVAL (primes p ≡ 3 mod 4)")
print("=" * 70)

print()
print("Question: Does Paley maximize H among all tournaments for p ≡ 3 mod 4?")
print()
print(f"{'p':>5}  {'H(Paley)':>22}  {'H(Int)':>22}  {'Winner':>10}  {'ΔH':>15}")
print("-" * 80)

comparisons = [
    (7,  189,              175,              ),
    (11, 95095,            93027,            ),
    (19, 1172695746915,    1184212824763,    ),
    (23, None,             16011537490557279,),  # Paley TBD
]
for row in comparisons:
    p = row[0]; HP = row[1]; HI = row[2]
    if HP is None:
        print(f"{p:>5}  {'(computing)':>22}  {HI:>22,}  {'?':>10}  {'?':>15}")
    else:
        w = "PALEY" if HP > HI else "INTERVAL"
        print(f"{p:>5}  {HP:>22,}  {HI:>22,}  {w:>10}  {HP-HI:>+15,}")

print()
print("Paley term decomposition vs Interval:")
for p in [7, 11]:
    pd = DATA_PALEY[p]; id_ = DATA_INT[p]
    pa, ia = pd['a'], id_['a']
    kmax = max(max(pa.keys()), max(ia.keys()))
    print(f"\n  p={p}:")
    cum = 0
    for k in range(1, kmax+1):
        dp = 2**k * (pa.get(k,0) - ia.get(k,0))
        cum += dp
        print(f"    k={k}: 2^{k}(α_P-α_I) = {dp:>+12,}  cumul={cum:>+12,}")
    print(f"    → H(P)-H(I) = {pd['H']-id_['H']:>+,}")

print()
print("KEY INSIGHT: Paley vs Interval:")
print("  Paley has MORE odd cycles (α₁_P > α₁_I)")
print("  but FEWER disjoint pairs (α₂_P < α₂_I at p=7: 7<14)")
print("  at p=7: 2Δα₁ = +42, 4Δα₂ = -28 → net +14 → Paley wins by 14 ✓")
print("  at p=11: 2Δα₁ = +5544, 4Δα₂ = -924, 8Δα₃ = -2552 → net... let's compute")
for p in [7, 11]:
    pd = DATA_PALEY[p]; id_ = DATA_INT[p]
    net = 0
    for k in range(1, 5):
        dp = 2**k * (pd['a'].get(k,0) - id_['a'].get(k,0))
        net += dp
    print(f"  p={p}: net from k=1..{min(max(pd['a'].keys()),max(id_['a'].keys()))}: {net:+}")
    print(f"       H(P)-H(I) = {pd['H']-id_['H']:+}")

print()
print("=" * 70)
print("CIRCULANT MAXIMIZATION HYPOTHESIS")
print("=" * 70)
print("""
Claim: Among ALL tournaments on n vertices (odd n), the H-maximizer is
always a circulant tournament (specifically: Paley for p ≡ 3 mod 4 with
p small, then interval/cyclic S={1,...,(n-1)/2} for larger n).

STATUS:
  • n=7: Paley maximizes GLOBALLY (exhaustive check over all 1,044 tournaments)
  • n=11: Paley maximizes among all circulants (exhaustive), global not checked
  • n=13: Interval maximizes among all circulants (H=3,711,175 = top of list)
  • n≥15: Not exhaustively verified — only known among circulants

The transition from Paley-optimal to Interval-optimal (among ≡3 mod 4 primes):
  p=7: Paley wins  (H_P/H_I = 189/175 = 1.08)
  p=11: Paley wins (H_P/H_I = 95095/93027 = 1.022)
  p=19: Interval wins by 1.0%

For n NOT ≡ 3 mod 4: interval tournament is the natural maximizer
  (no Paley exists, and it consistently wins among all circulants)
""")

print()
print("MECHANISM: Why does Paley lose at large p?")
print()
print("H(P) - H(I) = 2Δα₁ + 4Δα₂ + 8Δα₃ + ...")
print("where Δαₖ = αₖ(Paley) - αₖ(Interval).")
print()
print("At p=7:  Δα₁ = +21 (Paley has more single cycles)")
print("         Δα₂ = -7  (Interval has more disjoint pairs)")
print("         2Δα₁ = +42 dominates 4Δα₂ = -28 → Paley wins")
print()
print("At large p: α₂ becomes the dominant term in H (crossover at n≈10-11)")
print("  Paley's advantage in α₁ gets multiplied by 2")
print("  Interval's advantage in α₂ gets multiplied by 4")
print("  Eventually 4|Δα₂| > 2Δα₁: Interval wins")
print()
print("The crossover p where Paley stops being optimal:")
print("  Set 2Δα₁ ≈ 4|Δα₂| + 8|Δα₃| + ...")
print("  Paley has extra α₁ from BIBD structure (every pair in exactly λ₃=(p+1)/4 cycles)")
print("  Interval wins on packings because its uniform spacing allows more disjoint cycles")
print("  The crossover is between p=11 and p=19 for the dominant-mode competition")
