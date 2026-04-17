#!/usr/bin/env python3
"""
alpha_ratio_trends.py — Analytic crossover analysis for H = 1 + 2α₁ + 4α₂ + 8α₃ + ...

Uses all known exact data to:
1. Map the crossover between dominant terms in H
2. Fit growth models and predict when α₃ overtakes α₂
3. Analyze Paley vs Interval structure for ≡3 mod 4 primes
4. Compute exact small-n alphas using fast DP on Omega(T)

Author: opus-2026-04-16-S1
"""

import numpy as np
from itertools import combinations
import time

# ─── Known exact alpha data (interval/cyclic tournament) ─────────────────────
#     n: H, α₁, α₂, α₃, α₄, α₅, α₆, α₇

INTERVAL_DATA = {
    # From exhaustive computation
    3: dict(H=3, a={1:1}),
    5: dict(H=15, a={1:6}),  # Only 3-cycles; two disjoint 3-cycles need 6 vertices
    7: dict(H=175, a={1:59, 2:14}),
    9: dict(H=3267, a={}),  # Need to compute
    11: dict(H=93027, a={1:18397, 2:11110, 3:1474}),
    13: dict(H=3711175, a={}),  # Need to compute
    15: dict(H=198464295, a={}),  # Need to compute
    17: dict(H=13689269499, a={1:1651334601, 2:1482234998, 3:458011858,
                                4:45997104, 5:1800368}),
    19: dict(H=1184212824763, a={1:126443605257, 2:122111579294, 3:42960731622,
                                  4:5521030944, 5:331078344, 6:4100656}),
    21: dict(H=125547534942879, a={1:12030499746751, 2:12330182836208, 3:4796354751404,
                                    4:738531326288, 5:58868297768, 6:1454221328, 7:12571712}),
    23: dict(H=16011537490557279, a={1:1391602826199187, 2:1499656616321278, 3:632921002322216,
                                       4:111796734828336, 5:10945293151712, 6:412282843184,
                                       7:7454017376}),
}

PALEY_DATA = {
    7: dict(H=189, a={1:80, 2:7}),
    11: dict(H=95095, a={1:21169, 2:10879, 3:1155}),
    19: dict(H=1172695746915, a={1:130965270477}),  # only α₁ known
    # 23: TO BE COMPUTED (Paley SSC in progress)
}

# ─── Compute missing small-n alphas via exact DP ─────────────────────────────

def make_interval_adj(n):
    m = (n - 1) // 2
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for k in range(1, m + 1):
            j = (i + k) % n
            A[i][j] = 1
    return A

def make_paley_adj(p):
    QR = set(pow(a, 2, p) for a in range(1, p))
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in QR:
                A[i][j] = 1
    return A

def find_odd_cycles_fast(A):
    """Find all directed odd cycles by DP. Returns list of frozensets."""
    n = len(A)
    cycles = set()

    # DP: dp[S][v] = set of paths from min(S) using exactly vertex set S ending at v
    # S is a bitmask
    # Only track minimal starting vertex to avoid duplicates

    for start in range(n):
        # BFS/DFS paths starting from 'start', only using vertices >= start
        # A path is a sequence of distinct vertices starting at 'start'
        stack = [(start, (start,), 1 << start)]
        while stack:
            v, path, vmask = stack.pop()
            for u in range(n):
                if A[v][u]:
                    if u == start:
                        L = len(path)
                        if L >= 3 and L % 2 == 1:
                            cycles.add(frozenset(path))
                    elif u > start and not (vmask & (1 << u)):
                        stack.append((u, path + (u,), vmask | (1 << u)))
    return list(cycles)

def compute_alpha_decomp(A):
    """Compute {k: alpha_k} and H via Omega independence numbers."""
    n = len(A)
    cycles = find_odd_cycles_fast(A)
    cycle_masks = [sum(1 << v for v in cyc) for cyc in cycles]
    m = len(cycle_masks)
    kmax = n // 3

    alphas = {}
    for k in range(1, kmax + 1):
        count = 0
        for combo in combinations(range(m), k):
            union = 0
            ok = True
            for idx in combo:
                if union & cycle_masks[idx]:
                    ok = False
                    break
                union |= cycle_masks[idx]
            if ok:
                count += 1
        if count == 0:
            break
        alphas[k] = count

    H = 1 + sum(2**k * v for k, v in alphas.items())
    return alphas, H

# Compute missing entries
print("Computing missing small-n alpha decompositions...")
for n in [3, 5, 9, 13]:
    if not INTERVAL_DATA[n]['a']:
        print(f"  n={n} interval...")
        t0 = time.time()
        A = make_interval_adj(n)
        alphas, H_comp = compute_alpha_decomp(A)
        INTERVAL_DATA[n]['a'] = alphas
        if H_comp != INTERVAL_DATA[n]['H'] and INTERVAL_DATA[n]['H'] != 0:
            print(f"    WARNING: H mismatch! computed={H_comp}, expected={INTERVAL_DATA[n]['H']}")
        elif INTERVAL_DATA[n]['H'] == 0:
            INTERVAL_DATA[n]['H'] = H_comp
        print(f"    H={H_comp}, alphas={dict(sorted(alphas.items()))}, t={time.time()-t0:.1f}s")

for n in [7]:
    if not PALEY_DATA[n]['a']:
        print(f"  n={n} Paley...")
        A = make_paley_adj(n)
        alphas, H_comp = compute_alpha_decomp(A)
        PALEY_DATA[n]['a'] = alphas

print()

# ─── Main crossover table ─────────────────────────────────────────────────────
print("=" * 96)
print("COMPLETE ALPHA DECOMPOSITION TABLE (interval/cyclic tournament)")
print("=" * 96)
print(f"{'n':>4} {'H':>24} {'2α₁%':>7} {'4α₂%':>7} {'8α₃%':>7} {'dom_k':>6} "
      f"{'α₁/2α₂':>9} {'α₃/α₂':>8} {'α₄/α₃':>8}")
print("-" * 96)

for n in sorted(INTERVAL_DATA.keys()):
    d = INTERVAL_DATA[n]
    H = d['H']
    if H == 0 or not d['a']:
        continue
    a = d['a']
    a1 = a.get(1, 0)
    a2 = a.get(2, 0)
    a3 = a.get(3, 0)
    a4 = a.get(4, 0)

    # Term percentages
    terms = {k: 2**k * a.get(k, 0) for k in range(1, 8)}
    pct = {k: 100*v/H for k, v in terms.items() if v > 0}
    dom_k = max(pct, key=pct.get) if pct else 0

    r12 = a1/(2*a2) if a2 > 0 else float('inf')
    r32 = a3/a2 if a2 > 0 else 0
    r43 = a4/a3 if a3 > 0 else 0

    p12 = pct.get(1, 0)
    p22 = pct.get(2, 0)
    p32 = pct.get(3, 0)

    print(f"{n:>4} {H:>24,} {p12:>6.1f}% {p22:>6.1f}% {p32:>6.1f}% {dom_k:>6} "
          f"{r12:>9.4f} {r32:>8.4f} {r43:>8.4f}")

print()
print("NOTE: dom_k = k with largest 2^k·αₖ term in H")
print("      α₁/2α₂ = 0.5 marks crossover where 4α₂ first exceeds 2α₁")
print("      α₃/α₂  = 0.5 marks crossover where 8α₃ first exceeds 4α₂")

# ─── Crossover analysis ───────────────────────────────────────────────────────
print()
print("=" * 70)
print("CROSSOVER 1: When does 4α₂ overtake 2α₁? (ratio α₁/2α₂ = 0.5)")
print("=" * 70)

# Collect data points where both a1 and a2 are known
r12_pts = []
for n in sorted(INTERVAL_DATA.keys()):
    d = INTERVAL_DATA[n]
    a = d['a']
    if a.get(1, 0) > 0 and a.get(2, 0) > 0:
        r12 = a[1] / (2 * a[2])
        r12_pts.append((n, r12))
        flag = ""
        if abs(r12 - 0.5) < 0.06:
            flag = " ← NEAR CROSSOVER"
        elif r12 < 0.5:
            flag = " ← α₂ dominant"
        print(f"  n={n:>2}: α₁/(2α₂) = {r12:.6f}{flag}")

# Find bracketing pair
cross1_n = None
for i in range(len(r12_pts)-1):
    n1, r1 = r12_pts[i]
    n2, r2 = r12_pts[i+1]
    if r1 > 0.5 and r2 < 0.5:
        cross1_n = n1 + (n2-n1)*(r1-0.5)/(r1-r2)
        print(f"\n  → Crossover between n={n1} and n={n2}")
        print(f"  → Interpolated crossover: n ≈ {cross1_n:.2f}")

print()
print("=" * 70)
print("CROSSOVER 2: When does 8α₃ overtake 4α₂? (ratio α₃/α₂ = 0.5)")
print("=" * 70)

r32_pts = []
for n in sorted(INTERVAL_DATA.keys()):
    d = INTERVAL_DATA[n]
    a = d['a']
    if a.get(2, 0) > 0 and a.get(3, 0) > 0:
        r32 = a[3] / a[2]
        r32_pts.append((n, r32))
        flag = ""
        if r32 >= 0.5:
            flag = " ← α₃ dominant"
        print(f"  n={n:>2}: α₃/α₂ = {r32:.6f}{flag}")

# Check if we've crossed yet
if r32_pts and r32_pts[-1][1] < 0.5:
    print(f"\n  → NOT YET crossed (last known: n={r32_pts[-1][0]}, ratio={r32_pts[-1][1]:.4f})")

    # Fit exponential approach to asymptote
    ns = np.array([x[0] for x in r32_pts])
    rs = np.array([x[1] for x in r32_pts])

    from scipy.optimize import curve_fit

    def model_exp(n, A, B, c):
        return A - B * np.exp(-c * n)

    def model_power(n, A, B, c):
        return A - B * n**(-c)

    def model_sqrt(n, A, B):
        return A - B / np.sqrt(n)

    print("\n  Fitting asymptotic model r(n) = A - B·exp(-c·n):")
    try:
        popt, pcov = curve_fit(model_exp, ns, rs, p0=[0.65, 10.0, 0.07], maxfev=100000)
        A, B, c = popt
        perr = np.sqrt(np.diag(pcov))
        print(f"    A={A:.4f}±{perr[0]:.4f}, B={B:.4f}±{perr[1]:.4f}, c={c:.5f}±{perr[2]:.5f}")
        print(f"    Asymptote r(∞) = {A:.4f}")
        if A > 0.5:
            n_cross = -np.log((A-0.5)/B) / c
            print(f"    → Crossover predicted at n ≈ {n_cross:.1f}")
            for n_pred in [25, 27, 29, 31, 33, 35, 41, 51]:
                r_pred = model_exp(n_pred, A, B, c)
                print(f"    n={n_pred}: α₃/α₂ ≈ {r_pred:.4f}")
        else:
            print(f"    → Asymptote {A:.4f} < 0.5: crossover MAY NOT OCCUR (model extrapolation)")
    except Exception as e:
        print(f"    Fit failed: {e}")

    print("\n  Fitting power model r(n) = A - B·n^(-c):")
    try:
        popt2, _ = curve_fit(model_power, ns, rs, p0=[0.65, 5.0, 0.5], maxfev=100000)
        A2, B2, c2 = popt2
        print(f"    A={A2:.4f}, B={B2:.4f}, c={c2:.4f}")
        print(f"    Asymptote r(∞) = {A2:.4f}")
        if A2 > 0.5:
            n_cross2 = (B2/(A2-0.5))**(1/c2)
            print(f"    → Crossover predicted at n ≈ {n_cross2:.1f}")
    except Exception as e:
        print(f"    Fit failed: {e}")

print()
print("=" * 70)
print("CROSSOVER 3: When does 16α₄ overtake 8α₃? (ratio α₄/α₃ = 0.5)")
print("=" * 70)

r43_pts = []
for n in sorted(INTERVAL_DATA.keys()):
    d = INTERVAL_DATA[n]
    a = d['a']
    if a.get(3, 0) > 0 and a.get(4, 0) > 0:
        r43 = a[4] / a[3]
        r43_pts.append((n, r43))
        print(f"  n={n:>2}: α₄/α₃ = {r43:.6f}")

if r43_pts:
    ns = np.array([x[0] for x in r43_pts])
    rs = np.array([x[1] for x in r43_pts])
    try:
        popt, pcov = curve_fit(model_exp, ns, rs, p0=[0.4, 5.0, 0.1], maxfev=100000)
        A, B, c = popt
        print(f"\n  Fitted: A={A:.4f}, B={B:.4f}, c={c:.5f}")
        print(f"  Asymptote r(∞) = {A:.4f}")
        if A > 0.5:
            n_cross = -np.log((A-0.5)/B) / c
            print(f"  → Crossover at n ≈ {n_cross:.1f}")
        else:
            print(f"  → Asymptote {A:.4f} < 0.5: 16α₄ may NEVER overtake 8α₃")
    except Exception as e:
        print(f"  Fit failed: {e}")

# ─── General ratio pattern: α_{k+1}/α_k ────────────────────────────────────
print()
print("=" * 70)
print("RATIO PATTERN: α_{k+1}/α_k across k at fixed n")
print("=" * 70)
print("(NB: the dominant term shifts as n grows)")
for n in [17, 19, 21, 23]:
    d = INTERVAL_DATA[n]
    a = d['a']
    H = d['H']
    print(f"\n  n={n}, H={H:,}")
    kmax_k = max(a.keys()) if a else 0
    for k in range(1, kmax_k):
        if a.get(k, 0) > 0 and a.get(k+1, 0) > 0:
            r = a[k+1] / a[k]
            # Also show contribution fraction
            c_k = 100 * 2**k * a[k] / H
            c_k1 = 100 * 2**(k+1) * a[k+1] / H
            print(f"    α_{k+1}/α_{k} = {r:.5f}  "
                  f"(2^{k}α_{k}: {c_k:.1f}%, 2^{k+1}α_{k+1}: {c_k1:.1f}%)")

# ─── Paley vs Interval: deep comparison ──────────────────────────────────────
print()
print("=" * 70)
print("PALEY vs INTERVAL: detailed comparison")
print("=" * 70)
print(f"{'p':>4} {'H_P':>20} {'H_I':>20} {'Winner':>10} {'ΔH':>15}")
print("-" * 72)

p_series_HPI = [
    (7, 189, 175),
    (11, 95095, 93027),
    (19, 1172695746915, 1184212824763),
]
for p, HP, HI in p_series_HPI:
    w = "Paley" if HP > HI else "Interval"
    print(f"{p:>4} {HP:>20,} {HI:>20,} {w:>10} {HP-HI:>+15,}")

print(f"{'23':>4} {'(computing)':>20} {16011537490557279:>20,} {'?':>10} {'?':>15}")

print()
print("Paley vs Interval term-by-term breakdown:")
for p, pd, id_ in [(7, PALEY_DATA[7], INTERVAL_DATA[7]),
                    (11, PALEY_DATA[11], INTERVAL_DATA[11])]:
    HP = pd['H']
    HI = id_['H']
    pa, ia = pd['a'], id_['a']
    print(f"\n  p={p}:")
    kmax = max(max(pa.keys() if pa else [0]), max(ia.keys() if ia else [0]))
    for k in range(1, kmax+1):
        pk = pa.get(k, 0)
        ik = ia.get(k, 0)
        delta = 2**k * (pk - ik)
        print(f"    k={k}: α_P={pk:>8}, α_I={ik:>8}, "
              f"2^{k}·Δα = {delta:>+10,}  "
              f"({'+' if delta>=0 else '-'}Paley)")
    print(f"    Total: H(P)-H(I) = {HP-HI:+,}")

# Key insight: why does Paley lose at p=19?
print()
print("KEY QUESTION: Why does Paley lose at p=19?")
print("  α₁(P₁₉) - α₁(I₁₉) = {:,}  (Paley has MORE 3,5,7-cycles)".format(
    130965270477 - 126443605257))
print("  But 2·(α₁(P)-α₁(I)) = {:,}".format(2*(130965270477-126443605257)))
print("  ΔH(P-I) = {:,}  (Interval wins by {:,})".format(
    1172695746915 - 1184212824763,
    1184212824763 - 1172695746915))
print("  So higher αₖ (k≥2) must compensate: 4·Δα₂ + 8·Δα₃ + ... ≥ {:,}".format(
    abs((1172695746915 - 1184212824763) - 2*(130965270477-126443605257))))

# ─── Dominant-term shift and H-maximization principle ─────────────────────────
print()
print("=" * 70)
print("MAXIMIZATION PRINCIPLE")
print("=" * 70)
print("""
H = 1 + 2α₁ + 4α₂ + 8α₃ + ... = Σ_k 2^k α_k

To maximize H at large n, the dominant term is 2^k_dom * α_{k_dom}.
Strategy: maximize whichever αₖ is the current dominant contribution.

Observations:
  • Small n (≤ 11): α₁ dominates → Paley wins (maximizes c₃,c₅,c₇ via BIBD)
  • Medium n (≈ 13-19): α₁ and α₂ comparable → transition region
  • Large n (≥ 21): α₂ dominates (currently, at n=21-23) → Interval wins
    because interval maximizes 'disjoint cycle packings' better than Paley

As n grows:
  • The dominant term keeps shifting to larger k
  • k_dom(n) → floor(n/3) for very large n (all vertices packed into triangles)
  • Maximum is achieved by tournaments that maximize packings of disjoint cycles

For EVEN more concisely:
  H ≈ (# ways to tile all vertices with disjoint odd cycles) × 2^{n/3}
  for large n. This is the MAXIMUM PACKING (α_{floor(n/3)}) contribution.
""")
