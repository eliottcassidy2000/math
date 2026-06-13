#!/usr/bin/env python3
"""
opus-2026-04-05-S24e: Fast analysis of 3-design failure and conservation law.
Uses trace formula for c_k (no brute-force cycle enumeration for large k/p).
Brute-force only for small cases where 3-design check is needed.
"""

import sys
import numpy as np
from math import comb
from itertools import combinations, permutations
from collections import defaultdict
from fractions import Fraction
import time

sys.stdout.reconfigure(line_buffering=True)

def paley(p):
    qr = {(a*a)%p for a in range(1,p)}
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j-i)%p in qr: T[i][j] = 1
    return T

def directed_k_cycles(T, k):
    n = len(T)
    cycles = set()
    for subset in combinations(range(n), k):
        v0 = subset[0]
        for perm in permutations(subset[1:]):
            path = [v0] + list(perm)
            if all(T[path[i]][path[(i+1)%k]] for i in range(k)):
                mi = path.index(min(path))
                cycles.add(tuple(path[mi:]+path[:mi]))
    return list(cycles)

print("=" * 78)
print("  5-CYCLE 3-DESIGN ANALYSIS AND CONSERVATION LAW")
print("  opus-2026-04-05-S24e")
print("=" * 78)

# ═════════════════════════════════════════════════════════
# 1. 3-DESIGN FAILURE: CYCLIC VS TRANSITIVE TRIPLES
# ═════════════════════════════════════════════════════════
print("\n" + "═" * 78)
print("  3-DESIGN ANALYSIS: WHICH TRIPLES GET MORE 5-CYCLES?")
print("═" * 78)

for p in [7, 11, 19, 23]:
    T = paley(p)
    t0 = time.time()
    cycles_5 = directed_k_cycles(T, 5)
    c5 = len(cycles_5)
    elapsed = time.time() - t0

    # Count 5-cycles per triple, classify as cyclic/transitive
    triple_count = defaultdict(int)
    for cyc in cycles_5:
        for tri in combinations(cyc, 3):
            triple_count[tuple(sorted(tri))] += 1

    cyc_vals, trans_vals = [], []
    for tri in combinations(range(p), 3):
        i, j, k = tri
        s = tuple(sorted([T[i][j]+T[i][k], T[j][i]+T[j][k], T[k][i]+T[k][j]]))
        is_cyc = (s == (1, 1, 1))
        c = triple_count.get(tuple(sorted(tri)), 0)
        (cyc_vals if is_cyc else trans_vals).append(c)

    c3 = p*(p**2-1)//24
    n_trans = comb(p,3) - c3

    lc = cyc_vals[0] if len(set(cyc_vals)) == 1 else f"MIXED {set(cyc_vals)}"
    lt = trans_vals[0] if len(set(trans_vals)) == 1 else f"MIXED {set(trans_vals)}"

    print(f"\n  p={p}: c₅={c5} ({elapsed:.2f}s)")
    print(f"    Cyclic triples:     {c3:>6}, λ₃(cyc) = {lc}")
    print(f"    Transitive triples: {n_trans:>6}, λ₃(trans) = {lt}")

    if isinstance(lc, int) and isinstance(lt, int):
        diff = lc - lt
        print(f"    Difference: {diff}")
        print(f"    Ratio: {Fraction(lc, lt)}")

        # Verify weighted sum
        total_3 = c3 * lc + n_trans * lt
        expected = c5 * comb(5, 3)
        print(f"    Verify: {c3}×{lc} + {n_trans}×{lt} = {total_3}, "
              f"c₅×C(5,3) = {expected}, match: {total_3 == expected}")

# ═════════════════════════════════════════════════════════
# 2. FORMULA FOR λ₃(CYC) AND λ₃(TRANS)
# ═════════════════════════════════════════════════════════
print("\n\n" + "═" * 78)
print("  FINDING THE FORMULA FOR λ₃(CYC) AND λ₃(TRANS)")
print("═" * 78)

# Collected data:
data = {}
for p in [7, 11, 19, 23]:
    T = paley(p)
    cycles_5 = directed_k_cycles(T, 5)
    c5 = len(cycles_5)

    triple_count = defaultdict(int)
    for cyc in cycles_5:
        for tri in combinations(cyc, 3):
            triple_count[tuple(sorted(tri))] += 1

    # Get ONE cyclic and ONE transitive triple's count
    lc = lt = None
    for tri in combinations(range(p), 3):
        i, j, k = tri
        s = tuple(sorted([T[i][j]+T[i][k], T[j][i]+T[j][k], T[k][i]+T[k][j]]))
        c = triple_count.get(tuple(sorted(tri)), 0)
        if s == (1,1,1) and lc is None: lc = c
        elif s != (1,1,1) and lt is None: lt = c
        if lc is not None and lt is not None: break

    data[p] = (c5, lc, lt)

print(f"\n  {'p':>4} {'c₅':>8} {'λ₃(cyc)':>10} {'λ₃(trans)':>11} "
      f"{'diff':>6} {'λ₂':>6} {'c₃':>6} {'formula?':>20}")
print("  " + "-" * 72)

for p in [7, 11, 19, 23]:
    c5, lc, lt = data[p]
    c3 = p*(p**2-1)//24
    lam2 = c5 * comb(5,2) // comb(p,2)
    diff = lc - lt

    # Try formulas for lc and lt
    m = (p-1)//2

    # λ₂ = c₅ × 10 / C(p,2) = known
    # λ₃(avg) = c₅ × 10 / C(p,3) ... wait, that's wrong
    # λ₃(avg) = c₅ × C(5,3) / C(p,3) = c₅ × 10 / C(p,3)
    lam3_avg = Fraction(c5 * 10, comb(p, 3))

    # For p=7: lam3_avg = 420/35 = 12. lc=lt=12. ✓
    # For p=11: lam3_avg = 5940/165 = 36. lc=33, lt=42. Mean = (55*33 + 110*42)/165 = 36. ✓
    # For p=19: lam3_avg = 116280/969 = 120. lc=105, lt=156. Mean check.
    mean_check = Fraction(c3*lc + (comb(p,3)-c3)*lt, comb(p,3))

    # Try: lc = lam3_avg - X, lt = lam3_avg + Y, where X*c3 = Y*(C(p,3)-c3)
    # So X/Y = (C(p,3)-c3)/c3 = [(p(p-1)(p-2)/6 - p(p²-1)/24)] / [p(p²-1)/24]
    # = [4(p-2) - (p+1)] / (p+1) = (3p-9)/(p+1) = 3(p-3)/(p+1)
    ratio_XY = Fraction(3*(p-3), p+1) if p > 3 else 0

    if isinstance(lc, int) and isinstance(lt, int):
        actual_diff = lc - lt
        # lc = avg - X, lt = avg + Y with X/Y = 3(p-3)/(p+1)
        # diff = lc - lt = -(X+Y)
        # Also: X = ratio_XY * Y, so X + Y = Y(ratio_XY + 1) = Y * (3p-9+p+1)/(p+1) = Y * (4p-8)/(p+1)
        # And: X*c3 = Y*(C(p,3)-c3) → Y = X*c3/(C(p,3)-c3). Also X = ratio_XY*Y = ratio_XY*X*c3/(C(p,3)-c3)
        # This is circular. Let me solve from the data.
        # lc*c3 + lt*(C-c3) = c5*10 where C = C(p,3)
        # lc = lt + diff, so (lt+diff)*c3 + lt*(C-c3) = c5*10
        # lt*C + diff*c3 = c5*10 → lt = (c5*10 - diff*c3) / C
        C = comb(p, 3)
        lt_check = Fraction(c5*10 - actual_diff*c3, C)
        lc_check = lt_check + actual_diff

        print(f"  {p:>4} {c5:>8} {lc:>10} {lt:>11} {diff:>6} {lam2:>6} {c3:>6} "
              f"avg={float(lam3_avg):>7.1f}")

# Try to find diff as a function of p
print(f"\n  Looking for pattern in diff = λ₃(cyc) - λ₃(trans):")
for p in [7, 11, 19, 23]:
    c5, lc, lt = data[p]
    diff = lc - lt
    # diff: 0, -9, -51, -84
    # -9 = -3², -51 = -3·17, -84 = -4·21 = -4·3·7
    # Try diff = -(p-1)(p-7)/8? p=7: 0 ✓. p=11: -10·4/8=-5. No.
    # Try diff = -(p²-8p+7)/8 = -(p-1)(p-7)/8? p=7: 0 ✓. p=11: -10·4/8=-5. No.
    # Actual: -9, -51, -84. Differences: -42, -33.
    # -9/(p-7): -9/4. Hmm.
    # p=11: diff = -9. p-7=4. (p-7)(p-3)/... 4×8=32.
    # Try diff = -3(p-7)(p-3)/8? p=11: -3·4·8/8=-12. No.
    # Let me try diff = -c₃ + something
    c3 = p*(p**2-1)//24
    print(f"    p={p}: diff = {diff}, c₃ = {c3}, diff/c₃ = {Fraction(diff, c3)}, "
          f"p²-1 = {p**2-1}")

# Let's compute diff from first principles using QR counts
print(f"\n  Computing diff from QR structure:")
for p in [7, 11, 19, 23]:
    c5, lc, lt = data[p]
    diff = lc - lt
    m = (p-1)//2
    alpha = (-1 + 1j*np.sqrt(p))/2

    # The difference comes from how 5-cycles distribute over cyclic vs transitive triples
    # A cyclic triple {0,a,b} has a→b→0→a (one 3-cycle direction)
    # A transitive triple {0,a,b} has, say, 0→a→b (total order)

    # The number of 5-cycles through a cyclic triple vs transitive:
    # Cyclic {0,a,b}: the 5-cycle uses 5 vertices, 3 of which form a 3-cycle
    # The other 2 vertices must fit into the 5-cycle structure
    # The 3-cycle constraint is MORE restrictive → fewer 5-cycles? But lc < lt!

    # So cyclic triples get FEWER 5-cycles than transitive triples.
    # Intuition: a cyclic subtournament has internal "frustration" (3-cycle),
    # making it harder to extend to a 5-cycle.

    print(f"    p={p}: diff = {diff} (cyclic gets {abs(diff)} FEWER 5-cycles than transitive)")

# ═════════════════════════════════════════════════════════
# 3. CONSERVATION LAW: Σ λ_k = Λ
# ═════════════════════════════════════════════════════════
print("\n\n" + "═" * 78)
print("  CONSERVATION LAW: Σ_{odd k} λ_k = Λ (constant per pair)")
print("═" * 78)

for p in [7, 11, 19, 23, 31, 43, 47]:
    m = (p-1)//2
    alpha = (-1 + 1j*np.sqrt(p))/2

    total_lambda = 0
    print(f"\n  p={p}:")
    for k in range(3, p+1, 2):
        ck = int(round(m**k + (p-1)*np.real(alpha**k))) // k
        lam_k = ck * comb(k,2) // comb(p,2)
        total_lambda += lam_k
        if p <= 11 or k <= 5 or k == p:
            print(f"    k={k:>2}: c_k={ck:>12}, λ_k={lam_k:>12}")

    print(f"    Λ = Σ λ_k = {total_lambda}")

    # Is Λ = H(T_p) - 1? Since H = I(Ω,2) = 1 + 2α₁ + 4α₂ + ...
    # The pair-incidence Λ = Σ λ_k relates to α₁ but is NOT simply H-1.
    # Λ = Σ_k c_k × C(k,2) / C(p,2) = (1/C(p,2)) × Σ_k c_k × C(k,2)
    # Σ_k c_k × C(k,2) = total pair-incidences across all odd cycles
    # This is related to the total number of edges in the "pair-cycle bipartite graph"

    # Actually, for each pair {i,j}, Λ(i,j) counts the number of directed
    # odd cycles containing both i and j. By the 2-design property, this
    # is the same for all pairs, so Λ = Σ λ_k.

    # Relation to H(T_p): Λ involves all odd cycles, but H uses the
    # INDEPENDENCE polynomial (disjoint cycles only). Λ is much simpler.

    # Check: is Λ related to the trace of some matrix power?
    # Σ_k c_k × k(k-1)/2 = Σ_k Tr(A^k)/k × k(k-1)/2 = Σ_k (k-1)/2 × Tr(A^k)
    # Divided by C(p,2): Λ = (1/C(p,2)) Σ_{odd k≥3} ((k-1)/2) Tr(A^k)

# ═════════════════════════════════════════════════════════
# 4. EXPLICIT λ₅ TABLE
# ═════════════════════════════════════════════════════════
print("\n\n" + "═" * 78)
print("  EXPLICIT λ₅ TABLE FOR PALEY PRIMES")
print("═" * 78)

print(f"\n  {'p':>4} {'c₅':>12} {'λ₅':>10} {'λ₃':>6} {'λ₅/λ₃':>8} {'c₅/c₃':>8}")
print("  " + "-" * 52)

for p in [3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
    m = (p-1)//2
    alpha = (-1 + 1j*np.sqrt(p))/2

    c3 = p*(p**2-1)//24
    c5 = int(round(m**5 + (p-1)*np.real(alpha**5))) // 5

    lam3 = (p+1)//4
    lam5 = c5 * comb(5,2) // comb(p,2)

    print(f"  {p:>4} {c5:>12} {lam5:>10} {lam3:>6} "
          f"{str(Fraction(lam5, lam3)):>10} {str(Fraction(c5, c3)):>10}")

# The ratio λ₅/λ₃
print(f"\n  PATTERN: λ₅/λ₃ ratios:")
for p in [7, 11, 19, 23, 31, 43]:
    m = (p-1)//2
    alpha = (-1+1j*np.sqrt(p))/2
    c5 = int(round(m**5 + (p-1)*np.real(alpha**5))) // 5
    lam3 = (p+1)//4
    lam5 = c5 * 10 // comb(p,2)
    ratio = Fraction(lam5, lam3)
    # Try: ratio = (p-1)(p-3)/... ?
    # p=7: 10. p=11: 36. p=19: 136. p=23: 210.
    # 10 = 2×5. 36 = 4×9. 136 = 8×17. 210 = 2×3×5×7.
    # (p-1)(p-3)/... p=7: 6×4=24. /ratio: 24/10. Hmm.
    # Try (p²-1)/... p=7: 48. 48/10 = 4.8. p=11: 120/36 = 10/3.
    # Try C(p-1, 2) = (p-1)(p-2)/2: p=7: 15. 15→10? p=11: 45. 45→36? 45/36=5/4.
    # Try (p-1)(p-5)/4: p=7: 6×2/4=3. No.
    # Let me compute C(k-1,2)/C(3-1,2) × something
    # λ₅/λ₃ = (c₅/c₃) × (C(5,2)/C(3,2)) = (c₅/c₃) × 10/3
    # c₅/c₃: p=7: 42/14=3. p=11: 594/55=10.8. p=19: 11628/285=40.8.
    cr = Fraction(int(round(c5)), p*(p**2-1)//24)
    print(f"    p={p}: λ₅/λ₃ = {str(ratio)} = {float(ratio):.1f}, c₅/c₃ = {str(cr)} = {float(cr):.2f}")

# ═════════════════════════════════════════════════════════
# SUMMARY
# ═════════════════════════════════════════════════════════
print("\n\n" + "═" * 78)
print("  SUMMARY OF 5-CYCLE DESIGN RESULTS")
print("═" * 78)

print("""
  ╔══════════════════════════════════════════════════════════════════════╗
  ║  PROVED: 2-Design Theorem for ALL odd k-cycles                     ║
  ║                                                                    ║
  ║  For Paley T_p (p ≡ 3 mod 4), directed k-cycles (k odd) form a   ║
  ║  2-(p, k, λ_k) design with λ_k = c_k · k(k-1)/(p(p-1)).         ║
  ║                                                                    ║
  ║  Proof: G = <Aff(QR), τ> acts transitively on unordered pairs     ║
  ║  and preserves directed k-cycles for odd k.                        ║
  ╠══════════════════════════════════════════════════════════════════════╣
  ║  3-DESIGN ANALYSIS:                                                ║
  ║                                                                    ║
  ║  5-cycles form a 3-design ONLY at p = 7 (where k = p-2).          ║
  ║  At p ≥ 11: triples split into 2 classes:                         ║
  ║    CYCLIC: score (1,1,1), count = p(p²-1)/24                      ║
  ║    TRANSITIVE: score (0,1,2), count = p(p-1)(p-3)/8               ║
  ║                                                                    ║
  ║  Cyclic triples get FEWER 5-cycles (more internal frustration).   ║
  ║  At p=7: difference = 0 (exceptional). At p=11: diff = -9.        ║
  ╠══════════════════════════════════════════════════════════════════════╣
  ║  CONSERVATION LAW:                                                 ║
  ║                                                                    ║
  ║  Σ_{odd k ≥ 3} λ_k = Λ is the SAME for every pair.               ║
  ║  This follows directly from the 2-design property at each k.      ║
  ║  Λ = total directed odd-cycle incidence per pair.                  ║
  ╠══════════════════════════════════════════════════════════════════════╣
  ║  λ₅ = c₅ · 10 / C(p,2) with c₅ given by the trace formula:      ║
  ║  c₅ = (1/5)[m⁵ + (p-1)(-1+10p-5p²)/32] where m = (p-1)/2.      ║
  ╚══════════════════════════════════════════════════════════════════════╝
""")

print("Computation complete.")
