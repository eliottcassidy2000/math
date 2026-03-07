#!/usr/bin/env python3
"""
WHY does bc_coeff = 2 * t5_coeff in ALL moments?

The key insight: in the sigma decomposition, t5 and bc ONLY enter through
sub-tournament Hamiltonian path counts H(T[S]) for |S| >= 5.

For |S| = 5: H = 1 + 2*c3(S) + 2*c5(S)  (no bc term)
For |S| = 6: H = 1 + 2*c3(S) + 2*c5(S) + 4*bc(S)
For |S| = 3+3 (disjoint): H(S1)*H(S2) = (1+2*cyc1)(1+2*cyc2) -> bc from cross term

The t5 contribution comes from:
  - Single 5-sets: sigma((4,)) with coefficient 2 per 5-set in OCF
  - Composite patterns involving a 5-set component

The bc contribution comes from:
  - Single 6-sets: sigma((5,)) with coefficient 4 per 6-set in OCF
  - Disjoint 3+3 patterns: sigma((2,2)) with coefficient 8 per pair
  - Composite patterns: sigma((3,2)) with coefficient 8*(n-6) per pair

CLAIM: The ratio bc/t5 = 2 arises from the OCF weight function I(Omega, 2).
Specifically, t5 appears with weight 2 (single odd cycle), while bc appears
with weight 4 = 2^2 (pair of odd cycles). The factor of 2 is 2^(#components-1)
where bc has 2 components vs 1 for t5.

But wait -- t5 and bc enter through DIFFERENT sigma patterns at different
SIGMA_k levels. The 2:1 ratio must hold at each SIGMA_k level separately.

Let me verify: at each k, what are the t5 and bc contributions?

opus-2026-03-06-S28
"""
from itertools import combinations
from math import factorial, comb
from collections import Counter

def count_position_patterns(n, k):
    if k == 0: return {(): 1}
    patterns = Counter()
    for S in combinations(range(n-1), k):
        pos = sorted(S)
        comps, comp = [], [pos[0]]
        for i in range(1, len(pos)):
            if pos[i] == comp[-1] + 1:
                comp.append(pos[i])
            else:
                comps.append(len(comp))
                comp = [pos[i]]
        comps.append(len(comp))
        patterns[tuple(sorted(comps, reverse=True))] += 1
    return dict(patterns)

def sigma_full(pattern, n):
    sizes = list(pattern)
    total_verts = sum(s + 1 for s in sizes)
    free = n - total_verts
    if free < 0: return (0, 0, 0, 0)
    num_size1 = sizes.count(1)
    big_sizes = sorted([s for s in sizes if s > 1], reverse=True)
    F = factorial(free)
    def pair_prod(start):
        p, r = 1, start
        for i in range(num_size1):
            p *= comb(r, 2); r -= 2
        return p
    if len(big_sizes) == 0: return (factorial(n) // (2**len(sizes)), 0, 0, 0)
    if big_sizes == [2]:
        pp = pair_prod(n-3); return (F*comb(n,3)*pp, F*2*pp, 0, 0)
    if big_sizes == [3]:
        pp = pair_prod(n-4); return (F*comb(n,4)*pp, F*2*(n-3)*pp, 0, 0)
    if big_sizes == [4]:
        pp = pair_prod(n-5); return (F*comb(n,5)*pp, F*2*comb(n-3,2)*pp, F*2*pp, 0)
    if big_sizes == [5]:
        pp = pair_prod(n-6)
        return (F*comb(n,6)*pp, F*2*comb(n-3,3)*pp, F*2*(n-5)*pp, F*4*pp)
    if big_sizes == [2, 2]:
        pp = pair_prod(n-6)
        return (F*comb(n,3)*comb(n-3,3)*pp, F*4*comb(n-3,3)*pp, 0, F*8*pp)
    if big_sizes == [3, 2]:
        pp = pair_prod(n-7)
        return (F*comb(n,4)*comb(n-4,3)*pp,
                F*(2*(n-3)*comb(n-4,3)+2*comb(n-3,4))*pp,
                0, F*8*(n-6)*pp)
    return None

# =====================================================================
# Check bc/t5 ratio at each SIGMA_k level
# =====================================================================
print("=" * 70)
print("bc/t5 ratio per SIGMA_k level")
print("=" * 70)

for n in [7, 9, 11, 13]:
    print(f"\n  n={n}:")
    for k in range(2, 8):
        pats = count_position_patterns(n, k)
        t5_total = 0
        bc_total = 0
        for pat, cnt in pats.items():
            r = sigma_full(pat, n)
            if r is None:
                if any(s >= 4 for s in pat) or pat.count(2) >= 2 or (3 in pat and 2 in pat):
                    pass  # might have t5 or bc
                continue
            t5_total += cnt * r[2]
            bc_total += cnt * r[3]
        if t5_total > 0 or bc_total > 0:
            ratio = bc_total / t5_total if t5_total > 0 else "inf"
            print(f"    SIGMA_{k}: t5={t5_total}, bc={bc_total}, bc/t5={ratio}")

# =====================================================================
# Algebraic proof of bc/t5 = 2 at SIGMA_4
# =====================================================================
print(f"\n{'='*70}")
print("Algebraic proof: bc/t5 = 2 at SIGMA_4")
print(f"{'='*70}")
print("""
At k=4, the patterns contributing t5 or bc are:

  (4,):  count = n-4
         sigma_t5 = (n-5)! * 2
         sigma_bc = 0
         Contribution to SIGMA_4: t5 = 2*(n-4)*(n-5)! = 2*(n-4)!

  (2,2): count = C(n-4, 2)
         sigma_t5 = 0
         sigma_bc = (n-6)! * 8
         Contribution to SIGMA_4: bc = 8*C(n-4,2)*(n-6)! = 4*(n-4)!

  Ratio: bc/t5 = 4*(n-4)! / (2*(n-4)!) = 2. QED.

This is NOT a coincidence! It traces to:
  - Each 5-cycle contributes 2 to H(5-set) -> one factor of 2
  - Each disjoint 3-cycle pair contributes 4 to the product -> two factors of 2
  - The pattern (2,2) has C(n-4,2) positions while (4,) has n-4 positions
  - C(n-4,2) / (n-4) = (n-5)/2
  - 8*(n-5)/2 / (2*1) = 2*(n-5) ... hmm, need to include (n-5)! and (n-6)!

More precisely:
  SIGMA_4_t5 = (n-4) * 2 * (n-5)! = 2*(n-4)!
  SIGMA_4_bc = C(n-4,2) * 8 * (n-6)! = (n-4)(n-5)/2 * 8/(n-5) * (n-5)!
             = 4*(n-4) * (n-5)! = 4*(n-4)!
  Ratio = 4/2 = 2.
""")

# =====================================================================
# Verify at SIGMA_5
# =====================================================================
print(f"{'='*70}")
print("SIGMA_5 analysis: bc/t5 ratio")
print(f"{'='*70}")

for n in [7, 9, 11, 13]:
    pats5 = count_position_patterns(n, 5)
    t5_total = bc_total = 0
    details = []
    for pat, cnt in sorted(pats5.items()):
        r = sigma_full(pat, n)
        if r is None:
            details.append(f"    {pat} x {cnt}: UNHANDLED")
            continue
        if r[2] > 0 or r[3] > 0:
            details.append(f"    {pat} x {cnt}: t5={cnt*r[2]}, bc={cnt*r[3]}")
        t5_total += cnt * r[2]
        bc_total += cnt * r[3]

    ratio = bc_total / t5_total if t5_total > 0 else "N/A"
    print(f"\n  n={n}: SIGMA_5 t5={t5_total}, bc={bc_total}, bc/t5={ratio}")
    for d in details:
        print(d)

# =====================================================================
# Algebraic proof for SIGMA_5
# =====================================================================
print(f"\n{'='*70}")
print("Algebraic proof: bc/t5 = 2 at SIGMA_5")
print(f"{'='*70}")

print("""
Patterns contributing t5 or bc at k=5:

  (4,1): one 5-set + one isolated pair
    count(4,1) = (n-4) choices for block start, minus boundary overlaps...
    sigma_t5((4,1)) = (n-7)! * C(n-5,2) * 2  [H(5-set) has +2*t5, pair universal]
    sigma_bc((4,1)) = 0

  (5,): one 6-set
    count(5,) = n-5
    sigma_t5((5,)) = (n-6)! * 2*(n-5)  -> only contributes when n>=6
    sigma_bc((5,)) = (n-6)! * 4

  (2,2,1): two 3-sets + one isolated pair
    sigma_t5 = 0
    sigma_bc((2,2,1)) = (n-8)! * C(n-6,2) * 8

  (3,2): one 4-set + one 3-set
    sigma_t5 = 0
    sigma_bc((3,2)) = (n-7)! * 8*(n-6)
""")

# Compute algebraically
for n in [7, 9, 11, 13]:
    pats5 = count_position_patterns(n, 5)

    # (4,1)
    cnt_41 = pats5.get((4,1), 0)
    t5_41 = cnt_41 * factorial(max(n-7,0)) * comb(n-5, 2) * 2
    bc_41 = 0

    # (5,)
    cnt_5 = pats5.get((5,), 0)
    t5_5 = cnt_5 * factorial(max(n-6,0)) * 2 * (n-5) if n >= 6 else 0
    bc_5 = cnt_5 * factorial(max(n-6,0)) * 4

    # (2,2,1)
    cnt_221 = pats5.get((2,2,1), 0)
    t5_221 = 0
    bc_221 = cnt_221 * factorial(max(n-8,0)) * comb(n-6, 2) * 8 if n >= 8 else 0

    # (3,2)
    cnt_32 = pats5.get((3,2), 0)
    t5_32 = 0
    bc_32 = cnt_32 * factorial(max(n-7,0)) * 8 * (n-6) if n >= 7 else 0

    t5_tot = t5_41 + t5_5
    bc_tot = bc_41 + bc_5 + bc_221 + bc_32
    ratio = bc_tot / t5_tot if t5_tot > 0 else "N/A"

    print(f"  n={n}: t5={t5_tot}, bc={bc_tot}, bc/t5={ratio}")
    print(f"    (4,1) x {cnt_41}: t5={t5_41}")
    print(f"    (5,)  x {cnt_5}: t5={t5_5}, bc={bc_5}")
    print(f"    (2,2,1) x {cnt_221}: bc={bc_221}")
    print(f"    (3,2) x {cnt_32}: bc={bc_32}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
