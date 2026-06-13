#!/usr/bin/env python3
"""
CLOSED-FORM moment formulas at general n via the General Sigma Formula.

Key insight: sigma((c1,...,cm)) = (n-sum_ci-m)! * sum prod H(T[Gj])
over ordered partitions into groups of sizes (c1+1,...,cm+1).

Since H(T[2-set]) = 1 always, patterns with only size-1 components are universal.
Patterns with a single size-2 component bring in t3.
Patterns with size-3 or two size-2 components bring in t5 and bc.

opus-2026-03-06-S28
"""
from itertools import permutations, combinations
from math import factorial, comb
from collections import Counter
import random

# =====================================================================
# POSITION PATTERN COUNTING
# =====================================================================

def count_position_patterns(n, k):
    """Count patterns of size-k subsets of positions {0,...,n-2}."""
    if k == 0:
        return {(): 1}
    patterns = Counter()
    for S in combinations(range(n-1), k):
        pos = sorted(S)
        comps = []
        comp = [pos[0]]
        for i in range(1, len(pos)):
            if pos[i] == comp[-1] + 1:
                comp.append(pos[i])
            else:
                comps.append(len(comp))
                comp = [pos[i]]
        comps.append(len(comp))
        pat = tuple(sorted(comps, reverse=True))
        patterns[pat] += 1
    return dict(patterns)

# =====================================================================
# SIGMA FORMULAS from General Sigma Formula
# =====================================================================

def sigma_const_and_t3(pattern, n):
    """
    Return (const, t3_coeff) for sigma(pattern) at general n.
    Only works for patterns where the only non-trivial invariant is t3.
    That is: all components are size 1 or exactly one component is size 2.
    """
    sizes = list(pattern)  # component sizes
    total_verts = sum(s + 1 for s in sizes)
    free = n - total_verts

    if free < 0:
        return (0, 0)

    # Count how many size-1 components and size-2+ components
    num_size1 = sizes.count(1)
    big_sizes = [s for s in sizes if s > 1]

    if len(big_sizes) == 0:
        # All size-1: sigma = (n-2k)! * n!/((n-2k)! * 2^k) = n!/2^k
        k = len(sizes)
        const = factorial(n) // (2**k)
        return (const, 0)

    if len(big_sizes) == 1 and big_sizes[0] == 2:
        # One size-2 component (3 vertices) + num_size1 pairs
        # sigma = (free)! * C(n-3, 2)^... no, need to be more careful
        # G1 = 3-vertex group (for size-2), remaining groups are 2-vertex
        # H(3-set) = 1 + 2*cyc, H(2-set) = 1
        # sum over ordered (G1 of 3, G2 of 2, ..., Gm of 2, all disjoint from [n])
        # = [sum_{G1 in C(n,3)} H(T[G1])] * [C(n-3,2) * C(n-5,2) * ... ]
        # = (C(n,3) + 2*t3) * prod_{i=0}^{num_size1-1} C(n-3-2i, 2)

        pair_product = 1
        remaining = n - 3
        for i in range(num_size1):
            pair_product *= comb(remaining, 2)
            remaining -= 2

        const = factorial(free) * comb(n, 3) * pair_product
        t3_coeff = factorial(free) * 2 * pair_product
        return (const, t3_coeff)

    if len(big_sizes) == 1 and big_sizes[0] == 3:
        # One size-3 component (4 vertices): H(4-set) = 1 + 2*c3(4-set)
        # sum_{4-subsets} H = C(n,4) + 2*(n-3)*t3
        # (each 3-cycle is in (n-3) 4-subsets)
        pair_product = 1
        remaining = n - 4
        for i in range(num_size1):
            pair_product *= comb(remaining, 2)
            remaining -= 2

        const = factorial(free) * comb(n, 4) * pair_product
        t3_coeff = factorial(free) * 2 * (n - 3) * pair_product
        return (const, t3_coeff)

    # For other patterns, return None to indicate we need more invariants
    return None

# =====================================================================
# m2 CLOSED FORM
# =====================================================================
print("=" * 70)
print("m2 = S(2,1)*1!*SIGMA_1 + S(2,2)*2!*SIGMA_2")
print("   = SIGMA_1 + 2*SIGMA_2")
print("=" * 70)

for n in [3, 5, 7, 9, 11, 13, 15]:
    # SIGMA_1: all size-1 position subsets of size 1
    # Only pattern: (1,), count = n-1
    SIGMA_1 = (n - 1) * factorial(n) // 2

    # SIGMA_2: patterns of size 2
    pats2 = count_position_patterns(n, 2)
    SIGMA_2_const = 0
    SIGMA_2_t3 = 0
    for pat, cnt in pats2.items():
        c, t = sigma_const_and_t3(pat, n)
        SIGMA_2_const += cnt * c
        SIGMA_2_t3 += cnt * t

    m2_const = SIGMA_1 + 2 * SIGMA_2_const
    m2_t3 = 2 * SIGMA_2_t3

    # Closed form prediction
    pred_const = factorial(n) * (3*n*n - 5*n + 4) // 12
    pred_t3 = 4 * factorial(n - 2)

    print(f"  n={n:2d}: m2 = {m2_const} + {m2_t3}*t3")
    print(f"         closed: {pred_const} + {pred_t3}*t3  "
          f"{'OK' if m2_const == pred_const and m2_t3 == pred_t3 else 'FAIL'}")

# =====================================================================
# m3 CLOSED FORM
# =====================================================================
print(f"\n{'='*70}")
print("m3 = SIGMA_1 + 6*SIGMA_2 + 6*SIGMA_3")
print(f"{'='*70}")

for n in [3, 5, 7, 9, 11, 13]:
    SIGMA_1 = (n - 1) * factorial(n) // 2

    pats2 = count_position_patterns(n, 2)
    SIGMA_2_const = 0
    SIGMA_2_t3 = 0
    for pat, cnt in pats2.items():
        c, t = sigma_const_and_t3(pat, n)
        SIGMA_2_const += cnt * c
        SIGMA_2_t3 += cnt * t

    pats3 = count_position_patterns(n, 3)
    SIGMA_3_const = 0
    SIGMA_3_t3 = 0
    for pat, cnt in pats3.items():
        result = sigma_const_and_t3(pat, n)
        if result is None:
            print(f"    WARNING: pattern {pat} at n={n} needs more invariants")
            continue
        c, t = result
        SIGMA_3_const += cnt * c
        SIGMA_3_t3 += cnt * t

    m3_const = SIGMA_1 + 6 * SIGMA_2_const + 6 * SIGMA_3_const
    m3_t3 = 6 * SIGMA_2_t3 + 6 * SIGMA_3_t3

    print(f"  n={n:2d}: m3 = {m3_const} + {m3_t3}*t3")

    # Try to find closed form: m3_const/n! and m3_t3/(n-2)!
    ratio_const = m3_const / factorial(n)
    ratio_t3 = m3_t3 / factorial(n - 2) if n >= 2 else 0
    print(f"         const/n! = {ratio_const:.6f}, t3_coeff/(n-2)! = {ratio_t3:.6f}")

# =====================================================================
# m3 closed form derivation
# =====================================================================
print(f"\n{'='*70}")
print("m3 closed form: fitting polynomials")
print(f"{'='*70}")

# Collect data
m3_consts = []
m3_t3s = []
ns = [3, 5, 7, 9, 11, 13]

for n in ns:
    SIGMA_1 = (n - 1) * factorial(n) // 2

    pats2 = count_position_patterns(n, 2)
    SIGMA_2_const = SIGMA_2_t3 = 0
    for pat, cnt in pats2.items():
        c, t = sigma_const_and_t3(pat, n)
        SIGMA_2_const += cnt * c
        SIGMA_2_t3 += cnt * t

    pats3 = count_position_patterns(n, 3)
    SIGMA_3_const = SIGMA_3_t3 = 0
    for pat, cnt in pats3.items():
        result = sigma_const_and_t3(pat, n)
        if result is None: continue
        c, t = result
        SIGMA_3_const += cnt * c
        SIGMA_3_t3 += cnt * t

    m3_const = SIGMA_1 + 6 * SIGMA_2_const + 6 * SIGMA_3_const
    m3_t3 = 6 * SIGMA_2_t3 + 6 * SIGMA_3_t3
    m3_consts.append(m3_const)
    m3_t3s.append(m3_t3)

# Fit: m3_t3 / (n-2)!
print("\n  t3 coefficient analysis:")
for i, n in enumerate(ns):
    r = m3_t3s[i] / factorial(n-2)
    print(f"    n={n}: m3_t3/(n-2)! = {r}")

# Fit: 12 * m3_const / n!
print("\n  constant term analysis:")
for i, n in enumerate(ns):
    r = m3_consts[i] / factorial(n)
    # Try polynomial fit
    print(f"    n={n}: m3_const/n! = {r:.10f}")

# Try: m3_const = n! * P(n) / D for some polynomial P(n)
# 12*m3_const/n!:
print("\n  Try 12*m3_const/n!:")
for i, n in enumerate(ns):
    r = 12 * m3_consts[i] / factorial(n)
    print(f"    n={n}: {r:.6f}")

# Try polynomial of degree 4 in n for m3_const/n!
print("\n  Try m3_const * 120 / n!:")
for i, n in enumerate(ns):
    r = 120 * m3_consts[i] / factorial(n)
    print(f"    n={n}: {r:.6f}")

# =====================================================================
# m4: when do t5 and bc enter?
# =====================================================================
print(f"\n{'='*70}")
print("m4 analysis: m4 = S(4,1)*SIGMA_1 + S(4,2)*2*SIGMA_2 + S(4,3)*6*SIGMA_3 + S(4,4)*24*SIGMA_4")
print("S(4,1)=1, S(4,2)=7, S(4,3)=6, S(4,4)=1")
print("m4 = SIGMA_1 + 14*SIGMA_2 + 36*SIGMA_3 + 24*SIGMA_4")
print(f"{'='*70}")

for n in [5, 7, 9, 11]:
    SIGMA_1 = (n - 1) * factorial(n) // 2

    SIGMA_2_const = SIGMA_2_t3 = 0
    for pat, cnt in count_position_patterns(n, 2).items():
        c, t = sigma_const_and_t3(pat, n)
        SIGMA_2_const += cnt * c
        SIGMA_2_t3 += cnt * t

    SIGMA_3_const = SIGMA_3_t3 = 0
    for pat, cnt in count_position_patterns(n, 3).items():
        result = sigma_const_and_t3(pat, n)
        if result is None: continue
        c, t = result
        SIGMA_3_const += cnt * c
        SIGMA_3_t3 += cnt * t

    # SIGMA_4 patterns
    pats4 = count_position_patterns(n, 4)
    SIGMA_4_const = SIGMA_4_t3 = 0
    needs_more = []
    for pat, cnt in sorted(pats4.items()):
        result = sigma_const_and_t3(pat, n)
        if result is None:
            needs_more.append((pat, cnt))
        else:
            c, t = result
            SIGMA_4_const += cnt * c
            SIGMA_4_t3 += cnt * t

    m4_const_partial = SIGMA_1 + 14*SIGMA_2_const + 36*SIGMA_3_const + 24*SIGMA_4_const
    m4_t3_partial = 14*SIGMA_2_t3 + 36*SIGMA_3_t3 + 24*SIGMA_4_t3

    print(f"\n  n={n}: m4 (t3-only part) = {m4_const_partial} + {m4_t3_partial}*t3")
    if needs_more:
        print(f"    Patterns needing more invariants: {needs_more}")
    else:
        print(f"    ALL patterns handled by t3 alone!")

# =====================================================================
# Extend sigma formulas for size-4 component and two size-2 components
# =====================================================================
print(f"\n{'='*70}")
print("Extended sigma formulas for patterns requiring t5, bc")
print(f"{'='*70}")

def sigma_extended(pattern, n):
    """
    Return (const, t3_coeff, t5_coeff, bc_coeff) for sigma(pattern).
    """
    sizes = list(pattern)
    total_verts = sum(s + 1 for s in sizes)
    free = n - total_verts

    if free < 0:
        return (0, 0, 0, 0)

    num_size1 = sizes.count(1)
    big_sizes = sorted([s for s in sizes if s > 1], reverse=True)

    if len(big_sizes) == 0:
        # All isolated: n!/2^k, tournament-independent
        k = len(sizes)
        return (factorial(n) // (2**k), 0, 0, 0)

    if big_sizes == [2]:
        # One size-2 component: depends only on t3
        pair_product = 1
        remaining = n - 3
        for i in range(num_size1):
            pair_product *= comb(remaining, 2)
            remaining -= 2
        const = factorial(free) * comb(n, 3) * pair_product
        t3_coeff = factorial(free) * 2 * pair_product
        return (const, t3_coeff, 0, 0)

    if big_sizes == [3]:
        # One size-3 component: H(4-set) = 1 + 2*c3(4-set)
        # sum_{4-subsets} H = C(n,4) + 2*(n-3)*t3
        pair_product = 1
        remaining = n - 4
        for i in range(num_size1):
            pair_product *= comb(remaining, 2)
            remaining -= 2
        const = factorial(free) * comb(n, 4) * pair_product
        t3_coeff = factorial(free) * 2 * (n - 3) * pair_product
        return (const, t3_coeff, 0, 0)

    if big_sizes == [4]:
        # One size-4 component: H(5-set) = 1 + 2*c3(5-set) + 2*c5(5-set)
        # At n=5: H = 1 + 2*(c3+c5), no bc term (OCF at n=5)
        # sum_{5-subsets S} H(T[S]) = C(n,5) + 2*(n-4 choose 2)*t3 + 2*t5
        # Wait: each 3-cycle in [n] appears in C(n-3,2) 5-subsets
        # each 5-cycle is in exactly 1 5-subset
        pair_product = 1
        remaining = n - 5
        for i in range(num_size1):
            pair_product *= comb(remaining, 2)
            remaining -= 2
        const = factorial(free) * comb(n, 5) * pair_product
        t3_coeff = factorial(free) * 2 * comb(n - 3, 2) * pair_product
        t5_coeff = factorial(free) * 2 * pair_product
        return (const, t3_coeff, t5_coeff, 0)

    if big_sizes == [2, 2]:
        # Two size-2 components: H(3-set)*H(3-set) on disjoint triples
        # sum_{ordered (G1,G2) disjoint 3-sets} dp2(G1)*dp2(G2)
        # = sum_{6-subsets S} [sum of ordered pair partitions into 3+3]
        # Each 6-subset contributes: C(6,3) ordered partitions? No...
        # C(6,3)/2 unordered, but since G1 and G2 are assigned to different
        # position components, they ARE ordered.
        # Actually: # ordered (G1,G2) pairs of disjoint 3-subsets from [n]
        # = C(n,3)*C(n-3,3)
        # For each such pair: H(G1)*H(G2) = (1+2*cyc1)(1+2*cyc2)
        # = 1 + 2*cyc1 + 2*cyc2 + 4*cyc1*cyc2

        # sum over ordered disjoint pairs of 3-sets:
        # sum 1 = C(n,3)*C(n-3,3)
        # sum cyc1 = t3 * C(n-3,3) [for each cyclic G1, choose any G2]
        #          + t3 * C(n-3,3)... wait, need to be more careful

        # sum_{(G1,G2)} cyc(G1) = sum_{G1} cyc(G1) * C(n-3, 3)
        #   = t3 * C(n-3, 3)
        # By symmetry, sum cyc(G2) = t3 * C(n-3, 3)
        # sum_{(G1,G2)} cyc(G1)*cyc(G2)
        #   = sum_{disjoint (G1,G2)} cyc(G1)*cyc(G2)
        #   = [all pairs] - [intersecting pairs]
        #   = t3^2 - sum_{G1 cyclic} sum_{G2 cyclic, G1 cap G2 != empty} 1

        # Hmm this is getting complicated. Let me use the bc definition:
        # bc = sum over 6-subsets of #{disjoint cyclic 3-cycle pairs}
        # So sum_{unordered disjoint cyclic pairs} = bc
        # sum_{ordered disjoint cyclic pairs} = 2*bc

        # Wait, at n=7 the formula was sigma((2,2)) = 140 + 16*t3 + 8*bc
        # Let me derive it for general n.

        # With num_size1 additional isolated pairs:
        pair_product = 1
        remaining = n - 6
        for i in range(num_size1):
            pair_product *= comb(remaining, 2)
            remaining -= 2

        # For the two 3-groups:
        # sum_{ordered (G1,G2) disjoint 3-subsets} (1+2*cyc1)(1+2*cyc2)
        # = C(n,3)*C(n-3,3) + 2*t3*C(n-3,3) + 2*t3*C(n-3,3) + 4*(2*bc + correction?)
        # Hmm, the 2*bc counts ordered pairs of disjoint cyclic triples.
        # But we're summing over ordered (G1,G2) where G1 is assigned to first
        # position component and G2 to second. So:
        # sum cyc1*cyc2 = sum_{ordered disjoint cyclic pairs} 1 = 2*bc? No.
        #
        # bc = #{unordered pairs of vertex-disjoint 3-cycles among ALL 6-vertex subsets}
        # At n=7: bc = sum_v bc(T-v) where bc(T-v) = #{disjoint 3-cycle pairs in T-v}
        # Actually bc in THM-056 is sum over 6-vertex subsets of #{pairs of disjoint 3-cycles}
        #
        # For our purpose: ordered disjoint cyclic triple pairs from [n]
        # = 2 * (sum over 6-subsets of #{unordered disjoint cyclic pairs in that 6-subset})
        #   ... no, a 6-subset determines both triples uniquely (complement)
        # A 6-vertex subset S partitions into two complementary 3-subsets in C(6,3)/2 = 10 ways.
        # For each such partition (G1,G2), cyc(G1)*cyc(G2) is 1 iff both are cyclic.
        #
        # Hmm, this is the wrong framing. The ordered pairs (G1,G2) of disjoint 3-subsets
        # from [n] DON'T restrict to 6-subsets in a simple way when n>6.
        #
        # Let me just note that sum_{ordered (G1,G2) disjoint from [n]} cyc(G1)*cyc(G2)
        # counts ordered pairs of disjoint directed 3-cycles.
        # Unordered count = bc_total (the total over all of [n]).
        # So ordered = 2*bc_total.
        #
        # But bc in THM-056 was defined as sum over 6-subsets... which IS bc_total
        # since each pair of disjoint 3-cycles lives in exactly one 6-subset.
        # So bc_total = bc. And ordered = 2*bc.

        base_const = comb(n, 3) * comb(n-3, 3)
        base_t3 = 2 * comb(n-3, 3) + 2 * comb(n-3, 3)  # = 4*C(n-3,3)
        base_bc = 4 * 2  # 4 * (ordered pairs) = 4*2*bc = 8*bc

        const = factorial(free) * pair_product * base_const
        t3_coeff = factorial(free) * pair_product * (2 * comb(n-3, 3) * 2)
        bc_coeff = factorial(free) * pair_product * 8
        return (const, t3_coeff, 0, bc_coeff)

    if big_sizes == [3, 2]:
        # Size-3 comp (4 verts) + size-2 comp (3 verts)
        # H(4-set) * H(3-set) on disjoint groups
        # sum = sum_{G1 of 4, G2 of 3, disjoint} H(G1)*H(G2)
        # = sum_{G1} sum_{G2 from [n]-G1} (1+2*c3(G1)) * (1+2*cyc(G2))
        pair_product = 1
        remaining = n - 7
        for i in range(num_size1):
            pair_product *= comb(remaining, 2)
            remaining -= 2

        # sum_{G1 of 4} (1+2*c3(G1)) * sum_{G2 of 3 from [n]-G1} (1+2*cyc(G2))
        # The inner sum: C(n-4,3) + 2*[# cyclic triples in [n]-G1]
        # = C(n-4,3) + 2*(t3 - c3(G1)*(correction for overlap))
        # Hmm, # cyclic triples in [n]\G1 = t3 - #{cyclic triples intersecting G1}
        # This doesn't factor nicely.

        # Alternative: expand the product
        # sum_{disjoint (G1 of 4, G2 of 3)} 1 = C(n,4)*C(n-4,3)
        # sum c3(G1) = (n-3)*t3 * ... no, c3(G1) is the number of directed 3-cycles in G1.
        # For a 4-vertex subset, c3(G1) is 0-4 depending on the tournament.
        # sum_{G1 of 4} c3(G1) = (n-3)*t3 (each 3-cycle is in (n-3) 4-subsets)

        # sum_{disjoint} c3(G1) = sum_{G1} c3(G1) * C(n-4,3) = (n-3)*t3 * C(n-4,3)
        # ... no, this counts with G2 as well. The correct:
        # sum_{disjoint (G1,G2)} c3(G1) = sum_{G1 of 4} c3(G1) * C(n-4, 3)
        # = (n-3)*t3 * C(n-4,3)

        # sum_{disjoint} cyc(G2) = sum_{G2 of 3} cyc(G2) * C(n-3, 4)
        # ... no. G1 is chosen from [n]\G2, a set of n-3 elements.
        # sum_{disjoint} cyc(G2) = sum_{G2 of 3} cyc(G2) * C(n-3, 4)
        # = t3 * C(n-3, 4)

        # sum_{disjoint} c3(G1)*cyc(G2): this is harder. G1 and G2 are disjoint.
        # = sum_{G2 cyclic} sum_{G1 of 4 from [n]\G2} c3(G1)
        # = sum_{G2 cyclic} (n-3-3)*[t3 restricted to [n]\G2]... ugh
        # = sum_{G2 cyclic} [sum over 3-cycles disjoint from G2] * (n-6)... no
        # Actually: sum_{G1 of 4 from S} c3(G1) = |S|-3 choose 1 * (# 3-cycles in S)
        # where S = [n]\G2 has n-3 elements.
        # 3-cycles in [n]\G2: t3 - #{3-cycles intersecting G2}
        # #{3-cycles sharing >=1 vertex with G2}: not easy to express simply.

        # This is getting complex. For now, just verify at n=7.
        # At n=7: sigma((3,2)) = 35 + 10*t3 + 8*bc (from THM-056)
        # I'll compute numerically and try to find the pattern.
        return None  # TODO

    if big_sizes == [5]:
        # H(6-set): uses OCF at n=6
        # H(T[S]) = 1 + 2*c3(S) + 2*c5(S) + 4*bc(S)
        # sum_{6-subsets S} H = C(n,6) + 2*C(n-3,3)*t3 + 2*C(n-5,1)*t5 + 4*bc
        pair_product = 1
        remaining = n - 6
        for i in range(num_size1):
            pair_product *= comb(remaining, 2)
            remaining -= 2
        const = factorial(free) * comb(n, 6) * pair_product
        t3_coeff = factorial(free) * 2 * comb(n-3, 3) * pair_product
        t5_coeff = factorial(free) * 2 * (n - 5) * pair_product
        bc_coeff = factorial(free) * 4 * pair_product
        return (const, t3_coeff, t5_coeff, bc_coeff)

    return None  # Unknown pattern

# =====================================================================
# Verify extended formulas at n=7
# =====================================================================
print(f"\n{'='*70}")
print("Verify extended sigma formulas at n=7")
print(f"{'='*70}")

n = 7
known_n7 = {
    (): (5040, 0, 0, 0),
    (1,): (2520, 0, 0, 0),
    (1,1): (1260, 0, 0, 0),
    (1,1,1): (630, 0, 0, 0),
    (2,): (840, 48, 0, 0),
    (2,1): (420, 24, 0, 0),
    (2,1,1): (210, 12, 0, 0),
    (3,): (210, 48, 0, 0),
    (3,1): (105, 24, 0, 0),
    (4,): (42, 24, 4, 0),
    (4,1): (21, 12, 2, 0),
    (2,2): (140, 16, 0, 8),
    (5,): (7, 8, 4, 4),
}

for pat, known in known_n7.items():
    result = sigma_extended(pat, 7)
    if result is None:
        print(f"  {pat}: not implemented")
        continue
    match = result == known
    if not match:
        print(f"  {pat}: FAIL — got {result}, expected {known}")
    else:
        print(f"  {pat}: OK")

# =====================================================================
# Verify m2, m3 formulas with brute force at small n
# =====================================================================
print(f"\n{'='*70}")
print("Brute-force verification at n=5")
print(f"{'='*70}")

for n in [5]:
    for trial in range(5):
        random.seed(n*1000 + trial)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        # Count t3
        t3 = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if A[i][j]*A[j][k]*A[k][i]: t3 += 1
                    if A[i][k]*A[k][j]*A[j][i]: t3 += 1

        # Brute force moments
        moments = [0] * 5
        for p in permutations(range(n)):
            f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]])
            for j in range(5):
                moments[j] += f**j

        # Predicted m2
        m2_pred = factorial(n) * (3*n*n - 5*n + 4) // 12 + 4 * factorial(n-2) * t3
        print(f"  Trial {trial}: t3={t3}, m2_actual={moments[2]}, m2_pred={m2_pred}, "
              f"{'OK' if moments[2] == m2_pred else 'FAIL'}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
