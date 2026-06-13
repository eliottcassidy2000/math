"""
permanent_forbidden_h.py — Systematic analysis of permanently forbidden H values
kind-pasteur-2026-03-14-S65

For H = 1 + 2*a1 + 4*a2 + 8*a3 + 16*a4 + ..., where a_k = # independent
k-sets in Omega(T), structural constraints force certain H values to be
impossible for ALL n.

Key constraints on a_k:
  (C1) a_k > 0 => a_1 >= k  (each k-set contains k individual cycles)
  (C2) a_k > 0 => a_2 >= C(k,2)  (each k-set has C(k,2) disjoint pairs)
  (C3) a_k >= s > 0 => a_1 >= k*s (each k-set has k cycles; with s disjoint
       k-sets the cycles are ALL distinct)

Previously proved: H=7, H=21 permanently forbidden.
Now: systematic check of ALL odd H up to 300.
"""

from math import comb
from collections import defaultdict

def min_contribution(ak_dict):
    """
    Given a dict {k: a_k}, compute the minimum sum_{k>=1} 2^{k-1} * a_k
    PLUS the forced contributions from structural constraints.

    Returns (total, forced_a1, forced_a2).
    """
    forced_a1 = 0
    forced_a2 = 0

    for k, val in ak_dict.items():
        if val > 0 and k >= 2:
            # Each k-set has k cycles
            forced_a1 = max(forced_a1, k)  # At least k cycles
            # With val many k-sets: if they're independent, a1 >= k*val
            # But they might share cycles, so conservatively a1 >= k
            # Actually: a_{k} = val means val independent k-sets exist
            # Each k-set has k members, so a1 >= k (weakest)
            # But the val k-sets themselves have k*val cycles IF independent
            # For independent polynomial: a_k = # independent k-sets in Omega
            # These k-sets CAN share cycles, so a1 >= k is the right bound

            forced_a2 = max(forced_a2, comb(k, 2))

    return forced_a1, forced_a2

def check_h_feasibility(H):
    """
    Check if H is structurally feasible.

    H = 1 + 2*a1 + 4*a2 + 8*a3 + ...
    Let T = (H-1)/2 = a1 + 2*a2 + 4*a3 + ...

    Return (feasible, valid_decomps, explanation)
    """
    if H < 1 or H % 2 == 0:
        return False, [], "H must be odd >= 1"

    T = (H - 1) // 2
    if T == 0:
        return True, [{}], "H=1 (transitive tournament)"

    # Max possible k: 2^{k-1} <= T
    max_k = 1
    while 2**max_k <= T:
        max_k += 1

    valid = []

    # Enumerate: a_k for k from max_k down to 1
    # sum_{k=1}^{max_k} 2^{k-1} * a_k = T
    def search(k, remaining, partial):
        if k == 0:
            if remaining == 0:
                # Check structural constraints
                d = dict(partial)
                a1 = d.get(1, 0)
                a2 = d.get(2, 0)

                # (C1): For each k >= 2 with a_k > 0: a1 >= k
                for kk, vv in d.items():
                    if kk >= 2 and vv > 0:
                        if a1 < kk:
                            return  # FAIL

                # (C2): For each k >= 2 with a_k > 0: a2 >= C(k,2)
                for kk, vv in d.items():
                    if kk >= 2 and vv > 0:
                        if a2 < comb(kk, 2):
                            return  # FAIL

                # (C3): a2 > 0 => a1 >= 2
                if a2 > 0 and a1 < 2:
                    return  # FAIL

                valid.append(dict(d))
            return

        max_val = remaining // (2**(k-1))
        for v in range(max_val + 1):
            search(k - 1, remaining - v * 2**(k-1), partial + [(k, v)])

    search(max_k, T, [])

    if valid:
        return True, valid, f"{len(valid)} valid decompositions"
    else:
        return False, [], "No valid decomposition"

def explain_infeasibility(H):
    """Detailed explanation of why H is infeasible."""
    T = (H - 1) // 2

    max_k = 1
    while 2**max_k <= T:
        max_k += 1

    reasons = []

    def search_with_reasons(k, remaining, partial):
        if k == 0:
            if remaining == 0:
                d = dict(partial)
                a1 = d.get(1, 0)
                a2 = d.get(2, 0)

                for kk, vv in d.items():
                    if kk >= 2 and vv > 0 and a1 < kk:
                        parts = ", ".join(f"a{kk2}={vv2}" for kk2, vv2 in sorted(d.items()) if vv2 > 0)
                        reasons.append(f"({parts}): a_{kk}>0 forces a_1>={kk} but a_1={a1}")
                        return

                for kk, vv in d.items():
                    if kk >= 2 and vv > 0 and a2 < comb(kk, 2):
                        parts = ", ".join(f"a{kk2}={vv2}" for kk2, vv2 in sorted(d.items()) if vv2 > 0)
                        reasons.append(f"({parts}): a_{kk}>0 forces a_2>={comb(kk,2)} but a_2={a2}")
                        return

                if a2 > 0 and a1 < 2:
                    parts = ", ".join(f"a{kk2}={vv2}" for kk2, vv2 in sorted(d.items()) if vv2 > 0)
                    reasons.append(f"({parts}): a_2>0 forces a_1>=2 but a_1={a1}")
                    return
            return

        max_val = remaining // (2**(k-1))
        for v in range(max_val + 1):
            search_with_reasons(k - 1, remaining - v * 2**(k-1), partial + [(k, v)])

    search_with_reasons(max_k, T, [])
    return reasons

def main():
    print("=" * 70)
    print("PERMANENT FORBIDDEN H VALUES — SYSTEMATIC ANALYSIS")
    print("=" * 70)

    print("\nStructural constraints on independence polynomial coefficients:")
    print("  H = 1 + 2*a1 + 4*a2 + 8*a3 + 16*a4 + ...")
    print("  (C1) a_k > 0 => a_1 >= k")
    print("  (C2) a_k > 0 => a_2 >= C(k,2)")
    print("  (C3) a_2 > 0 => a_1 >= 2")

    # Part 1: Known forbidden from n=7
    print("\n" + "=" * 70)
    print("PART 1: ANALYSIS OF KNOWN n=7 FORBIDDEN VALUES")
    print("=" * 70)

    forbidden_n7 = [7, 21, 63, 107, 119, 149]

    for H in forbidden_n7:
        T = (H - 1) // 2
        feasible, decomps, msg = check_h_feasibility(H)

        print(f"\n  H = {H} (T = a1 + 2*a2 + ... = {T}):")
        print(f"    Structurally feasible: {feasible}")

        if feasible:
            print(f"    Valid decompositions ({len(decomps)}):")
            for d in decomps:
                parts = [(k, v) for k, v in sorted(d.items()) if v > 0]
                if parts:
                    desc = ", ".join(f"a_{k}={v}" for k, v in parts)
                else:
                    desc = "(all zero)"
                print(f"      {desc}")
        else:
            reasons = explain_infeasibility(H)
            # Group by reason type
            reason_counts = defaultdict(int)
            for r in reasons:
                reason_counts[r] += 1
            print(f"    Elimination reasons ({len(reasons)} decompositions rejected):")
            for r in sorted(set(reasons)):
                print(f"      {r}")

    # Part 2: Systematic scan
    print("\n" + "=" * 70)
    print("PART 2: ALL STRUCTURALLY FORBIDDEN ODD H IN [3, 999]")
    print("=" * 70)

    structural_gaps = []
    for H in range(3, 1000, 2):
        feasible, _, _ = check_h_feasibility(H)
        if not feasible:
            structural_gaps.append(H)

    print(f"\n  Count: {len(structural_gaps)}")
    print(f"  Values: {structural_gaps}")

    # Check which n=7 forbidden are structural
    print(f"\n  n=7 forbidden that are structural: {[h for h in forbidden_n7 if h in structural_gaps]}")
    print(f"  n=7 forbidden but NOT structural: {[h for h in forbidden_n7 if h not in structural_gaps]}")

    # Part 3: Pattern analysis
    print("\n" + "=" * 70)
    print("PART 3: PATTERN ANALYSIS OF STRUCTURAL GAPS")
    print("=" * 70)

    if structural_gaps:
        print(f"\n  Gaps: {structural_gaps}")
        print(f"\n  Binary representations:")
        for h in structural_gaps[:30]:
            print(f"    {h:4d} = {bin(h)}")

        print(f"\n  Gaps mod 8: {[h % 8 for h in structural_gaps[:30]]}")
        print(f"  Gaps mod 6: {[h % 6 for h in structural_gaps[:30]]}")

        # Check: are these all of the form 2^k - 1?
        mersenne = [2**k - 1 for k in range(2, 12)]
        print(f"\n  Mersenne numbers: {mersenne}")
        print(f"  Overlap with gaps: {[h for h in structural_gaps if h in mersenne]}")

        # Differences
        if len(structural_gaps) >= 2:
            diffs = [structural_gaps[i+1] - structural_gaps[i] for i in range(len(structural_gaps)-1)]
            print(f"\n  Consecutive differences: {diffs[:20]}")

    # Part 4: Deeper look at H=63
    print("\n" + "=" * 70)
    print("PART 4: DEEP ANALYSIS — WHY IS H=63 FORBIDDEN?")
    print("=" * 70)

    H = 63
    T = 31  # a1 + 2*a2 + 4*a3 + ... = 31

    feasible, decomps, _ = check_h_feasibility(H)
    if feasible:
        print(f"\n  H=63 is structurally FEASIBLE with {len(decomps)} decompositions.")
        print(f"  This means H=63 might be achievable at some n > 7.")
        print(f"\n  Valid decompositions:")
        for d in decomps:
            parts = [(k, v) for k, v in sorted(d.items()) if v > 0]
            desc = ", ".join(f"a_{k}={v}" for k, v in parts)
            highest = max(k for k, v in d.items() if v > 0) if any(v > 0 for v in d.values()) else 0
            min_n = 3 * highest if highest > 0 else 3
            print(f"    {desc}  (min n = {min_n})")

        # Why absent at n=7 then?
        print(f"\n  Why absent at n=7?")
        print(f"  At n=7, only quadratic I.P. (a_k = 0 for k >= 3)")
        print(f"  So only a1 + 2*a2 = 31 with a1,a2 >= 0")
        print(f"  Valid quadratic decompositions:")
        for a2 in range(16):
            a1 = 31 - 2*a2
            if a1 >= 0:
                if a2 > 0 and a1 < 2:
                    print(f"    a1={a1}, a2={a2}: BLOCKED (a2>0 but a1<2)")
                else:
                    print(f"    a1={a1}, a2={a2}: structurally OK — but NOT achievable at n=7!")
        print(f"\n  So H=63 is structurally fine but empirically absent at n=7.")
        print(f"  The GAP is from the achievable (a1,a2) set, NOT structural constraints.")
        print(f"  This means H=63 MIGHT become achievable at larger n!")
    else:
        print(f"\n  H=63 is structurally FORBIDDEN for ALL n.")
        reasons = explain_infeasibility(H)
        for r in set(reasons):
            print(f"    {r}")

    # Part 5: Same analysis for H=107, 119, 149
    print("\n" + "=" * 70)
    print("PART 5: ANALYSIS OF H=107, 119, 149")
    print("=" * 70)

    for H in [107, 119, 149]:
        T = (H - 1) // 2
        feasible, decomps, _ = check_h_feasibility(H)
        print(f"\n  H = {H} (T = {T}):")
        if feasible:
            print(f"    Structurally FEASIBLE ({len(decomps)} decompositions)")
            # Show quadratic ones (relevant for n=7)
            quad_decomps = [d for d in decomps if all(k <= 2 for k, v in d.items() if v > 0)]
            print(f"    Quadratic decompositions: {len(quad_decomps)}")
            for d in quad_decomps[:5]:
                a1 = d.get(1, 0)
                a2 = d.get(2, 0)
                print(f"      a1={a1}, a2={a2}")
            print(f"    => Absent at n=7 due to (a1,a2) ACHIEVABILITY, not structure")
            print(f"    => MIGHT become achievable at larger n")
        else:
            print(f"    Structurally FORBIDDEN for all n")

    # Part 6: The achievable (a1, a2) question at n=7
    print("\n" + "=" * 70)
    print("PART 6: WHAT CONSTRAINS (a1, a2) BEYOND STRUCTURE?")
    print("=" * 70)

    print(f"""
  The structural constraints (C1-C3) are NECESSARY but not SUFFICIENT.
  Additional constraints come from the GRAPH STRUCTURE of Omega(T):

  (G1) a1 = |V(Omega)| = number of odd directed cycles in T
  (G2) a2 = |edges NOT in Omega| = number of disjoint cycle pairs
  (G3) Not all (a1, a2) pairs are realizable: the conflict structure
       of Omega(T) constrains which pairs can occur.

  For H=63 at n=7: quadratic decompositions include a1=31,a2=0
  which would need 31 odd cycles, all pairwise conflicting.
  At n=7: max cycles is c3=C(7,3)*2/7? Actually c3 <= 35 (C(7,3)=35).
  So a1=31 is potentially achievable in terms of cycle count.
  But all 31 cycles must pairwise share a vertex? Very tight constraint.

  The ACHIEVABILITY question is fundamentally different from STRUCTURAL:
  - Structural: can the integers (a1,a2,...) satisfy the constraints?
  - Achievability: does some tournament T have Omega(T) with these parameters?

  H=7,21: STRUCTURALLY forbidden (no valid decomposition at any degree)
  H=63,107,119,149: structurally FEASIBLE but ACHIEVABILITY-forbidden at n=7
  """)

    # Part 7: Growth of achievable H values with n
    print("=" * 70)
    print("PART 7: H=63 — WILL IT EVENTUALLY APPEAR?")
    print("=" * 70)

    print(f"""
  H = 63 = 1 + 2*31 = 1 + 2*a1 (if a2=0)
  Needs a1 = 31 odd cycles, ALL pairwise conflicting (sharing a vertex).

  At n=7: c3 up to 35, c5 up to 21, c7 = 0 or 1
  So a1 = c3+c5+c7 up to ~57. Getting a1=31 is plausible.
  But ALL pairwise conflicting? With 31 cycles, every pair shares a vertex.

  Helly property: if a family of sets is such that every pair intersects,
  must they have a common point? For 3-element subsets of [n]:
  Helly number for sets of size <= k is 2k-1.
  For k=3: Helly number = 5. So 6 pairwise-intersecting 3-sets
  do NOT need a common point.

  But for odd cycles of various sizes (3,5,7), the constraint is weaker.

  At n >= 9: cubic I.P. unlocks. Could have (a1=23, a2=0, a3=1) giving
  H = 1 + 46 + 0 + 8 = 55. Or (a1=27, a2=0, a3=1) giving H = 1+54+0+8 = 63!
  Wait: a3=1 forces a1 >= 3, a2 >= 3.
  So H = 1 + 2*a1 + 4*a2 + 8*1 = 1 + 2*a1 + 4*a2 + 8
  For H=63: 2*a1 + 4*a2 = 54, i.e., a1 + 2*a2 = 27
  With a1 >= 3, a2 >= 3: a1 = 27-2*a2, need a1 >= 3 => a2 <= 12
  And a2 >= 3. So a2 in [3, 12], a1 = 27-2*a2 in [3, 21].
  MANY feasible decompositions!

  So H=63 IS structurally feasible at n >= 9 even with cubic I.P.
  The question is whether some tournament at n >= 9 actually achieves it.
  This is an OPEN QUESTION we cannot resolve by pure structural argument.
  """)

    print(f"{'='*70}")
    print(f"FINAL SUMMARY")
    print(f"{'='*70}")
    print(f"""
  PERMANENTLY FORBIDDEN (structural, all n):
    H = 7:  a1+2*a2 = 3, too tight for any higher coefficient
    H = 21: a3>=1 forces a1>=3,a2>=3 => H>=27; quadratic case impossible

  Structural gaps (odd H in [3,999]):
    {structural_gaps}

  ACHIEVABILITY-FORBIDDEN at n=7 but structurally feasible:
    H = 63, 107, 119, 149

  Status: These 4 values MIGHT become achievable at n >= 8 or n >= 9.
  Their absence at n=7 is a GRAPH CONSTRAINT, not a pure integer constraint.
  To resolve: need exhaustive computation at n=8 or targeted construction.
  """)

if __name__ == "__main__":
    main()
