"""
beta2_n5_gap_check.py — opus-2026-04-04-S1
Check whether Case 2 of THM-109 (free cycles exist) is vacuous at n=5.
If so, the entire beta_2=0 proof is algebraic with NO exhaustive verification needed.
"""
import numpy as np
from beta2_n5_case2_analysis import (tournament_from_bits, find_3cycles,
                                      is_dominated, is_strongly_connected,
                                      vertex_connectivity, compute_b1,
                                      share_directed_edge)

def full_analysis():
    n = 5
    num_arcs = n*(n-1)//2  # 10

    counts = {
        'total': 0,
        'b1_0': 0,
        'b1_0_SC': 0,
        'b1_0_SC_kappa1': 0,
        'b1_0_SC_kappa2+': 0,
        'b1_0_notSC': 0,
        'b1_0_has_free': 0,
        'b1_0_SC_has_free': 0,
        'b1_0_SC_kappa2_has_free': 0,
        'b1_0_all_dominated': 0,
        'b1_0_no_cycles': 0,
        'b1_1': 0,
    }

    for bits in range(1 << num_arcs):
        T = tournament_from_bits(n, bits)
        counts['total'] += 1

        b1 = compute_b1(T, n)
        if b1 > 0:
            counts['b1_1'] += 1
            continue

        counts['b1_0'] += 1

        cycles = find_3cycles(T, n)
        if not cycles:
            counts['b1_0_no_cycles'] += 1
            continue

        has_free = any(len(is_dominated(T, n, c)) == 0 for c in cycles)
        all_dom = all(len(is_dominated(T, n, c)) > 0 for c in cycles)

        if has_free:
            counts['b1_0_has_free'] += 1
        if all_dom:
            counts['b1_0_all_dominated'] += 1

        sc = is_strongly_connected(T, n)
        if sc:
            counts['b1_0_SC'] += 1
            kappa = vertex_connectivity(T, n)
            if kappa == 1:
                counts['b1_0_SC_kappa1'] += 1
            else:
                counts['b1_0_SC_kappa2+'] += 1

            if has_free:
                counts['b1_0_SC_has_free'] += 1
                if kappa >= 2:
                    counts['b1_0_SC_kappa2_has_free'] += 1
        else:
            counts['b1_0_notSC'] += 1

    print("=== COMPLETE n=5 TOURNAMENT CLASSIFICATION ===")
    print(f"Total tournaments: {counts['total']}")
    print(f"  b1=1: {counts['b1_1']}")
    print(f"  b1=0: {counts['b1_0']}")
    print(f"    No 3-cycles (transitive): {counts['b1_0_no_cycles']}")
    print(f"    Has free cycles: {counts['b1_0_has_free']}")
    print(f"    All cycles dominated: {counts['b1_0_all_dominated']}")
    print(f"    Strongly connected: {counts['b1_0_SC']}")
    print(f"      kappa=1: {counts['b1_0_SC_kappa1']}")
    print(f"      kappa>=2: {counts['b1_0_SC_kappa2+']}")
    print(f"      SC + has free: {counts['b1_0_SC_has_free']}")
    print(f"      SC + kappa>=2 + has free: {counts['b1_0_SC_kappa2_has_free']}")
    print(f"    Not SC: {counts['b1_0_notSC']}")

    # THE KEY QUESTION: Is Case 2 (has free cycles) vacuous at n=5 for b1=0?
    print(f"\n=== THE KEY FINDING ===")
    if counts['b1_0_has_free'] == 0:
        print("*** Case 2 of THM-109 (free cycles exist, b1=0) is VACUOUS at n=5! ***")
        print("*** Therefore the ENTIRE beta_2=0 proof is algebraic! ***")
        print("*** No exhaustive verification needed at n=5! ***")
    else:
        print(f"Case 2 applies to {counts['b1_0_has_free']} tournaments at n=5")
        print("Need to check if these all have good vertices algebraically")

    # Also check: does EVERY b1=0 tournament with cycles have all dominated?
    print(f"\n=== VERIFICATION ===")
    print(f"b1=0 with cycles: {counts['b1_0'] - counts['b1_0_no_cycles']}")
    print(f"b1=0 all dominated: {counts['b1_0_all_dominated']}")
    if counts['b1_0'] - counts['b1_0_no_cycles'] == counts['b1_0_all_dominated']:
        print("CONFIRMED: b1=0 WITH cycles implies ALL cycles dominated at n=5!")
        print("This means: b1(T)=0 iff every 3-cycle is dominated (for n=5).")

    # Additional: verify this at n=4 too
    print(f"\n=== n=4 CHECK ===")
    n4 = 4
    n4_arcs = n4*(n4-1)//2  # 6
    n4_b1_0 = 0
    n4_b1_0_has_free = 0
    n4_b1_0_all_dom = 0
    for bits in range(1 << n4_arcs):
        T = tournament_from_bits(n4, bits)
        b1 = compute_b1(T, n4)
        if b1 == 0:
            n4_b1_0 += 1
            cycles = find_3cycles(T, n4)
            if cycles:
                has_free = any(len(is_dominated(T, n4, c)) == 0 for c in cycles)
                all_dom = all(len(is_dominated(T, n4, c)) > 0 for c in cycles)
                if has_free:
                    n4_b1_0_has_free += 1
                if all_dom:
                    n4_b1_0_all_dom += 1
    print(f"n=4: b1=0 count={n4_b1_0}, has_free={n4_b1_0_has_free}, all_dom={n4_b1_0_all_dom}")

    # Check n=6 — this is where the gap matters
    print(f"\n=== n=6 CHECK (sampled) ===")
    n6 = 6
    n6_arcs = n6*(n6-1)//2  # 15
    import random
    random.seed(42)
    n6_b1_0 = 0
    n6_b1_0_has_free = 0
    for _ in range(5000):
        bits = random.randint(0, (1 << n6_arcs) - 1)
        T = tournament_from_bits(n6, bits)
        b1 = compute_b1(T, n6)
        if b1 == 0:
            n6_b1_0 += 1
            cycles = find_3cycles(T, n6)
            if cycles:
                has_free = any(len(is_dominated(T, n6, c)) == 0 for c in cycles)
                if has_free:
                    n6_b1_0_has_free += 1
    print(f"n=6 (5000 sampled): b1=0 count={n6_b1_0}, has_free={n6_b1_0_has_free}")
    if n6_b1_0_has_free > 0:
        print(f"  Free cycle rate among b1=0: {n6_b1_0_has_free/n6_b1_0:.1%}")
        print("  This confirms Case 2 is NON-vacuous at n=6 (as expected)")

if __name__ == "__main__":
    full_analysis()
