#!/usr/bin/env python3
"""
verify_factorization_s263.py — Verify the Burnside exponent factorization
opus-2026-03-23-S263

THEOREM: For odd-cycle-type σ with k cycles:
    Fix_even(σ) = 2^{arc_orbits(σ) - k + 1}

This gives V_even as a CLOSED-FORM Burnside sum (no constraint matrix needed):
    V_even = (1/n!) Σ_{σ odd} ccs(σ) × 2^{arc_orbits(σ) - k(σ) + 1}

Also explore: what does this imply for the edge formula E(G_n)?
"""

from math import comb, factorial, gcd
from collections import Counter

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield []; return
    for f in range(min(n, mx), 0, -1):
        for r in partitions(n - f, f): yield [f] + r

def ccs(n, ct):
    c = Counter(ct); r = factorial(n)
    for l, k in c.items(): r //= (l ** k) * factorial(k)
    return r

def arc_orbits(ct):
    total = 0
    for i in range(len(ct)):
        total += (ct[i] - 1) // 2
        for j in range(i + 1, len(ct)):
            total += gcd(ct[i], ct[j])
    return total

def V_tournament(n):
    """Standard Burnside for A000568."""
    nf = factorial(n)
    total = 0
    for ct in partitions(n):
        ct_list = list(ct)
        if any(c % 2 == 0 for c in ct_list): continue
        total += ccs(n, ct_list) * (2 ** arc_orbits(ct_list))
    return total // nf

def V_even_new(n):
    """NEW formula using factorization: Fix_even = 2^{arc_orbits - k + 1}."""
    nf = factorial(n)
    total = 0
    for ct in partitions(n):
        ct_list = list(ct)
        # ALL cycle types contribute to even graphs, not just odd ones!
        # We need to compute Fix_even for ALL permutations, not just odd-cycle-type.
        # The factorization Fix_even = 2^{arc_orbits - k + 1} was derived for
        # ODD cycle types only. For even cycle types, we need the constraint matrix.
        # But wait: even cycle types ALSO fix some even graphs.
        # Let me compute both ways and compare.
        pass

    # Method 1: Use factorization for ALL cycle types (conjecture)
    total_all = 0
    for ct in partitions(n):
        ct_list = list(ct)
        k = len(ct_list)  # number of cycles
        ao = arc_orbits_general(ct_list)  # need general formula
        # For even cycle types, arc_orbits formula is different
        # Actually, the edge_orbits formula works for all cycle types:
        eo = edge_orbits(ct_list)
        # Conjecture: Fix_even = 2^{eo - rank(A)} for all cycle types
        # And the factorization arc_orbits = eo - rank(A) + (k-1) only holds for odd types

    # Method 2: Use factorization only for odd cycle types,
    # compute the rest via constraint matrix
    # Actually, EVERY permutation can fix some even graph.
    # The formula Fix_even(σ) = 2^{eo - rank(A_σ)} is general.
    # My theorem says: for ODD cycle types, rank(A_σ) = eo - arc_orbits + (k-1)
    # For EVEN cycle types, rank(A_σ) is something else.

    # Let me just compute V_even using the closed-form for odd types
    # and the constraint matrix for even types.
    pass

def edge_orbits(ct):
    total = 0
    for i in range(len(ct)):
        total += ct[i] // 2
        for j in range(i + 1, len(ct)):
            total += gcd(ct[i], ct[j])
    return total

def arc_orbits_general(ct):
    """For ALL cycle types (not just odd).
    For odd c: (c-1)/2 self-pairing orbits
    For even c: c/2 self-pairing orbits (but some are orientation-constrained)
    Between cycles: gcd(c_i, c_j) orbits
    """
    # This is the same formula as edge_orbits for general cycle types
    # because edge orbits and arc orbits are related by:
    # arc_orbits = edge_orbits for odd cycle types (each edge orbit splits into 2 arc orbits)
    # For even cycle types, some edge orbits have self-pairing (arc can't be independently oriented)
    # So Fix_tournament = 0 for even cycle types.
    return edge_orbits(ct)

if __name__ == '__main__':
    print("=" * 70)
    print("  VERIFICATION: V_even via factorization")
    print("  opus-2026-03-23-S263")
    print("=" * 70)

    # Test: compute V_even using ONLY the closed form for odd cycle types
    # V_even_odd = (1/n!) Σ_{σ odd} ccs(σ) × 2^{arc_orbits - k + 1}
    # Compare with V_even computed via constraint matrix (from previous session)

    print(f"\n{'n':>3} {'V_tourn':>12} {'V_even_odd':>12} {'V_even_all':>12} {'match':>6}")

    for n in range(3, 16):
        nf = factorial(n)

        # V_tournament
        vt = V_tournament(n)

        # V_even via factorization (odd cycle types only)
        v_even_odd = 0
        for ct in partitions(n):
            ct_list = list(ct)
            if any(c % 2 == 0 for c in ct_list): continue
            k = len(ct_list)
            ao = arc_orbits(ct_list)
            v_even_odd += ccs(n, ct_list) * (2 ** (ao - k + 1))
        v_even_odd //= nf

        # V_even via ALL cycle types (need constraint matrix for even types)
        # For now, use the formula that V_even_all includes even-type contributions
        # These are computed via the constraint matrix, but for comparison,
        # let me compute using the general edge_orbits formula
        v_even_all_approx = 0
        for ct in partitions(n):
            ct_list = list(ct)
            k = len(ct_list)
            eo = edge_orbits(ct_list)
            # Conjecture: Fix_even = 2^{eo - k + 1} for ALL cycle types?
            v_even_all_approx += ccs(n, ct_list) * (2 ** (eo - k + 1))
        v_even_all_approx //= nf

        # Check if v_even_odd matches the known V_even (A002854)
        # A002854: 2, 3, 7, 16, 54, 243, 2038, 33120, 1182004, 87723296
        a002854 = {3:2, 4:3, 5:7, 6:16, 7:54, 8:243, 9:2038, 10:33120,
                   11:1182004, 12:87723296}

        match_odd = "✓" if n in a002854 and v_even_odd == a002854[n] else "?"
        match_all = "✓" if n in a002854 and v_even_all_approx == a002854[n] else "?"

        print(f"{n:>3} {vt:>12} {v_even_odd:>12} {v_even_all_approx:>12} "
              f"odd={match_odd} all={match_all}")

    # Key relationships
    print(f"\n{'='*70}")
    print("  KEY RELATIONSHIPS")
    print(f"{'='*70}")

    for n in range(3, 10):
        nf = factorial(n)
        vt = V_tournament(n)

        # Compute contribution from identity permutation
        identity_ao = comb(n, 2)  # arc_orbits for identity = C(n,2)
        identity_k = n  # k cycles (all fixed points)
        identity_fix_even = 2 ** (identity_ao - n + 1)  # = 2^{C(n-1,2)}

        # Total labeled even graphs = 2^{C(n-1,2)}
        total_labeled_even = 2 ** comb(n-1, 2)

        print(f"  n={n}: 2^{{C(n-1,2)}} = {total_labeled_even}, "
              f"identity_fix_even = {identity_fix_even}, "
              f"match = {total_labeled_even == identity_fix_even}")

    # The factorization theorem
    print(f"\n{'='*70}")
    print("  FACTORIZATION THEOREM (final form)")
    print(f"{'='*70}")
    print("""
  For all-odd cycle type σ with k cycles:

    Fix_tournament(σ) = 2^{arc_orbits(σ)}
    Fix_even_graph(σ) = 2^{arc_orbits(σ) - k + 1}

  Therefore:
    Fix_tournament(σ) = Fix_even_graph(σ) × 2^{k-1}

  The 2^{k-1} factor counts the CUT SPACE free bits.

  The Burnside sum for V_even:
    V_even = (1/n!) Σ_{all-odd σ} ccs(σ) × 2^{arc_orbits(σ) - k(σ) + 1}
           + (1/n!) Σ_{some-even σ} ccs(σ) × Fix_even(σ)

  The first sum (odd-cycle types only) misses even-cycle contributions.
  For even cycle types, Fix_even ≠ 0 in general.
  The FULL formula requires the constraint matrix for even types.
""")

    print("DONE.")
