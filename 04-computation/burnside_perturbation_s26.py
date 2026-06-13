#!/usr/bin/env python3
"""
opus-2026-04-05-S26: Burnside Perturbation Theory — a general principle.

THE ABSTRACTION:
In ANY Burnside counting problem a(n) = (1/|G|) Σ_{g∈G} |Fix(g)|,
the identity element dominates. Non-identity elements contribute
"perturbation corrections" that decay exponentially.

This gives: a(n) = base_term × (1 + Σ_k ε_k(n))

where ε_k → 0 exponentially. For large n, only O(1) terms matter.

APPLICATION 1: A000568 (tournaments) — we proved this above
APPLICATION 2: A000088 (simple graphs) — same principle!
APPLICATION 3: A000273 (digraphs)
APPLICATION 4: A000665 (3-uniform hypergraphs)

In each case, the "base term" is base^{edges}/n! and the corrections
come from non-identity permutations grouped by cycle type.

CREATIVE TWIST: View this as a RENORMALIZATION GROUP.
The "bare" count (identity term) is renormalized by "quantum corrections"
(non-identity permutations). The RG flow is n → n+1.
Each "defect type" (cycle type) has a "running coupling constant"
that decreases with n. The theory is "asymptotically free" —
at large n, the bare count dominates.
"""

from math import factorial, gcd, comb
import time

def z_lambda(parts_mults):
    z = 1
    for k, m in parts_mults:
        z *= k**m * factorial(m)
    return z

def enumerate_all_partitions(n, odd_only=False):
    """Enumerate all partitions (or odd-part partitions) of n."""
    results = []
    if odd_only:
        parts = list(range(n if n % 2 else n-1, 0, -2))
    else:
        parts = list(range(n, 0, -1))

    def recurse(idx, remaining, current):
        if remaining == 0:
            results.append(list(current))
            return
        if idx >= len(parts):
            return
        k = parts[idx]
        if k > remaining:
            recurse(idx + 1, remaining, current)
            return
        recurse(idx + 1, remaining, current)
        for m in range(1, remaining // k + 1):
            recurse(idx + 1, remaining - m*k, current + [(k, m)])

    recurse(0, n, [])
    return results

def edge_count_tournament(parts_mults):
    """Edge count for tournaments (A000568): gcd_sum + self_floor."""
    total = 0
    for i, (ki, mi) in enumerate(parts_mults):
        total += mi * (ki - 1) // 2
        total += mi * (mi - 1) // 2 * ki
        for j in range(i+1, len(parts_mults)):
            kj, mj = parts_mults[j]
            total += mi * mj * gcd(ki, kj)
    return total

def edge_count_graph(parts_mults):
    """Edge count for simple graphs (A000088): same formula as tournaments."""
    return edge_count_tournament(parts_mults)

def edge_count_digraph(parts_mults):
    """Edge count for digraphs (A000273): 2×gcd_sum + sum_v - num_v."""
    gcd_sum = 0
    sum_v = 0
    num_v = 0
    for i, (ki, mi) in enumerate(parts_mults):
        sum_v += mi * ki
        num_v += mi
        gcd_sum += mi * (mi - 1) // 2 * ki
        for j in range(i+1, len(parts_mults)):
            kj, mj = parts_mults[j]
            gcd_sum += mi * mj * gcd(ki, kj)
    return 2 * gcd_sum + sum_v - num_v

def burnside_count(n, edge_func, base, odd_only=False):
    """General Burnside counting: a(n) = (1/n!) Σ_λ (n!/z_λ) × base^{e(λ)}."""
    if n <= 1:
        return 1
    nfact = factorial(n)
    parts = enumerate_all_partitions(n, odd_only=odd_only)
    total = 0
    for p in parts:
        z = z_lambda(p)
        e = edge_func(p)
        total += (nfact // z) * (base ** e)
    return total // nfact

def perturbation_analysis(n, edge_func, base, odd_only=False):
    """
    Analyze the perturbation structure: how much does each
    "defect type" contribute relative to the identity?
    """
    nfact = factorial(n)
    parts = enumerate_all_partitions(n, odd_only=odd_only)

    # Identity partition
    if odd_only:
        id_part = [(1, n)]
    else:
        id_part = [(1, n)]
    e_id = edge_func(id_part)
    w_id = (nfact // z_lambda(id_part)) * (base ** e_id)

    results = []
    for p in parts:
        z = z_lambda(p)
        e = edge_func(p)
        w = (nfact // z) * (base ** e)
        ratio = w / w_id if w_id > 0 else 0
        results.append((p, e, w, ratio))

    results.sort(key=lambda x: -x[3])
    return results, w_id

# ═══════════════════════════════════════════════════════
# TESTS
# ═══════════════════════════════════════════════════════

# Known values
A000568 = {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536, 10:9733056}
A000088 = {1:1, 2:2, 3:4, 4:11, 5:34, 6:156, 7:1044, 8:12346, 9:274668, 10:12005168}
A000273 = {1:1, 2:3, 3:16, 4:218, 5:9608, 6:1540944, 7:882033440, 8:1793359192848}

if __name__ == "__main__":
    print("="*70)
    print("BURNSIDE PERTURBATION THEORY — GENERAL PRINCIPLE")
    print("="*70)

    # Verify A000568
    print("\n── A000568 (Tournaments, base=2, odd parts) ──")
    for n in range(3, 11):
        result = burnside_count(n, edge_count_tournament, 2, odd_only=True)
        expected = A000568.get(n)
        match = "✓" if result == expected else f"✗ (got {result})"
        print(f"  a({n}) = {result:>12d}  {match}")

    # Verify A000088
    print("\n── A000088 (Simple graphs, base=2, all parts) ──")
    for n in range(1, 11):
        result = burnside_count(n, edge_count_graph, 2, odd_only=False)
        expected = A000088.get(n)
        match = "✓" if result == expected else f"✗ (got {result})"
        print(f"  a({n}) = {result:>12d}  {match}")

    # Verify A000273
    print("\n── A000273 (Digraphs, base=2, all parts) ──")
    for n in range(1, 8):
        result = burnside_count(n, edge_count_digraph, 2, odd_only=False)
        expected = A000273.get(n)
        match = "✓" if result == expected else f"✗ (got {result})"
        print(f"  a({n}) = {result:>15d}  {match}")

    # ═══════════════════════════════════════════════════════
    # PERTURBATION ANALYSIS
    # ═══════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PERTURBATION STRUCTURE — comparing sequences")
    print("="*70)

    sequences = [
        ("A000568 (Tournaments)", edge_count_tournament, 2, True),
        ("A000088 (Graphs)", edge_count_graph, 2, False),
        ("A000273 (Digraphs)", edge_count_digraph, 2, False),
    ]

    for name, efunc, base, odd_only in sequences:
        print(f"\n── {name} at n=10 ──")
        results, w_id = perturbation_analysis(10, efunc, base, odd_only)

        # Show top 5 contributions
        total = sum(w for _, _, w, _ in results)
        cumulative = 0
        for i, (p, e, w, ratio) in enumerate(results[:8]):
            cumulative += w
            frac = cumulative / total
            # Simplify partition display
            pstr = "+".join(f"{k}^{m}" for k, m in p)
            print(f"  #{i+1}: λ={pstr:20s} e={e:4d} ratio={ratio:12.8f} cum={frac:.8f}")

        identity_frac = results[0][2] / total
        print(f"  Identity fraction: {identity_frac:.6f}")
        print(f"  Top-3 cumulative: {sum(r[2] for r in results[:3])/total:.8f}")

    # ═══════════════════════════════════════════════════════
    # THE RENORMALIZATION GROUP PICTURE
    # ═══════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("RG FLOW: Identity fraction vs n")
    print("="*70)

    for name, efunc, base, odd_only in sequences:
        print(f"\n── {name} ──")
        for n in range(4, 16):
            results, w_id = perturbation_analysis(n, efunc, base, odd_only)
            total = sum(w for _, _, w, _ in results)
            id_frac = results[0][2] / total
            n_parts = len(results)

            # "Coupling constant" = 1 - identity fraction
            coupling = 1 - id_frac
            print(f"  n={n:2d}: id_frac={id_frac:.8f}, coupling={coupling:.6e}, "
                  f"#partitions={n_parts}")

    # ═══════════════════════════════════════════════════════
    # PREDICT LARGE-n VALUES USING IDENTITY + CORRECTIONS
    # ═══════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PRACTICAL: Approximate values for large n")
    print("="*70)

    for name, efunc, base, odd_only in sequences:
        print(f"\n── {name} ──")
        for n in [20, 50, 100]:
            nfact = factorial(n)
            # Identity partition
            id_part = [(1, n)]
            e_id = efunc(id_part)
            identity_val = (base ** e_id) // nfact

            # Identity + 1st correction (single transposition for graphs/digraphs,
            # single 3-cycle for tournaments)
            if odd_only:
                corr_part = [(3, 1), (1, n-3)]
            else:
                corr_part = [(2, 1), (1, n-2)]
            e_corr = efunc(corr_part)
            z_corr = z_lambda(corr_part)
            corr_val = ((nfact // z_corr) * (base ** e_corr) + (base ** e_id)) // nfact

            ndigits = len(str(identity_val))
            ratio = (corr_val - identity_val) / identity_val if identity_val > 0 else 0
            print(f"  n={n:3d}: ~{ndigits} digits, 1st correction = {ratio:.4e} of identity")

    print("\n\nDONE.")
