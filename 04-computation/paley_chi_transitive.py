#!/usr/bin/env python3
"""
paley_chi_transitive.py — opus-2026-03-13-S71b

Compute transitive sub-tournament counts t_k for Paley tournaments P_p
(p ≡ 3 mod 4 only).

Compare with GLMY Omega_d to test whether t_{d+1} = Omega_d.
If not, compare the alternating sums (Euler characteristics).

NOTE: GLMY Omega_d uses a WEAKER condition than "all faces allowed" — it only
requires ∂σ ∈ A_{d-1} (the alternating sum, not each face individually).
So Omega_d ≥ t_{d+1} in general.

The question is whether chi_transitive = chi_GLMY despite different Omega.
"""

from itertools import combinations
import time

def paley_tournament(p):
    """Return adjacency: adj[v] = set of vertices v beats.
    Only valid for p ≡ 3 mod 4."""
    assert p % 4 == 3, f"Paley tournament only exists for p ≡ 3 (mod 4), got p={p}"
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    adj = [set() for _ in range(p)]
    for a in range(p):
        for b in range(p):
            if a != b and (b - a) % p in qr:
                adj[a].add(b)
    return adj, qr

def count_transitive(adj, p, k):
    """Count transitive k-vertex sub-tournaments."""
    count = 0
    for subset in combinations(range(p), k):
        vertices = list(subset)
        scores = {}
        for v in vertices:
            scores[v] = sum(1 for w in vertices if w != v and w in adj[v])
        score_vals = sorted(scores.values())
        if score_vals == list(range(k)):
            count += 1
    return count

def main():
    print("=" * 70)
    print("TRANSITIVE SUB-TOURNAMENT COUNTS FOR PALEY PRIMES (p ≡ 3 mod 4)")
    print("=" * 70)

    # Known GLMY Omega/p values for comparison
    known_omega = {
        7: [1, 3, 6, 9, 9, 6, 3],   # From circulant_homology
        11: [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30],
    }

    for p in [3, 7, 11, 19, 23]:
        adj, qr = paley_tournament(p)
        m = (p - 1) // 2
        t0 = time.time()

        print(f"\nP_{p} (m={m}, QR={sorted(qr)}):")

        t_vals = []
        chi_t = 0
        for k in range(1, p + 1):
            t_k = count_transitive(adj, p, k)
            t_vals.append(t_k)
            d = k - 1
            chi_t += (-1) ** d * t_k
            if t_k > 0 or k <= 4:
                print(f"  t_{k:2d} = {t_k:>12d}  (d={d})")
            if p >= 23 and k > 12 and t_k == 0:
                # Skip printing zeros for large p
                break

        elapsed = time.time() - t0
        print(f"  t_k/p: {[t // p for t in t_vals if t > 0]}")
        print(f"  chi_transitive = {chi_t}")
        print(f"  chi_t / p = {chi_t / p}")

        if p in known_omega:
            omega = known_omega[p]
            chi_glmy = sum((-1)**d * o * p for d, o in enumerate(omega))
            print(f"  GLMY Omega/p = {omega}")
            print(f"  chi_GLMY = {chi_glmy}")
            print(f"  chi match: {chi_t == chi_glmy}")

            # Show diff at each degree
            print(f"  Omega_d/p vs t_{'{d+1}'}/p:")
            for d in range(min(len(omega), len(t_vals))):
                o_per_p = omega[d] if d < len(omega) else 0
                t_per_p = t_vals[d] // p if d < len(t_vals) else 0
                diff = o_per_p - t_per_p
                if o_per_p > 0 or t_per_p > 0:
                    print(f"    d={d}: Omega/p={o_per_p}, t/p={t_per_p}, diff={diff}")

        print(f"  Time: {elapsed:.1f}s")

    # Pattern analysis
    print(f"\n{'='*70}")
    print("PATTERN ANALYSIS")
    print("=" * 70)
    print("For Paley primes p ≡ 3 mod 4:")
    print("  p=3:  chi_t = ? (computed above)")
    print("  p=7:  chi_t = ? (computed above)")
    print("  p=11: chi_t = ? (computed above)")
    print("  Conjecture: chi(P_p) = p for p ≡ 3 mod 4, p >= 7")

    print("\nDONE.")

if __name__ == "__main__":
    main()
