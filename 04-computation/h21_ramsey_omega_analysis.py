"""
h21_ramsey_omega_analysis.py — Ramsey-theoretic analysis of Omega(T) for H=21 gap.

Key observation: H=21 requires w=10 with alpha_3=0.
alpha_3=0 means independence number of Omega(T) <= 2.
So complement of Omega is triangle-free (K_3-free).

By Ramsey R(3,3)=6: any graph on >=6 vertices has either a clique of size 3
or an independent set of size 3.
If ind(Omega) <= 2, then for alpha_1 >= 6, Omega must contain K_3
(3 pairwise conflicting odd cycles).

Key question: what does K_3 in Omega mean for tournaments?
Three pairwise vertex-intersecting odd cycles. How constraining is this?

Author: opus-2026-03-07-S43
"""
from math import comb
from itertools import combinations

def analyze_w_decompositions(w, max_alpha1=None):
    """Analyze all possible (alpha_1, alpha_2) decompositions for a given w."""
    if max_alpha1 is None:
        max_alpha1 = w + 1

    results = []
    for a1 in range(max_alpha1 + 1):
        rem = w - a1
        if rem < 0 or rem % 2 != 0:
            continue
        a2 = rem // 2

        # Basic feasibility: a2 <= C(a1, 2)
        if a2 > comb(a1, 2):
            status = "INFEASIBLE (pairs exceed vertices)"
        else:
            # Turán-type bound from alpha_3=0:
            # ind(Omega) <= 2 => complement is K3-free
            # => complement has <= floor(a1^2/4) edges (Turán)
            # => a2 (= non-edges of Omega = edges of complement) <= floor(a1^2/4)
            turan_bound = a1 * a1 // 4
            if a2 > turan_bound:
                status = f"INFEASIBLE by Turán (a2={a2} > a1²/4={turan_bound})"
            else:
                # Check: edges of Omega = C(a1,2) - a2
                edges_omega = comb(a1, 2) - a2
                # Min degree of Omega: by complement being K3-free (Turán),
                # complement has max degree <= ...
                # Actually, complement being bipartite means it's 2-chromatic
                # So Omega's complement can be partitioned into 2 independent sets
                # <=> Omega can be covered by 2 cliques

                status = f"feasible (Omega: {a1} vertices, {edges_omega} edges)"

        results.append((a1, a2, status))

    return results

print("=" * 70)
print("RAMSEY-TURÁN ANALYSIS OF w-DECOMPOSITIONS")
print("=" * 70)

for w in [3, 10]:
    print(f"\n{'='*50}")
    print(f"w = {w} (H = {2*w+1})")
    print(f"{'='*50}")
    decomps = analyze_w_decompositions(w, max_alpha1=w)
    feasible_count = 0
    for a1, a2, status in decomps:
        print(f"  (α₁={a1}, α₂={a2}): {status}")
        if "feasible" in status.lower() and "INFEASIBLE" not in status:
            feasible_count += 1
    print(f"  Graph-theoretically feasible: {feasible_count}")

# Now analyze what tournament constraints add beyond graph theory
print("\n" + "=" * 70)
print("TOURNAMENT-SPECIFIC OBSTRUCTIONS NEEDED")
print("=" * 70)

print("""
For w=3 (H=7), graph-theoretically feasible: (3,0)
  (3,0): 3 cycles, all pairwise conflicting (K_3 in Omega).
  Tournament obstruction: THM-029 shows 3 mutually conflicting 3-cycles
  force alpha_1 >= 4. So (3,0) is impossible.

For w=10 (H=21), graph-theoretically feasible: (4,3), (6,2), (8,1), (10,0)
  (4,3): Part D — 4 cycles with 3 disjoint pairs creates a specific graph
         structure that forces additional cycles.
  (6,2): Part F — 6 cycles with only 2 disjoint pairs, independence <= 2.
         Very dense Omega with 13 edges on 6 vertices (of 15 possible).
  (8,1): Part N — 8 cycles with only 1 disjoint pair. Almost complete Omega
         (27 of 28 edges). One missing edge = one disjoint pair.
  (10,0): Part M — 10 mutually conflicting cycles. Complete graph K_10 in Omega.
          Each pair shares a vertex. Very constrained star-like structure.
""")

# New angle: can Turán + Ramsey eliminate more decompositions?
print("=" * 70)
print("ENHANCED BOUNDS: Turán + chromatic number")
print("=" * 70)

for w in [3, 5, 7, 10, 13, 21, 31]:
    decomps = analyze_w_decompositions(w, max_alpha1=min(w, 20))
    feasible = [(a1, a2) for a1, a2, s in decomps if "INFEASIBLE" not in s]
    graph_feasible = len(feasible)

    # Additional: for alpha_3=0, complement of Omega is bipartite
    # So Omega is co-bipartite (union of 2 cliques)
    # A co-bipartite graph on a1 vertices with 2 cliques of sizes p, q (p+q=a1)
    # has exactly p*q non-edges (between the two cliques)
    # So alpha_2 = p*q for some partition p+q=a1

    cobip_feasible = []
    for a1, a2 in feasible:
        # Check if a2 = p*q for some p+q=a1
        achievable = False
        for p in range(a1 + 1):
            q = a1 - p
            if p * q == a2:
                achievable = True
                break
        if achievable:
            cobip_feasible.append((a1, a2))

    print(f"w={w} (H={2*w+1}): {graph_feasible} graph-feasible, "
          f"{len(cobip_feasible)} co-bipartite feasible: {cobip_feasible}")

# KEY INSIGHT: For alpha_3=0, Omega is co-bipartite (union of 2 cliques).
# So alpha_2 = p*q where p+q=alpha_1 (cross-edges between the 2 cliques
# are the non-edges of Omega = independent pairs).
# This means alpha_2 must equal p*(alpha_1-p) for some p.
# Given alpha_1 + 2*alpha_2 = w:
# w = alpha_1 + 2*p*(alpha_1-p) = alpha_1 + 2p*alpha_1 - 2p^2

print("\n" + "=" * 70)
print("CO-BIPARTITE CONSTRAINT: alpha_2 = p*(alpha_1-p)")
print("=" * 70)

for w in [3, 10]:
    print(f"\nw={w}:")
    for a1 in range(w + 1):
        rem = w - a1
        if rem < 0 or rem % 2 != 0:
            continue
        a2 = rem // 2
        # Check a2 = p*(a1-p) for some integer p
        solutions = []
        for p in range(a1 + 1):
            if p * (a1 - p) == a2:
                solutions.append((p, a1 - p))
        if solutions:
            print(f"  (α₁={a1}, α₂={a2}): co-bipartite split = {solutions}")
        elif a2 <= comb(a1, 2):
            print(f"  (α₁={a1}, α₂={a2}): NOT co-bipartite realizable!")

# This is powerful: if alpha_3=0 forces co-bipartite Omega, then
# alpha_2 = p*(alpha_1-p) is an additional constraint that may eliminate
# some decompositions purely from graph theory!
