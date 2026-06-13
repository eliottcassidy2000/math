"""
h21_cobipartite_verify.py — Verify whether alpha_3=0 forces co-bipartite Omega.

Claim: If independence_number(G) <= 2, then G is co-bipartite (union of 2 cliques).

Proof attempt: If ind(G) <= 2, then the complement G_bar has clique_number <= 2.
A graph with clique number <= 2 is triangle-free.
By König... no. Triangle-free != bipartite in general.

COUNTEREXAMPLE: C_5 (5-cycle) is triangle-free but NOT bipartite (odd cycle).
So complement being triangle-free does NOT mean complement is bipartite.
So ind(G) <= 2 does NOT imply G is co-bipartite!

The correct statement: ind(G) <= 2 => complement G_bar is triangle-free.
Triangle-free graphs include bipartite graphs AND odd-cycle-containing graphs.

So the co-bipartite constraint was TOO STRONG. Let me re-analyze.

What IS true: by Ramsey R(3,3)=6, if G has >= 6 vertices and ind(G) <= 2,
then G has a triangle (K_3). But G need not be co-bipartite.

For the alpha_2 bound: alpha_2 = number of independent pairs in G.
For G with ind(G) <= 2:
alpha_2 = |non-edges of G| = C(n,2) - |E(G)|.

Turán's theorem for complement: G_bar is triangle-free => |E(G_bar)| <= n^2/4.
So alpha_2 = |E(G_bar)| <= n^2/4 = alpha_1^2/4.

This bound is tight (achieved by complete bipartite K_{n/2,n/2} as complement).
But the complement need not actually BE bipartite.

Author: opus-2026-03-07-S43
"""

# Verify: for small graphs, enumerate all graphs with ind(G) <= 2
# and check their alpha_2 values

from itertools import combinations

def independence_number(n, edges):
    """Compute independence number of graph on n vertices with given edges."""
    edge_set = set(edges)
    for k in range(n, 0, -1):
        for S in combinations(range(n), k):
            independent = True
            for i in range(len(S)):
                for j in range(i+1, len(S)):
                    if (S[i], S[j]) in edge_set or (S[j], S[i]) in edge_set:
                        independent = False
                        break
                if not independent:
                    break
            if independent:
                return k
    return 0

def is_cobipartite(n, edges):
    """Check if graph is co-bipartite (union of 2 cliques)."""
    edge_set = set(edges)
    # Try all 2-partitions
    for mask in range(1 << n):
        clique_A = [v for v in range(n) if mask & (1 << v)]
        clique_B = [v for v in range(n) if not (mask & (1 << v))]
        # Check both are cliques in G
        is_valid = True
        for group in [clique_A, clique_B]:
            for i in range(len(group)):
                for j in range(i+1, len(group)):
                    if (group[i], group[j]) not in edge_set and (group[j], group[i]) not in edge_set:
                        is_valid = False
                        break
                if not is_valid:
                    break
        if is_valid:
            return True
    return False

print("Checking: does ind(G) <= 2 imply co-bipartite?")
print()

# Enumerate graphs on n=5 vertices
n = 5
all_possible_edges = list(combinations(range(n), 2))
num_edges = len(all_possible_edges)

counterexamples = 0
total_ind2 = 0

for mask in range(1 << num_edges):
    edges = [all_possible_edges[i] for i in range(num_edges) if mask & (1 << i)]
    ind = independence_number(n, edges)
    if ind <= 2:
        total_ind2 += 1
        if not is_cobipartite(n, edges):
            counterexamples += 1
            alpha2 = num_edges - len(edges)  # non-edges
            if counterexamples <= 5:
                print(f"  Counterexample: n={n}, edges={edges}, α₂={alpha2}")

print(f"\nn={n}: {total_ind2} graphs with ind<=2, {counterexamples} NOT co-bipartite")

# So the co-bipartite claim is FALSE!
# But what alpha_2 values are achievable for ind(G) <= 2?

print(f"\nAchievable α₂ values for graphs on n={n} with ind(G) <= 2:")
a2_values = set()
for mask in range(1 << num_edges):
    edges = [all_possible_edges[i] for i in range(num_edges) if mask & (1 << i)]
    ind = independence_number(n, edges)
    if ind <= 2:
        a2 = num_edges - len(edges)
        a2_values.add(a2)

print(f"  α₂ ∈ {sorted(a2_values)}")
print(f"  Turán bound: α₂ <= {n*n//4} = {n*n//4}")
print(f"  Co-bipartite values: {sorted(p*(n-p) for p in range(n+1))}")

# Now do n=6
n = 6
all_possible_edges = list(combinations(range(n), 2))
num_edges = len(all_possible_edges)
print(f"\nn={n}: Checking achievable α₂ values for ind(G) <= 2...")

a2_values = set()
total_ind2 = 0
not_cobip = 0

for mask in range(1 << num_edges):
    edges = [all_possible_edges[i] for i in range(num_edges) if mask & (1 << i)]
    ind = independence_number(n, edges)
    if ind <= 2:
        total_ind2 += 1
        a2 = num_edges - len(edges)
        a2_values.add(a2)
        if not is_cobipartite(n, edges):
            not_cobip += 1

print(f"  Total with ind<=2: {total_ind2}")
print(f"  NOT co-bipartite: {not_cobip}")
print(f"  α₂ ∈ {sorted(a2_values)}")
print(f"  Turán bound: α₂ <= {n*n//4}")
cobip_vals = sorted(set(p*(n-p) for p in range(n+1)))
print(f"  Co-bipartite α₂: {cobip_vals}")
extra = sorted(a2_values - set(cobip_vals))
if extra:
    print(f"  EXTRA α₂ values (non-co-bipartite): {extra}")
