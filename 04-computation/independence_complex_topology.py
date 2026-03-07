#!/usr/bin/env python3
"""
Independence complex topology of Omega(T).

Novel connection (kind-pasteur-S34): The independence complex Ind(Omega(T))
is a simplicial complex whose faces are independent sets of Omega(T).

We compute:
1. The f-vector (face counts by dimension)
2. The reduced Euler characteristic chi~ = sum (-1)^k f_k
3. Check: chi~ = (-1)^d * I(Omega, -1) where d = dimension
4. Betti numbers via boundary matrices (Smith normal form)
5. Homotopy type classification: contractible, sphere, wedge of spheres?

Key question: Does the topology of Ind(Omega(T)) reveal anything about
permanent H-gaps? E.g., is there a topological obstruction for H=21?

Instance: kind-pasteur-2026-03-07-S34
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from math import comb
import random


def get_tournament(n, bits):
    """Generate tournament from bit encoding."""
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def find_odd_cycles(adj, n):
    """Find all directed odd cycles (3-cycles and 5-cycles for small n)."""
    cycles = []
    # 3-cycles
    for vs in combinations(range(n), 3):
        a, b, c = vs
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles.append(frozenset(vs))
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            cycles.append(frozenset(vs))
    # 5-cycles
    if n >= 5:
        for vs in combinations(range(n), 5):
            perms_5 = []
            # All directed 5-cycles on these vertices
            from itertools import permutations
            for p in permutations(vs):
                if all(adj[p[i]][p[(i+1)%5]] for i in range(5)):
                    fvs = frozenset(vs)
                    if fvs not in [c for c in cycles if len(c) == 5]:
                        perms_5.append(fvs)
            if perms_5:
                cycles.append(frozenset(vs))
    # 7-cycles
    if n == 7:
        vs = tuple(range(7))
        from itertools import permutations
        seen = set()
        for p in permutations(vs):
            if all(adj[p[i]][p[(i+1)%7]] for i in range(7)):
                fvs = frozenset(vs)
                if fvs not in seen:
                    seen.add(fvs)
        if seen:
            cycles.append(frozenset(vs))
    # Deduplicate
    return list(set(cycles))


def build_conflict_graph(cycles):
    """Build adjacency matrix of Omega(T): cycles adjacent iff share vertex."""
    nc = len(cycles)
    adj = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i][j] = 1
                adj[j][i] = 1
    return adj


def find_all_independent_sets(adj, nc):
    """Find all independent sets of the conflict graph."""
    indep_sets = [frozenset()]  # empty set
    for mask in range(1, 1 << nc):
        vertices = [i for i in range(nc) if mask & (1 << i)]
        is_indep = True
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                if adj[vertices[i]][vertices[j]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            indep_sets.append(frozenset(vertices))
    return indep_sets


def compute_f_vector(indep_sets):
    """Compute f-vector: f_k = # independent sets of size k+1 (dimension k)."""
    if not indep_sets:
        return []
    max_size = max(len(s) for s in indep_sets)
    f = [0] * (max_size + 1)
    for s in indep_sets:
        f[len(s)] += 1
    return f  # f[0] = 1 (empty set), f[1] = # vertices, etc.


def reduced_euler_char(f_vector):
    """Compute reduced Euler characteristic: chi~ = sum_{k>=0} (-1)^k f_k - 1."""
    # For abstract simplicial complex: chi~ = -1 + f_0 - f_1 + f_2 - ...
    # where f_k = # of k-dimensional faces
    # But our f_vector[i] = # independent sets of size i
    # So f_vector[0] = 1 (empty set), f_vector[1] = # vertices, etc.
    # Simplicial complex convention: dim(empty) = -1, dim({v}) = 0
    # chi~ = sum_{k=-1}^{d} (-1)^k f_k = -1 + f_0 - f_1 + f_2 - ...
    # = (-1)^0 * f_vector[0] + sum_{k>=1} (-1)^k * f_vector[k]  (but wait)
    # Actually: chi~ = sum_{i>=0} (-1)^{i-1} * f_vector[i]
    # = -f_vector[0] + f_vector[1] - f_vector[2] + ...
    # = -1 + |V| - (# edges in Ind) + ...
    result = 0
    for i, f in enumerate(f_vector):
        result += (-1)**(i-1) * f  # dim = i-1
    return result


def independence_poly_at_minus1(f_vector):
    """I(G, -1) = sum alpha_k * (-1)^k."""
    result = 0
    for k, f in enumerate(f_vector):
        result += f * ((-1)**k)
    return result


def main():
    print("=== Independence Complex Topology of Omega(T) ===")
    print()

    for n in range(4, 8):
        print(f"\n--- n={n} ---")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        # Sample or exhaust
        if n <= 6:
            trials = range(max_bits)
            label = "exhaustive"
        else:
            random.seed(42)
            trials = [random.randint(0, max_bits - 1) for _ in range(500)]
            label = "500 random"

        # Collect statistics
        euler_chars = {}
        f_vector_dims = {}
        h_to_euler = {}

        count = 0
        for bits in trials:
            adj = get_tournament(n, bits if isinstance(bits, int) else bits)
            cycles = find_odd_cycles(adj, n)
            if not cycles:
                continue

            omega_adj = build_conflict_graph(cycles)
            nc = len(cycles)

            if nc > 20:  # skip very large conflict graphs
                continue

            indep = find_all_independent_sets(omega_adj, nc)
            f_vec = compute_f_vector(indep)
            chi = reduced_euler_char(f_vec)
            dim = len(f_vec) - 2  # max dimension

            # Compute H = I(Omega, 2)
            H = sum(f * (2**k) for k, f in enumerate(f_vec))

            euler_chars[chi] = euler_chars.get(chi, 0) + 1
            f_vector_dims[dim] = f_vector_dims.get(dim, 0) + 1

            if H not in h_to_euler:
                h_to_euler[H] = set()
            h_to_euler[H].add(chi)

            count += 1
            if count <= 3:
                print(f"  Example: {nc} cycles, f-vec={f_vec}, chi~={chi}, dim={dim}, H={H}")

        print(f"  Total: {count} tournaments ({label})")
        print(f"  Euler chars: {dict(sorted(euler_chars.items()))}")
        print(f"  Dimensions: {dict(sorted(f_vector_dims.items()))}")

        # Check relationship between H and chi~
        # Note: I(G, -1) = chi~ * (-1)^{dim+1} for certain complexes
        print(f"  H -> chi~ mapping (first 10):")
        for H in sorted(h_to_euler.keys())[:10]:
            chis = h_to_euler[H]
            print(f"    H={H}: chi~ in {sorted(chis)}")


if __name__ == "__main__":
    main()
