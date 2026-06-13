#!/usr/bin/env python3
"""
INV-028: Hard-core lattice gas connection to OCF.

H(T) = I(Omega(T), 2) means H(T) is the partition function of the hard-core
lattice gas on Omega(T) at fugacity lambda=2.

Key questions:
1. What statistical mechanics properties does Omega(T) have?
2. Is lambda=2 special (phase transition threshold, etc.)?
3. Does the Lovász Local Lemma or cluster expansion converge at lambda=2?
4. Can correlation decay give structural information?

The hard-core model partition function is:
  Z(G, lambda) = sum_{independent sets I} lambda^|I|

For lambda=2: Z(G,2) = I(G,2) = H(T) by OCF.

The cluster expansion converges when lambda < lambda_c(G) where lambda_c
depends on the max degree Delta: lambda_c >= 1/(Delta-1) (Shearer's bound)
or more refined: lambda < (Delta-1)^{Delta-1} / Delta^Delta (Lovász Local Lemma).

For Omega(T), what is the max degree?

Instance: opus-2026-03-05-S4b
"""
import sys, random
from itertools import combinations, permutations

def random_tournament(n, seed=None):
    if seed is not None:
        random.seed(seed)
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

def find_odd_cycles(T):
    n = len(T)
    cycles = []
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            start = combo[0]
            others = list(combo[1:])
            mo = len(others)
            if mo == 0:
                continue
            full = (1 << mo) - 1
            dp = [[0]*mo for _ in range(1 << mo)]
            for k, o in enumerate(others):
                if T[start][o]:
                    dp[1<<k][k] = 1
            for mask in range(1, full+1):
                for li in range(mo):
                    c = dp[mask][li]
                    if c == 0: continue
                    for ni in range(mo):
                        if mask & (1<<ni): continue
                        if T[others[li]][others[ni]]:
                            dp[mask|(1<<ni)][ni] += c
            cnt = sum(dp[full][li] for li in range(mo) if T[others[li]][start])
            if cnt > 0:
                cycles.append(frozenset(combo))
    return cycles

def build_conflict_graph(cycles):
    """Build adjacency list for conflict graph."""
    n = len(cycles)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if cycles[i] & cycles[j]:  # share a vertex
                adj[i].append(j)
                adj[j].append(i)
    return adj

def ham_count(T):
    n = len(T)
    FULL = (1 << n) - 1
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for mask in range(1, 1<<n):
        for last in range(n):
            c = dp[mask][last]
            if c == 0: continue
            for nxt in range(n):
                if mask & (1<<nxt): continue
                if T[last][nxt]:
                    dp[mask|(1<<nxt)][nxt] += c
    return sum(dp[FULL])

def ind_poly_at_x(cycles, adj, x):
    """Compute I(Omega, x) = sum_{indep sets} x^|set|."""
    n = len(cycles)
    total = 1  # empty set
    # Enumerate all independent sets by bitmask
    for mask in range(1, 1 << n):
        # Check independence
        is_indep = True
        verts = [i for i in range(n) if mask & (1 << i)]
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if verts[j] in adj[verts[i]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            total += x ** len(verts)
    return total

print("=" * 60)
print("HARD-CORE LATTICE GAS ANALYSIS OF Omega(T)")
print("=" * 60)

for n in [5, 6, 7]:
    print(f"\n--- n = {n} ---")
    random.seed(42)
    trials = 200 if n <= 6 else 50

    max_degrees = []
    num_cycles_list = []
    densities = []
    chromatic_numbers = []
    lambda_c_bounds = []

    for _ in range(trials):
        T = random_tournament(n)
        cycles = find_odd_cycles(T)
        if len(cycles) == 0:
            continue
        adj = build_conflict_graph(cycles)
        nc = len(cycles)
        num_cycles_list.append(nc)

        # Max degree of Omega(T)
        degrees = [len(adj[i]) for i in range(nc)]
        max_deg = max(degrees) if degrees else 0
        avg_deg = sum(degrees) / nc if nc > 0 else 0
        max_degrees.append(max_deg)

        # Edge count and density
        edges = sum(len(adj[i]) for i in range(nc)) // 2
        max_edges = nc * (nc - 1) // 2 if nc > 1 else 1
        density = edges / max_edges if max_edges > 0 else 0
        densities.append(density)

        # Shearer bound: lambda_c >= 1/(Delta-1) for Delta >= 2
        if max_deg >= 2:
            shearer = 1.0 / (max_deg - 1)
            # LLL bound: (Delta-1)^{Delta-1} / Delta^Delta
            lll = ((max_deg - 1) ** (max_deg - 1)) / (max_deg ** max_deg) if max_deg > 0 else float('inf')
            lambda_c_bounds.append((shearer, lll, max_deg))

    if num_cycles_list:
        print(f"  Avg #cycles: {sum(num_cycles_list)/len(num_cycles_list):.1f}")
        print(f"  Max degree range: {min(max_degrees)}-{max(max_degrees)}")
        print(f"  Avg max degree: {sum(max_degrees)/len(max_degrees):.1f}")
        print(f"  Avg density: {sum(densities)/len(densities):.3f}")

    if lambda_c_bounds:
        shearer_vals = [s for s, l, d in lambda_c_bounds]
        lll_vals = [l for s, l, d in lambda_c_bounds]
        print(f"  Shearer bound range: {min(shearer_vals):.3f}-{max(shearer_vals):.3f}")
        print(f"  LLL bound range: {min(lll_vals):.4f}-{max(lll_vals):.4f}")
        print(f"  lambda=2 > Shearer bound in {sum(1 for s,_,_ in lambda_c_bounds if 2 > s)}/{len(lambda_c_bounds)} cases")
        print(f"  lambda=2 > LLL bound in {sum(1 for _,l,_ in lambda_c_bounds if 2 > l)}/{len(lambda_c_bounds)} cases")

print("\n" + "=" * 60)
print("VERIFICATION: I(Omega, lambda) at various fugacities")
print("=" * 60)

for n in [5, 6, 7]:
    print(f"\n--- n = {n} ---")
    random.seed(42)
    trials = 50 if n <= 6 else 20

    for trial in range(min(5, trials)):
        T = random_tournament(n)
        H = ham_count(T)
        cycles = find_odd_cycles(T)
        adj = build_conflict_graph(cycles)

        I2 = ind_poly_at_x(cycles, adj, 2)
        I1 = ind_poly_at_x(cycles, adj, 1)  # total independent sets
        I3 = ind_poly_at_x(cycles, adj, 3)

        print(f"  T#{trial}: H={H}, I(2)={I2}, I(1)={I1}, I(3)={I3}, "
              f"#cyc={len(cycles)}, match={'YES' if H==I2 else 'NO'}")

print("\n" + "=" * 60)
print("CLUSTER EXPANSION ANALYSIS")
print("=" * 60)
print("The cluster expansion for log Z(G,lambda) converges when")
print("lambda * e * (Delta+1) < 1 (Kotecky-Preiss condition)")
print("or lambda < lambda_c(tree) = (Delta-1)^{Delta-1}/Delta^Delta")
print()
print("For lambda=2:")
for Delta in range(2, 15):
    kp = 1.0 / (2.718281828 * (Delta + 1))
    tree_bound = ((Delta-1)**(Delta-1)) / (Delta**Delta) if Delta > 0 else 0
    shearer = 1.0 / (Delta - 1) if Delta > 1 else float('inf')
    status_kp = "converges" if 2 < kp else "diverges"
    status_tree = "below" if 2 < tree_bound else "ABOVE"
    print(f"  Delta={Delta:2d}: KP bound={kp:.4f}, tree bound={tree_bound:.4f} ({status_tree}), Shearer={shearer:.4f}")

print("\nConclusion: lambda=2 is ALWAYS above the cluster expansion")
print("convergence radius for any Delta >= 2. This means the hard-core")
print("model at lambda=2 is in the 'non-perturbative' regime.")
print("Standard cluster expansion / polymer expansion methods do NOT apply.")
print("This suggests OCF requires genuinely non-perturbative methods.")

print("\n" + "=" * 60)
print("INDEPENDENCE NUMBER AND CLIQUE COVER")
print("=" * 60)

for n in [5, 6, 7, 8]:
    print(f"\n--- n = {n} ---")
    random.seed(42)
    trials = 100 if n <= 6 else (30 if n == 7 else 10)

    alpha_max = 0
    clique_sizes = []

    for _ in range(trials):
        T = random_tournament(n)
        cycles = find_odd_cycles(T)
        adj = build_conflict_graph(cycles)
        nc = len(cycles)

        # Independence number (max independent set)
        best = 0
        for mask in range(1 << min(nc, 20)):  # cap at 20 cycles
            verts = [i for i in range(min(nc, 20)) if mask & (1 << i)]
            is_indep = True
            for i in range(len(verts)):
                for j in range(i+1, len(verts)):
                    if verts[j] in adj[verts[i]]:
                        is_indep = False
                        break
                if not is_indep:
                    break
            if is_indep:
                best = max(best, len(verts))
        alpha_max = max(alpha_max, best)

        # Max clique size (greedy)
        max_clique = 0
        for i in range(nc):
            clique = {i}
            for j in range(nc):
                if j == i: continue
                if all(j in adj[c] for c in clique):
                    clique.add(j)
            max_clique = max(max_clique, len(clique))
        if nc > 0:
            clique_sizes.append(max_clique)

    print(f"  Max independence number seen: {alpha_max}")
    if clique_sizes:
        print(f"  Max clique size range: {min(clique_sizes)}-{max(clique_sizes)}")
        print(f"  Avg max clique: {sum(clique_sizes)/len(clique_sizes):.1f}")
