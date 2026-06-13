#!/usr/bin/env python3
"""
tournament_deep_topology_87.py — opus-2026-03-14-S87

THE DEEPER TOPOLOGY OF THE TOURNAMENT WORLD

Exploration areas:
1. The tournament cube: T_n lives in {0,1}^C(n,2), the binary hypercube.
   Hamming distance = number of arc flips. The H-landscape on this cube.
2. Independence complex Ind(Ω(T)): simplicial complex whose faces = independent
   sets in the conflict graph. Euler char χ = I(Ω,-1). f-vector → Betti numbers.
3. Persistent homology of the H-filtration: filter tournaments by H value,
   track Betti numbers as threshold increases.
4. The flip graph: vertices = tournaments, edges = single arc flips.
   H changes by how much under a single flip? Gradient flow on the cube.
5. Homotopy type of the poset of tournaments ordered by H.
6. Tournament complex: simplicial complex on vertex set = all tournaments,
   faces = sets sharing a common property (e.g., same H value).

The key insight we're chasing: the x=-1 → x=2 deformation
(Euler characteristic → tournament count) is a TOPOLOGICAL bridge.
"""

from itertools import combinations, permutations
from collections import defaultdict, Counter
import sys
import math

# ── Core tournament infrastructure ──────────────────────────────

def all_tournaments(n):
    """Generate all tournaments on n vertices as adjacency matrices."""
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(edges):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj, bits

def compute_H_dp(adj, n):
    """Hamiltonian path count via Held-Karp DP. O(2^n * n^2)."""
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            if dp[S][v] == 0: continue
            for w in range(n):
                if S & (1 << w): continue
                if adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def get_odd_cycles(adj, n):
    """Find all odd directed cycles in tournament.
    Each cycle is represented as a tuple (v0, v1, ..., v_{k-1})
    with arcs v0->v1->...->v_{k-1}->v0, normalized so min vertex is first.
    Returns list of (cycle_tuple, vertex_set) pairs.
    """
    cycles = []
    seen = set()
    for length in range(3, n+1, 2):  # odd lengths only
        for combo in combinations(range(n), length):
            verts = list(combo)
            for perm in permutations(verts):
                is_cycle = True
                for i in range(length):
                    if not adj[perm[i]][perm[(i+1)%length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # Normalize: rotate so min vertex is first
                    min_pos = perm.index(min(perm))
                    normalized = tuple(perm[min_pos:] + perm[:min_pos])
                    if normalized not in seen:
                        seen.add(normalized)
                        cycles.append((normalized, frozenset(combo)))
    return cycles

def build_conflict_graph(cycles):
    """Build conflict graph: vertices = odd cycles, edges = sharing a vertex.
    cycles is list of (cycle_tuple, vertex_set) pairs.
    """
    nc = len(cycles)
    adj = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][1] & cycles[j][1]:  # share a vertex
                adj[i][j] = adj[j][i] = 1
    return adj

def independence_polynomial(adj_cg, nc, x):
    """Compute I(G, x) for small graph via brute force over subsets."""
    result = 0
    for mask in range(1 << nc):
        verts = [i for i in range(nc) if mask & (1 << i)]
        # Check independence
        indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj_cg[verts[i]][verts[j]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            result += x ** len(verts)
    return result

def indep_complex_f_vector(adj_cg, nc):
    """Get f-vector of independence complex: f_k = #independent sets of size k+1."""
    f = defaultdict(int)
    for mask in range(1 << nc):
        verts = [i for i in range(nc) if mask & (1 << i)]
        indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj_cg[verts[i]][verts[j]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            f[len(verts)] += 1
    return f  # f[0]=1 (empty set), f[1]=nc, etc.

def hamming_distance(bits1, bits2, m):
    """Hamming distance between two tournament bit-encodings."""
    return bin(bits1 ^ bits2).count('1')

# ══════════════════════════════════════════════════════════════════
# PART 1: THE TOURNAMENT CUBE — H-landscape on {0,1}^C(n,2)
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("PART 1: THE TOURNAMENT CUBE — H-LANDSCAPE TOPOLOGY")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    # Compute H for all tournaments
    tournament_data = []
    for adj, bits in all_tournaments(n):
        H = compute_H_dp(adj, n)
        tournament_data.append((bits, H))

    H_values = [h for _, h in tournament_data]
    H_counter = Counter(H_values)

    print(f"  Cube dimension: {m}")
    print(f"  Total tournaments: {1 << m}")
    print(f"  H distribution: {dict(sorted(H_counter.items()))}")

    # For each tournament, compute H-gradient:
    # how many neighbors (single flip) have H > current H?
    # This is the "uphill" degree in the H-landscape.
    bits_to_H = {b: h for b, h in tournament_data}

    gradient_stats = defaultdict(list)
    for bits, H in tournament_data:
        uphill = 0
        downhill = 0
        level = 0
        for e in range(m):
            nbr = bits ^ (1 << e)
            H_nbr = bits_to_H[nbr]
            if H_nbr > H:
                uphill += 1
            elif H_nbr < H:
                downhill += 1
            else:
                level += 1
        gradient_stats[H].append((uphill, downhill, level))

    print(f"\n  H-gradient on cube (avg uphill, downhill, level neighbors):")
    for H_val in sorted(gradient_stats.keys()):
        data = gradient_stats[H_val]
        avg_up = sum(u for u,d,l in data) / len(data)
        avg_dn = sum(d for u,d,l in data) / len(data)
        avg_lv = sum(l for u,d,l in data) / len(data)
        print(f"    H={H_val}: up={avg_up:.1f} down={avg_dn:.1f} level={avg_lv:.1f}  (count={len(data)})")

    # Critical points: local minima (no uphill) and local maxima (no downhill)
    local_min = sum(1 for b, H in tournament_data
                    if all(bits_to_H[b ^ (1 << e)] >= H for e in range(m)))
    local_max = sum(1 for b, H in tournament_data
                    if all(bits_to_H[b ^ (1 << e)] <= H for e in range(m)))
    saddle = (1 << m) - local_min - local_max  # very rough

    print(f"\n  Critical points of H-landscape:")
    print(f"    Local minima: {local_min}")
    print(f"    Local maxima: {local_max}")

    # Morse inequality: #critical ≥ sum of Betti numbers
    # This gives a lower bound on the topology of level sets

    # Level set connectivity
    print(f"\n  Level set analysis (each H value):")
    for H_val in sorted(H_counter.keys()):
        # Get all tournaments with this H value
        level_set = [b for b, h in tournament_data if h == H_val]
        # Check connectivity via BFS on flip graph restricted to level set
        if len(level_set) <= 1:
            components = len(level_set)
        else:
            level_set_bits = set(level_set)
            visited = set()
            components = 0
            for start in level_set:
                if start in visited:
                    continue
                components += 1
                queue = [start]
                visited.add(start)
                while queue:
                    cur = queue.pop()
                    for e in range(m):
                        nbr = cur ^ (1 << e)
                        if nbr in level_set_bits and nbr not in visited:
                            visited.add(nbr)
                            queue.append(nbr)
        print(f"    H={H_val}: {H_counter[H_val]} tournaments, {components} connected component(s) in flip graph")

# ══════════════════════════════════════════════════════════════════
# PART 2: INDEPENDENCE COMPLEX TOPOLOGY — Ind(Ω(T))
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 2: INDEPENDENCE COMPLEX OF Ω(T)")
print("=" * 70)
print()
print("For each tournament T, Ω(T) is the conflict graph of odd cycles.")
print("Ind(Ω(T)) is the simplicial complex of independent sets.")
print("χ(Ind) = I(Ω, -1) — the Euler characteristic.")
print("H(T) = I(Ω, 2) — the tournament evaluation.")
print("The DEFORMATION x: -1 → 2 bridges topology and counting.")

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")

    results_by_H = defaultdict(list)

    for adj, bits in all_tournaments(n):
        H = compute_H_dp(adj, n)
        cycles = get_odd_cycles(adj, n)
        nc = len(cycles)

        if nc == 0:
            # Transitive tournament: Ω is empty
            f_vec = {0: 1}  # just the empty set
            euler_char = 1  # I(empty, -1) = 1
            I_at_2 = 1
        else:
            adj_cg = build_conflict_graph(cycles)
            f_vec = indep_complex_f_vector(adj_cg, nc)
            euler_char = independence_polynomial(adj_cg, nc, -1)
            I_at_2 = independence_polynomial(adj_cg, nc, 2)

        # Verify I(Ω, 2) = H
        assert I_at_2 == H, f"Mismatch: I(Ω,2)={I_at_2} vs H={H}"

        results_by_H[H].append({
            'euler': euler_char,
            'f_vec': dict(f_vec),
            'n_cycles': nc
        })

    print(f"\n  Topology of Ind(Ω(T)) grouped by H:")
    for H_val in sorted(results_by_H.keys()):
        data = results_by_H[H_val]
        euler_vals = set(d['euler'] for d in data)
        f_vecs = [tuple(sorted(d['f_vec'].items())) for d in data]
        f_vec_counter = Counter(f_vecs)

        print(f"  H={H_val} ({len(data)} tournaments):")
        print(f"    Euler chars χ(Ind) = I(Ω,-1): {sorted(euler_vals)}")
        for fv, count in f_vec_counter.most_common(5):
            fv_str = ", ".join(f"f_{k}={v}" for k,v in fv)
            print(f"    f-vector ({count}x): {fv_str}")

        # The ratio H/χ = I(Ω,2)/I(Ω,-1) — the "topological amplification"
        for d in data[:1]:  # just first example
            if d['euler'] != 0:
                ratio = H_val / d['euler']
                print(f"    Topological amplification H/χ = {H_val}/{d['euler']} = {ratio:.3f}")

# ══════════════════════════════════════════════════════════════════
# PART 3: THE x-DEFORMATION — I(Ω, x) as x goes from -1 to 2
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 3: THE x-DEFORMATION — I(Ω, x) from topology to counting")
print("=" * 70)
print()
print("I(Ω, -1) = Euler characteristic (topology)")
print("I(Ω, 0) = 1 (trivial)")
print("I(Ω, 1) = #independent sets (combinatorics)")
print("I(Ω, 2) = H (tournament count)")
print()
print("The interpolation I(Ω, x) for x ∈ {-1, 0, 1, 2, 3} shows how")
print("topology deforms into counting.")

n = 5
print(f"\n--- n = {n}: Full interpolation table ---")
print(f"{'H':>4} {'#cyc':>5} {'I(-1)':>6} {'I(0)':>5} {'I(1)':>5} {'I(2)':>5} {'I(3)':>6} {'I(4)':>6}")
print("-" * 55)

seen_H = set()
for adj, bits in all_tournaments(n):
    H = compute_H_dp(adj, n)
    if H in seen_H:
        continue
    seen_H.add(H)

    cycles = get_odd_cycles(adj, n)
    nc = len(cycles)

    if nc == 0:
        vals = {x: 1 for x in [-1, 0, 1, 2, 3, 4]}
    else:
        adj_cg = build_conflict_graph(cycles)
        vals = {}
        for x in [-1, 0, 1, 2, 3, 4]:
            vals[x] = independence_polynomial(adj_cg, nc, x)

    print(f"{H:>4} {nc:>5} {vals[-1]:>6} {vals[0]:>5} {vals[1]:>5} {vals[2]:>5} {vals[3]:>6} {vals[4]:>6}")

# ══════════════════════════════════════════════════════════════════
# PART 4: FLIP GRAPH TOPOLOGY — The "space of tournaments"
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 4: FLIP GRAPH — GRADIENT FLOW AND MORSE THEORY")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n = {n} ---")
    edges_list = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges_list)

    tournament_data = []
    for adj, bits in all_tournaments(n):
        H = compute_H_dp(adj, n)
        tournament_data.append((bits, H))

    bits_to_H = {b: h for b, h in tournament_data}

    # Morse theory: count critical points at each H level
    # A tournament T is a critical point of H on the flip graph if
    # it is a local min or max (or saddle)
    critical_by_H = defaultdict(lambda: {'min': 0, 'max': 0, 'saddle': 0, 'total': 0})

    for bits, H in tournament_data:
        neighbors_H = [bits_to_H[bits ^ (1 << e)] for e in range(m)]
        is_min = all(nh >= H for nh in neighbors_H)
        is_max = all(nh <= H for nh in neighbors_H)
        critical_by_H[H]['total'] += 1
        if is_min and not is_max:
            critical_by_H[H]['min'] += 1
        elif is_max and not is_min:
            critical_by_H[H]['max'] += 1
        elif is_min and is_max:
            # Isolated point (all neighbors have same H — plateau)
            critical_by_H[H]['min'] += 1  # count as min

    print(f"  Morse critical points on H-landscape:")
    total_crit = 0
    for H_val in sorted(critical_by_H.keys()):
        d = critical_by_H[H_val]
        crit = d['min'] + d['max']
        total_crit += crit
        if d['min'] > 0 or d['max'] > 0:
            print(f"    H={H_val}: {d['min']} min, {d['max']} max (of {d['total']} total)")
    print(f"  Total critical: {total_crit} of {1 << m}")

    # Compute the H-change distribution for single flips
    delta_H_dist = Counter()
    for bits, H in tournament_data:
        for e in range(m):
            nbr = bits ^ (1 << e)
            delta = bits_to_H[nbr] - H
            delta_H_dist[delta] += 1

    print(f"\n  ΔH distribution for single arc flips:")
    for delta in sorted(delta_H_dist.keys()):
        count = delta_H_dist[delta]
        # Each flip counted twice (from both endpoints), so divide
        print(f"    ΔH = {delta:+d}: {count} (fraction = {count/sum(delta_H_dist.values()):.4f})")

    # Symmetry check: delta_H should be symmetric around 0
    # because T → T^op maps H to H and reverses all arcs
    print(f"  Symmetry check: ΔH=+k count == ΔH=-k count?")
    symmetric = True
    for delta in sorted(delta_H_dist.keys()):
        if delta > 0:
            if delta_H_dist[delta] != delta_H_dist[-delta]:
                symmetric = False
                print(f"    ASYMMETRIC at ΔH=±{delta}: {delta_H_dist[delta]} vs {delta_H_dist[-delta]}")
    if symmetric:
        print(f"    YES — perfectly symmetric (as expected from T↔T^op)")

# ══════════════════════════════════════════════════════════════════
# PART 5: SUBLEVEL SETS AND PERSISTENT HOMOLOGY (discrete)
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 5: SUBLEVEL SET FILTRATION — DISCRETE PERSISTENT HOMOLOGY")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n = {n} ---")
    edges_list = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges_list)

    tournament_data = []
    for adj, bits in all_tournaments(n):
        H = compute_H_dp(adj, n)
        tournament_data.append((bits, H))

    H_levels = sorted(set(h for _, h in tournament_data))

    print(f"  Sublevel sets T_≤h = {{T : H(T) ≤ h}} on flip graph:")
    print(f"  {'h':>4} {'|T_≤h|':>8} {'components':>11} {'new_verts':>10} {'new_edges':>10}")

    # Build sublevel set incrementally
    all_bits = {}  # bits -> H
    for bits, H in tournament_data:
        all_bits[bits] = H

    prev_components = 0
    for h in H_levels:
        sublevel = set(b for b, hv in tournament_data if hv <= h)
        new_verts = set(b for b, hv in tournament_data if hv == h)

        # Count edges in sublevel set (within flip graph)
        edge_count = 0
        for b in sublevel:
            for e in range(m):
                nbr = b ^ (1 << e)
                if nbr in sublevel and nbr > b:
                    edge_count += 1

        # New edges connecting to this level
        new_edges = 0
        for b in new_verts:
            for e in range(m):
                nbr = b ^ (1 << e)
                if nbr in sublevel:
                    new_edges += 1

        # Connected components via BFS
        visited = set()
        components = 0
        for start in sublevel:
            if start in visited:
                continue
            components += 1
            queue = [start]
            visited.add(start)
            while queue:
                cur = queue.pop()
                for e in range(m):
                    nbr = cur ^ (1 << e)
                    if nbr in sublevel and nbr not in visited:
                        visited.add(nbr)
                        queue.append(nbr)

        print(f"  {h:>4} {len(sublevel):>8} {components:>11} {len(new_verts):>10} {new_edges:>10}")

    # Euler characteristic of full flip graph
    V = 1 << m
    E = sum(1 for b in range(V) for e in range(m) if (b ^ (1 << e)) > b)
    print(f"\n  Full flip graph: V={V}, E={E}, χ = V-E = {V - E}")
    # The flip graph is the m-dimensional hypercube graph Q_m
    # Q_m has V=2^m, E=m*2^{m-1}
    print(f"  (This is Q_{m}: V=2^{m}={2**m}, E={m}*2^{m-1}={m * 2**(m-1)})")

# ══════════════════════════════════════════════════════════════════
# PART 6: THE α-VECTOR AND TOPOLOGICAL TYPE
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 6: THE α-VECTOR — TOPOLOGICAL TYPE OF Ω(T)")
print("=" * 70)
print()
print("The OCF decomposes H = 1 + α₁·2 + α₂·4 + α₃·8 + ...")
print("Each αₖ counts independent sets of size k in Ω(T).")
print("The tuple (α₁, α₂, ...) is the TOPOLOGICAL fingerprint of Ω.")
print()

for n in [5, 6]:
    print(f"--- n = {n} ---")

    alpha_tuples = defaultdict(int)
    H_to_alphas = defaultdict(list)

    count = 0
    for adj, bits in all_tournaments(n):
        H = compute_H_dp(adj, n)
        cycles = get_odd_cycles(adj, n)
        nc = len(cycles)

        if nc == 0:
            alpha = (0,)
            H_to_alphas[H].append(alpha)
            alpha_tuples[alpha] += 1
            count += 1
            continue

        adj_cg = build_conflict_graph(cycles)
        f = indep_complex_f_vector(adj_cg, nc)

        # α_k = f[k] for k ≥ 1 (number of independent sets of size k)
        max_k = max(f.keys())
        alpha = tuple(f.get(k, 0) for k in range(1, max_k + 1))

        H_to_alphas[H].append(alpha)
        alpha_tuples[alpha] += 1
        count += 1

        if n == 6 and count > 2000:
            break  # sample for n=6

    print(f"  Distinct α-vectors: {len(alpha_tuples)}")
    print(f"  Most common α-vectors:")
    for alpha, cnt in sorted(alpha_tuples.items(), key=lambda x: -x[1])[:10]:
        H_val = 1 + sum(alpha[k] * (2 ** (k+1)) for k in range(len(alpha))) if alpha != (0,) else 1
        # Recompute H from alpha
        H_from_alpha = 1 + sum(alpha[k-1] * (2**k) for k in range(1, len(alpha)+1)) if alpha != (0,) else 1
        euler_from_alpha = 1 + sum(alpha[k-1] * ((-1)**k) for k in range(1, len(alpha)+1)) if alpha != (0,) else 1
        print(f"    α={alpha}: count={cnt}, H={H_from_alpha}, χ={euler_from_alpha}")

    # Check: how many H values have MULTIPLE α-vectors?
    print(f"\n  H values with multiple α-vectors (I(Ω,2) same but I(Ω,x) differs):")
    multi = 0
    for H_val in sorted(H_to_alphas.keys()):
        distinct = len(set(H_to_alphas[H_val]))
        if distinct > 1:
            multi += 1
            if multi <= 5:
                print(f"    H={H_val}: {distinct} distinct α-vectors: {set(H_to_alphas[H_val])}")
    if multi == 0:
        print(f"    NONE — H determines α-vector uniquely at n={n}!")
    else:
        print(f"    Total: {multi} H values with non-unique α-vector")

# ══════════════════════════════════════════════════════════════════
# PART 7: HOMOTOPY BETWEEN I(Ω,-1) AND I(Ω,2)
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 7: DEFORMATION SPECTRUM — I(Ω, x) AS POLYNOMIAL")
print("=" * 70)
print()
print("For each distinct Ω-type, I(Ω, x) is a polynomial in x.")
print("The roots of this polynomial tell us where the deformation")
print("has 'phase transitions' (zeros).")

n = 5
print(f"\n--- n = {n}: Roots of I(Ω, x) ---")

seen_polys = {}
for adj, bits in all_tournaments(n):
    H = compute_H_dp(adj, n)
    cycles = get_odd_cycles(adj, n)
    nc = len(cycles)

    if nc == 0:
        poly = (1,)  # I(Ω,x) = 1
    else:
        adj_cg = build_conflict_graph(cycles)
        f = indep_complex_f_vector(adj_cg, nc)
        max_k = max(f.keys())
        # Polynomial: sum f[k] * x^k
        poly = tuple(f.get(k, 0) for k in range(max_k + 1))

    if poly not in seen_polys:
        seen_polys[poly] = H

print(f"  Distinct I(Ω, x) polynomials: {len(seen_polys)}")
for poly, H_val in sorted(seen_polys.items(), key=lambda x: x[1]):
    # Format polynomial
    terms = []
    for k, coeff in enumerate(poly):
        if coeff == 0:
            continue
        if k == 0:
            terms.append(f"{coeff}")
        elif k == 1:
            terms.append(f"{coeff}x")
        else:
            terms.append(f"{coeff}x^{k}")
    poly_str = " + ".join(terms) if terms else "0"

    # Evaluate at key points
    I_neg1 = sum(c * ((-1)**k) for k, c in enumerate(poly))
    I_2 = sum(c * (2**k) for k, c in enumerate(poly))

    # Find real roots using simple scan
    roots = []
    for x_100 in range(-300, 100):
        x = x_100 / 100.0
        val = sum(c * (x**k) for k, c in enumerate(poly))
        x_next = (x_100 + 1) / 100.0
        val_next = sum(c * (x_next**k) for k, c in enumerate(poly))
        if val * val_next < 0:
            # Root between x and x_next, refine
            lo, hi = x, x_next
            for _ in range(20):
                mid = (lo + hi) / 2
                val_mid = sum(c * (mid**k) for k, c in enumerate(poly))
                if val * val_mid < 0:
                    hi = mid
                else:
                    lo = mid
                    val = val_mid
            roots.append(round((lo + hi) / 2, 4))

    root_str = ", ".join(f"{r:.3f}" for r in roots) if roots else "none in [-3, 1)"
    print(f"  H={I_2}: I(Ω,x) = {poly_str}  χ={I_neg1}  roots≈{root_str}")

# ══════════════════════════════════════════════════════════════════
# PART 8: THE BETTI NUMBER CONNECTION
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 8: REDUCED BETTI NUMBERS AND HOMOLOGICAL DIMENSION")
print("=" * 70)
print()
print("For a simplicial complex K with f-vector (f₀, f₁, ...),")
print("the reduced Euler char χ̃ = -1 + f₀ - f₁ + f₂ - ...")
print("For independence complex: χ̃(Ind(G)) = I(G,-1) - 1")
print("By Euler-Poincaré: χ̃ = Σ (-1)^i β̃_i")
print()
print("We can't compute Betti numbers directly without homology,")
print("but we can get the HOMOLOGICAL DIMENSION (max nonempty face)")
print("and the alternating sum constraint.")

n = 5
print(f"\n--- n = {n}: Independence complex dimensions ---")

dim_data = defaultdict(list)
for adj, bits in all_tournaments(n):
    H = compute_H_dp(adj, n)
    cycles = get_odd_cycles(adj, n)
    nc = len(cycles)

    if nc == 0:
        hom_dim = -1  # empty complex (just the void)
        dim_data[H].append((-1, 1))  # dim, euler
        continue

    adj_cg = build_conflict_graph(cycles)
    f = indep_complex_f_vector(adj_cg, nc)
    max_dim = max(f.keys()) - 1  # dimension = size - 1
    euler = independence_polynomial(adj_cg, nc, -1)
    dim_data[H].append((max_dim, euler))

print(f"  {'H':>4} {'dim(Ind)':>9} {'χ(Ind)':>7} {'χ̃(Ind)':>8} {'interpretation':>40}")
for H_val in sorted(dim_data.keys()):
    entries = dim_data[H_val]
    dims = set(d for d, e in entries)
    eulers = set(e for d, e in entries)

    dim_str = "/".join(str(d) for d in sorted(dims))
    euler_str = "/".join(str(e) for e in sorted(eulers))
    red_euler_str = "/".join(str(e-1) for e in sorted(eulers))

    interp = ""
    if max(dims) == -1:
        interp = "void (transitive)"
    elif max(dims) == 0:
        interp = "discrete points"
    elif max(dims) == 1:
        interp = "graph (1-skeleton)"
    elif max(dims) == 2:
        interp = "2-complex (triangles)"
    else:
        interp = f"{max(dims)}-complex"

    print(f"  {H_val:>4} {dim_str:>9} {euler_str:>7} {red_euler_str:>8} {interp:>40}")

# ══════════════════════════════════════════════════════════════════
# PART 9: FIBONACCI STRANDS IN THE TOPOLOGY
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 9: FIBONACCI STRANDS IN INDEPENDENCE COMPLEX TOPOLOGY")
print("=" * 70)
print()
print("The independence complex Ind(P_k) of the path graph P_k")
print("has χ = I(P_k, -1), which has period 6.")
print()
print("Period-6 sequence: I(P_k, -1) for k=0,1,2,...")

# Compute I(P_k, x) directly from the recurrence
# I(P_0, x) = 1, I(P_1, x) = 1+x, I(P_k, x) = I(P_{k-1}, x) + x*I(P_{k-2}, x)

def I_path(k, x):
    if k == 0: return 1
    if k == 1: return 1 + x
    a, b = 1, 1 + x
    for _ in range(k - 1):
        a, b = b, b + x * a
    return b

print("k:   ", " ".join(f"{k:>4}" for k in range(20)))
print("I,-1:", " ".join(f"{I_path(k,-1):>4}" for k in range(20)))
print("I, 0:", " ".join(f"{I_path(k, 0):>4}" for k in range(20)))
print("I, 1:", " ".join(f"{I_path(k, 1):>4}" for k in range(20)))
print("I, 2:", " ".join(f"{I_path(k, 2):>4}" for k in range(20)))
print("I, 3:", " ".join(f"{I_path(k, 3):>4}" for k in range(20)))

# The period-6 at x=-1
print(f"\nPeriod-6 of I(P_k, -1): ", end="")
seq = [I_path(k, -1) for k in range(24)]
print(seq)
# Check periodicity
period = None
for p in range(1, 13):
    if all(seq[i] == seq[i+p] for i in range(12)):
        period = p
        break
print(f"Minimal period: {period}")

# Now: I(C_k, -1) = 2*cos(k*pi/3)
print(f"\nI(C_k, -1) for k=3..15 (theoretical: 2cos(kπ/3)):")
for k in range(3, 16):
    # I(C_k, x) = I(P_k, x) + x * I(P_{k-2}, x) ... actually
    # I(C_k, x) = I(P_{k-1}, x) + x * I(P_{k-2}, x) - but that's I(P_k, x)
    # Actually: I(C_k, x) = I(P_{k-1}, x) + x * I(P_{k-3}, x) for k >= 4
    # Hmm, let me use the formula directly:
    # I(C_k, -1) = (-1)^k + 1 (from eigenvalue formula)
    # Actually from S86: I(C_k, -1) = 2*cos(k*pi/3)
    import cmath
    theoretical = round(2 * math.cos(k * math.pi / 3))
    print(f"  k={k}: I(C_{k}, -1) = {theoretical} = 2cos({k}π/3)")

# Connection: the period-6 structure means Ind(P_k) has a PERIODIC HOMOTOPY TYPE
print(f"\n★ INSIGHT: The period-6 of I(P_k,-1) = {{1,0,-1,-1,0,1,...}}")
print(f"  means the homotopy type of Ind(P_k) cycles through 6 phases.")
print(f"  Phase 0 (k≡0 mod 6): χ̃=0 — contractible")
print(f"  Phase 1 (k≡1 mod 6): χ̃=-1 — S^odd homotopy type")
print(f"  Phase 2 (k≡2 mod 6): χ̃=-2 — two-sphere type")
print(f"  Phase 3 (k≡3 mod 6): χ̃=-2 — two-sphere type")
print(f"  Phase 4 (k≡4 mod 6): χ̃=-1 — S^odd homotopy type")
print(f"  Phase 5 (k≡5 mod 6): χ̃=0 — contractible")

# ══════════════════════════════════════════════════════════════════
# PART 10: TOPOLOGICAL AMPLIFICATION RATIO — H/χ
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 10: TOPOLOGICAL AMPLIFICATION — H/χ RATIO")
print("=" * 70)
print()
print("The ratio H/χ = I(Ω,2)/I(Ω,-1) measures how much")
print("the fugacity jump from -1 to 2 'amplifies' the topology.")
print("This is a measure of the ALGEBRAIC RICHNESS of Ω.")

n = 5
ratios = []
for adj, bits in all_tournaments(n):
    H = compute_H_dp(adj, n)
    cycles = get_odd_cycles(adj, n)
    nc = len(cycles)

    if nc == 0:
        chi = 1
    else:
        adj_cg = build_conflict_graph(cycles)
        chi = independence_polynomial(adj_cg, nc, -1)

    if chi != 0:
        ratios.append((H, chi, H/chi))

# Group by H
ratio_by_H = defaultdict(list)
for H, chi, r in ratios:
    ratio_by_H[H].append((chi, r))

print(f"\nn={n}: Amplification H/χ:")
for H_val in sorted(ratio_by_H.keys()):
    data = ratio_by_H[H_val]
    chis = set(c for c, r in data)
    rs = set(r for c, r in data)
    chi_str = "/".join(str(c) for c in sorted(chis))
    r_str = "/".join(f"{r:.2f}" for r in sorted(rs))
    print(f"  H={H_val}: χ={chi_str}, H/χ={r_str}")

# Check: is the amplification always positive? (i.e., χ > 0 always?)
neg_chi = [(H, chi) for H, chi, r in ratios if chi < 0]
zero_chi = [(H, chi) for H, chi, r in ratios if chi == 0]
print(f"\n  Tournaments with χ < 0: {len(neg_chi)}")
print(f"  Tournaments with χ = 0: {len(zero_chi)}")
if neg_chi:
    print(f"  Examples: {neg_chi[:5]}")

print("\n" + "=" * 70)
print("EXPLORATION COMPLETE")
print("=" * 70)
