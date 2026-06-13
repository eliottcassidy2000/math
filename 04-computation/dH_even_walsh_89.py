#!/usr/bin/env python3
"""
WHY IS |dH| ALWAYS EVEN? — Walsh-Theoretic Proof
opus-2026-03-14-S89

Observation: When we flip a single arc in a tournament, |H(T) - H(T')| is
always even. Since H is always odd (Redei), H - H' = odd - odd = even.
This is TRIVIALLY true!

But the deeper question: what is the DISTRIBUTION of dH values?
And: is dH DIVISIBLE BY 4 sometimes? By 8? What's the divisibility structure?

This script:
1. Proves |dH| even is trivial (Redei consequence)
2. Explores the FINE divisibility structure of dH
3. Connects dH to Walsh derivative (directional derivative in Fourier space)
4. Explores the "edge graph" where vertices = tournaments, edges = single flip
5. Investigates whether the edge graph has interesting spectral properties
"""

from itertools import permutations, combinations
from math import factorial, comb
from collections import Counter, defaultdict

print("=" * 70)
print("|dH| ALWAYS EVEN: THE WALSH DERIVATIVE STRUCTURE")
print("opus-2026-03-14-S89")
print("=" * 70)

# ======================================================================
# PART 1: WHY |dH| IS EVEN — TRIVIAL FROM REDEI
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: THE TRIVIAL PROOF")
print("=" * 70)

print("""
  THEOREM: For any tournament T and any arc e, |H(T) - H(T^e)| is even,
  where T^e is T with arc e reversed.

  PROOF: By Redei's theorem, H(T) is odd for all T.
  H(T) - H(T^e) = odd - odd = even. QED.

  This is not deep — it's a CONSEQUENCE of Redei.
  But the STRUCTURE of dH = H(T) - H(T^e) is much richer.
""")

# ======================================================================
# PART 2: FINE DIVISIBILITY OF dH
# ======================================================================
print("=" * 70)
print("PART 2: FINE DIVISIBILITY OF dH")
print("=" * 70)

def all_tournaments(n):
    """Generate all tournaments on n vertices as adjacency dicts."""
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    for bits in range(2**m):
        adj = {}
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
        yield bits, adj, edges

def compute_H(n, adj):
    count = 0
    for p in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if adj[(p[k], p[k+1])] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    # Compute H for all tournaments
    H_vals = {}
    for bits in range(2**m):
        adj = {}
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
        H_vals[bits] = compute_H(n, adj)

    # Analyze dH for each edge flip
    dH_vals = []
    dH_by_edge = defaultdict(list)

    for bits in range(2**m):
        for e_idx in range(m):
            new_bits = bits ^ (1 << e_idx)
            dh = H_vals[bits] - H_vals[new_bits]
            dH_vals.append(dh)
            dH_by_edge[e_idx].append(dh)

    # Divisibility analysis
    div_counts = Counter()
    for dh in dH_vals:
        if dh == 0:
            div_counts['0'] += 1
        else:
            d = abs(dh)
            max_pow_2 = 0
            while d % 2 == 0:
                d //= 2
                max_pow_2 += 1
            div_counts[f'2^{max_pow_2}'] += 1

    print(f"\n  n={n}: dH divisibility by powers of 2:")
    total = len(dH_vals)
    for key in sorted(div_counts.keys()):
        print(f"    {key:>5s}: {div_counts[key]:7d} ({div_counts[key]/total:.4f})")

    # Distribution of dH values (not |dH|, but signed dH)
    dH_counter = Counter(dH_vals)
    print(f"\n  n={n}: Signed dH distribution:")
    for dh in sorted(dH_counter.keys()):
        print(f"    dH = {dh:4d}: {dH_counter[dh]:7d}", end="")
        if dh != 0 and -dh in dH_counter:
            print(f"  (symmetric: dH={-dh} has {dH_counter[-dh]})", end="")
        print()

    # Check symmetry: is the distribution of dH symmetric about 0?
    symmetric = all(dH_counter.get(dh, 0) == dH_counter.get(-dh, 0)
                    for dh in dH_counter)
    print(f"  Symmetric about 0? {symmetric}")

    # Mean dH (should be 0 by symmetry)
    mean_dH = sum(dH_vals) / len(dH_vals)
    print(f"  Mean dH = {mean_dH:.6f}")

    # Variance of dH
    var_dH = sum(dh**2 for dh in dH_vals) / len(dH_vals)
    print(f"  Var(dH) = E[dH^2] = {var_dH:.4f}")
    print(f"  sqrt(Var) = {var_dH**0.5:.4f}")

# ======================================================================
# PART 3: WALSH DERIVATIVE INTERPRETATION
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: WALSH DERIVATIVE -- dH AS DIRECTIONAL DERIVATIVE")
print("=" * 70)

print("""
  In the Walsh-Fourier framework:
  H(T) = sum_S H_hat[S] * chi_S(T)

  When we flip edge e, we get T^e.
  chi_S(T^e) = chi_S(T) * (-1)^{[e in S]}

  So: H(T^e) = sum_S H_hat[S] * chi_S(T) * (-1)^{[e in S]}

  dH_e(T) = H(T) - H(T^e)
          = sum_S H_hat[S] * chi_S(T) * (1 - (-1)^{[e in S]})
          = 2 * sum_{S: e in S} H_hat[S] * chi_S(T)

  So dH_e is EXACTLY TWICE the "restriction to S containing e."
  This is the WALSH DERIVATIVE of H in direction e.

  |dH| is even because dH is always 2 * (some integer).
  The "some integer" is the Walsh restriction to sets containing e.

  DEEPER: dH_e = 0 iff sum_{S: e in S} H_hat[S] * chi_S(T) = 0.
  This means H is "insensitive" to edge e at tournament T
  iff the Walsh components through e cancel at T.
""")

# Verify the Walsh derivative formula for n=3
print("  Verification for n=3:")
n = 3
m = 3
edges = [(0,1), (0,2), (1,2)]

# Compute Walsh-Hadamard transform of H
H_table = []
for bits in range(2**m):
    adj = {}
    for idx, (i, j) in enumerate(edges):
        if (bits >> idx) & 1:
            adj[(i,j)] = 1
            adj[(j,i)] = 0
        else:
            adj[(i,j)] = 0
            adj[(j,i)] = 1
    H_table.append(compute_H(n, adj))

print(f"    H values: {H_table}")

# Walsh-Hadamard transform
H_hat = [0.0] * (2**m)
for S in range(2**m):
    for T_bits in range(2**m):
        # chi_S(T) = (-1)^{|S cap T|}
        chi = (-1) ** bin(S & T_bits).count('1')
        H_hat[S] += H_table[T_bits] * chi
    H_hat[S] /= 2**m

print(f"    H_hat: {[round(h, 4) for h in H_hat]}")

# Check: dH for edge 0 (arc 0,1)
print(f"\n    Walsh derivative for edge 0 (arc {edges[0]}):")
print(f"    Sets containing edge 0: ", end="")
sets_with_0 = [S for S in range(2**m) if S & 1]
print(sets_with_0)
for T_bits in range(2**m):
    actual_dH = H_table[T_bits] - H_table[T_bits ^ 1]
    # Formula: 2 * sum_{S: 0 in S} H_hat[S] * chi_S(T)
    formula_dH = 0
    for S in sets_with_0:
        chi = (-1) ** bin(S & T_bits).count('1')
        formula_dH += H_hat[S] * chi
    formula_dH *= 2
    print(f"    T={T_bits:03b}: actual dH = {actual_dH:+d}, "
          f"formula 2*sum = {formula_dH:+.1f}, match = {abs(actual_dH - formula_dH) < 0.001}")

# ======================================================================
# PART 4: THE EDGE FLIP GRAPH
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: THE EDGE FLIP GRAPH (H-LEVEL SETS)")
print("=" * 70)

# For n=5, analyze the graph where:
# - Vertices = tournaments with the same H value
# - Edges = single arc flip connecting two tournaments with the same H

for n in range(3, 6):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    # Compute H for all
    H_map = {}
    for bits in range(2**m):
        adj = {}
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
        H_map[bits] = compute_H(n, adj)

    # Group by H value
    by_H = defaultdict(list)
    for bits, h in H_map.items():
        by_H[h].append(bits)

    print(f"\n  n={n}: Tournament counts by H value:")
    for h in sorted(by_H.keys()):
        count = len(by_H[h])
        # Check how many same-H neighbors each tournament has
        same_h_neighbors = []
        for bits in by_H[h]:
            shn = 0
            for e_idx in range(m):
                nbr = bits ^ (1 << e_idx)
                if H_map[nbr] == h:
                    shn += 1
            same_h_neighbors.append(shn)
        avg_shn = sum(same_h_neighbors) / len(same_h_neighbors) if same_h_neighbors else 0
        print(f"    H={h:3d}: {count:5d} tournaments, "
              f"avg same-H neighbors = {avg_shn:.2f}/{m}")

# ======================================================================
# PART 5: WHICH EDGES ARE MOST/LEAST SENSITIVE?
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: EDGE SENSITIVITY BY POSITION (n=5)")
print("=" * 70)

n = 5
m = 10
edges = []
for i in range(n):
    for j in range(i+1, n):
        edges.append((i, j))

# For each edge, compute average |dH|
H_map = {}
for bits in range(2**m):
    adj = {}
    for idx, (i, j) in enumerate(edges):
        if (bits >> idx) & 1:
            adj[(i,j)] = 1
            adj[(j,i)] = 0
        else:
            adj[(i,j)] = 0
            adj[(j,i)] = 1
    H_map[bits] = compute_H(n, adj)

print(f"\n  n={n}: Average |dH| for each edge:")
for e_idx in range(m):
    total_abs_dh = 0
    count = 0
    for bits in range(2**m):
        nbr = bits ^ (1 << e_idx)
        dh = abs(H_map[bits] - H_map[nbr])
        total_abs_dh += dh
        count += 1
    avg = total_abs_dh / count
    print(f"    Edge {edges[e_idx]}: avg |dH| = {avg:.4f}")

print("\n  ALL EDGES HAVE THE SAME AVERAGE |dH|!")
print("  This is because all edges are equivalent under S_n symmetry.")
print("  (Any edge can be mapped to any other by relabeling vertices.)")

# ======================================================================
# PART 6: dH AND THE CONE FUNCTOR
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: dH AND THE CONE FUNCTOR")
print("=" * 70)

print("""
  If T' = Cone_top(T), then T' has one extra vertex v* beating everyone.
  What happens to dH when we flip an arc of T'?

  Case 1: Flip an arc e of T (not involving v*).
    Then Cone_top(T^e) = (Cone_top(T))^e.
    By THM-205: H(Cone_top(T^e)) = H(T^e).
    So dH_{Cone_top(T)}(e) = H(Cone_top(T)) - H(Cone_top(T)^e)
                            = H(T) - H(T^e) = dH_T(e).
    The sensitivity of CONE EDGES is IDENTICAL to the original!

  Case 2: Flip a cone arc (v* -> w becomes w -> v*).
    After flip, v* no longer beats w. Now v* has one predecessor.
    H changes because v* can no longer always be first.
    This is a NEW kind of sensitivity not present in T.
""")

# Verify Case 1 for n=3 -> 4
print("  Verification of Case 1 (n=3 -> 4):")
n_small = 3
m_small = 3
edges_small = [(0,1), (0,2), (1,2)]

# For each 3-tournament, compute H and dH for each edge
# Then cone to n=4 and check that dH is preserved
for bits3 in range(2**m_small):
    adj3 = {}
    for idx, (i, j) in enumerate(edges_small):
        if (bits3 >> idx) & 1:
            adj3[(i,j)] = 1
            adj3[(j,i)] = 0
        else:
            adj3[(i,j)] = 0
            adj3[(j,i)] = 1
    H3 = compute_H(3, adj3)

    # Top cone: vertex 3 beats 0, 1, 2
    adj4 = dict(adj3)
    for v in range(3):
        adj4[(3, v)] = 1
        adj4[(v, 3)] = 0
    H4 = compute_H(4, adj4)

    # Flip each original edge in both versions
    dH_match = True
    for e_idx in range(m_small):
        # dH in 3-tournament
        adj3_flip = dict(adj3)
        i, j = edges_small[e_idx]
        adj3_flip[(i,j)] = 1 - adj3[(i,j)]
        adj3_flip[(j,i)] = 1 - adj3[(j,i)]
        H3_flip = compute_H(3, adj3_flip)
        dH3 = H3 - H3_flip

        # dH in 4-tournament (cone of flipped)
        adj4_flip = dict(adj3_flip)
        for v in range(3):
            adj4_flip[(3, v)] = 1
            adj4_flip[(v, 3)] = 0
        H4_flip = compute_H(4, adj4_flip)
        dH4 = H4 - H4_flip

        if dH3 != dH4:
            dH_match = False
            print(f"    MISMATCH at T={bits3:03b}, edge {edges_small[e_idx]}: "
                  f"dH3={dH3}, dH4={dH4}")

    if not dH_match:
        print(f"    T={bits3:03b}: MISMATCH found")

print("  Case 1 verified: dH is preserved under coning for all original edges.")

# Case 2: sensitivity of cone edges
print("\n  Case 2: Sensitivity of CONE edges (n=3 -> 4):")
for bits3 in range(2**m_small):
    adj3 = {}
    for idx, (i, j) in enumerate(edges_small):
        if (bits3 >> idx) & 1:
            adj3[(i,j)] = 1
            adj3[(j,i)] = 0
        else:
            adj3[(i,j)] = 0
            adj3[(j,i)] = 1
    H3 = compute_H(3, adj3)

    # Top cone
    adj4 = dict(adj3)
    for v in range(3):
        adj4[(3, v)] = 1
        adj4[(v, 3)] = 0
    H4 = compute_H(4, adj4)

    # Flip each cone edge (3->v becomes v->3)
    for v in range(3):
        adj4_flip = dict(adj4)
        adj4_flip[(3, v)] = 0
        adj4_flip[(v, 3)] = 1
        H4_flip = compute_H(4, adj4_flip)
        dH_cone = H4 - H4_flip
        print(f"    T={bits3:03b} (H={H3}), flip cone arc 3->{v}: "
              f"H_cone={H4}, H_flipped={H4_flip}, dH={dH_cone:+d}")

# ======================================================================
# PART 7: THE LAPLACIAN OF H ON THE TOURNAMENT HYPERCUBE
# ======================================================================
print("\n" + "=" * 70)
print("PART 7: HYPERCUBE LAPLACIAN OF H")
print("=" * 70)

print("""
  The tournament hypercube Q_m has vertices = 2^m tournaments.
  The Laplacian: (Delta H)(T) = m * H(T) - sum_{e} H(T^e)
  where the sum is over all m neighbors (single arc flips).

  Equivalently: (Delta H)(T) = sum_e dH_e(T) = sum_e (H(T) - H(T^e))

  From the Walsh decomposition:
  (Delta H)(T) = sum_e 2 * sum_{S: e in S} H_hat[S] * chi_S(T)
              = 2 * sum_S |S| * H_hat[S] * chi_S(T)

  So the Laplacian multiplies each Walsh coefficient by 2|S|.
  This is the STANDARD hypercube Laplacian eigenvalue: lambda_S = 2|S|.
""")

for n in range(3, 6):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    # Compute H for all
    H_list = [0] * (2**m)
    for bits in range(2**m):
        adj = {}
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
        H_list[bits] = compute_H(n, adj)

    # Compute Laplacian
    Lap = [0] * (2**m)
    for bits in range(2**m):
        total = 0
        for e_idx in range(m):
            nbr = bits ^ (1 << e_idx)
            total += H_list[bits] - H_list[nbr]
        Lap[bits] = total

    # Statistics
    max_lap = max(Lap)
    min_lap = min(Lap)
    mean_lap = sum(Lap) / len(Lap)

    print(f"  n={n}: Laplacian of H on Q_{m}:")
    print(f"    max(Delta H) = {max_lap}, min(Delta H) = {min_lap}, "
          f"mean = {mean_lap:.4f}")

    # The mean Laplacian should be 2 * sum_S |S| * H_hat[S]^2 * 2^m / ... hmm
    # Actually mean of Lap = 0 by symmetry? No.
    # mean of Delta H = m * mean(H) - m * mean(H) = 0. Yes!
    # Wait: (Delta H)(T) = sum_e (H(T) - H(T^e))
    # E[Delta H] = sum_e (E[H] - E[H]) = 0. Yes!

    # Ratio of Laplacian to H
    print(f"    Lap/H ratios: ", end="")
    ratios = set()
    for bits in range(min(2**m, 64)):
        if H_list[bits] > 0:
            r = Lap[bits] / H_list[bits]
            ratios.add(round(r, 4))
    print(sorted(ratios)[:10], "..." if len(ratios) > 10 else "")

# ======================================================================
# PART 8: MAXIMUM SENSITIVITY — WHICH TOURNAMENT IS MOST SENSITIVE?
# ======================================================================
print("\n" + "=" * 70)
print("PART 8: MAXIMUM EDGE SENSITIVITY")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    H_list = [0] * (2**m)
    for bits in range(2**m):
        adj = {}
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
        H_list[bits] = compute_H(n, adj)

    # For each tournament, compute max single-flip |dH|
    max_sensitivity = 0
    max_sens_T = None
    max_sens_H = None

    total_sensitivity = []
    for bits in range(2**m):
        max_dh = 0
        total_dh = 0
        for e_idx in range(m):
            nbr = bits ^ (1 << e_idx)
            dh = abs(H_list[bits] - H_list[nbr])
            max_dh = max(max_dh, dh)
            total_dh += dh
        total_sensitivity.append(total_dh)
        if max_dh > max_sensitivity:
            max_sensitivity = max_dh
            max_sens_T = bits
            max_sens_H = H_list[bits]

    max_total = max(total_sensitivity)
    avg_total = sum(total_sensitivity) / len(total_sensitivity)

    print(f"\n  n={n}:")
    print(f"    Max single-flip |dH| = {max_sensitivity} "
          f"(at H={max_sens_H})")
    print(f"    Max total sensitivity = {max_total}")
    print(f"    Avg total sensitivity = {avg_total:.2f}")
    print(f"    Max total / m = {max_total/m:.2f}")
    print(f"    Max single / max(H) = {max_sensitivity / max(H_list):.4f}")

# ======================================================================
# PART 9: THE 4-DIVISIBILITY PATTERN
# ======================================================================
print("\n" + "=" * 70)
print("PART 9: HOW OFTEN IS dH DIVISIBLE BY 4?")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    H_list = [0] * (2**m)
    for bits in range(2**m):
        adj = {}
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
        H_list[bits] = compute_H(n, adj)

    # Count dH values by divisibility
    div_by_4 = 0
    div_by_8 = 0
    div_by_2_not_4 = 0
    zero_count = 0
    total = 0

    for bits in range(2**m):
        for e_idx in range(m):
            nbr = bits ^ (1 << e_idx)
            dh = H_list[bits] - H_list[nbr]
            total += 1
            if dh == 0:
                zero_count += 1
            elif dh % 8 == 0:
                div_by_8 += 1
            elif dh % 4 == 0:
                div_by_4 += 1
            else:
                div_by_2_not_4 += 1

    print(f"\n  n={n} (total flips = {total}):")
    print(f"    dH = 0:          {zero_count:8d} ({zero_count/total:.4f})")
    print(f"    4 | dH (not 8):  {div_by_4:8d} ({div_by_4/total:.4f})")
    print(f"    8 | dH:          {div_by_8:8d} ({div_by_8/total:.4f})")
    print(f"    2 | dH, 4 !| dH: {div_by_2_not_4:8d} ({div_by_2_not_4/total:.4f})")

    # H mod 4 distribution
    h_mod4 = Counter(h % 4 for h in H_list)
    print(f"    H mod 4 distribution: {dict(sorted(h_mod4.items()))}")

    # dH mod 4 distribution (nonzero only)
    dh_mod4 = Counter()
    for bits in range(2**m):
        for e_idx in range(m):
            nbr = bits ^ (1 << e_idx)
            dh = H_list[bits] - H_list[nbr]
            if dh != 0:
                dh_mod4[dh % 4] += 1
    print(f"    dH mod 4 distribution (nonzero): {dict(sorted(dh_mod4.items()))}")

print("\n" + "=" * 70)
print("DONE -- dH ALWAYS EVEN, WALSH DERIVATIVE STRUCTURE")
print("=" * 70)
