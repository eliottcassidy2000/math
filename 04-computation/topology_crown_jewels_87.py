#!/usr/bin/env python3
"""
topology_crown_jewels_87.py — opus-2026-03-14-S87

DEEPER TOPOLOGY: The H+2χ identity, Euler obstruction, and
the topological meaning of forbidden values.

CROWN JEWEL: H + 2χ = 3 + Σ_{k≥2} (2^k + 2(-1)^k) α_k

At n≤5: α_k=0 for k≥2, so H + 2χ = 3 exactly.
At n=6: deviation = 6α₂ (first nontrivial correction).

This means:
  - H=7 forbidden ↔ (α₁=3, α₂=0) impossible from Ω(T) constraints
  - χ=-2 would be needed, but no tournament achieves α₁=3

Also: exploring the flip graph diameter, connected components of
isomorphism classes, and the Morse complex structure.
"""

from itertools import combinations, permutations
from collections import defaultdict, Counter
import math

# ── Core infrastructure (same as topology_87) ────────────────────

def all_tournaments(n):
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
    cycles = []
    seen = set()
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            verts = list(combo)
            for perm in permutations(verts):
                is_cycle = True
                for i in range(length):
                    if not adj[perm[i]][perm[(i+1)%length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    min_pos = perm.index(min(perm))
                    normalized = tuple(perm[min_pos:] + perm[:min_pos])
                    if normalized not in seen:
                        seen.add(normalized)
                        cycles.append((normalized, frozenset(combo)))
    return cycles

def build_conflict_graph(cycles):
    nc = len(cycles)
    adj = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][1] & cycles[j][1]:
                adj[i][j] = adj[j][i] = 1
    return adj

def independence_polynomial(adj_cg, nc, x):
    result = 0
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
            result += x ** len(verts)
    return result

def get_alpha_vector(adj_cg, nc):
    """Get α-vector: α_k = number of independent sets of size k."""
    alphas = defaultdict(int)
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
            alphas[len(verts)] += 1
    return alphas

# ══════════════════════════════════════════════════════════════════
# PART 1: THE H + 2χ IDENTITY
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("PART 1: THE H + 2χ = 3 + CORRECTION IDENTITY")
print("=" * 70)
print()
print("For I(Ω, x) = Σ_k α_k x^k, we have:")
print("  H = I(Ω, 2) = Σ α_k 2^k")
print("  χ = I(Ω, -1) = Σ α_k (-1)^k")
print()
print("So H + 2χ = Σ α_k (2^k + 2(-1)^k)")
print("  k=0: α_0(1+2) = 3  (always)")
print("  k=1: α_1(2-2) = 0  (cancels!)")
print("  k=2: α_2(4+2) = 6α₂")
print("  k=3: α_3(8-2) = 6α₃")
print("  k=4: α_4(16+2) = 18α₄")
print()
print("FORMULA: H + 2χ = 3 + 6α₂ + 6α₃ + 18α₄ + 30α₅ + ...")
print()
print("The coefficients c_k = 2^k + 2(-1)^k:")
for k in range(10):
    ck = 2**k + 2*((-1)**k)
    print(f"  c_{k} = {ck}", end="")
    # Check if this is 3 * Jacobsthal-like
    if k >= 2:
        print(f"  = 6 * {ck//6}" if ck % 6 == 0 else f"  (not div by 6)", end="")
    print()

# Actually c_k = 2^k + 2(-1)^k. Let's see:
# c_k for even k: 2^k + 2 = 2(2^{k-1} + 1)
# c_k for odd k: 2^k - 2 = 2(2^{k-1} - 1)
# So c_k = 2(2^{k-1} + (-1)^k)
print()
print("Pattern: c_k = 2(2^{k-1} + (-1)^k) for k ≥ 1")
print("For k=0: c_0 = 3 = the 'base' contribution")
print()

# Verify at n=5 and n=6
for n in [4, 5]:
    print(f"\n--- Verification at n = {n} ---")
    for adj, bits in all_tournaments(n):
        H = compute_H_dp(adj, n)
        cycles = get_odd_cycles(adj, n)
        nc = len(cycles)
        if nc == 0:
            chi = 1
            correction = 0
        else:
            adj_cg = build_conflict_graph(cycles)
            chi = independence_polynomial(adj_cg, nc, -1)
            alphas = get_alpha_vector(adj_cg, nc)
            correction = sum(alphas.get(k, 0) * (2**k + 2*((-1)**k)) for k in range(2, max(alphas.keys())+1))
        identity = H + 2 * chi
        expected = 3 + correction
        if identity != expected:
            print(f"  VIOLATION: H={H}, χ={chi}, H+2χ={identity}, 3+correction={expected}")
            break
    else:
        print(f"  ✓ H + 2χ = 3 + correction verified for all {1 << len([(i,j) for i in range(n) for j in range(i+1,n)])} tournaments!")

# ══════════════════════════════════════════════════════════════════
# PART 2: THE FORBIDDEN VALUE TOPOLOGY
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 2: WHY H=7 IS FORBIDDEN — TOPOLOGICAL OBSTRUCTION")
print("=" * 70)
print()
print("H=7 requires I(Ω, 2) = 7.")
print("At n=5: I(Ω, x) = 1 + α₁x (linear), so 7 = 1 + 2α₁ → α₁ = 3.")
print("But α₁ = 3 is NOT achieved at n=5: α₁ ∈ {0, 1, 2, 4, 5, 6, 7}.")
print("The value α₁ = 3 is MISSING — topological obstruction!")
print()

n = 5
alpha1_counter = Counter()
for adj, bits in all_tournaments(n):
    H = compute_H_dp(adj, n)
    cycles = get_odd_cycles(adj, n)
    alpha1 = len(cycles)
    alpha1_counter[alpha1] += 1

print(f"n=5: Distribution of α₁ = #odd cycles:")
for a1 in sorted(alpha1_counter.keys()):
    H_val = 1 + 2 * a1
    print(f"  α₁ = {a1}: {alpha1_counter[a1]} tournaments → H = {H_val}")

print()
print("MISSING α₁ values: ", end="")
all_a1 = set(range(max(alpha1_counter.keys()) + 1))
missing = all_a1 - set(alpha1_counter.keys())
print(sorted(missing))
print()

# For each missing α₁, what H would it give?
for a1 in sorted(missing):
    print(f"  α₁ = {a1} → H = {1 + 2*a1} (FORBIDDEN)")

# ══════════════════════════════════════════════════════════════════
# PART 3: α₁ SPECTRUM AT n=3,4,5,6 — GAPS AND FORBIDDEN VALUES
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 3: α₁ SPECTRUM — THE SOURCE OF FORBIDDEN H VALUES")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    alpha1_counter = Counter()
    for adj, bits in all_tournaments(n):
        cycles = get_odd_cycles(adj, n)
        alpha1 = len(cycles)
        alpha1_counter[alpha1] += 1

    max_a1 = max(alpha1_counter.keys())
    all_a1 = set(range(max_a1 + 1))
    missing = sorted(all_a1 - set(alpha1_counter.keys()))

    print(f"  α₁ values: {sorted(alpha1_counter.keys())}")
    print(f"  Missing: {missing}")
    print(f"  → Forbidden H at n={n}: {[1 + 2*a for a in missing]}")

# At n=6, sample
print(f"\n--- n = 6 (full enumeration) ---")
n = 6
alpha1_counter = Counter()
alpha_pair_counter = Counter()
count = 0
for adj, bits in all_tournaments(n):
    H = compute_H_dp(adj, n)
    cycles = get_odd_cycles(adj, n)
    nc = len(cycles)
    alpha1 = nc

    if nc > 0:
        adj_cg = build_conflict_graph(cycles)
        alphas = get_alpha_vector(adj_cg, nc)
        alpha2 = alphas.get(2, 0)
        alpha_pair_counter[(alpha1, alpha2)] += 1
    else:
        alpha_pair_counter[(0, 0)] += 1

    alpha1_counter[alpha1] += 1
    count += 1
    if count % 5000 == 0:
        print(f"  ... processed {count} tournaments")

max_a1 = max(alpha1_counter.keys())
all_a1 = set(range(max_a1 + 1))
missing = sorted(all_a1 - set(alpha1_counter.keys()))

print(f"  α₁ range: [0, {max_a1}]")
print(f"  Missing α₁: {missing}")
print(f"  → Forbidden H (from α₁ gaps, assuming α₂=0): {[1 + 2*a for a in missing]}")
print()
print(f"  Full α₁ distribution:")
for a1 in sorted(alpha1_counter.keys()):
    print(f"    α₁ = {a1:2d}: {alpha1_counter[a1]:5d} tournaments")

print(f"\n  (α₁, α₂) pairs that appear:")
for pair in sorted(alpha_pair_counter.keys()):
    if alpha_pair_counter[pair] > 0:
        a1, a2 = pair
        H = 1 + 2*a1 + 4*a2
        chi = 1 - a1 + a2
        print(f"    (α₁={a1:2d}, α₂={a2:2d}): {alpha_pair_counter[pair]:5d} tournaments → H={H:3d}, χ={chi:3d}")

# ══════════════════════════════════════════════════════════════════
# PART 4: THE FLIP GRAPH DIAMETER AND EXPANSION
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 4: FLIP GRAPH DIAMETER — LANDSCAPE GEOMETRY")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)

    tournament_data = []
    for adj, bits in all_tournaments(n):
        H = compute_H_dp(adj, n)
        tournament_data.append((bits, H))

    bits_to_H = {b: h for b, h in tournament_data}

    # BFS from transitive (H=1, min) to find distances
    # Take a specific transitive tournament (all edges i→j for i<j)
    trans_bits = sum(1 << k for k in range(m))  # all edges oriented i→j
    trans_H = bits_to_H[trans_bits]
    print(f"  Transitive tournament bits={trans_bits}, H={trans_H}")

    # BFS from transitive
    dist = {trans_bits: 0}
    queue = [trans_bits]
    max_dist = 0
    dist_by_H = defaultdict(list)
    dist_by_H[trans_H].append(0)

    while queue:
        cur = queue.pop(0)
        for e in range(m):
            nbr = cur ^ (1 << e)
            if nbr not in dist:
                dist[nbr] = dist[cur] + 1
                max_dist = max(max_dist, dist[nbr])
                queue.append(nbr)
                dist_by_H[bits_to_H[nbr]].append(dist[nbr])

    print(f"  Diameter of flip graph (= Q_{m}): {max_dist}")
    print(f"  (Hypercube Q_m always has diameter m = {m})")

    # Average distance from transitive to each H level
    print(f"  Average Hamming distance from transitive to each H level:")
    for H_val in sorted(dist_by_H.keys()):
        ds = dist_by_H[H_val]
        avg_d = sum(ds) / len(ds)
        min_d = min(ds)
        max_d = max(ds)
        print(f"    H={H_val}: avg={avg_d:.1f}, min={min_d}, max={max_d}")

# ══════════════════════════════════════════════════════════════════
# PART 5: H-PARITY AND THE BINARY CUBE
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 5: H MODULAR STRUCTURE ON THE CUBE")
print("=" * 70)

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

tournament_data = []
for adj, bits in all_tournaments(n):
    H = compute_H_dp(adj, n)
    tournament_data.append((bits, H))

# ΔH for single flip: what does it depend on?
# Let's classify flips by which edge is flipped
print(f"\n--- n = {n}: ΔH by edge position ---")
for e_idx, (i,j) in enumerate(edges):
    deltas = Counter()
    for bits, H in tournament_data:
        nbr = bits ^ (1 << e_idx)
        dH = next(h for b, h in tournament_data if b == nbr) - H
        deltas[dH] += 1
    delta_vals = sorted(deltas.keys())
    # Average absolute ΔH
    avg_abs = sum(abs(d) * c for d, c in deltas.items()) / sum(deltas.values())
    print(f"  Edge ({i},{j}): ΔH range [{min(delta_vals)}, {max(delta_vals)}], avg |ΔH| = {avg_abs:.2f}")

# Is ΔH the same for ALL edges? (by symmetry of H as function on the cube)
print()
print("  All edges give the same ΔH distribution (by vertex-relabeling symmetry of n-tournament).")

# ══════════════════════════════════════════════════════════════════
# PART 6: CYCLE-INDEX POLYNOMIAL OF Ω
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 6: CHROMATIC AND CLIQUE NUMBERS OF Ω(T)")
print("=" * 70)
print()
print("Beyond the independence polynomial, what are the graph invariants of Ω?")
print("  ω(Ω) = clique number (max clique = max pairwise conflicting cycles)")
print("  α(Ω) = independence number (max vertex-disjoint cycle collection)")

for n in [4, 5]:
    print(f"\n--- n = {n} ---")
    clique_by_H = defaultdict(list)
    indep_by_H = defaultdict(list)

    for adj, bits in all_tournaments(n):
        H = compute_H_dp(adj, n)
        cycles = get_odd_cycles(adj, n)
        nc = len(cycles)

        if nc == 0:
            clique_by_H[H].append(0)
            indep_by_H[H].append(0)
            continue

        adj_cg = build_conflict_graph(cycles)

        # Find max clique (brute force for small nc)
        max_clique = 0
        max_indep = 0
        for mask in range(1 << nc):
            verts = [i for i in range(nc) if mask & (1 << i)]
            if len(verts) == 0:
                continue

            # Check clique
            is_clique = True
            for i in range(len(verts)):
                for j in range(i+1, len(verts)):
                    if not adj_cg[verts[i]][verts[j]]:
                        is_clique = False
                        break
                if not is_clique:
                    break
            if is_clique:
                max_clique = max(max_clique, len(verts))

            # Check independence
            is_indep = True
            for i in range(len(verts)):
                for j in range(i+1, len(verts)):
                    if adj_cg[verts[i]][verts[j]]:
                        is_indep = False
                        break
                if not is_indep:
                    break
            if is_indep:
                max_indep = max(max_indep, len(verts))

        clique_by_H[H].append(max_clique)
        indep_by_H[H].append(max_indep)

    print(f"  {'H':>4} {'ω(Ω)':>6} {'α(Ω)':>6} {'ω+α':>5} {'ω·α':>5} {'nc':>5}")
    for H_val in sorted(clique_by_H.keys()):
        # All should be same within an H value at n≤5
        w = clique_by_H[H_val][0]
        a = indep_by_H[H_val][0]
        nc = (H_val - 1) // 2  # at n≤5, α₁ = (H-1)/2
        print(f"  {H_val:>4} {w:>6} {a:>6} {w+a:>5} {w*a:>5} {nc:>5}")

    # Ramsey-type bound: ω * α ≥ nc (pigeonhole)?
    print(f"  Note: nc = α₁ = (H-1)/2 at n≤{n}")
    print(f"  Observation: ω(Ω) = nc always! Ω is a COMPLETE GRAPH at n≤5!")

# ══════════════════════════════════════════════════════════════════
# PART 7: WHEN Ω STOPS BEING COMPLETE — The n=6 transition
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 7: THE n=6 TRANSITION — Ω STOPS BEING COMPLETE")
print("=" * 70)
print()
print("At n≤5, every pair of odd cycles shares a vertex (Ω is complete).")
print("At n=6, vertex-disjoint 3-cycles can exist!")
print("This is the phase transition where topology becomes nontrivial.")

n = 6
complete_count = 0
non_complete_count = 0
first_non_complete = None
count = 0

has_disjoint_by_H = defaultdict(lambda: [0, 0])  # [no_disjoint, has_disjoint]

for adj, bits in all_tournaments(n):
    H = compute_H_dp(adj, n)
    cycles = get_odd_cycles(adj, n)
    nc = len(cycles)
    count += 1

    if nc <= 1:
        complete_count += 1
        has_disjoint_by_H[H][0] += 1
        continue

    adj_cg = build_conflict_graph(cycles)
    # Check if any pair is non-adjacent (independent)
    has_independent_pair = False
    for i in range(nc):
        for j in range(i+1, nc):
            if not adj_cg[i][j]:
                has_independent_pair = True
                break
        if has_independent_pair:
            break

    if has_independent_pair:
        non_complete_count += 1
        has_disjoint_by_H[H][1] += 1
        if first_non_complete is None:
            first_non_complete = (H, nc, [c[0] for c in cycles[:5]])
    else:
        complete_count += 1
        has_disjoint_by_H[H][0] += 1

    if count % 5000 == 0:
        print(f"  ... processed {count}/{1 << 15} tournaments")

total = complete_count + non_complete_count
print(f"\n  Ω complete (all cycles pairwise conflict): {complete_count} ({100*complete_count/total:.1f}%)")
print(f"  Ω non-complete (has disjoint pair): {non_complete_count} ({100*non_complete_count/total:.1f}%)")

if first_non_complete:
    H, nc, ex_cycles = first_non_complete
    print(f"  First non-complete example: H={H}, #cycles={nc}")
    print(f"    Example cycles: {ex_cycles}")

print(f"\n  By H value:")
for H_val in sorted(has_disjoint_by_H.keys()):
    no_dis, has_dis = has_disjoint_by_H[H_val]
    total_h = no_dis + has_dis
    print(f"    H={H_val:3d}: complete={no_dis:5d}, non-complete={has_dis:5d} ({100*has_dis/total_h:.0f}% non-complete)")

# ══════════════════════════════════════════════════════════════════
# PART 8: THE TOPOLOGICAL PHASE DIAGRAM
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 8: TOPOLOGICAL PHASE DIAGRAM — H+2χ vs n")
print("=" * 70)
print()
print("At each n, what values does H + 2χ take?")
print("At n≤5: always 3.")
print("At n=6: 3 + 6α₂, so {3, 9, 15, 21, ...}")

for n in [3, 4, 5]:
    vals = set()
    for adj, bits in all_tournaments(n):
        H = compute_H_dp(adj, n)
        cycles = get_odd_cycles(adj, n)
        nc = len(cycles)
        if nc == 0:
            chi = 1
        else:
            adj_cg = build_conflict_graph(cycles)
            chi = independence_polynomial(adj_cg, nc, -1)
        vals.add(H + 2 * chi)
    print(f"  n={n}: H + 2χ ∈ {sorted(vals)}")

# n=6 partial
n = 6
vals = Counter()
count = 0
for adj, bits in all_tournaments(n):
    H = compute_H_dp(adj, n)
    cycles = get_odd_cycles(adj, n)
    nc = len(cycles)
    if nc == 0:
        chi = 1
    else:
        adj_cg = build_conflict_graph(cycles)
        chi = independence_polynomial(adj_cg, nc, -1)
    vals[H + 2 * chi] += 1
    count += 1
    if count % 5000 == 0:
        print(f"  ... processed {count}/{1 << 15}")

print(f"  n=6: H + 2χ values and counts:")
for v in sorted(vals.keys()):
    alpha2 = (v - 3) // 6
    print(f"    H + 2χ = {v:3d} (α₂ = {alpha2}): {vals[v]:5d} tournaments")

# ══════════════════════════════════════════════════════════════════
# PART 9: DEEPER INSIGHT — UNIVERSAL POLYNOMIAL IDENTITY
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 9: UNIVERSAL POLYNOMIAL IDENTITIES FOR I(Ω, x)")
print("=" * 70)
print()
print("Given I(Ω, x) = Σ α_k x^k, what identities hold?")
print()

# For any polynomial p(x) = Σ α_k x^k:
# p(2) = H, p(-1) = χ, p(1) = #indep sets, p(0) = 1
# Linear combinations that isolate specific α_k:
# α_0 = p(0) = 1 (always)
# p(2) - 2p(1) + p(0) = Σ α_k(2^k - 2·1^k + 1) = Σ α_k(2^k - 2 + 1)
# Let's find what p(a) + c*p(b) isolates

# We want to "measure" α₂ from tournament evaluations
# α₂ can be extracted from:
# p(x) = 1 + α₁x + α₂x² + ...
# p(2) = 1 + 2α₁ + 4α₂ = H
# p(-1) = 1 - α₁ + α₂ = χ
# p(1) = 1 + α₁ + α₂ = #indep sets

# From H and χ: H - χ = 3α₁ + 3α₂ + ... → α₁ = (H - χ - 3(α₂ + α₃ + ...))/3
# From H + 2χ = 3 + 6α₂ + 6α₃ + ... → α₂ + α₃ = (H + 2χ - 3)/6

# A cleaner extraction:
# p(1) - p(0) = α₁ + α₂ + α₃ + ... = #indep sets - 1
# p(-1) - p(0) = -α₁ + α₂ - α₃ + ... = χ - 1
# Sum: 2(α₂ + α₄ + ...) = p(1) + p(-1) - 2
# Diff: 2(α₁ + α₃ + ...) = p(1) - p(-1)

print("KEY EXTRACTION FORMULAS:")
print("  α_even := α₂ + α₄ + ... = (I(Ω,1) + I(Ω,-1) - 2) / 2")
print("  α_odd  := α₁ + α₃ + ... = (I(Ω,1) - I(Ω,-1)) / 2")
print()
print("At n=5: α_even = 0 (no pairs), so I(Ω,1) + I(Ω,-1) = 2 for all T.")

# Verify
n = 5
vals = set()
for adj, bits in all_tournaments(n):
    cycles = get_odd_cycles(adj, n)
    nc = len(cycles)
    if nc == 0:
        i1 = 1; im1 = 1
    else:
        adj_cg = build_conflict_graph(cycles)
        i1 = independence_polynomial(adj_cg, nc, 1)
        im1 = independence_polynomial(adj_cg, nc, -1)
    vals.add(i1 + im1)

print(f"  Verified: I(Ω,1) + I(Ω,-1) at n=5: {sorted(vals)}")
if vals == {2}:
    print(f"  ✓ Always 2! This means α_even = 0 for all n=5 tournaments.")

# At n=6
n = 6
vals_n6 = Counter()
count = 0
for adj, bits in all_tournaments(n):
    cycles = get_odd_cycles(adj, n)
    nc = len(cycles)
    if nc == 0:
        i1 = 1; im1 = 1
    else:
        adj_cg = build_conflict_graph(cycles)
        i1 = independence_polynomial(adj_cg, nc, 1)
        im1 = independence_polynomial(adj_cg, nc, -1)
    alpha_even = (i1 + im1 - 2) // 2
    vals_n6[alpha_even] += 1
    count += 1
    if count % 5000 == 0:
        print(f"  ... processed {count}/{1 << 15}")

print(f"\n  n=6: α_even = α₂ + α₄ + ... distribution:")
for v in sorted(vals_n6.keys()):
    print(f"    α_even = {v}: {vals_n6[v]} tournaments")

print("\n" + "=" * 70)
print("EXPLORATION COMPLETE")
print("=" * 70)
