#!/usr/bin/env python3
"""
CATEGORICAL AND TOPOLOGICAL ANALYSIS OF IP(G, 2) = H
opus-2026-03-14-S89b

The Independence Polynomial formula H = IP(G(T), 2) has deep structure:

1. DELETION-CONTRACTION: IP(G,x) = IP(G-v,x) + x*IP(G-N[v],x)
   This means H satisfies a RECURSIVE formula by removing cycles.

2. CLUSTER EXPANSION: H = exp(sum over connected clusters)
   The log of the independence polynomial gives the "free energy"

3. SIMPLICIAL COMPLEX: The faces of the independence complex
   form a simplicial complex whose f-vector gives H.

4. HOMOTOPY TYPE: Independence complexes of well-covered graphs
   are homotopy equivalent to wedges of spheres.

5. BAER SUBPLANE: The odd-cycle disjointness graph G(T) for the
   "regular tournament" might have Baer subplane structure.

6. CONE FUNCTOR: How does G(Cone(T)) relate to G(T)?
   Cone adds a new vertex that beats everyone, creating new 3-cycles.

Let's compute and analyze these structures.
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict
from math import factorial, comb, log, exp

def compute_H(n, adj):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            val = dp[(mask, v)]
            for u in range(n):
                if mask & (1 << u): continue
                if adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def tournament_adj_list(n, bits):
    adj = [[False]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
            idx += 1
    return adj

def get_all_odd_cycles(n, adj):
    """Get all directed odd cycles with their vertex sets."""
    cycles = []  # list of (frozenset(vertices), cycle_length)

    # 3-cycles
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    cycles.append((frozenset([i,j,k]), 3))

    # 5-cycles
    for combo in combinations(range(n), 5):
        count = 0
        for perm in permutations(combo):
            ok = True
            for idx in range(5):
                if not adj[perm[idx]][perm[(idx+1)%5]]:
                    ok = False
                    break
            if ok:
                count += 1
        num_5 = count // 5
        for _ in range(num_5):
            cycles.append((frozenset(combo), 5))

    return cycles

def build_disjointness_graph(cycles):
    """Build adjacency list of the cycle disjointness graph.
    Edge = cycles SHARE a vertex (i.e., NOT disjoint)."""
    n = len(cycles)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if len(cycles[i][0] & cycles[j][0]) > 0:
                adj[i].append(j)
                adj[j].append(i)
    return adj

def enumerate_independent_sets(adj_list, n):
    """Enumerate all independent sets in a graph. Return f-vector."""
    f_vec = Counter()  # f_vec[k] = number of independent sets of size k

    for mask in range(1 << n):
        selected = [i for i in range(n) if mask & (1 << i)]
        k = len(selected)
        # Check independence
        ok = True
        for a in range(k):
            for b in range(a+1, k):
                if selected[b] in adj_list[selected[a]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            f_vec[k] += 1

    return f_vec

print("="*70)
print("CATEGORICAL & TOPOLOGICAL ANALYSIS OF IP(G,2) = H")
print("opus-2026-03-14-S89b")
print("="*70)

# ===== PART 1: f-VECTORS AND EULER CHARACTERISTICS =====
print("\n" + "="*70)
print("PART 1: f-VECTORS OF THE INDEPENDENCE COMPLEX")
print("="*70)

for n in [4, 5]:
    m = n*(n-1)//2
    print(f"\nn={n}: f-vector statistics")
    fvec_counter = Counter()

    for bits in range(1 << m):
        adj = tournament_adj_list(n, bits)
        cycles = get_all_odd_cycles(n, adj)
        if len(cycles) > 20:
            continue  # skip if too many cycles
        cycle_adj = build_disjointness_graph(cycles)
        fvec = enumerate_independent_sets(cycle_adj, len(cycles))
        fvec_tuple = tuple(fvec[k] for k in range(max(fvec.keys())+1)) if fvec else (1,)
        fvec_counter[fvec_tuple] += 1

        # Compute H from f-vector
        H_from_fvec = sum(fvec[k] * (2**k) for k in fvec)
        H = compute_H(n, adj)
        assert H_from_fvec == H, f"f-vector mismatch: {H_from_fvec} != {H}"

    print(f"  Distinct f-vectors: {len(fvec_counter)}")
    for fv, count in sorted(fvec_counter.items(), key=lambda x: -x[1])[:10]:
        euler = sum((-1)**k * (fv[k] if k < len(fv) else 0) for k in range(len(fv)))
        H_val = sum((2**k) * (fv[k] if k < len(fv) else 0) for k in range(len(fv)))
        print(f"  f-vec={fv}: count={count}, H={H_val}, chi={euler}")

# ===== PART 2: DELETION-CONTRACTION FOR H =====
print("\n" + "="*70)
print("PART 2: DELETION-CONTRACTION RECURSION")
print("="*70)

print("""
IP(G, x) = IP(G-v, x) + x * IP(G-N[v], x)

where G-v removes vertex v, and G-N[v] removes v and all its neighbors.

Applied to H(T): pick any odd cycle C in T.
  H(T) = (H restricted to configs without C) + 2*(H restricted to configs
           using C and no cycle sharing a vertex with C)

This gives a RECURSIVE computation of H through the cycle structure!
""")

# Demonstrate on n=5
n = 5
# Pick a tournament with interesting cycle structure
bits = 0b1111100000  # some tournament
adj = tournament_adj_list(n, bits)
H = compute_H(n, adj)
cycles = get_all_odd_cycles(n, adj)
print(f"Example: n=5, bits={bits}, H={H}")
print(f"  Odd cycles: {len(cycles)}")
for i, (vs, length) in enumerate(cycles):
    print(f"    C{i}: {sorted(vs)} (length {length})")

if len(cycles) > 0:
    cycle_adj = build_disjointness_graph(cycles)
    fvec = enumerate_independent_sets(cycle_adj, len(cycles))
    print(f"  f-vector: {dict(fvec)}")
    print(f"  H = sum 2^k * f_k = {sum(2**k * fvec[k] for k in fvec)}")

    # Deletion-contraction on first cycle
    v = 0  # remove cycle 0
    # G - v: remove cycle 0
    remaining = list(range(1, len(cycles)))
    if remaining:
        sub_adj = [[] for _ in range(len(remaining))]
        for i, ci in enumerate(remaining):
            for j, cj in enumerate(remaining):
                if j > i and cj in cycle_adj[ci]:
                    sub_adj[i].append(j)
                    sub_adj[j].append(i)
        fvec_minus = enumerate_independent_sets(sub_adj, len(remaining))
        ip_minus = sum(2**k * fvec_minus[k] for k in fvec_minus)
    else:
        ip_minus = 1

    # G - N[v]: remove cycle 0 and all its neighbors
    neighbors = set(cycle_adj[0])
    independent = [c for c in range(1, len(cycles)) if c not in neighbors]
    if independent:
        sub_adj2 = [[] for _ in range(len(independent))]
        for i, ci in enumerate(independent):
            for j, cj in enumerate(independent):
                if j > i and cj in cycle_adj[ci]:
                    sub_adj2[i].append(j)
                    sub_adj2[j].append(i)
        fvec_closed = enumerate_independent_sets(sub_adj2, len(independent))
        ip_closed = sum(2**k * fvec_closed[k] for k in fvec_closed)
    else:
        ip_closed = 1

    print(f"\n  Deletion-contraction on C0 = {sorted(cycles[0][0])}:")
    print(f"    IP(G-C0, 2) = {ip_minus}")
    print(f"    IP(G-N[C0], 2) = {ip_closed}")
    print(f"    IP(G-C0, 2) + 2*IP(G-N[C0], 2) = {ip_minus + 2*ip_closed}")
    print(f"    H = {H}")
    print(f"    Match: {ip_minus + 2*ip_closed == H}")

# ===== PART 3: CONE FUNCTOR =====
print("\n" + "="*70)
print("PART 3: CONE FUNCTOR — How does coning affect the IP?")
print("="*70)

print("""
Cone(T) adds vertex n that beats all others.
This creates n-1 new 3-cycles (for each arc (i,j) where j->i in T,
  we get the 3-cycle n->i->j->... wait, let's think:

In T, for each edge i->j, in Cone(T):
  n beats i and j. So n->i->j exists.
  If j->i in T, then j->i->... no, i->j means i beats j.
  New 3-cycles: for each pair (i,j) with j->i in T (i.e., j beats i),
  we get the 3-cycle n->i, i<-j, j<-n: n->i<-j<-n? No...
  n->i (n beats i), j->i (j beats i), n->j (n beats j).
  3-cycle: need a,b,c with a->b->c->a.
  n->j->i->n: need i->n, but n beats everyone, so n->i. No 3-cycle.

Actually: n beats everyone, so no one beats n. Therefore n CANNOT
be in any 3-cycle (a 3-cycle a->b->c->a requires each vertex to
lose to one other vertex in the cycle, but n doesn't lose to anyone).

Wait — that means Cone(T) has EXACTLY THE SAME odd cycles as T!
(Plus any new cycles among the original vertices, but the adjacency
among original vertices is unchanged.)

So G(Cone(T)) = G(T), and H(Cone(T)) should satisfy:
  H(Cone(T)) = 1 + 2*(same cycles) + 4*(same pairs) + ...

But we KNOW that H(Cone(T)) = (n+1)*H(T) from THM-205!
""")

# Verify: H(Cone(T)) vs (n+1)*H(T)
n = 4
m = n*(n-1)//2
print(f"Testing cone functor at n={n}:")

mismatches = 0
for bits in range(1 << m):
    adj = tournament_adj_list(n, bits)
    H_T = compute_H(n, adj)

    # Build Cone(T): new vertex n beats everyone
    n_cone = n + 1
    adj_cone = [[False]*n_cone for _ in range(n_cone)]
    for i in range(n):
        for j in range(n):
            adj_cone[i][j] = adj[i][j]
    for i in range(n):
        adj_cone[n][i] = True  # n beats i

    H_cone = compute_H(n_cone, adj_cone)

    # Also check: do the odd cycles change?
    cycles_T = get_all_odd_cycles(n, adj)
    cycles_cone = get_all_odd_cycles(n_cone, adj_cone)

    # Filter cone cycles: those NOT involving vertex n
    cycles_cone_no_n = [(vs, l) for vs, l in cycles_cone if n not in vs]
    # Cycles involving n?
    cycles_cone_with_n = [(vs, l) for vs, l in cycles_cone if n in vs]

    if len(cycles_cone_with_n) > 0 and bits < 5:
        print(f"  bits={bits}: Cone has cycles involving new vertex! {cycles_cone_with_n}")

    expected = (n+1) * H_T
    if H_cone != expected:
        mismatches += 1
        if mismatches <= 3:
            print(f"  bits={bits}: H(Cone)={H_cone}, (n+1)*H(T)={expected}, diff={H_cone-expected}")

    if bits < 3:
        print(f"  bits={bits}: H(T)={H_T}, H(Cone)={H_cone}, (n+1)*H={expected}, cycles_T={len(cycles_T)}, cycles_cone_no_n={len(cycles_cone_no_n)}")

if mismatches == 0:
    print(f"  H(Cone(T)) = (n+1)*H(T) verified for all {1<<m} tournaments!")
else:
    print(f"  {mismatches} mismatches!")

# ===== PART 4: THE PARADOX — CONE PRESERVES CYCLES BUT MULTIPLIES H =====
print("\n" + "="*70)
print("PART 4: THE CONE PARADOX")
print("="*70)

print("""
PARADOX:
  1. Cone(T) has the SAME odd cycles as T (new vertex can't be in cycles)
  2. So G(Cone(T)) = G(T) as a graph
  3. Therefore IP(G(Cone(T)), 2) = IP(G(T), 2) = H(T)
  4. But H(Cone(T)) = (n+1)*H(T) ≠ H(T)

THIS MEANS THE FORMULA H = IP(G, 2) CANNOT HOLD FOR ALL n!

Wait — let me check. Does Cone(T) really have the same odd cycles?
The induced tournament on the original n vertices is IDENTICAL.
But vertex n beats everyone, so:
  - No 3-cycle can include n (n beats both others, so no return arc to n)
  - BUT: new 5-cycles could include n!
    Example: n->a->b->c->d->n. This needs d->n, but n beats d. IMPOSSIBLE.
  - ALL cycle types need a return arc TO n, which is impossible.
  So YES, Cone(T) has exactly the same directed cycles as T.

So the formula H = IP(G, 2) MUST FAIL at some n ≥ 7 where the cone
creates a discrepancy!

Actually wait — at n=4, Cone gives n=5.
  H(T) at n=4: determined by IP formula
  H(Cone(T)) at n=5: should be 5*H(T)
  But IP(G(Cone(T)), 2) = IP(G(T), 2) = H(T)
  So we need 5*H(T) = H(T)? That's only true for H(T)=0, impossible.

THIS MEANS H ≠ IP(G, 2) IN GENERAL!

Unless... the cone introduces new cycles I'm not counting.
Let me check very carefully.
""")

# Very careful check: enumerate ALL odd cycles in Cone(T)
n = 4
bits = 5  # arbitrary tournament on 4 vertices
adj = tournament_adj_list(n, bits)
H_T = compute_H(n, adj)

# Build Cone
n_cone = 5
adj_cone = [[False]*n_cone for _ in range(n_cone)]
for i in range(n):
    for j in range(n):
        adj_cone[i][j] = adj[i][j]
for i in range(n):
    adj_cone[n][i] = True

H_cone = compute_H(n_cone, adj_cone)
print(f"T on 4 vertices (bits={bits}): H(T) = {H_T}")
print(f"Cone(T) on 5 vertices: H(Cone(T)) = {H_cone}, 5*H(T) = {5*H_T}")

cycles_T = get_all_odd_cycles(n, adj)
cycles_cone = get_all_odd_cycles(n_cone, adj_cone)
print(f"\nCycles in T: {len(cycles_T)}")
for vs, l in cycles_T:
    print(f"  {sorted(vs)} (length {l})")
print(f"\nCycles in Cone(T): {len(cycles_cone)}")
for vs, l in cycles_cone:
    marker = " *** NEW" if n in vs else ""
    print(f"  {sorted(vs)} (length {l}){marker}")

# Check IP
ip_T = sum(2**k * enumerate_independent_sets(
    build_disjointness_graph(cycles_T), len(cycles_T))[k]
    for k in enumerate_independent_sets(
        build_disjointness_graph(cycles_T), len(cycles_T)))
ip_cone = sum(2**k * enumerate_independent_sets(
    build_disjointness_graph(cycles_cone), len(cycles_cone))[k]
    for k in enumerate_independent_sets(
        build_disjointness_graph(cycles_cone), len(cycles_cone)))

print(f"\nIP(G(T), 2) = {ip_T}")
print(f"IP(G(Cone(T)), 2) = {ip_cone}")
print(f"H(T) = {H_T}")
print(f"H(Cone(T)) = {H_cone}")

if ip_T == H_T and ip_cone == H_cone:
    print("BOTH match! Cone DID introduce new cycles.")
elif ip_T == H_T and ip_cone != H_cone:
    print("T matches but Cone doesn't! Formula fails at n=5 for cone tournaments.")
    print(f"Discrepancy: IP={ip_cone}, H={H_cone}, diff={H_cone-ip_cone}")

# ===== PART 5: WHAT CYCLES DOES THE CONE CREATE? =====
print("\n" + "="*70)
print("PART 5: UNDERSTANDING THE CONE'S CYCLE CREATION")
print("="*70)

# Check if vertex 4 (the cone point) appears in any cycle
for vs, l in cycles_cone:
    if n in vs:
        print(f"  Cone vertex {n} IS in cycle {sorted(vs)} length {l}")
        # Show the actual directed cycle
        verts = sorted(vs)
        for perm in permutations(verts):
            ok = True
            for idx in range(l):
                if not adj_cone[perm[idx]][perm[(idx+1)%l]]:
                    ok = False
                    break
            if ok:
                print(f"    Directed: {' -> '.join(map(str, perm))} -> {perm[0]}")
                break

# Check more tournaments
print(f"\nSystematic check: how many new cycles does cone create?")
for bits in range(1 << m):
    adj = tournament_adj_list(n, bits)
    adj_cone = [[False]*n_cone for _ in range(n_cone)]
    for i in range(n):
        for j in range(n):
            adj_cone[i][j] = adj[i][j]
    for i in range(n):
        adj_cone[n][i] = True

    cycles_T = get_all_odd_cycles(n, adj)
    cycles_cone = get_all_odd_cycles(n_cone, adj_cone)
    new_cycles = [(vs, l) for vs, l in cycles_cone if n in vs]

    if len(new_cycles) > 0 and bits < 10:
        print(f"  bits={bits}: {len(new_cycles)} new cycles involving cone vertex")
        for vs, l in new_cycles:
            print(f"    {sorted(vs)} (len {l})")

print("\n" + "="*70)
print("DONE — Categorical analysis")
print("="*70)
