#!/usr/bin/env python3
"""
simplicial_H_s112.py — Deep investigation of simplicial Hamiltonian paths
kind-pasteur-2026-03-15-S112

KEY FINDING from preliminary computation:
  At n=4: sim_H in {0, 1} while H in {1, 3, 5}.
  sim_H = 0 when transitive triples create contradictory constraints.
  This is a NEW INVARIANT distinct from H(T).

Questions:
1. What is the distribution of sim_H at n=4,5,6?
2. Is sim_H always 0 or 1? Or can it be larger?
3. Does sim_H = |linear extensions of the TRANSITIVE CLOSURE|?
4. What is the "simplicial OCF"? Is sim_H related to an independence polynomial?
5. How does this connect to the Cayley-Delannoy formula?
"""

from itertools import permutations, combinations
from collections import defaultdict, Counter
from math import factorial, comb

def analyze_tournament(T, n):
    """Full analysis of a tournament."""
    # Hamiltonian paths
    ham_paths = []
    for perm in permutations(range(n)):
        ok = True
        for p in range(n-1):
            if T[perm[p]][perm[p+1]] != 1:
                ok = False
                break
        if ok:
            ham_paths.append(perm)
    H = len(ham_paths)

    # 3-cycles and transitive triples
    cycles_3 = 0
    trans_triples = []  # (a,b,c) with a>b>c in T
    for i, j, k in combinations(range(n), 3):
        s = T[i][j] + T[j][k] + T[k][i]
        if s == 3 or s == 0:  # 3-cycle (one direction or the other)
            cycles_3 += 1
        # Find the transitive ordering
        for a, b, c in permutations([i,j,k]):
            if T[a][b] and T[b][c] and T[a][c]:
                trans_triples.append((a,b,c))
                break

    # Simplicial H: paths consistent with ALL transitive triples
    sim_paths = []
    for perm in ham_paths:
        pos = {perm[p]: p for p in range(n)}
        ok = True
        for a, b, c in trans_triples:
            if not (pos[a] < pos[b] < pos[c]):
                ok = False
                break
        if ok:
            sim_paths.append(perm)
    sim_H = len(sim_paths)

    # Transitive closure analysis
    # The transitive triples define a partial order (if consistent)
    # Actually: they define constraints. Are they always consistent?

    return H, sim_H, cycles_3, len(trans_triples), trans_triples

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    m = n * (n-1) // 2
    for bits in range(2**m):
        T = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    T[i][j] = 1
                else:
                    T[j][i] = 1
                idx += 1
        yield T, bits

# ============================================================
# PART 1: Exhaustive analysis at n=4
# ============================================================
print("="*70)
print("PART 1: EXHAUSTIVE ANALYSIS AT n=4")
print("="*70)

stats_4 = defaultdict(int)
for T, bits in all_tournaments(4):
    H, sH, c3, nt, _ = analyze_tournament(T, 4)
    stats_4[(H, sH, c3)] += 1

print("(H, sim_H, #3cycles): count")
for key in sorted(stats_4.keys()):
    print(f"  {key}: {stats_4[key]}")

# ============================================================
# PART 2: Exhaustive analysis at n=5
# ============================================================
print()
print("="*70)
print("PART 2: EXHAUSTIVE ANALYSIS AT n=5")
print("="*70)

stats_5 = defaultdict(int)
sH_dist_5 = Counter()
for T, bits in all_tournaments(5):
    H, sH, c3, nt, _ = analyze_tournament(T, 5)
    stats_5[(H, sH)] += 1
    sH_dist_5[sH] += 1

print("(H, sim_H): count")
for key in sorted(stats_5.keys()):
    print(f"  ({key[0]:2d}, {key[1]:2d}): {stats_5[key]:5d}")

print()
print("sim_H distribution:")
for sH in sorted(sH_dist_5.keys()):
    print(f"  sim_H={sH}: {sH_dist_5[sH]} ({sH_dist_5[sH]/1024*100:.1f}%)")

# ============================================================
# PART 3: What IS sim_H?
# ============================================================
print()
print("="*70)
print("PART 3: WHAT IS sim_H?")
print("="*70)

# Hypothesis: sim_H = number of linear extensions of the
# INTERSECTION partial order (the partial order where a > b
# iff EVERY transitive triple containing a,b has a before b)

# Actually, simpler: the transitive triples ARE the partial order.
# a > b in the partial order iff T[a][b] = 1 AND there's no path
# b -> ... -> a through transitive triples.

# Wait: the transitive triples (a,b,c) with a>b>c in T give us
# a>b, b>c, a>c as constraints. The union of all such constraints
# IS the edge set of the tournament MINUS the 3-cycle edges.

# More precisely: an edge a>b of T is "transitive-confirmed" if
# there exists some c such that (a,b,c) or (c,a,b) or (a,c,b)
# is a transitive triple containing a>b.

# Let's check: is sim_H = #linear extensions of the "transitive core"?

print("Checking: sim_H = #linear extensions of transitive core?")
print()

def transitive_core(T, n):
    """Return the partial order of 'transitively confirmed' edges."""
    # An edge a>b is confirmed if it appears in some transitive triple
    confirmed = [[False]*n for _ in range(n)]
    for i, j, k in combinations(range(n), 3):
        for a, b, c in permutations([i,j,k]):
            if T[a][b] and T[b][c] and T[a][c]:
                confirmed[a][b] = True
                confirmed[b][c] = True
                confirmed[a][c] = True
                break
    return confirmed

def count_linear_extensions(P, n):
    """Count linear extensions of partial order P via brute force."""
    count = 0
    for perm in permutations(range(n)):
        pos = {perm[i]: i for i in range(n)}
        ok = True
        for i in range(n):
            for j in range(n):
                if P[i][j] and pos[i] > pos[j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            count += 1
    return count

# Test at n=4
mismatches = 0
for T, bits in all_tournaments(4):
    H, sH, c3, nt, _ = analyze_tournament(T, 4)
    P = transitive_core(T, 4)
    le = count_linear_extensions(P, 4)
    if sH != le:
        mismatches += 1
        if mismatches <= 3:
            print(f"  MISMATCH at bits={bits}: sim_H={sH}, LE(core)={le}")

if mismatches == 0:
    print("  n=4: sim_H = #LE(transitive core) for ALL 64 tournaments!")
else:
    print(f"  n=4: {mismatches} mismatches")

# Test at n=5
mismatches_5 = 0
for T, bits in all_tournaments(5):
    H, sH, c3, nt, _ = analyze_tournament(T, 5)
    P = transitive_core(T, 5)
    le = count_linear_extensions(P, 5)
    if sH != le:
        mismatches_5 += 1

print(f"  n=5: {mismatches_5} mismatches out of 1024 tournaments")
if mismatches_5 == 0:
    print("  sim_H = #LE(transitive core) for ALL 1024 tournaments!")

# ============================================================
# PART 4: Is sim_H always 0 or 1? Parity?
# ============================================================
print()
print("="*70)
print("PART 4: PARITY OF sim_H")
print("="*70)

# At n=4: sim_H in {0, 1}
# At n=5: check
sH_vals_4 = sorted(set(key[1] for key in stats_4.keys()))
print(f"  n=4: sim_H values = {sH_vals_4}")
print(f"  n=5: sim_H values = {sorted(sH_dist_5.keys())}")

# Is sim_H always even? Always odd? Or mixed?
# At n=5: sim_H can be 0, 1, or possibly more
print()
print("  Parity check at n=5:")
for sH in sorted(sH_dist_5.keys()):
    print(f"    sim_H={sH}: parity={'even' if sH%2==0 else 'odd'}, count={sH_dist_5[sH]}")

# ============================================================
# PART 5: sim_H as independence polynomial?
# ============================================================
print()
print("="*70)
print("PART 5: SIMPLICIAL CONFLICT GRAPH")
print("="*70)
print("""
For tournament H(T), the OCF gives H = I(Omega(T), 2) where Omega is
the conflict graph of odd cycles.

Question: is there a 'simplicial OCF' where sim_H = I(Omega_simp, 2)
for some simplicial conflict structure?

The simplicial constraints are: for each transitive triple (a,b,c),
the ordering must have a before b before c. Two triples CONFLICT
if they impose contradictory orderings on shared vertices.

This is exactly a CONSTRAINT SATISFACTION problem, and the number
of solutions is sim_H.
""")

# Build the constraint conflict graph for a specific tournament
# and check if sim_H relates to its independence polynomial

def constraint_conflict_graph(trans_triples, n):
    """Build graph where vertices = transitive triples,
    edges = pairs that impose contradictory constraints."""
    # Two triples (a1,b1,c1) and (a2,b2,c2) conflict if
    # they share vertices and impose contradictory orderings.
    # Specifically: if the union of their constraints is infeasible.

    nt = len(trans_triples)
    conflicts = [[False]*nt for _ in range(nt)]

    for i in range(nt):
        a1, b1, c1 = trans_triples[i]
        for j in range(i+1, nt):
            a2, b2, c2 = trans_triples[j]
            # Check if constraints are contradictory
            # Constraints from i: a1<b1, b1<c1, a1<c1 (in position)
            # Constraints from j: a2<b2, b2<c2, a2<c2
            # Contradictory if there's a cycle in the union
            edges = {}
            for x, y in [(a1,b1),(b1,c1),(a1,c1),(a2,b2),(b2,c2),(a2,c2)]:
                if x not in edges:
                    edges[x] = set()
                edges[x].add(y)

            # Check for cycle via DFS
            visited = set()
            def has_cycle(node, path):
                if node in path:
                    return True
                if node in visited:
                    return False
                visited.add(node)
                path.add(node)
                for next_node in edges.get(node, []):
                    if has_cycle(next_node, path):
                        return True
                path.remove(node)
                return False

            cycle = False
            for start in edges:
                visited.clear()
                if has_cycle(start, set()):
                    cycle = True
                    break

            conflicts[i][j] = conflicts[j][i] = cycle

    return conflicts

# Analyze for a specific n=5 tournament with interesting sim_H
print("Example analysis for n=5 tournaments:")
for T, bits in all_tournaments(5):
    H, sH, c3, nt, tris = analyze_tournament(T, 5)
    if H == 5 and sH == 1:
        print(f"\n  bits={bits}: H={H}, sim_H={sH}, c3={c3}, #trans_triples={nt}")
        # Show the transitive triples
        print(f"  Transitive triples: {tris[:10]}")
        # Build conflict graph
        G = constraint_conflict_graph(tris, 5)
        n_conflicts = sum(1 for i in range(nt) for j in range(i+1, nt) if G[i][j])
        print(f"  Conflicts between triples: {n_conflicts}")
        break

# ============================================================
# PART 6: The key invariant: number of ACYCLIC orientations
# ============================================================
print()
print("="*70)
print("PART 6: sim_H = LINEAR EXTENSIONS OF TRANSITIVE REDUCTION")
print("="*70)

# Since sim_H = #LE(transitive core), and the transitive core
# is the tournament minus 3-cycle edges, let's understand this better.

# The transitive core of a tournament T is the DAG obtained by
# removing edges that participate in 3-cycles.

# Actually, let's check: is the transitive core always a DAG?
non_dag = 0
for T, bits in all_tournaments(5):
    P = transitive_core(T, 5)
    # Check for cycles in P
    visited = [0]*5  # 0=unvisited, 1=in_stack, 2=done
    def has_cycle_dag(v):
        visited[v] = 1
        for u in range(5):
            if P[v][u]:
                if visited[u] == 1:
                    return True
                if visited[u] == 0 and has_cycle_dag(u):
                    return True
        visited[v] = 2
        return False

    is_dag = True
    for v in range(5):
        visited = [0]*5
        if has_cycle_dag(v):
            is_dag = False
            break

    if not is_dag:
        non_dag += 1

print(f"  n=5: transitive core is DAG for all {1024-non_dag}/1024 tournaments")
if non_dag > 0:
    print(f"  WARNING: {non_dag} tournaments have CYCLIC transitive cores!")
else:
    print(f"  The transitive core is ALWAYS a DAG!")

# ============================================================
# PART 7: Connection to Redei and OCF
# ============================================================
print()
print("="*70)
print("PART 7: SIMPLICIAL REDEI THEOREM?")
print("="*70)

# Redei: H(T) is always ODD.
# Question: is sim_H always of a specific parity?

print("Parity of sim_H at n=4:")
for sH_val in sorted(set(key[1] for key in stats_4.keys())):
    total = sum(v for key, v in stats_4.items() if key[1] == sH_val)
    print(f"  sim_H={sH_val} ({'odd' if sH_val%2==1 else 'even'}): {total} tournaments")

print()
print("Parity of sim_H at n=5:")
for sH in sorted(sH_dist_5.keys()):
    p = 'odd' if sH % 2 == 1 else 'even'
    print(f"  sim_H={sH} ({p}): {sH_dist_5[sH]} tournaments")

# Check: is sim_H odd when nonzero? Or some other parity pattern?
print()
nonzero_parity = Counter()
for T, bits in all_tournaments(5):
    H, sH, c3, nt, _ = analyze_tournament(T, 5)
    if sH > 0:
        nonzero_parity[sH % 2] += 1

print(f"When sim_H > 0 at n=5:")
print(f"  sim_H odd: {nonzero_parity[1]}")
print(f"  sim_H even: {nonzero_parity[0]}")

# ============================================================
# PART 8: Statistics and expected value
# ============================================================
print()
print("="*70)
print("PART 8: STATISTICS OF sim_H")
print("="*70)

for n in [4, 5]:
    print(f"\nn={n}:")
    total_T = 2**(n*(n-1)//2)

    sum_H = 0
    sum_sH = 0
    sum_H2 = 0
    sum_sH2 = 0

    for T, bits in all_tournaments(n):
        H, sH, _, _, _ = analyze_tournament(T, n)
        sum_H += H
        sum_sH += sH
        sum_H2 += H*H
        sum_sH2 += sH*sH

    E_H = Fraction(sum_H, total_T)
    E_sH = Fraction(sum_sH, total_T)
    E_H2 = Fraction(sum_H2, total_T)
    E_sH2 = Fraction(sum_sH2, total_T)

    Var_H = E_H2 - E_H * E_H
    Var_sH = E_sH2 - E_sH * E_sH

    print(f"  E[H]    = {E_H} = {float(E_H):.4f}")
    print(f"  E[sH]   = {E_sH} = {float(E_sH):.4f}")
    print(f"  Var[H]  = {float(Var_H):.4f}")
    print(f"  Var[sH] = {float(Var_sH):.4f}")
    print(f"  E[sH]/E[H] = {float(E_sH/E_H):.4f}")
    print(f"  Ratio of expectations: {float(E_sH/E_H)*100:.1f}%")
    print(f"  Fraction with sH=0: {sum(1 for T,_ in all_tournaments(n) if analyze_tournament(T,n)[1]==0)/total_T*100:.1f}%")

print("\nDone!")
