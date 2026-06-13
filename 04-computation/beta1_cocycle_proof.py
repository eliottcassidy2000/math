#!/usr/bin/env python3
"""
β₁ ≤ 1 via cocycle analysis.

Key discovery: H¹(T) = cocycles/coboundaries has dim ≤ 1.
A cocycle w assigns weights to edges such that w(a,b)+w(b,c)=w(a,c)
for all transitive triples.

The question: why is dim(H¹) ≤ 1? Can we find an explicit cocycle
formula that works for ALL tournaments with β₁=1?

opus-2026-03-08
"""
import numpy as np
from collections import defaultdict, Counter
from math import comb

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def edges_list(A, n):
    return [(i,j) for i in range(n) for j in range(n) if i != j and A[i][j]]

def transitive_triples(A, n):
    tt = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    tt.append((a,b,c))
    return tt

def compute_cocycle_space(A, n):
    """Compute H¹ = cocycles/coboundaries."""
    el = edges_list(A, n)
    edge_idx = {e: i for i, e in enumerate(el)}
    ne = len(el)
    tt = transitive_triples(A, n)

    # Cocycle condition: w(a,b) + w(b,c) - w(a,c) = 0 for each transitive (a,b,c)
    # Build constraint matrix C: each row is a constraint
    C = np.zeros((len(tt), ne))
    for j, (a,b,c) in enumerate(tt):
        C[j, edge_idx[(a,b)]] = 1
        C[j, edge_idx[(b,c)]] = 1
        C[j, edge_idx[(a,c)]] = -1

    # Cocycles = ker(C)
    if C.shape[0] > 0:
        U, S, Vt = np.linalg.svd(C, full_matrices=True)
        rank_C = sum(s > 1e-10 for s in S)
        cocycles = Vt[rank_C:].T  # columns are cocycle basis
    else:
        cocycles = np.eye(ne)

    # Coboundaries = im(δ₀) where δ₀(f)(a→b) = f(b) - f(a)
    delta0 = np.zeros((ne, n))
    for j, (a,b) in enumerate(el):
        delta0[j, b] = 1
        delta0[j, a] = -1

    # Project coboundaries into cocycle space
    cobdy_proj = cocycles.T @ delta0

    rank_cobdy = np.linalg.matrix_rank(cobdy_proj, tol=1e-10)
    dim_H1 = cocycles.shape[1] - rank_cobdy

    # Extract H¹ representative
    representative = None
    if dim_H1 > 0:
        U_c, S_c, Vt_c = np.linalg.svd(cobdy_proj, full_matrices=True)
        noncobdy = U_c[:, rank_cobdy:]
        rep = cocycles @ noncobdy[:, 0]
        rep[np.abs(rep) < 1e-10] = 0
        nz = rep[rep != 0]
        if len(nz) > 0:
            scale = min(abs(nz))
            rep = rep / scale
            representative = np.round(rep).astype(int)

    return {
        'dim_H1': dim_H1,
        'dim_cocycles': cocycles.shape[1],
        'rank_cobdy': rank_cobdy,
        'representative': representative,
        'edges': el,
        'edge_idx': edge_idx,
    }


print("="*72)
print("COCYCLE ANALYSIS FOR β₁ ≤ 1")
print("="*72)

# =================================================================
# COMPUTE COCYCLES FOR ALL TOURNAMENTS
# =================================================================
for n in [3, 4, 5]:
    print(f"\n{'='*72}")
    print(f"n = {n}")
    print(f"{'='*72}")

    max_dim = 0
    cocycle_examples = []

    for A in all_tournaments(n):
        result = compute_cocycle_space(A, n)
        if result['dim_H1'] > max_dim:
            max_dim = result['dim_H1']
        if result['dim_H1'] == 1 and len(cocycle_examples) < 5:
            cocycle_examples.append((A, result))

    print(f"  max dim(H¹) = {max_dim}")

    if max_dim <= 1:
        print(f"  CONFIRMED: H¹ ≤ 1 for all n={n} tournaments")

    for idx, (A, result) in enumerate(cocycle_examples[:3]):
        print(f"\n  Example {idx+1} (dim H¹ = {result['dim_H1']}):")
        w = result['representative']
        el = result['edges']
        if w is not None:
            for i in range(len(w)):
                if w[i] != 0:
                    print(f"    w({el[i][0]}→{el[i][1]}) = {w[i]}")

            # Evaluate on all 3-cycles
            cycles3 = []
            for i in range(n):
                for j in range(n):
                    if j == i or not A[i][j]: continue
                    for k in range(n):
                        if k == i or k == j: continue
                        if A[j][k] and A[k][i]:
                            cycles3.append((i,j,k))

            for c in cycles3:
                edge_idx = result['edge_idx']
                val = sum(w[edge_idx[(c[j], c[(j+1)%3])]] for j in range(3))
                print(f"    w(3-cycle {c}) = {val}")


# =================================================================
# PATTERN IN COCYCLE VALUES
# =================================================================
print(f"\n\n{'='*72}")
print("PATTERN IN COCYCLE VALUES")
print("="*72)

n = 5
cocycle_patterns = Counter()
cocycle_3cycle_vals = Counter()

for A in all_tournaments(n):
    result = compute_cocycle_space(A, n)
    if result['dim_H1'] != 1:
        continue

    w = result['representative']
    el = result['edges']
    edge_idx = result['edge_idx']

    if w is None:
        continue

    # Values of w
    vals = tuple(sorted(set(w[w != 0].tolist())))
    cocycle_patterns[vals] += 1

    # Value on 3-cycles
    for i in range(n):
        for j in range(n):
            if j == i or not A[i][j]: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] and A[k][i]:
                    val = w[edge_idx[(i,j)]] + w[edge_idx[(j,k)]] + w[edge_idx[(k,i)]]
                    cocycle_3cycle_vals[int(val)] += 1

print(f"\nn=5: Cocycle value patterns (distinct nonzero values):")
for pattern, count in sorted(cocycle_patterns.items()):
    print(f"  {pattern}: {count} tournaments")

print(f"\nCocycle evaluated on 3-cycles:")
for val, count in sorted(cocycle_3cycle_vals.items()):
    print(f"  value {val}: {count} occurrences")


# =================================================================
# DEEPER ANALYSIS: what determines the cocycle?
# =================================================================
print(f"\n\n{'='*72}")
print("COCYCLE STRUCTURE: What determines w?")
print("="*72)

print("""
For a transitive triple (a,b,c) with a→b→c, a→c:
  w(a→b) + w(b→c) = w(a→c)

This means w is "additive along transitive paths":
  the weight of a shortcut = sum of weights of the two-step path.

Starting from any vertex v, the transitive triples from v determine
w for all edges from v, given w on one initial edge.

This is like a POTENTIAL function, but only along transitive paths.
On 3-cycles, the potential is not well-defined (it accumulates a
nonzero winding number).

The cocycle condition means: w is determined up to coboundaries
(vertex potentials) and one global parameter (the winding number
on 3-cycles).

CLAIM: For any tournament T, the space of cocycles modulo coboundaries
has dimension at most 1. This is because:

1. Fix a spanning tree of T (viewed as undirected K_n).
2. A coboundary δf is determined by f: V → R (n-dimensional).
   Modulo a global constant, coboundaries have n-1 parameters.
3. A cocycle w is determined by its values on non-tree edges
   (= fundamental cycle generators) BUT subject to the transitive
   triple constraints.
4. There are C(n,2)-n+1 = C(n-1,2) non-tree edges.
5. The transitive constraints propagate w values along transitive paths.
   How many degrees of freedom survive?
""")

# Count the effective degrees of freedom
for n in [4, 5, 6]:
    print(f"\nn={n}: C(n,2)={comb(n,2)} edges, n-1={n-1} tree edges, "
          f"C(n-1,2)={comb(n-1,2)} cycle edges")

    if n > 6:
        continue

    dof_counts = Counter()

    for A in all_tournaments(n):
        result = compute_cocycle_space(A, n)
        dof = result['dim_H1']
        dof_counts[dof] += 1

    print(f"  dim(H¹) distribution: {dict(sorted(dof_counts.items()))}")


# =================================================================
# THE KEY OBSERVATION: Transitive closure and connectivity
# =================================================================
print(f"\n\n{'='*72}")
print("TRANSITIVE CLOSURE AND THE COCYCLE CONSTRAINT GRAPH")
print("="*72)

print("""
Consider the constraint graph G_T:
  - Vertices = directed edges of T (the C(n,2) edges)
  - For each transitive triple (a,b,c), add a CONSTRAINT:
    w(a,b) + w(b,c) = w(a,c)

This means: knowing ANY TWO of {w(a,b), w(b,c), w(a,c)} determines
the third. So the constraint graph has "propagation paths."

Coboundaries use n-1 parameters (vertex potentials mod constant).
Cocycles use dim(cocycles) parameters.
H¹ = cocycles/coboundaries uses dim(H¹) parameters.

The transitive triple constraint connects THREE edge-variables.
The coboundary space provides n-1 free variables.

For dim(H¹) ≤ 1, we need:
  dim(cocycles) ≤ n = (n-1) + 1

Let's verify: dim(cocycles) ∈ {n-1, n} for all tournaments.
""")

for n in [3, 4, 5, 6]:
    print(f"\nn={n}:")
    dim_cocycles_dist = Counter()

    for A in all_tournaments(n):
        result = compute_cocycle_space(A, n)
        dim_cocycles_dist[result['dim_cocycles']] += 1

    print(f"  dim(cocycles) distribution: {dict(sorted(dim_cocycles_dist.items()))}")
    print(f"  n-1={n-1}, n={n}")


# =================================================================
# PROOF VIA COCYCLE DIMENSION
# =================================================================
print(f"\n\n{'='*72}")
print("PROOF APPROACH: dim(cocycles) ≤ n")
print("="*72)

print("""
THEOREM: For any tournament T on n vertices, dim(cocycles) ≤ n.

Since dim(coboundaries) = n-1 (n vertex potentials minus 1 constant),
this gives H¹ = dim(cocycles) - dim(coboundaries) ≤ n - (n-1) = 1.

PROOF SKETCH:
  The cocycle condition w(a,b)+w(b,c)=w(a,c) for transitive (a,b,c)
  can be rewritten as: for any vertex v and any two out-neighbors a,b
  of v with a→b: w(v,a) + w(a,b) = w(v,b).

  So: w(v,b) - w(v,a) = w(a,b) whenever v→a→b and v→b (transitive).

  This means: for vertex v, the values {w(v,a) : v→a} are constrained
  by the tournament structure among v's out-neighbors.

  Similarly for in-neighbors.

  KEY: For a vertex v with out-degree d⁺, the values w(v,a₁),...,w(v,a_{d⁺})
  are linked by the tournament on {a₁,...,a_{d⁺}}. How many independent
  parameters?

  If the tournament on out-neighbors is TRANSITIVE: w(v,a_i) are
  completely determined by ONE value (and the w(a_i,a_j) values).
  But the w(a_i,a_j) values are themselves part of the cocycle.

  The cocycle has C(n,2) variables and C(n,3)-t₃ constraints.
  dim(cocycles) = C(n,2) - rank(constraints).

  We need: rank(constraints) ≥ C(n,2) - n = C(n-1,2) - 1 + (n-1)
          Wait: C(n,2) - n = n(n-1)/2 - n = n(n-3)/2 for n≥3.
  Hmm, that doesn't simplify nicely.

  Actually: dim(cocycles) = C(n,2) - rank(constraint matrix).
  dim(coboundaries) = n-1.
  H¹ = dim(cocycles) - (n-1).
  H¹ ≤ 1 ⟺ dim(cocycles) ≤ n ⟺ rank(constraints) ≥ C(n,2)-n.

  Number of constraints = C(n,3) - t₃ (one per transitive triple).
  We need C(n,3) - t₃ constraints to have rank ≥ C(n,2) - n.

  For this to always hold: need C(n,3) - t₃ ≥ C(n,2) - n, i.e.,
  t₃ ≤ C(n,3) - C(n,2) + n = n(n-1)(n-2)/6 - n(n-1)/2 + n
     = n[(n-1)(n-2)/6 - (n-1)/2 + 1]
     = n[(n-1)(n-2) - 3(n-1) + 6] / 6
     = n[n²-3n+2-3n+3+6] / 6
     = n[n²-6n+11] / 6

  For n=5: 5[25-30+11]/6 = 5*6/6 = 5. Max t₃ at n=5 is 5. TIGHT!
  For n=7: 7[49-42+11]/6 = 7*18/6 = 21. Max t₃ at n=7 is 7*6*5/6-7*6/2+7=35-21+7=21. TIGHT!

  Hmm, is this always tight? The max t₃ for tournaments is C(n,3) - C(n,2) + n?
  No: max t₃ = C(n,3)/4 * (something)... Actually by the formula
  t₃ = C(n,3) - (1/2)Σs_i(n-1-s_i), where s_i are scores.
  Max t₃ occurs for regular tournaments: s_i=(n-1)/2 for all i.
  Then t₃ = C(n,3) - n*(n-1)/2*(n-1)/2/2 = C(n,3) - n(n-1)²/8.

  Wait, this is a different expression. Let me just check:
  n=5: C(5,3)=10, C(5,2)=10, n=5. C(n,3)-C(n,2)+n = 10-10+5=5. And max t₃=5. ✓
  n=3: C(3,3)=1, C(3,2)=3, n=3. 1-3+3=1. Max t₃=1. ✓
  n=7: C(7,3)=35, C(7,2)=21, n=7. 35-21+7=21. Max t₃? For regular n=7,
    t₃=C(7,3)-7*C(3,2)/4=35-7*3/... let me compute properly.
    s_i=3 for all i. Σs_i(6-s_i) = 7*3*3=63. t₃=35-63/2=35-31.5=3.5?
    No that can't be right. Let me use t₃ = C(n,3) - Σs_i*(s_i-1)/2/... actually
    t₃ = [C(n,2)*(n-2)/3 - Σ C(s_i,2)] / 2? No.

  OK let me just verify computationally.
""")

# VERIFY: dim(cocycles) ≤ n always
print("\n--- Verifying dim(cocycles) ≤ n ---")

for n in [3, 4, 5, 6]:
    max_dim = 0
    for A in all_tournaments(n):
        result = compute_cocycle_space(A, n)
        if result['dim_cocycles'] > max_dim:
            max_dim = result['dim_cocycles']

    print(f"n={n}: max dim(cocycles) = {max_dim}, n = {n}, "
          f"holds: {'YES' if max_dim <= n else 'NO'}")


# =================================================================
# CONSTRAINT RANK ANALYSIS
# =================================================================
print(f"\n\n{'='*72}")
print("CONSTRAINT RANK ANALYSIS")
print("="*72)

for n in [3, 4, 5, 6]:
    print(f"\nn={n}: C(n,2)={comb(n,2)}")

    rank_dist = Counter()
    for A in all_tournaments(n):
        tt = transitive_triples(A, n)
        el = edges_list(A, n)
        edge_idx = {e: i for i, e in enumerate(el)}
        ne = len(el)

        # Constraint matrix
        C = np.zeros((len(tt), ne))
        for j, (a,b,c) in enumerate(tt):
            C[j, edge_idx[(a,b)]] = 1
            C[j, edge_idx[(b,c)]] = 1
            C[j, edge_idx[(a,c)]] = -1

        if C.shape[0] > 0:
            r = np.linalg.matrix_rank(C, tol=1e-10)
        else:
            r = 0

        rank_dist[r] += 1

    target = comb(n,2) - n
    print(f"  Target rank ≥ {target} (= C(n,2)-n)")
    print(f"  Constraint rank distribution: {dict(sorted(rank_dist.items()))}")

    min_rank = min(rank_dist.keys())
    print(f"  Min rank = {min_rank}, target = {target}: "
          f"{'MET' if min_rank >= target else 'NOT MET'}")


# =================================================================
# EXPLICIT COCYCLE CONSTRUCTION
# =================================================================
print(f"\n\n{'='*72}")
print("EXPLICIT COCYCLE FOR β₁=1 TOURNAMENTS")
print("="*72)

print("""
The cocycle w: edges → Z satisfies:
  w(a,b) + w(b,c) = w(a,c) for transitive triples.

For the regular tournament on 5 vertices (C₅ with adjacency 1,2):
  w values were: some edges get -1, some get -2.

Let me check: is w related to the "winding number" around 3-cycles?
Define w(a→b) = (number of directed 3-cycles containing edge a→b).

Check cocycle condition: for transitive (a,b,c):
  #3cycles through (a,b) + #3cycles through (b,c) = #3cycles through (a,c)?

This would be amazing if true...
""")

for n in [5]:
    for A in all_tournaments(n):
        t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))
        if t3 != 5:
            continue

        el = edges_list(A, n)
        edge_idx = {e: i for i, e in enumerate(el)}
        ne = len(el)
        tt = transitive_triples(A, n)

        # Count 3-cycles through each edge
        cycle_count = np.zeros(ne)
        for i in range(n):
            for j in range(n):
                if j == i or not A[i][j]: continue
                for k in range(n):
                    if k == i or k == j: continue
                    if A[j][k] and A[k][i]:
                        cycle_count[edge_idx[(i,j)]] += 1
                        cycle_count[edge_idx[(j,k)]] += 1
                        cycle_count[edge_idx[(k,i)]] += 1

        # Each cycle counted 3 times (once per edge), but each edge
        # counted once per cycle through it. Divide by 3? No.
        # cycle_count[e] = number of directed 3-cycles containing edge e.

        # Check cocycle condition
        is_cocycle = True
        for (a,b,c) in tt:
            lhs = cycle_count[edge_idx[(a,b)]] + cycle_count[edge_idx[(b,c)]]
            rhs = cycle_count[edge_idx[(a,c)]]
            if abs(lhs - rhs) > 1e-10:
                is_cocycle = False
                break

        print(f"\nRegular tournament n=5 (t₃=5):")
        print(f"  3-cycle count per edge: {cycle_count.astype(int).tolist()}")
        print(f"  Is 3-cycle count a cocycle? {is_cocycle}")

        if not is_cocycle:
            # Try: w(a→b) = #3cycles through (a→b) minus average
            avg = np.mean(cycle_count)
            w_centered = cycle_count - avg

            is_cocycle2 = True
            for (a,b,c) in tt:
                lhs = w_centered[edge_idx[(a,b)]] + w_centered[edge_idx[(b,c)]]
                rhs = w_centered[edge_idx[(a,c)]]
                if abs(lhs - rhs) > 1e-10:
                    is_cocycle2 = False
                    break
            print(f"  Is centered 3-cycle count a cocycle? {is_cocycle2}")

            # Try: w(a→b) = -(out_degree(a) + out_degree(b))
            scores = [sum(A[v]) for v in range(n)]
            w_score = np.array([-(scores[a]+scores[b]) for a,b in el], dtype=float)

            is_cocycle3 = True
            for (a,b,c) in tt:
                lhs = w_score[edge_idx[(a,b)]] + w_score[edge_idx[(b,c)]]
                rhs = w_score[edge_idx[(a,c)]]
                if abs(lhs - rhs) > 1e-10:
                    is_cocycle3 = False
                    break
            print(f"  Is -(s_a+s_b) a cocycle? {is_cocycle3}")

        # Print the actual H¹ cocycle
        result = compute_cocycle_space(A, n)
        w = result['representative']
        print(f"\n  Actual H¹ cocycle:")
        for i in range(ne):
            a, b = el[i]
            s_a = scores[a] if 'scores' in dir() else sum(A[a])
            s_b = scores[b] if 'scores' in dir() else sum(A[b])
            cc = int(cycle_count[i])
            print(f"    ({a}→{b}): w={w[i]:+d}, s_a={s_a}, s_b={s_b}, #3cyc={cc}")

        break


# =================================================================
# DOES THE COCYCLE HAVE A UNIVERSAL FORMULA?
# =================================================================
print(f"\n\n{'='*72}")
print("UNIVERSAL COCYCLE FORMULA SEARCH")
print("="*72)

n = 5
print(f"\nSearching for universal cocycle formula at n={n}...")

# For each β₁=1 tournament, compute the H¹ cocycle (up to scaling)
# and try to express it in terms of tournament invariants

formulas_tested = {
    "s_a": lambda A, n, a, b: sum(A[a]),
    "s_b": lambda A, n, a, b: sum(A[b]),
    "s_a-s_b": lambda A, n, a, b: sum(A[a]) - sum(A[b]),
    "s_a+s_b": lambda A, n, a, b: sum(A[a]) + sum(A[b]),
    "s_a*s_b": lambda A, n, a, b: sum(A[a]) * sum(A[b]),
    "#common_out": lambda A, n, a, b: sum(1 for k in range(n) if k!=a and k!=b and A[a][k] and A[b][k]),
    "#common_in": lambda A, n, a, b: sum(1 for k in range(n) if k!=a and k!=b and A[k][a] and A[k][b]),
    "s_a^2": lambda A, n, a, b: sum(A[a])**2,
    "s_b^2": lambda A, n, a, b: sum(A[b])**2,
    "s_a^2-s_b^2": lambda A, n, a, b: sum(A[a])**2 - sum(A[b])**2,
}

for name, f in formulas_tested.items():
    # Check if f defines a cocycle for ALL β₁=1 tournaments
    cocycle_count = 0
    exact_match = 0
    total_b1 = 0

    for A in all_tournaments(n):
        result = compute_cocycle_space(A, n)
        if result['dim_H1'] != 1:
            continue
        total_b1 += 1

        el = result['edges']
        edge_idx = result['edge_idx']
        ne = len(el)
        tt = transitive_triples(A, n)

        # Compute f-based weight
        w_f = np.array([f(A, n, a, b) for a, b in el], dtype=float)

        # Check cocycle
        is_coc = True
        for (a,b,c) in tt:
            if abs(w_f[edge_idx[(a,b)]] + w_f[edge_idx[(b,c)]] - w_f[edge_idx[(a,c)]]) > 1e-10:
                is_coc = False
                break

        if is_coc:
            cocycle_count += 1

            # Check if it's proportional to the H¹ representative
            w_true = result['representative'].astype(float)
            # Check proportionality
            nz_true = np.where(w_true != 0)[0]
            if len(nz_true) > 0:
                ratio = w_f[nz_true[0]] / w_true[nz_true[0]] if w_true[nz_true[0]] != 0 else None
                if ratio is not None and all(abs(w_f[i] - ratio * w_true[i]) < 1e-10 for i in range(ne)):
                    exact_match += 1

    if cocycle_count > 0:
        print(f"  {name}: cocycle for {cocycle_count}/{total_b1} β₁=1 tournaments, "
              f"exact match: {exact_match}/{total_b1}")


# =================================================================
# TRY A RECURSIVE/STRUCTURAL COCYCLE
# =================================================================
print(f"\n\n{'='*72}")
print("STRUCTURAL ANALYSIS: The 'flow resistance' cocycle")
print("="*72)

print("""
The cocycle w satisfies: w(a→b) + w(b→c) = w(a→c) on transitive triples.
This is the SAME as saying: w restricted to any transitive sub-tournament
is a coboundary (potential difference).

On 3-cycles, w has a nonzero "winding number" = w(a→b)+w(b→c)+w(c→a).
This winding number is the same for ALL 3-cycles (since H¹ = 1-dim).

So the cocycle measures: "how much does a path wind around the
tournament's cycle structure?"

For a tournament-theoretic interpretation: w is a "potential"
that is additive on transitive paths but accumulates ±1 winding
on each 3-cycle.

CLAIM: w(a→b) = (# directed paths from a to b using only vertices
that beat both a and b) - (some correction term).

Or more simply: w might be related to the graph distance in the
"dominance graph" of T.
""")

# Compute the actual cocycle for several tournaments and look for pattern
print("\n--- Cocycle values vs tournament structure (n=5) ---")

n = 5
for idx_A, A in enumerate(all_tournaments(n)):
    result = compute_cocycle_space(A, n)
    if result['dim_H1'] != 1:
        continue

    w = result['representative']
    el = result['edges']
    ne = len(el)

    # Print adjacency matrix and cocycle
    scores = [sum(A[i]) for i in range(n)]
    t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))

    # For each edge, compute: shared outneighborhood size
    print(f"\nTournament #{idx_A} (scores={scores}, t₃={t3}):")

    # Create a matrix of w values
    W = np.zeros((n,n), dtype=int)
    for i in range(ne):
        a, b = el[i]
        W[a,b] = w[i]

    print(f"  Cocycle matrix W[a][b] (w(a→b), 0 if b→a):")
    for i in range(n):
        row = [f"{W[i,j]:+d}" if A[i][j] else "  ." for j in range(n)]
        print(f"    {i}: {' '.join(row)}")

    # Check: does W[a][b] = -(some function of position)?
    # In the dominance order: vertex with score s is at "level" s.
    # If a→b: the "cost" of the arc might depend on relative scores.

    # Check: is W[a][b] = -(s_a + s_b - (n-1))?
    # For regular: s_a=s_b=2, so -(2+2-4)=0. But we see nonzero values.

    # Check if W is antisymmetric: W[a][b] + W[b][a] = const?
    print(f"  W[a][b] + W[b][a] for each edge pair:")
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j]:
                sums = W[i,j] + W[j,i] if A[j][i] else W[i,j]
                # Only one direction exists
                pass

    # Actually in tournament, exactly one of (i,j) or (j,i) is an edge.
    # So W[i][j] and W[j][i] can't both be nonzero (one is always 0).

    # What about W[i][j] for i→j vs |W[j][i]| = 0: the cocycle is defined
    # only on existing edges. So no antisymmetry to check.

    # Try: define f(v) = sum of w(v→u) for all out-neighbors u of v
    for v in range(n):
        out_sum = sum(W[v,u] for u in range(n) if A[v][u])
        in_sum = sum(W[u,v] for u in range(n) if A[u][v])
        print(f"    v={v}: Σw(v→·)={out_sum}, Σw(·→v)={in_sum}, total={out_sum+in_sum}")

    break  # Just first example


# =================================================================
# THE DEFINITIVE APPROACH: dim(cocycles) = C(n,2) - rank(C)
# =================================================================
print(f"\n\n{'='*72}")
print("DEFINITIVE: Why rank(constraints) ≥ C(n,2) - n?")
print("="*72)

print("""
The constraint matrix C has:
  - One row per transitive triple (a,b,c): [1 at (a,b), 1 at (b,c), -1 at (a,c)]
  - C(n,2) columns (one per directed edge)

We need rank(C) ≥ C(n,2) - n.

OBSERVATION: Each constraint w(a,b)+w(b,c)=w(a,c) says that knowing
2 of {w(a,b), w(b,c), w(a,c)} determines the third.

Think of the CONSTRAINT GRAPH on edges:
  - Vertices = directed edges of T
  - Connect (a,b)—(b,c) and (a,b)—(a,c) and (b,c)—(a,c) for each
    transitive triple.

Connected components of this graph determine the rank of C.
If the graph has k components, then rank(C) = C(n,2) - k.
So dim(cocycles) = k.
And H¹ = k - (n-1).

For H¹ ≤ 1: need k ≤ n.

The constraint graph on edges has at most n connected components
because... let's verify.
""")

# Build constraint graph and count components
for n in [3, 4, 5, 6]:
    print(f"\nn={n}:")
    comp_dist = Counter()

    for A in all_tournaments(n):
        el = edges_list(A, n)
        ne = len(el)
        edge_idx = {e: i for i, e in enumerate(el)}
        tt = transitive_triples(A, n)

        # Build adjacency on edges
        adj = defaultdict(set)
        for (a,b,c) in tt:
            e1 = edge_idx[(a,b)]
            e2 = edge_idx[(b,c)]
            e3 = edge_idx[(a,c)]
            adj[e1].add(e2)
            adj[e1].add(e3)
            adj[e2].add(e1)
            adj[e2].add(e3)
            adj[e3].add(e1)
            adj[e3].add(e2)

        # BFS components
        visited = set()
        num_components = 0
        for start in range(ne):
            if start in visited:
                continue
            num_components += 1
            queue = [start]
            visited.add(start)
            while queue:
                node = queue.pop(0)
                for nbr in adj[node]:
                    if nbr not in visited:
                        visited.add(nbr)
                        queue.append(nbr)

        comp_dist[num_components] += 1

    print(f"  Components of constraint graph: {dict(sorted(comp_dist.items()))}")
    print(f"  n = {n}: H¹ = components - (n-1)")
    for k, cnt in sorted(comp_dist.items()):
        h1 = k - (n-1)
        print(f"    {k} components → H¹ = {h1}: {cnt} tournaments")


# =================================================================
# WAIT: Is constraint graph components = dim(cocycles)?
# =================================================================
print(f"\n\n{'='*72}")
print("VERIFYING: components = dim(cocycles)?")
print("="*72)

for n in [3, 4, 5]:
    fail = 0
    total = 0
    for A in all_tournaments(n):
        el = edges_list(A, n)
        ne = len(el)
        edge_idx = {e: i for i, e in enumerate(el)}
        tt = transitive_triples(A, n)

        # Components
        adj = defaultdict(set)
        for (a,b,c) in tt:
            e1, e2, e3 = edge_idx[(a,b)], edge_idx[(b,c)], edge_idx[(a,c)]
            for x in [e1,e2,e3]:
                for y in [e1,e2,e3]:
                    if x != y:
                        adj[x].add(y)
        visited = set()
        num_comp = 0
        for s in range(ne):
            if s in visited: continue
            num_comp += 1
            q = [s]
            visited.add(s)
            while q:
                node = q.pop(0)
                for nbr in adj[node]:
                    if nbr not in visited:
                        visited.add(nbr)
                        q.append(nbr)

        result = compute_cocycle_space(A, n)
        total += 1
        if num_comp != result['dim_cocycles']:
            fail += 1
            if fail <= 3:
                print(f"  n={n}: MISMATCH: {num_comp} components vs {result['dim_cocycles']} dim(cocycles)")

    print(f"n={n}: {total} tournaments, {fail} mismatches")


print(f"\n\n{'='*72}")
print("CONCLUSION")
print("="*72)

print("""
THEOREM: β₁(T) ≤ 1 for any tournament T on n vertices.

PROOF:

By GLMY duality, β₁ = dim(H₁) = dim(H¹) where H¹ = cocycles/coboundaries.

Step 1: dim(coboundaries) = n-1.
  A coboundary is δf(a→b) = f(b)-f(a) for some f: V→R.
  The space of coboundaries has dimension n-1 (n vertex potentials
  minus the constant kernel).

Step 2: dim(cocycles) ≤ n.
  A cocycle w satisfies w(a,b)+w(b,c)=w(a,c) for every transitive triple.
  Define the CONSTRAINT GRAPH G on the edge set:
    - Vertices of G = directed edges of T (there are C(n,2) of them)
    - Edges of G: connect (a,b), (b,c), (a,c) for each transitive triple
  dim(cocycles) = number of connected components of G.

  CLAIM: G has at most n connected components.

  PROOF OF CLAIM: For each vertex v of T, all outgoing edges from v
  are in the SAME component of G. This is because for any two out-neighbors
  a, b of v with a→b (which exists since T is a tournament): the triple
  (v,a,b) is transitive, linking edges (v,a), (a,b), (v,b).

  Similarly, all incoming edges to v are in the same component (for any
  two in-neighbors a, b of v with a→b: (a,b,v) is transitive if b→v and
  a→v, linking (a,b), (b,v), (a,v)).

  So: for each vertex v, all edges incident to v are in ≤ 2 components
  (one for outgoing, one for incoming). But an outgoing edge (v,u) is
  also an incoming edge of u, so these components overlap.

  Actually: for vertex v with out-neighbors {a₁,...,a_d}, the edges
  (v,a₁),...,(v,a_d) are all linked via transitive triples.

  Now consider vertex u with out-neighbors including v. Then (u,v) is
  an edge, and (u,v) is linked to (u,x) for any x ∈ out(u) with v→x
  (via transitive (u,v,x)). And (u,v) is linked to (v,a₁) if (u,v,a₁)
  is transitive (need u→a₁, which holds iff a₁ ∈ out(u)).

  The full tournament T is strongly connected if it's not transitive.
  For transitive T: β₁=0 (no 3-cycles, dim(cocycles)=n-1, H¹=0).
  For non-transitive T: the constraint graph is more connected.

  CLAIM PROOF (cleaner): Consider any edge (a,b) where a→b.
  For any vertex c ≠ a,b, either:
    - a→c, b→c: transitive (a,b,c) links (a,b)↔(a,c)↔(b,c)
    - a→c, c→b: transitive (a,c,b)... no, need c→b and a→b, but (a,c,b)
      needs a→c, c→b, a→b. That's transitive! Links (a,c)↔(c,b)↔(a,b).
    - c→a, b→c: transitive (b,c,a)... need b→c, c→a, b→a.
      But a→b, so b→a is false. NOT transitive.
      Try (c,a,b): c→a, a→b, c→b? Need c→b. But b→c. NOT.
      Hmm: c→a, b→c, a→b. Triangle (a,b,c): a→b, b→c, c→a — 3-cycle!
    - c→a, c→b: transitive (c,a,b): c→a, a→b, c→b. YES!
      Links (c,a)↔(a,b)↔(c,b).

  So for edge (a,b), vertex c gives:
    - a→c, b→c: (a,b) linked to (a,c) and (b,c)
    - a→c, c→b: (a,b) linked to (a,c) and (c,b)
    - c→a, b→c: NO LINK (3-cycle!)
    - c→a, c→b: (a,b) linked to (c,a) and (c,b)

  In ALL cases except the 3-cycle case, (a,b) gets linked to 2 other edges.
  In the 3-cycle case, (a,b) is NOT linked to any new edge via vertex c.

  For (a,b) to be isolated: ALL other vertices c must form 3-cycles
  with {a,b}. But that's impossible for n ≥ 4: if c forms 3-cycle
  (a,b,c) and d forms 3-cycle (a,b,d), then consider {c,d}:
  c→a, b→c and d→a, b→d. So c,d ∈ in(a) ∩ out(b).
  Arc between c and d: either c→d or d→c.
  If c→d: is (b,c,d) transitive? b→c, c→d, b→d. Need b→d: YES.
    So (b,c,d) transitive, linking (b,c)↔(c,d)↔(b,d).
    Now (a,b,c) is a 3-cycle: edges (a,b),(b,c),(c,a).
    And (b,c) is linked to (a,b) via... wait, (a,b,c) is a 3-cycle,
    so NO transitive triple involving all three.
    But (b,c,d) IS transitive: links (b,c) to (b,d) to (c,d).
    Does (b,c) get linked to (a,b)? Not directly via 3-cycle {a,b,c}.
    But via d: (a,b,d) is also a 3-cycle (assumption). So no link there.

    We need another vertex e to create transitive triples involving (a,b).

  This is getting complex. Let me just verify the component count bound.
""")

# DEFINITIVE CHECK: Is #components ≤ n always?
print("\n--- DEFINITIVE: Is #components ≤ n always? ---")
for n in [3, 4, 5, 6, 7]:
    if n > 6:
        # Sample
        import random
        max_comp = 0
        for trial in range(1000):
            A = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1

            el = edges_list(A, n)
            ne = len(el)
            edge_idx = {e: i for i, e in enumerate(el)}
            tt = transitive_triples(A, n)
            adj = defaultdict(set)
            for (a,b,c) in tt:
                e1, e2, e3 = edge_idx[(a,b)], edge_idx[(b,c)], edge_idx[(a,c)]
                for x in [e1,e2,e3]:
                    for y in [e1,e2,e3]:
                        if x != y: adj[x].add(y)
            visited = set()
            num_comp = 0
            for s in range(ne):
                if s in visited: continue
                num_comp += 1
                q = [s]
                visited.add(s)
                while q:
                    node = q.pop(0)
                    for nbr in adj[node]:
                        if nbr not in visited:
                            visited.add(nbr)
                            q.append(nbr)
            max_comp = max(max_comp, num_comp)

        print(f"n={n} (1000 random): max components = {max_comp}, n = {n}: "
              f"{'YES' if max_comp <= n else 'NO'}")
    else:
        max_comp = 0
        for A in all_tournaments(n):
            el = edges_list(A, n)
            ne = len(el)
            edge_idx = {e: i for i, e in enumerate(el)}
            tt = transitive_triples(A, n)
            adj = defaultdict(set)
            for (a,b,c) in tt:
                e1, e2, e3 = edge_idx[(a,b)], edge_idx[(b,c)], edge_idx[(a,c)]
                for x in [e1,e2,e3]:
                    for y in [e1,e2,e3]:
                        if x != y: adj[x].add(y)
            visited = set()
            num_comp = 0
            for s in range(ne):
                if s in visited: continue
                num_comp += 1
                q = [s]
                visited.add(s)
                while q:
                    node = q.pop(0)
                    for nbr in adj[node]:
                        if nbr not in visited:
                            visited.add(nbr)
                            q.append(nbr)
            max_comp = max(max_comp, num_comp)

        print(f"n={n} (exhaustive): max components = {max_comp}, n = {n}: "
              f"{'YES' if max_comp <= n else 'NO'}")


print("\nDone.")
