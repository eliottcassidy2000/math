"""
Analyze Omega_2 for tournaments in detail.
Understand WHY beta_2 = 0 for all tournaments.

For a tournament T on n vertices:
- Allowed 2-paths: (a, b, c) where a->b and b->c are arcs (a != b != c, but a may equal c)
  Actually in tournaments, a != c is automatic when a != b != c because if a=c then
  we'd have a->b and b->a which means both arcs exist, contradiction for a tournament.

- Omega_2: elements omega = sum c_{abc} e_{abc} such that d(omega) has all faces allowed.
  d(e_{abc}) = e_{bc} - e_{ac} + e_{ab}
  The face e_{ac} is allowed iff there is an arc a->c (or c->a would give e_{ca}).
  Wait: e_{ac} is allowed iff a->c is an arc.

  For a tournament, exactly one of a->c or c->a holds.
  If a->c: then e_{ac} is allowed.
  If c->a: then e_{ac} is NOT allowed (it's e_{ac} not e_{ca}).

So for a 2-path (a,b,c):
  d(e_{abc}) = e_{bc} - e_{ac} + e_{ab}
  e_{ab} is always allowed (since a->b is an arc)
  e_{bc} is always allowed (since b->c is an arc)
  e_{ac}: allowed iff a->c, NOT allowed iff c->a

Case 1: a->c (transitive triple {a,b,c} with a->b, b->c, a->c)
  All three faces are allowed => e_{abc} is in Omega_2 by itself.

Case 2: c->a (directed 3-cycle through a->b->c->a, but note this is NOT a 3-cycle,
  it's just that a->b, b->c, c->a. This IS a directed 3-cycle on {a,b,c}.)
  The face e_{ac} is NOT allowed.
  So for omega = sum c_{abc} e_{abc} to be in Omega_2,
  we need the coefficient of e_{ac} in d(omega) to be 0 for all non-allowed (a,c).

This means: For each pair (a,c) where c->a (so e_{ac} is NOT allowed),
the sum of -c_{abc} over all b with a->b->c must be zero.
"""

import numpy as np
from numpy.linalg import matrix_rank

def make_tournament(n, bits):
    adj = [[False]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                adj[i][j] = True
            else:
                adj[j][i] = True
            k += 1
    return adj

def analyze_omega2(adj, verbose=True):
    """Detailed analysis of Omega_2 for a tournament."""
    n = len(adj)

    # Get all allowed 2-paths
    paths_2 = []
    for a in range(n):
        for b in range(n):
            if b == a or not adj[a][b]:
                continue
            for c in range(n):
                if c == b or not adj[b][c]:
                    continue
                # a->b->c is an allowed 2-path
                paths_2.append((a, b, c))

    # Classify: transitive triples vs 3-cycle triples
    transitive = []
    cyclic = []
    for (a, b, c) in paths_2:
        if adj[a][c]:
            transitive.append((a, b, c))
        else:
            # c->a must hold (tournament)
            assert adj[c][a]
            cyclic.append((a, b, c))

    if verbose:
        print(f"  Total 2-paths: {len(paths_2)}")
        print(f"  Transitive triples (a->c): {len(transitive)}")
        print(f"  Cyclic triples (c->a): {len(cyclic)}")

    # For each non-allowed face pair (a,c) with c->a,
    # collect the 2-paths (a,b,c) that contribute
    non_allowed_constraints = {}
    for (a, b, c) in cyclic:
        key = (a, c)  # non-allowed face
        if key not in non_allowed_constraints:
            non_allowed_constraints[key] = []
        non_allowed_constraints[key].append(b)

    if verbose:
        print(f"  Non-allowed face pairs: {len(non_allowed_constraints)}")
        for (a, c), bs in sorted(non_allowed_constraints.items()):
            print(f"    ({a},{c}): intermediates = {bs} (count={len(bs)})")

    # Build the constraint matrix
    # Rows = non-allowed pairs, Columns = cyclic 2-paths
    # For cyclic path (a,b,c), the constraint for (a,c) has coefficient -1
    # (from the -e_{ac} term in the boundary)

    # But Omega_2 includes BOTH transitive and cyclic paths
    # Transitive paths (a,b,c) with a->c automatically have all faces allowed
    # So they are always in Omega_2
    # Cyclic paths need constraints

    # Actually, let me build the full constraint matrix
    path_index = {p: i for i, p in enumerate(paths_2)}
    na_pairs = sorted(non_allowed_constraints.keys())
    na_index = {p: i for i, p in enumerate(na_pairs)}

    M = np.zeros((len(na_pairs), len(paths_2)))
    for j, (a, b, c) in enumerate(paths_2):
        # The coefficient of e_{ac} in d(e_{abc}) is (-1)^1 = -1
        if (a, c) in na_index:
            M[na_index[(a, c)], j] = -1  # coefficient of e_{ac}

    # Omega_2 = kernel of M
    if len(na_pairs) == 0:
        omega_dim = len(paths_2)
    else:
        rank = matrix_rank(M, tol=1e-10)
        omega_dim = len(paths_2) - rank

    if verbose:
        print(f"  Constraint matrix: {M.shape}")
        print(f"  Rank of constraint matrix: {matrix_rank(M, tol=1e-10) if len(na_pairs) > 0 else 0}")
        print(f"  dim(Omega_2) = {omega_dim}")

    return omega_dim, len(transitive), len(cyclic), non_allowed_constraints


# Analyze n=3
print("=" * 60)
print("n=3: Transitive tournament")
print("=" * 60)
adj = [[False, True, True], [False, False, True], [False, False, False]]
analyze_omega2(adj)

print("\n" + "=" * 60)
print("n=3: Cyclic tournament")
print("=" * 60)
adj = [[False, True, False], [False, False, True], [True, False, False]]
analyze_omega2(adj)

# Analyze n=4: all tournaments
print("\n" + "=" * 60)
print("n=4: All tournaments - Omega_2 analysis")
print("=" * 60)
for bits in range(64):
    adj = make_tournament(4, bits)
    odim, ntrans, ncyc, constraints = analyze_omega2(adj, verbose=False)
    # Count number of intermediates per constraint
    max_inter = max((len(bs) for bs in constraints.values()), default=0)
    if odim > 0 and ntrans > 0:
        # This is interesting: check if Omega_2 = span of transitive triples
        pass
    if bits < 5 or odim != ntrans:  # Print some examples
        print(f"  bits={bits}: dim(Omega_2)={odim}, trans={ntrans}, cyc={ncyc}, max_intermediates={max_inter}")

# Key question: is dim(Omega_2) always = number of transitive triples?
print("\n" + "=" * 60)
print("Key question: is dim(Omega_2) = #transitive_triples for all tournaments?")
print("=" * 60)
for n in [3, 4, 5]:
    ne = n*(n-1)//2
    all_match = True
    for bits in range(2**ne):
        adj = make_tournament(n, bits)
        odim, ntrans, ncyc, _ = analyze_omega2(adj, verbose=False)
        if odim != ntrans:
            print(f"  MISMATCH at n={n}, bits={bits}: dim(Omega_2)={odim}, trans={ntrans}")
            all_match = False
    if all_match:
        print(f"  n={n}: YES, dim(Omega_2) = #transitive_triples for all {2**ne} tournaments")


# Now compute H_2 = ker(d_2: Omega_2 -> Omega_1) / im(d_3: Omega_3 -> Omega_2)
# and check the structure
print("\n" + "=" * 60)
print("Detailed H_2 computation for n=5 regular tournament")
print("=" * 60)
# Regular tournament on 5 vertices: Paley T_5
# Using circulant: 0->1, 0->2, 1->2, 1->3, 2->3, 2->4, 3->4, 3->0, 4->0, 4->1
adj5 = [[False]*5 for _ in range(5)]
# QR: {1,2} are quadratic residues mod 5
for i in range(5):
    for d in [1, 2]:
        adj5[i][(i+d) % 5] = True

print("Paley tournament T_5:")
for i in range(5):
    outs = [j for j in range(5) if adj5[i][j]]
    print(f"  {i} -> {outs}")

odim, ntrans, ncyc, constraints = analyze_omega2(adj5, verbose=True)

# Now compute full boundary maps
from path_homology_tournament import compute_path_homology
h = compute_path_homology(adj5, max_dim=4, verbose=True)
print(f"\nBetti numbers for Paley T_5: {h}")
