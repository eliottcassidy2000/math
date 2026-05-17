"""
Computational experiments for the beta1*beta3=0 seesaw conjecture.

Key formula (from THM-290, LES argument):
  When beta3(T\v) = 0 for all v, we have:
    H3(T) ≅ H2tilde(I(H_v))   for EVERY vertex v

where H_v is the bipartite graph:
  Left = N+(v), Right = N-(v)
  Edge {u,w}: u in N+(v), w in N-(v), u->w in T

Questions to test:
  Q1: Is H2tilde(I(H_v)) = 0 for all v when beta1(T)=1? (would prove seesaw)
  Q2: Does H2tilde(I(H_v)) = beta3(T) for all v? (validates the formula)
  Q3: When beta3(T)>0, is H2tilde(I(H_v)) > 0 for all v?
  Q4: Do cone points (isolated vertices in H_v) always exist when beta1(T)=1?
"""

import sys
import itertools
import numpy as np
from collections import defaultdict

# Try importing existing path homology module
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')

# --- Tournament representation ---

def all_tournaments(n):
    """Generate all tournaments on vertices 0..n-1 as adjacency sets."""
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(pairs)):
        arc_set = set()
        for k, (i,j) in enumerate(pairs):
            if bits & (1 << k):
                arc_set.add((i,j))
            else:
                arc_set.add((j,i))
        yield arc_set

def beats(arc_set, u, v):
    return (u,v) in arc_set

def out_nbrs(arc_set, n, v):
    return [u for u in range(n) if u != v and beats(arc_set, v, u)]

def in_nbrs(arc_set, n, v):
    return [u for u in range(n) if u != v and beats(arc_set, u, v)]

# --- Independence complex and its homology ---

def build_independence_complex(graph_edges, vertices):
    """
    Build the independence complex of an undirected graph.
    graph_edges: set of frozensets {u,v}
    vertices: list of vertices
    Returns: list of maximal independent sets (as frozensets)
    """
    adj = defaultdict(set)
    for e in graph_edges:
        u, v = sorted(e)
        adj[u].add(v)
        adj[v].add(u)

    def is_independent(S):
        S = list(S)
        for i in range(len(S)):
            for j in range(i+1, len(S)):
                if S[j] in adj[S[i]]:
                    return False
        return True

    # All faces (independent sets)
    faces = []
    for k in range(len(vertices)+1):
        for S in itertools.combinations(vertices, k):
            if is_independent(S):
                faces.append(frozenset(S))
    return faces

def compute_homology_of_complex(faces):
    """
    Compute reduced homology H_k(complex) using integer matrices.
    Returns dict: k -> rank of H_k
    """
    if not faces:
        return {}

    max_dim = max(len(f)-1 for f in faces) if faces else -1

    # Organize by dimension
    by_dim = defaultdict(list)
    for f in faces:
        by_dim[len(f)-1].append(tuple(sorted(f)))

    # Sort for consistent indexing
    for d in by_dim:
        by_dim[d].sort()

    def boundary_matrix(dim):
        """Build boundary matrix partial_dim: C_dim -> C_{dim-1}"""
        if dim <= 0 or dim-1 not in by_dim:
            return np.zeros((len(by_dim.get(dim-1,[])), len(by_dim.get(dim,[]))), dtype=int)

        rows = by_dim[dim-1]  # (dim-1)-faces
        cols = by_dim[dim]    # dim-faces

        row_idx = {f: i for i, f in enumerate(rows)}
        M = np.zeros((len(rows), len(cols)), dtype=int)

        for j, sigma in enumerate(cols):
            for k in range(len(sigma)):
                face = sigma[:k] + sigma[k+1:]
                if face in row_idx:
                    M[row_idx[face], j] = (-1)**k

        return M

    def rank_over_integers(M):
        """Compute rank of integer matrix via Smith normal form (approximate via QR)."""
        if M.size == 0:
            return 0
        # Use SVD for rank (works for integer matrices)
        try:
            _, s, _ = np.linalg.svd(M.astype(float))
            return int(np.sum(s > 1e-9))
        except:
            return np.linalg.matrix_rank(M.astype(float))

    result = {}
    dims = sorted(by_dim.keys())

    # For each dimension k, compute H_k = ker(d_k) / im(d_{k+1})
    # Note: for REDUCED homology, we augment with a (-1)-dim face

    for k in range(max_dim + 1):
        C_k = by_dim.get(k, [])
        C_km1 = by_dim.get(k-1, []) if k > 0 else ['empty']  # augmented
        C_kp1 = by_dim.get(k+1, [])

        nk = len(C_k)
        if nk == 0:
            result[k] = 0
            continue

        # Boundary map d_k: C_k -> C_{k-1}
        dk = boundary_matrix(k)
        # Boundary map d_{k+1}: C_{k+1} -> C_k
        dkp1 = boundary_matrix(k+1)

        rank_dk = rank_over_integers(dk)
        rank_dkp1 = rank_over_integers(dkp1)

        dim_ker_dk = nk - rank_dk
        dim_im_dkp1 = rank_dkp1

        h_k = dim_ker_dk - dim_im_dkp1
        result[k] = max(0, h_k)

    # Reduced H_{-1}: for non-empty complex, H_{-1}=0; for empty complex, H_{-1}=Z
    # H_0 reduced: H_0 - 1 (subtract 1 for the augmentation)
    if 0 in result:
        # count connected components via augmented boundary
        C_0 = by_dim.get(0, [])
        C_m1 = [()]  # one (-1)-face
        if C_0:
            # d_0: C_0 -> C_{-1}: every vertex maps to the (-1)-face
            rank_d0 = 1 if C_0 else 0  # rank is 1 (all map to the one (-1)-face)
            # ker(d_0) = kernel = 0 for non-empty (actually ker has dim n_0 - 1)
            n_0 = len(C_0)
            dim_ker_d0 = n_0 - rank_d0  # = n_0 - 1
            dim_im_d1 = rank_over_integers(boundary_matrix(1))
            result[-1] = 0  # reduced H_{-1} = 0 for non-empty complex
            result[0] = dim_ker_d0 - dim_im_d1  # = (components - 1)

    return result

# --- Path homology for tournaments ---

def compute_allowed_paths(arc_set, n, max_p=4):
    """
    Compute Omega_p for p=0,1,2,3 (allowed paths of length p).
    Omega_p = set of (v0,...,vp) with vi->vj for ALL i<j (transitive orderings).
    """
    allowed = defaultdict(set)
    # p=0: all vertices
    for v in range(n):
        allowed[0].add((v,))
    # p>=1: extend by one vertex at a time
    for p in range(1, max_p+1):
        for path in allowed[p-1]:
            last = path[-1]
            for v in range(n):
                if v in path:
                    continue
                # v must be beaten by all vertices in path
                if all(beats(arc_set, path[i], v) for i in range(len(path))):
                    allowed[p].add(path + (v,))
    return allowed

def build_boundary_matrix(allowed_p, allowed_pm1):
    """Build boundary matrix d_p: C_p -> C_{p-1}"""
    paths_p = sorted(allowed_p)
    paths_pm1 = sorted(allowed_pm1)

    if not paths_p or not paths_pm1:
        return np.zeros((len(paths_pm1), len(paths_p)), dtype=int), paths_p, paths_pm1

    row_idx = {path: i for i, path in enumerate(paths_pm1)}
    M = np.zeros((len(paths_pm1), len(paths_p)), dtype=int)

    for j, sigma in enumerate(paths_p):
        for k in range(len(sigma)):
            face = sigma[:k] + sigma[k+1:]
            if face in row_idx:
                M[row_idx[face], j] = (-1)**k

    return M, paths_p, paths_pm1

def compute_tournament_betti(arc_set, n, max_p=4):
    """Compute Betti numbers of tournament path homology."""
    allowed = compute_allowed_paths(arc_set, n, max_p)

    def rank(M):
        if M.size == 0:
            return 0
        try:
            _, s, _ = np.linalg.svd(M.astype(float))
            return int(np.sum(s > 1e-9))
        except:
            return np.linalg.matrix_rank(M.astype(float))

    betti = {}
    for p in range(max_p):
        C_p = sorted(allowed[p])
        C_pm1 = sorted(allowed.get(p-1, set())) if p > 0 else []
        C_pp1 = sorted(allowed.get(p+1, set()))

        n_p = len(C_p)
        if n_p == 0:
            betti[p] = 0
            continue

        M_dp, _, _ = build_boundary_matrix(allowed[p], allowed.get(p-1, set()))
        M_dpp1, _, _ = build_boundary_matrix(allowed.get(p+1, set()), allowed[p])

        rank_dp = rank(M_dp)
        rank_dpp1 = rank(M_dpp1)

        dim_ker = n_p - rank_dp
        dim_im = rank_dpp1
        betti[p] = dim_ker - dim_im

    return betti

# --- Bipartite graph H_v ---

def compute_H_v(arc_set, n, v):
    """
    Compute the bipartite graph H_v.
    H_v = (Left=N+(v), Right=N-(v), edges: u->w in T for u in Left, w in Right)
    Returns: (left_vertices, right_vertices, edge_set)
    """
    left = out_nbrs(arc_set, n, v)   # N+(v): v beats these
    right = in_nbrs(arc_set, n, v)   # N-(v): these beat v

    edges = set()
    for u in left:
        for w in right:
            if beats(arc_set, u, w):  # u->w creates 3-cycle v->u->w->v
                edges.add(frozenset([u, w]))

    return left, right, edges

def compute_H2tilde_I_Hv(arc_set, n, v):
    """
    Compute H2tilde(I(H_v)) - the second reduced homology of the independence complex of H_v.
    """
    left, right, edges = compute_H_v(arc_set, n, v)
    vertices = left + right

    if not vertices:
        return 0

    # Build independence complex
    faces = build_independence_complex(edges, vertices)

    # Compute homology
    hom = compute_homology_of_complex(faces)

    return hom.get(2, 0)

def has_cone_point(arc_set, n, v):
    """
    Check if H_v has a cone point (isolated vertex in H_v).
    An isolated vertex in H_v means:
    - u in N+(v) with no edges to any w in N-(v), i.e., all w in N-(v) beat u
    - w in N-(v) with no edges from any u in N+(v), i.e., w beats all u in N+(v)
    """
    left, right, edges = compute_H_v(arc_set, n, v)

    # Check for isolated left vertex
    for u in left:
        if not any(frozenset([u,w]) in edges for w in right):
            return True, ('left', u)

    # Check for isolated right vertex
    for w in right:
        if not any(frozenset([u,w]) in edges for u in left):
            return True, ('right', w)

    return False, None

# --- Main experiments ---

def find_free_cycle(arc_set, n):
    """Find a 3-cycle that is a generator of H1 (if beta1=1)."""
    for a,b,c in itertools.permutations(range(n), 3):
        if beats(arc_set,a,b) and beats(arc_set,b,c) and beats(arc_set,c,a):
            return (a,b,c)
    return None

def experiment_seesaw(n_max=7):
    """Main experiment: test the seesaw conjecture computationally."""
    print(f"=== Seesaw Experiment: β₁·β₃=0 via H̃₂(I(H_v)) ===\n")

    for n in range(3, n_max+1):
        print(f"--- n={n} ---")

        beta1_1_count = 0
        beta3_pos_count = 0
        seesaw_holds = True
        formula_holds = True
        cone_always_exists = True

        # Track cases where formulas fail
        failures = []

        for arc_set in all_tournaments(n):
            betti = compute_tournament_betti(arc_set, n, max_p=4)
            b1 = betti.get(1, 0)
            b3 = betti.get(3, 0)

            if b1 > 0:
                beta1_1_count += 1
            if b3 > 0:
                beta3_pos_count += 1

            # Seesaw check
            if b1 > 0 and b3 > 0:
                seesaw_holds = False
                failures.append(f"  SEESAW VIOLATION: beta1={b1}, beta3={b3}")

            # Formula check: beta3(T) = H2tilde(I(H_v)) for all v (when beta3(T\v)=0)
            # (Only reliable for n<=7 since beta3=0 for n<=6)
            if n <= 7:
                for v in range(n):
                    h2 = compute_H2tilde_I_Hv(arc_set, n, v)
                    if h2 != b3:
                        # This is expected to fail when beta3(T\v) != 0
                        # but for n<=6, beta3(T\v)=0 always, so n<=6 should be exact
                        if n <= 6:
                            formula_holds = False
                            failures.append(f"  FORMULA FAIL at v={v}: b3={b3}, H2tilde={h2}")
                        # For n=7, T\v is n=6 which has beta3=0, so formula should hold
                        elif n == 7:
                            formula_holds = False
                            failures.append(f"  FORMULA FAIL at n=7, v={v}: b3={b3}, H2tilde={h2}")

            # Cone point experiment: when beta1=1, does some v have a cone point?
            if b1 == 1:
                free_cyc = find_free_cycle(arc_set, n)
                if free_cyc is None:
                    failures.append(f"  beta1=1 but no 3-cycle found!")
                    continue

                a,b,c = free_cyc
                other_vertices = [v for v in range(n) if v not in (a,b,c)]

                found_cone = False
                for v in other_vertices:
                    has_cp, _ = has_cone_point(arc_set, n, v)
                    if has_cp:
                        found_cone = True
                        break

                if not found_cone and other_vertices:
                    cone_always_exists = False
                    # More detail
                    h2vals = [compute_H2tilde_I_Hv(arc_set, n, v) for v in other_vertices]
                    if any(h > 0 for h in h2vals):
                        failures.append(f"  No cone point AND H2tilde>0 for v not in free cycle! h2={h2vals}")

        total = 2**(n*(n-1)//2)
        print(f"  Total tournaments: {total}")
        print(f"  beta1=1: {beta1_1_count}")
        print(f"  beta3>0: {beta3_pos_count}")
        print(f"  Seesaw holds: {seesaw_holds}")
        if n <= 7:
            print(f"  Formula H3(T)=H2tilde(I(H_v)) holds for all v: {formula_holds}")
        print(f"  Cone point exists for some v∉{{a,b,c}} when beta1=1: {cone_always_exists}")

        if failures:
            print("  FAILURES:")
            for f in failures[:10]:
                print(f)
        print()

def experiment_bipartite_structure(n=6):
    """Examine the structure of H_v when beta1=1."""
    print(f"=== Bipartite Structure of H_v for beta1=1 tournaments (n={n}) ===\n")

    cases = []
    for arc_set in all_tournaments(n):
        betti = compute_tournament_betti(arc_set, n, max_p=4)
        if betti.get(1, 0) == 1:
            free_cyc = find_free_cycle(arc_set, n)
            if free_cyc is None:
                continue
            a, b, c = free_cyc

            for v in range(n):
                if v in (a,b,c):
                    continue
                left, right, edges = compute_H_v(arc_set, n, v)
                h2 = compute_H2tilde_I_Hv(arc_set, n, v)
                has_cp, cp_info = has_cone_point(arc_set, n, v)
                cases.append({
                    'n': n,
                    'v': v,
                    'free_cyc': (a,b,c),
                    'left_size': len(left),
                    'right_size': len(right),
                    'edges': len(edges),
                    'h2tilde': h2,
                    'has_cone': has_cp,
                    'cone_info': cp_info
                })

    # Summarize
    from collections import Counter
    h2_dist = Counter(c['h2tilde'] for c in cases)
    cone_dist = Counter(c['has_cone'] for c in cases)

    print(f"Cases (all v not in free cycle, beta1=1 tournaments):")
    print(f"  Total: {len(cases)}")
    print(f"  H2tilde distribution: {dict(h2_dist)}")
    print(f"  Has cone point: {dict(cone_dist)}")
    print()

    # When H2tilde=0, what's the structure?
    no_cone_no_h2 = [c for c in cases if not c['has_cone'] and c['h2tilde']==0]
    has_cone = [c for c in cases if c['has_cone']]
    h2_pos = [c for c in cases if c['h2tilde']>0]

    print(f"  Cases with cone point: {len(has_cone)}")
    print(f"  Cases without cone but H2tilde=0: {len(no_cone_no_h2)}")
    print(f"  Cases with H2tilde>0 (unexpected for beta1=1): {len(h2_pos)}")

    if no_cone_no_h2:
        print("\n  [No cone, H2tilde=0 cases] -- what makes H2tilde=0?")
        for c in no_cone_no_h2[:5]:
            print(f"    v={c['v']}, free_cyc={c['free_cyc']}, |L|={c['left_size']}, |R|={c['right_size']}, |edges|={c['edges']}")

    if h2_pos:
        print("\n  WARNING: H2tilde>0 cases (would violate seesaw!):")
        for c in h2_pos[:5]:
            print(f"    v={c['v']}, free_cyc={c['free_cyc']}, |L|={c['left_size']}, |R|={c['right_size']}, |edges|={c['edges']}, H2={c['h2tilde']}")

def experiment_formula_exact(n_max=6):
    """Verify: for all v, beta3(T) = H2tilde(I(H_v)) exactly (when n<=6 so beta3(T\v)=0)."""
    print(f"=== Exact Formula Verification: β₃(T) = H̃₂(I(H_v)) ===\n")

    for n in range(3, n_max+1):
        print(f"n={n}: ", end='', flush=True)
        mismatches = 0
        total = 0

        for arc_set in all_tournaments(n):
            betti = compute_tournament_betti(arc_set, n, max_p=4)
            b3 = betti.get(3, 0)

            for v in range(n):
                h2 = compute_H2tilde_I_Hv(arc_set, n, v)
                total += 1
                if h2 != b3:
                    mismatches += 1

        print(f"{total} (v,T) pairs, {mismatches} mismatches. {'OK' if mismatches==0 else 'FAIL'}")

if __name__ == '__main__':
    import sys

    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'

    if mode == 'formula':
        experiment_formula_exact(n_max=6)
    elif mode == 'bipartite':
        experiment_bipartite_structure(n=6)
    elif mode == 'seesaw':
        experiment_seesaw(n_max=6)
    else:
        # Run all experiments
        experiment_formula_exact(n_max=5)
        print()
        experiment_seesaw(n_max=6)
        print()
        experiment_bipartite_structure(n=5)
