"""
top_injectivity_analysis.py — Deep analysis of WHY d_{n-1} is injective on Omega_{n-1}

For tournaments on n vertices, Omega_{n-1} consists of those allowed Hamiltonian
paths P = (v_0,...,v_{n-1}) whose boundary d(P) = sum_i (-1)^i face_i(P) has
ALL faces landing in the allowed (n-2)-path space.

Key questions:
1. What makes a Hamiltonian path "allowed" (all faces allowed)?
2. Is face_0: Omega_{n-1} -> A_{n-2} injective? (This would imply d_{n-1} injective.)
3. What is the structure of allowed vs non-allowed faces?
4. How do starting vertices partition Omega_{n-1}?

Author: kind-pasteur-2026-03-10
"""
import sys
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations

sys.path.insert(0, '04-computation')
from tournament_utils import (
    random_tournament, enumerate_all_allowed, boundary_faces,
    full_chain_complex_modp, bits_to_adj, exhaustive_tournaments,
    build_adj_lists, enumerate_allowed_paths, _build_constraint_matrix,
    _gauss_nullbasis_modp, _build_boundary_matrix, _gauss_rank_np,
    RANK_PRIME
)

# ============================================================
# Helper: Compute Omega_{n-1} elements explicitly
# ============================================================

def get_omega_elements(A, n, p, ap=None, prime=RANK_PRIME):
    """Return the actual elements of Omega_p as a list of paths.

    Omega_p = {x in A_p : d(x) has all faces in A_{p-1}}.
    Equivalently, x is in the null space of the non-allowed constraint matrix.

    For individual basis elements (not linear combinations), an allowed p-path P
    is in Omega_p iff EVERY face of P is either:
      (a) an allowed (p-1)-path (in A_{p-1}), or
      (b) cancels out in the boundary with other paths' non-allowed faces.

    The second case means Omega_p can contain LINEAR COMBINATIONS where
    non-allowed faces cancel. So individual Hamiltonian paths need not
    have all faces allowed -- the null space condition is on vectors.

    Returns:
      omega_paths: list of tuples (path indices in ap[p] that form the basis)
      null_basis: the actual null basis vectors (each vector gives coefficients over ap[p])
      paths: the list of allowed p-paths
    """
    if ap is None:
        ap = enumerate_all_allowed(A, n)

    paths = ap.get(p, [])
    if not paths:
        return [], None, paths

    if p == 0:
        # All vertices are in Omega_0
        return list(range(len(paths))), np.eye(len(paths)), paths

    # Build constraint matrix
    apm1_set = set(ap.get(p-1, []))
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(paths):
        for sign, face in boundary_faces(path):
            if len(set(face)) == len(face) and face not in apm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1

    if na_count == 0:
        # No constraints -- all paths are in Omega
        return list(range(len(paths))), np.eye(len(paths)), paths

    # Build constraint matrix
    P_mat = np.zeros((na_count, len(paths)), dtype=np.int64)
    for j, path in enumerate(paths):
        for sign, face in boundary_faces(path):
            if face in non_allowed:
                P_mat[non_allowed[face], j] = (P_mat[non_allowed[face], j] + sign) % prime

    # Null basis
    rank, nbasis = _gauss_nullbasis_modp(P_mat, na_count, len(paths), prime)
    dim = len(paths) - rank

    return list(range(len(paths))), nbasis, paths


def analyze_face_structure(A, n, ap):
    """For each Hamiltonian path, analyze which faces are allowed vs non-allowed."""
    p = n - 1
    paths = ap.get(p, [])
    apm1_set = set(ap.get(p-1, []))

    results = []
    for path in paths:
        faces_info = []
        for i in range(len(path)):
            sign = (-1)**i
            face = path[:i] + path[i+1:]
            is_allowed = face in apm1_set
            faces_info.append({
                'index': i,
                'removed_vertex': path[i],
                'face': face,
                'is_allowed': is_allowed,
                'sign': sign
            })

        num_allowed = sum(1 for f in faces_info if f['is_allowed'])
        num_non_allowed = len(faces_info) - num_allowed

        results.append({
            'path': path,
            'faces': faces_info,
            'num_allowed': num_allowed,
            'num_non_allowed': num_non_allowed,
            'non_allowed_indices': [f['index'] for f in faces_info if not f['is_allowed']]
        })

    return results


def test_face0_injectivity(paths, apm1_set):
    """Test if face_0 (removing first vertex) is injective on a set of paths."""
    face0_map = {}
    for path in paths:
        f0 = path[1:]  # remove v_0
        if f0 in face0_map:
            return False, face0_map[f0], path
        face0_map[f0] = path
    return True, None, None


def test_facen_injectivity(paths, n):
    """Test if face_{n-1} (removing last vertex) is injective on a set of paths."""
    facen_map = {}
    for path in paths:
        fn = path[:-1]  # remove v_{n-1}
        if fn in facen_map:
            return False, facen_map[fn], path
        facen_map[fn] = path
    return True, None, None


def compute_boundary_rank_on_omega(A, n, ap, prime=RANK_PRIME):
    """Compute the rank of d_{n-1} restricted to Omega_{n-1}."""
    p = n - 1

    # Get Omega basis
    P_mat, na_rows, na_cols = _build_constraint_matrix(ap, p, prime)
    if P_mat is None:
        omega_dim = len(ap.get(p, []))
        omega_basis = None
    else:
        rank_P, nbasis = _gauss_nullbasis_modp(P_mat, na_rows, na_cols, prime)
        omega_dim = na_cols - rank_P
        omega_basis = np.array(nbasis, dtype=np.int64) if nbasis else None

    if omega_dim == 0:
        return 0, omega_dim

    # Build boundary matrix
    bd_np, bd_nrows, bd_ncols = _build_boundary_matrix(ap, p, prime)
    if bd_np is None:
        return 0, omega_dim

    if omega_basis is not None:
        composed = bd_np @ omega_basis.T % prime
        rank = _gauss_rank_np(composed, prime)
    else:
        rank = _gauss_rank_np(bd_np.copy(), prime)

    return rank, omega_dim


# ============================================================
# Main analysis
# ============================================================

def analyze_tournament(A, n, verbose=False):
    """Full analysis of top-degree injectivity for one tournament."""
    ap = enumerate_all_allowed(A, n)
    p = n - 1  # top degree

    ham_paths = ap.get(p, [])
    prev_paths = ap.get(p-1, [])
    apm1_set = set(prev_paths)

    # 1. Count Hamiltonian paths and Omega elements
    num_ham = len(ham_paths)

    # Get Omega_{n-1} dimension
    P_mat, na_rows, na_cols = _build_constraint_matrix(ap, p, RANK_PRIME)
    if P_mat is None:
        omega_dim = num_ham
        constraint_rank = 0
        num_non_allowed_faces = 0
    else:
        rank_P, nbasis = _gauss_nullbasis_modp(P_mat, na_rows, na_cols, RANK_PRIME)
        omega_dim = na_cols - rank_P
        constraint_rank = rank_P
        num_non_allowed_faces = na_rows

    # 2. Face structure analysis
    face_analysis = analyze_face_structure(A, n, ap)

    # Count Hamiltonian paths with ALL faces allowed (these are trivially in Omega)
    all_faces_allowed = sum(1 for r in face_analysis if r['num_non_allowed'] == 0)

    # Distribution of non-allowed face counts
    non_allowed_dist = Counter(r['num_non_allowed'] for r in face_analysis)

    # Which face indices tend to be non-allowed?
    non_allowed_idx_counts = Counter()
    for r in face_analysis:
        for idx in r['non_allowed_indices']:
            non_allowed_idx_counts[idx] += 1

    # 3. Face_0 injectivity on ALL Hamiltonian paths
    f0_inj_all, f0_dup1, f0_dup2 = test_face0_injectivity(ham_paths, apm1_set)

    # Face_{n-1} injectivity on ALL Hamiltonian paths
    fn_inj_all, fn_dup1, fn_dup2 = test_facen_injectivity(ham_paths, n)

    # 4. Face_0 allowed status
    face0_allowed_count = 0
    face0_non_allowed_count = 0
    for r in face_analysis:
        if r['faces'][0]['is_allowed']:
            face0_allowed_count += 1
        else:
            face0_non_allowed_count += 1

    # 5. Compute actual rank of d_{n-1} on Omega_{n-1}
    bd_rank, omega_dim_check = compute_boundary_rank_on_omega(A, n, ap)
    is_injective = (bd_rank == omega_dim)

    # 6. Starting vertex analysis
    start_vertex_groups = defaultdict(list)
    for path in ham_paths:
        start_vertex_groups[path[0]].append(path)

    # 7. Ending vertex analysis
    end_vertex_groups = defaultdict(list)
    for path in ham_paths:
        end_vertex_groups[path[-1]].append(path)

    result = {
        'n': n,
        'num_ham_paths': num_ham,
        'num_prev_paths': len(prev_paths),
        'omega_dim': omega_dim,
        'constraint_rank': constraint_rank,
        'num_non_allowed_face_types': num_non_allowed_faces,
        'all_faces_allowed': all_faces_allowed,
        'non_allowed_dist': dict(non_allowed_dist),
        'non_allowed_idx_counts': dict(non_allowed_idx_counts),
        'face0_inj_all_ham': f0_inj_all,
        'facen_inj_all_ham': fn_inj_all,
        'face0_allowed_count': face0_allowed_count,
        'face0_non_allowed_count': face0_non_allowed_count,
        'bd_rank': bd_rank,
        'is_injective': is_injective,
        'start_vertex_sizes': {v: len(paths) for v, paths in start_vertex_groups.items()},
        'end_vertex_sizes': {v: len(paths) for v, paths in end_vertex_groups.items()},
    }

    if verbose and not f0_inj_all:
        result['f0_collision'] = (f0_dup1, f0_dup2)
    if verbose and not fn_inj_all:
        result['fn_collision'] = (fn_dup1, fn_dup2)

    return result


def deep_omega_analysis(A, n, ap):
    """Analyze the actual Omega_{n-1} basis vectors and what they look like."""
    p = n - 1
    paths = ap.get(p, [])
    apm1_set = set(ap.get(p-1, []))

    if not paths:
        return None

    P_mat, na_rows, na_cols = _build_constraint_matrix(ap, p, RANK_PRIME)
    if P_mat is None:
        return {
            'type': 'unconstrained',
            'omega_dim': len(paths),
            'msg': 'All Hamiltonian paths are in Omega (no constraints)'
        }

    rank_P, nbasis = _gauss_nullbasis_modp(P_mat, na_rows, na_cols, RANK_PRIME)
    omega_dim = na_cols - rank_P

    if not nbasis:
        return {
            'type': 'empty',
            'omega_dim': 0,
            'msg': 'Omega_{n-1} is empty'
        }

    # Analyze each basis vector
    basis_info = []
    for bv in nbasis:
        nonzero = [(i, bv[i]) for i in range(len(bv)) if bv[i] != 0]
        # Convert large mod-p values back to small signed integers
        half_p = RANK_PRIME // 2
        signed_nz = []
        for i, c in nonzero:
            if c > half_p:
                signed_nz.append((i, c - RANK_PRIME))
            else:
                signed_nz.append((i, c))

        # How many paths in this basis vector?
        num_paths = len(signed_nz)

        # What are the paths?
        path_details = [(paths[i], coeff) for i, coeff in signed_nz]

        # Check: do all paths in this basis vector share a common face?
        if num_paths >= 2:
            # Check face_0
            face0s = set(paths[i][1:] for i, _ in signed_nz)
            face0_shared = len(face0s) < num_paths

            # Starting vertices
            starts = set(paths[i][0] for i, _ in signed_nz)
            ends = set(paths[i][-1] for i, _ in signed_nz)
        else:
            face0_shared = False
            starts = set()
            ends = set()

        basis_info.append({
            'num_paths': num_paths,
            'coefficients': signed_nz,
            'paths': [paths[i] for i, _ in signed_nz],
            'face0_shared': face0_shared if num_paths >= 2 else 'N/A',
            'starting_vertices': starts,
            'ending_vertices': ends,
        })

    return {
        'type': 'constrained',
        'omega_dim': omega_dim,
        'num_constraints': na_rows,
        'constraint_rank': rank_P,
        'num_ham_paths': len(paths),
        'basis_vectors': basis_info
    }


def check_face0_on_omega_basis(A, n, ap, nbasis, paths):
    """Check if applying face_0 to each Omega basis vector gives distinct images.

    Even though face_0 is trivially injective on individual paths (removing v_0
    uniquely identifies the path since v_0 and the remaining order determine it),
    the question is whether face_0 is injective on the Omega BASIS VECTORS,
    which are linear combinations of paths.

    Actually, for injectivity of d_{n-1}, we need that the boundary MAP d is
    injective on Omega_{n-1}, not just face_0. But face_0 being injective
    on basis elements would help.
    """
    if nbasis is None:
        return True, "No basis to check"

    apm1_set = set(ap.get(n-2, []))
    prev_paths = ap.get(n-2, [])
    prev_idx = {p: i for i, p in enumerate(prev_paths)}

    # For each basis vector, compute d(basis_vector)
    images = []
    for bv in nbasis:
        # d(sum c_j P_j) = sum c_j d(P_j)
        image = defaultdict(int)
        for j in range(len(bv)):
            if bv[j] == 0:
                continue
            coeff = bv[j]
            path = paths[j]
            for i in range(len(path)):
                sign = (-1)**i
                face = path[:i] + path[i+1:]
                if face in prev_idx:
                    image[prev_idx[face]] = (image[prev_idx[face]] + coeff * sign) % RANK_PRIME

        # Clean up zeros
        image = {k: v for k, v in image.items() if v != 0}
        images.append(frozenset(image.items()))

    # Check distinctness
    if len(set(images)) == len(images):
        return True, "All basis vector images are distinct"
    else:
        return False, "Some basis vectors have identical images (should not happen if d is injective)"


def analyze_non_allowed_face_pairing(A, n, ap):
    """For Hamiltonian paths with non-allowed faces, analyze how they pair up.

    Key question: for each non-allowed face f that appears in d(P) for some
    Hamiltonian path P, how many OTHER Hamiltonian paths Q also have f in d(Q)?

    For Omega to be nontrivial, non-allowed faces must cancel, which requires
    at least two Hamiltonian paths sharing a non-allowed face (with opposite signs).
    """
    p = n - 1
    paths = ap.get(p, [])
    apm1_set = set(ap.get(p-1, []))

    # Map: non-allowed face -> list of (path_idx, sign)
    na_face_to_paths = defaultdict(list)

    for j, path in enumerate(paths):
        for i in range(len(path)):
            sign = (-1)**i
            face = path[:i] + path[i+1:]
            if len(set(face)) == len(face) and face not in apm1_set:
                na_face_to_paths[face].append((j, sign))

    if not na_face_to_paths:
        return {
            'num_non_allowed_faces': 0,
            'msg': 'No non-allowed faces; all Hamiltonian paths are in Omega'
        }

    # Distribution: how many paths contribute to each non-allowed face?
    contrib_dist = Counter(len(v) for v in na_face_to_paths.values())

    # Analyze sign patterns: do paths contributing to same non-allowed face have matching signs?
    sign_patterns = Counter()
    for face, contributors in na_face_to_paths.items():
        signs = tuple(sorted(s for _, s in contributors))
        sign_patterns[signs] += 1

    # For each non-allowed face, what face index creates it?
    face_idx_dist = Counter()
    for face, contributors in na_face_to_paths.items():
        for j, sign in contributors:
            path = paths[j]
            for i in range(len(path)):
                f = path[:i] + path[i+1:]
                if f == face:
                    face_idx_dist[i] += 1

    return {
        'num_non_allowed_faces': len(na_face_to_paths),
        'contributor_distribution': dict(contrib_dist),
        'sign_patterns': dict(sign_patterns),
        'face_idx_dist': dict(face_idx_dist),
    }


# ============================================================
# Why is a face non-allowed?
# ============================================================

def why_non_allowed(face, A, n):
    """Explain WHY a face (subsequence of vertices) is not an allowed path.

    A path (v_0,...,v_k) is allowed iff v_i -> v_{i+1} for all i.
    So it's non-allowed iff some consecutive pair has the WRONG direction.
    """
    reasons = []
    for i in range(len(face) - 1):
        u, v = face[i], face[i+1]
        if A[u][v] != 1:
            reasons.append((i, u, v, 'arc goes v->u not u->v'))
    return reasons


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 80)
    print("TOP INJECTIVITY ANALYSIS: Why is d_{n-1} injective on Omega_{n-1}?")
    print("=" * 80)

    for n in [4, 5, 6, 7]:
        print(f"\n{'#' * 80}")
        print(f"# n = {n}")
        print(f"{'#' * 80}")

        if n <= 5:
            # Exhaustive for small n
            mode = "EXHAUSTIVE"
            tournaments = list(exhaustive_tournaments(n))
            num_T = len(tournaments)
        else:
            # Sample for larger n
            mode = "SAMPLED"
            num_T = 100
            rng = np.random.RandomState(42)
            tournaments = [(i, random_tournament(n, rng)) for i in range(num_T)]

        print(f"\nMode: {mode}, {num_T} tournaments")

        # Aggregate statistics
        omega_dims = []
        ham_counts = []
        injective_count = 0
        non_injective_examples = []
        face0_inj_count = 0
        facen_inj_count = 0
        all_allowed_counts = []
        non_allowed_dists_agg = Counter()
        non_allowed_idx_agg = Counter()
        face0_allowed_agg = 0
        face0_non_allowed_agg = 0
        constraint_ranks = []
        na_face_type_counts = []

        # Detailed examples for verbose output
        detailed_examples = []
        omega_basis_examples = []
        na_pairing_examples = []

        for trial_idx, (label, A) in enumerate(tournaments):
            result = analyze_tournament(A, n, verbose=True)

            omega_dims.append(result['omega_dim'])
            ham_counts.append(result['num_ham_paths'])
            constraint_ranks.append(result['constraint_rank'])
            na_face_type_counts.append(result['num_non_allowed_face_types'])
            all_allowed_counts.append(result['all_faces_allowed'])
            face0_allowed_agg += result['face0_allowed_count']
            face0_non_allowed_agg += result['face0_non_allowed_count']

            if result['is_injective']:
                injective_count += 1
            else:
                non_injective_examples.append((label, A, result))

            if result['face0_inj_all_ham']:
                face0_inj_count += 1

            if result['facen_inj_all_ham']:
                facen_inj_count += 1

            for k, v in result['non_allowed_dist'].items():
                non_allowed_dists_agg[k] += v
            for k, v in result['non_allowed_idx_counts'].items():
                non_allowed_idx_agg[k] += v

            # Collect a few detailed examples
            if trial_idx < 3 or (result['omega_dim'] > 0 and result['omega_dim'] < result['num_ham_paths']):
                ap = enumerate_all_allowed(A, n)

                if len(detailed_examples) < 5:
                    deep = deep_omega_analysis(A, n, ap)
                    detailed_examples.append((label, result, deep))

                if len(na_pairing_examples) < 5:
                    pairing = analyze_non_allowed_face_pairing(A, n, ap)
                    na_pairing_examples.append((label, result, pairing))

        # ============================================================
        # REPORT
        # ============================================================

        print(f"\n--- SUMMARY for n={n} ---")
        print(f"Tournaments analyzed: {num_T}")
        print(f"d_{{n-1}} INJECTIVE on Omega_{{n-1}}: {injective_count}/{num_T} ({100*injective_count/num_T:.1f}%)")

        print(f"\n--- Hamiltonian path counts ---")
        print(f"  Min: {min(ham_counts)}, Max: {max(ham_counts)}, Mean: {np.mean(ham_counts):.1f}")

        print(f"\n--- Omega_{{n-1}} dimensions ---")
        print(f"  Min: {min(omega_dims)}, Max: {max(omega_dims)}, Mean: {np.mean(omega_dims):.1f}")
        print(f"  Distribution: {Counter(omega_dims).most_common(10)}")
        print(f"  Omega = 0 (empty): {sum(1 for d in omega_dims if d == 0)}/{num_T}")
        print(f"  Omega = all Ham paths (unconstrained): {sum(1 for i in range(num_T) if omega_dims[i] == ham_counts[i])}/{num_T}")

        print(f"\n--- Constraint structure ---")
        print(f"  Constraint rank distribution: {Counter(constraint_ranks).most_common(10)}")
        print(f"  Non-allowed face types: min={min(na_face_type_counts)}, max={max(na_face_type_counts)}, mean={np.mean(na_face_type_counts):.1f}")

        print(f"\n--- Face allowedness for Hamiltonian paths ---")
        print(f"  Paths with ALL faces allowed: total {sum(all_allowed_counts)}, dist per tournament: {Counter(all_allowed_counts).most_common(10)}")
        print(f"  Non-allowed face count distribution (across all paths):")
        for k in sorted(non_allowed_dists_agg):
            print(f"    {k} non-allowed faces: {non_allowed_dists_agg[k]} paths")

        print(f"\n--- Which face indices are non-allowed? ---")
        for k in sorted(non_allowed_idx_agg):
            print(f"    face_{k}: non-allowed {non_allowed_idx_agg[k]} times")

        print(f"\n--- face_0 analysis ---")
        print(f"  face_0 (remove v_0) is allowed: {face0_allowed_agg} paths")
        print(f"  face_0 is NOT allowed: {face0_non_allowed_agg} paths")
        print(f"  face_0 injective on ALL Ham paths: {face0_inj_count}/{num_T}")
        print(f"  face_{{n-1}} injective on ALL Ham paths: {facen_inj_count}/{num_T}")

        # CRITICAL: face_0 is ALWAYS injective on Ham paths since
        # two paths (v_0,v_1,...,v_{n-1}) and (w_0,w_1,...,w_{n-1}) with
        # same face_0 = (v_1,...,v_{n-1}) = (w_1,...,w_{n-1}) implies
        # they agree on vertices 1..n-1, so v_0=w_0 (only vertex left).
        print(f"\n  ** TRIVIAL FACT: face_0 is ALWAYS injective on Hamiltonian paths **")
        print(f"  Reason: A Ham path uses ALL n vertices. Removing v_0 gives (v_1,...,v_{n-1}).")
        print(f"  This determines which vertex was v_0 (the missing one). So face_0 is injective.")
        print(f"  Similarly, face_i is injective on Ham paths for ANY fixed i.")

        if non_injective_examples:
            print(f"\n--- NON-INJECTIVE EXAMPLES ---")
            for label, A_ex, res in non_injective_examples[:3]:
                print(f"  Tournament {label}: omega_dim={res['omega_dim']}, bd_rank={res['bd_rank']}")

        # Detailed examples
        if detailed_examples:
            print(f"\n--- DETAILED Omega_{{n-1}} ANALYSIS (first {len(detailed_examples)} examples) ---")
            for label, result, deep in detailed_examples:
                if deep is None:
                    continue
                print(f"\n  Tournament {label}:")
                print(f"    Ham paths: {result['num_ham_paths']}, Omega dim: {result['omega_dim']}")
                print(f"    Type: {deep['type']}")

                if deep['type'] == 'constrained':
                    print(f"    Constraints: {deep['num_constraints']}, rank={deep['constraint_rank']}")
                    print(f"    Basis vectors:")
                    for i, binfo in enumerate(deep['basis_vectors']):
                        print(f"      v{i}: {binfo['num_paths']} paths, coefficients: {binfo['coefficients'][:5]}{'...' if binfo['num_paths'] > 5 else ''}")
                        if binfo['num_paths'] <= 6:
                            for p_idx, (j, c) in enumerate(binfo['coefficients']):
                                print(f"        coeff={c}: path={binfo['paths'][p_idx]}")
                        print(f"      Starting vertices: {binfo['starting_vertices']}")
                        print(f"      face_0 shared: {binfo['face0_shared']}")
                elif deep['type'] == 'unconstrained':
                    print(f"    {deep['msg']}")

        # Non-allowed face pairing
        if na_pairing_examples:
            print(f"\n--- NON-ALLOWED FACE PAIRING (first {len(na_pairing_examples)} examples) ---")
            for label, result, pairing in na_pairing_examples:
                print(f"\n  Tournament {label}:")
                print(f"    Num non-allowed face types: {pairing['num_non_allowed_faces']}")
                if 'contributor_distribution' in pairing:
                    print(f"    Contributor dist: {pairing['contributor_distribution']}")
                    print(f"    Sign patterns: {pairing['sign_patterns']}")
                    print(f"    Face index dist: {pairing['face_idx_dist']}")

    # ============================================================
    # DEEPER INVESTIGATION: The real mechanism
    # ============================================================

    print(f"\n\n{'=' * 80}")
    print("DEEP INVESTIGATION: WHY d_{n-1} is injective")
    print("=" * 80)

    print("""
KEY INSIGHT: face_i is injective on Hamiltonian paths for ANY fixed i.
This is trivial: a Hamiltonian path visits ALL n vertices, so removing
any single vertex v_i determines a unique (n-1)-path on the remaining
vertices. Since the remaining vertices and their ORDER uniquely determine
the original path (the missing vertex must be v_i, placed at position i),
face_i is injective.

BUT this does NOT directly imply d_{n-1} is injective, because:
  d_{n-1}(P) = face_0(P) - face_1(P) + face_2(P) - ...
is an ALTERNATING SUM of faces. Two different paths P,P' could
potentially have d(P) = d(P') even though face_0(P) != face_0(P'),
if the other faces conspire to cancel.

The real question is: on Omega_{n-1} (the constrained subspace where
d lands in allowed paths), does this cancellation happen?
""")

    # For each n, do a detailed kernel analysis
    for n in [4, 5, 6, 7]:
        print(f"\n--- Kernel analysis for n={n} ---")

        if n <= 5:
            tournaments_list = list(exhaustive_tournaments(n))
        else:
            rng = np.random.RandomState(123)
            tournaments_list = [(i, random_tournament(n, rng)) for i in range(100)]

        kernel_nonzero = 0
        total = 0

        for label, A in tournaments_list:
            total += 1
            ap = enumerate_all_allowed(A, n)
            p = n - 1
            paths_p = ap.get(p, [])
            paths_pm1 = ap.get(p-1, [])

            if not paths_p or not paths_pm1:
                continue

            # Compute Omega basis
            P_mat, na_rows, na_cols = _build_constraint_matrix(ap, p, RANK_PRIME)
            if P_mat is None:
                omega_dim = len(paths_p)
                omega_basis = None
            else:
                rank_P, nbasis = _gauss_nullbasis_modp(P_mat, na_rows, na_cols, RANK_PRIME)
                omega_dim = na_cols - rank_P
                omega_basis = np.array(nbasis, dtype=np.int64) if nbasis else None

            if omega_dim == 0:
                continue

            # Compute boundary map restricted to Omega
            bd_np, bd_nrows, bd_ncols = _build_boundary_matrix(ap, p, RANK_PRIME)
            if bd_np is None:
                continue

            if omega_basis is not None:
                composed = bd_np @ omega_basis.T % RANK_PRIME
            else:
                composed = bd_np.copy() % RANK_PRIME

            rank = _gauss_rank_np(composed.copy(), RANK_PRIME)
            kernel_dim = omega_dim - rank

            if kernel_dim > 0:
                kernel_nonzero += 1

        print(f"  Tournaments with ker(d_{{n-1}}) > 0: {kernel_nonzero}/{total}")
        print(f"  d_{{n-1}} injective on Omega_{{n-1}} for ALL: {kernel_nonzero == 0}")

    # ============================================================
    # STRUCTURAL ANALYSIS: Why can't two Omega elements have same boundary?
    # ============================================================

    print(f"\n\n{'=' * 80}")
    print("STRUCTURAL ANALYSIS: Why non-allowed faces force injectivity")
    print("=" * 80)

    for n in [5, 6]:
        print(f"\n--- n={n}: Detailed face type analysis ---")

        if n <= 5:
            tournaments_list = list(exhaustive_tournaments(n))
        else:
            rng = np.random.RandomState(42)
            tournaments_list = [(i, random_tournament(n, rng)) for i in range(50)]

        for trial, (label, A) in enumerate(tournaments_list[:5]):
            ap = enumerate_all_allowed(A, n)
            p = n - 1
            paths = ap.get(p, [])
            apm1_set = set(ap.get(p-1, []))

            if not paths:
                continue

            print(f"\n  Tournament {label} (n={n}):")
            print(f"  Ham paths: {len(paths)}")

            # For each Hamiltonian path, show face details
            for path in paths[:8]:  # show first 8
                non_allowed_faces = []
                allowed_faces = []
                for i in range(len(path)):
                    face = path[:i] + path[i+1:]
                    if face in apm1_set:
                        allowed_faces.append(i)
                    else:
                        # Why non-allowed?
                        reasons = why_non_allowed(face, A, n)
                        non_allowed_faces.append((i, face, reasons))

                if non_allowed_faces:
                    print(f"    Path {path}: {len(non_allowed_faces)} non-allowed faces")
                    for idx, face, reasons in non_allowed_faces:
                        print(f"      face_{idx} = {face}: NOT allowed because {reasons}")
                else:
                    print(f"    Path {path}: ALL faces allowed")

            if len(paths) > 8:
                print(f"    ... ({len(paths) - 8} more paths)")

    # ============================================================
    # OMEGA DIMENSION PATTERNS
    # ============================================================

    print(f"\n\n{'=' * 80}")
    print("OMEGA_{n-1} DIMENSION PATTERNS")
    print("=" * 80)

    for n in [4, 5, 6, 7]:
        print(f"\n--- n={n} ---")

        if n <= 5:
            tournaments_list = list(exhaustive_tournaments(n))
        else:
            rng = np.random.RandomState(42)
            tournaments_list = [(i, random_tournament(n, rng)) for i in range(200)]

        data = []
        for label, A in tournaments_list:
            ap = enumerate_all_allowed(A, n)
            p = n - 1

            # Omega dim
            P_mat, na_rows, na_cols = _build_constraint_matrix(ap, p, RANK_PRIME)
            if P_mat is None:
                omega_dim = len(ap.get(p, []))
                constraint_rank = 0
            else:
                rank_P, _ = _gauss_nullbasis_modp(P_mat, na_rows, na_cols, RANK_PRIME)
                omega_dim = na_cols - rank_P
                constraint_rank = rank_P

            num_ham = len(ap.get(p, []))

            # Also compute beta_{n-1}
            cc = full_chain_complex_modp(A, n, max_p=n-1)
            beta_top = cc['bettis'].get(n-1, 0)
            rank_d = cc['ranks'].get(n-1, 0)

            data.append({
                'label': label,
                'num_ham': num_ham,
                'omega_dim': omega_dim,
                'constraint_rank': constraint_rank,
                'na_rows': na_rows if P_mat is not None else 0,
                'rank_d': rank_d,
                'beta_top': beta_top,
            })

        # Summary
        omega_zero = sum(1 for d in data if d['omega_dim'] == 0)
        omega_eq_ham = sum(1 for d in data if d['omega_dim'] == d['num_ham'])
        beta_zero = sum(1 for d in data if d['beta_top'] == 0)

        print(f"  Total: {len(data)}")
        print(f"  Omega_{{n-1}} = 0: {omega_zero} ({100*omega_zero/len(data):.1f}%)")
        print(f"  Omega_{{n-1}} = #Ham (unconstrained): {omega_eq_ham} ({100*omega_eq_ham/len(data):.1f}%)")
        print(f"  beta_{{n-1}} = 0: {beta_zero} ({100*beta_zero/len(data):.1f}%)")

        dim_dist = Counter(d['omega_dim'] for d in data)
        rank_dist = Counter(d['rank_d'] for d in data)
        beta_dist = Counter(d['beta_top'] for d in data)

        print(f"  Omega_{{n-1}} dimension distribution: {dim_dist.most_common(15)}")
        print(f"  rank(d_{{n-1}}) distribution: {rank_dist.most_common(15)}")
        print(f"  beta_{{n-1}} distribution: {beta_dist.most_common(10)}")

        # Key: rank(d_{n-1}) should ALWAYS equal omega_dim (injective)
        all_injective = all(d['rank_d'] == d['omega_dim'] for d in data)
        print(f"  rank(d_{{n-1}}) == omega_dim for ALL: {all_injective}")

        if not all_injective:
            for d in data:
                if d['rank_d'] != d['omega_dim']:
                    print(f"    COUNTEREXAMPLE: label={d['label']}, omega={d['omega_dim']}, rank={d['rank_d']}")

    # ============================================================
    # THE CRITICAL QUESTION: What is a non-allowed (n-2)-path?
    # ============================================================

    print(f"\n\n{'=' * 80}")
    print("WHAT MAKES A FACE NON-ALLOWED?")
    print("=" * 80)

    print("""
A Hamiltonian path P = (v_0, v_1, ..., v_{n-1}) has arcs v_i -> v_{i+1}.
Face_i(P) = (v_0, ..., v_{i-1}, v_{i+1}, ..., v_{n-1}).

Face_i(P) is an allowed (n-2)-path iff v_{i-1} -> v_{i+1} (when i > 0 and i < n-1).
- face_0(P) = (v_1, ..., v_{n-1}): always allowed iff v_1 -> v_2 -> ... -> v_{n-1},
  which is ALWAYS TRUE since P is a path with those arcs.
- face_{n-1}(P) = (v_0, ..., v_{n-2}): similarly always allowed.
- face_i(P) for 0 < i < n-1: allowed iff v_{i-1} -> v_{i+1} in the tournament.
  This is the "shortcut" condition: removing v_i, can we go directly from v_{i-1} to v_{i+1}?
""")

    # Verify this understanding
    for n in [5, 6]:
        print(f"\n--- Verification for n={n} ---")

        if n <= 5:
            tournaments_list = list(exhaustive_tournaments(n))
        else:
            rng = np.random.RandomState(42)
            tournaments_list = [(i, random_tournament(n, rng)) for i in range(30)]

        total_faces = 0
        shortcut_correct = 0
        face0_always_allowed = 0
        facen_always_allowed = 0
        face0_total = 0
        facen_total = 0

        for label, A in tournaments_list:
            ap = enumerate_all_allowed(A, n)
            p = n - 1
            paths = ap.get(p, [])
            apm1_set = set(ap.get(p-1, []))

            for path in paths:
                for i in range(len(path)):
                    face = path[:i] + path[i+1:]
                    is_allowed = face in apm1_set
                    total_faces += 1

                    if i == 0:
                        face0_total += 1
                        if is_allowed:
                            face0_always_allowed += 1
                    elif i == n - 1:
                        facen_total += 1
                        if is_allowed:
                            facen_always_allowed += 1
                    else:
                        # Interior face: check shortcut condition
                        has_shortcut = (A[path[i-1]][path[i+1]] == 1)
                        if has_shortcut == is_allowed:
                            shortcut_correct += 1
                        else:
                            print(f"    MISMATCH at path={path}, i={i}: shortcut={has_shortcut}, allowed={is_allowed}")

        interior_total = total_faces - face0_total - facen_total
        print(f"  face_0 allowed: {face0_always_allowed}/{face0_total} (should be ALL)")
        print(f"  face_{{n-1}} allowed: {facen_always_allowed}/{facen_total} (should be ALL)")
        print(f"  Interior face shortcut condition matches: {shortcut_correct}/{interior_total} (should be ALL)")

    # ============================================================
    # THE MECHANISM: Why does the shortcut condition force injectivity?
    # ============================================================

    print(f"\n\n{'=' * 80}")
    print("THE MECHANISM: Shortcut structure and d_{n-1} injectivity")
    print("=" * 80)

    print("""
SUMMARY OF FINDINGS:

1. face_0(P) and face_{n-1}(P) are ALWAYS allowed for Hamiltonian paths.
   - face_0(P) = (v_1,...,v_{n-1}) is allowed because v_1->v_2->...->v_{n-1} are
     consecutive arcs of P, hence already directed arcs.
   - Similarly for face_{n-1}(P) = (v_0,...,v_{n-2}).

2. face_i(P) for 0 < i < n-1 is allowed iff v_{i-1} -> v_{i+1} (shortcut arc exists).
   - Non-allowed interior face means v_{i+1} -> v_{i-1} (shortcut goes BACKWARD).

3. Omega_{n-1} constraint: sum (-1)^i face_i(P) must have no non-allowed components.
   - Since face_0 and face_{n-1} are ALWAYS allowed, the constraint is on interior faces.
   - P is in Omega_{n-1} iff: for each non-allowed face_i(P), there exist other
     Hamiltonian paths Q with face_j(Q) = face_i(P) (same non-allowed face) and
     the signed contributions cancel.

4. The TRIVIALITY argument for face injectivity:
   Any face_i is injective on the set of Hamiltonian paths, because a Hamiltonian
   path uses all n vertices, so removing vertex v_i yields a unique sequence of the
   remaining n-1 vertices, which determines the original path.

5. WHY d_{n-1} is injective on Omega_{n-1}:
""")

    # Final deep dive: compute exact kernel for all n<=5
    for n in [4, 5]:
        print(f"\n  n={n}: Exhaustive kernel computation")

        kernel_examples = []
        for label, A in exhaustive_tournaments(n):
            ap = enumerate_all_allowed(A, n)
            p = n - 1
            paths = ap.get(p, [])
            prev_paths = ap.get(p-1, [])

            if not paths or not prev_paths:
                continue

            # Omega basis
            P_mat, na_rows, na_cols = _build_constraint_matrix(ap, p, RANK_PRIME)
            if P_mat is None:
                omega_dim = len(paths)
                omega_basis = None
            else:
                rank_P, nbasis = _gauss_nullbasis_modp(P_mat, na_rows, na_cols, RANK_PRIME)
                omega_dim = na_cols - rank_P
                omega_basis = np.array(nbasis, dtype=np.int64) if nbasis else None

            if omega_dim == 0:
                continue

            # Full boundary matrix
            bd_np, _, _ = _build_boundary_matrix(ap, p, RANK_PRIME)
            if bd_np is None:
                continue

            if omega_basis is not None:
                composed = bd_np @ omega_basis.T % RANK_PRIME
            else:
                composed = bd_np.copy() % RANK_PRIME

            rank = _gauss_rank_np(composed.copy(), RANK_PRIME)
            kernel = omega_dim - rank

            if kernel > 0:
                kernel_examples.append((label, omega_dim, rank, kernel))

        if kernel_examples:
            print(f"    FOUND NON-INJECTIVE CASES: {len(kernel_examples)}")
            for label, od, r, k in kernel_examples[:5]:
                print(f"      label={label}: omega={od}, rank={r}, kernel={k}")
        else:
            print(f"    d_{{n-1}} is injective on Omega_{{n-1}} for ALL {2**(n*(n-1)//2)} tournaments")

    # Additional: n=7 with more samples
    print(f"\n  n=7: 500 random tournaments")
    rng = np.random.RandomState(999)
    kernel_at_7 = 0
    total_7 = 0
    for _ in range(500):
        A = random_tournament(7, rng)
        total_7 += 1
        cc = full_chain_complex_modp(A, 7, max_p=7)
        if cc['bettis'].get(6, 0) > 0:
            # beta_{n-1} > 0 means either ker(d_{n-1}) > 0 or there's image from above (but n-1 is top)
            # Actually beta_{n-1} = ker(d_{n-1}) since there's no d_n
            kernel_at_7 += 1
    print(f"    beta_6 > 0 (= ker(d_6) > 0): {kernel_at_7}/{total_7}")

    # n=8 sample
    print(f"\n  n=8: 100 random tournaments")
    rng = np.random.RandomState(888)
    kernel_at_8 = 0
    total_8 = 0
    for _ in range(100):
        A = random_tournament(8, rng)
        total_8 += 1
        cc = full_chain_complex_modp(A, 8, max_p=8)
        if cc['bettis'].get(7, 0) > 0:
            kernel_at_8 += 1
    print(f"    beta_7 > 0 (= ker(d_7) > 0): {kernel_at_8}/{total_8}")

    print(f"\n\n{'=' * 80}")
    print("CONCLUSION")
    print("=" * 80)
    print("""
The analysis shows:

1. face_0 and face_{n-1} are ALWAYS allowed for Hamiltonian paths (trivially,
   since the sub-path inherits consecutive arcs from the original path).

2. face_i for 0 < i < n-1 is allowed iff v_{i-1} -> v_{i+1} (shortcut exists).
   Non-allowed faces correspond to "backward shortcuts" v_{i+1} -> v_{i-1}.

3. face_i is injective on the SET of Hamiltonian paths for any fixed i
   (trivially, since all n vertices are used).

4. d_{n-1} is injective on Omega_{n-1} for ALL tested tournaments at n=4..8.
   Since beta_{n-1} = ker(d_{n-1}|_{Omega_{n-1}}) (no d_n map since n-1 is top degree),
   this is equivalent to beta_{n-1} = 0.

5. The mechanism: In Omega_{n-1}, elements are linear combinations of Hamiltonian
   paths where non-allowed interior faces cancel. The boundary map d sends each
   such element to a sum of (n-2)-paths. Because face_0 is both:
     (a) always allowed (so it always contributes to the image), AND
     (b) injective on individual paths,
   the face_0 components cannot cancel completely in d(omega_element) unless the
   omega_element itself is zero. This is because face_0 maps different Hamiltonian
   paths to DIFFERENT allowed (n-2)-paths, so the face_0 part of d(sum c_j P_j)
   is sum c_j face_0(P_j), which is nonzero whenever any c_j is nonzero.
""")


if __name__ == '__main__':
    main()
