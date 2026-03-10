"""
hyp398_test_n8.py — Test HYP-398 at n=8 for beta_3=2 tournaments

HYP-398: "NEW BOUNDARIES TARGET ONLY NEW CYCLES"
  For beta_3=1 at n=7: im(d_4^new) ∩ span(old ker_d3) ⊂ im(d_4^old)
  New d_4 content (5-paths through v) never kills old cycles beyond what
  im(d_4^old) already kills.

Does HYP-398 FAIL at n=8 for beta_3=2 tournaments?
If yes, then new 5-paths through v DO kill old 3-cycles that weren't killed before,
which is the mechanism allowing beta_3=2.

Method:
1. Find a beta_3=2 tournament T at n=8
2. For each vertex v:
   a. Compute Z_3(T) = ker(d_3) in Omega_3(T)
   b. Identify "old" Z_3: cycles not involving v (= ker(d_3) in Omega_3(T\v), embedded)
   c. Identify "new" im(d_4): boundaries from 5-paths through v
   d. Test: does im(d_4^new) kill any old cycle that im(d_4^old) doesn't?

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import gc
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, enumerate_allowed_paths, build_adj_lists,
    boundary_faces, _gauss_rank_np, _gauss_nullbasis_modp,
    RANK_PRIME
)
from beta3_lean import fast_beta3_lean, _gauss_rank_i32, _gauss_nullbasis_i32

PRIME = 32749  # Small prime for int32


def build_omega_and_boundary(A, n, prime=PRIME):
    """Build Omega_3, Omega_4, d_3, d_4 matrices and bases."""
    adj = build_adj_lists(A, n)
    ap = {}
    for p in range(min(6, n)):
        ap[p] = enumerate_allowed_paths(A, n, p, adj)

    paths_sets = {p: set(ap.get(p, [])) for p in range(6)}

    def build_constraint(paths_p, paths_pm1_set):
        if not paths_p:
            return None, 0, 0
        non_allowed = {}
        na_count = 0
        entries = []
        for j, path in enumerate(paths_p):
            for sign, face in boundary_faces(path):
                if len(set(face)) == len(face) and face not in paths_pm1_set:
                    if face not in non_allowed:
                        non_allowed[face] = na_count
                        na_count += 1
                    entries.append((non_allowed[face], j, sign))
        if na_count == 0:
            return None, 0, len(paths_p)
        P = np.zeros((na_count, len(paths_p)), dtype=np.int32)
        for row, col, sign in entries:
            P[row, col] = (P[row, col] + sign) % prime
        return P, na_count, len(paths_p)

    # Get Omega_3 basis
    P3, na3, nc3 = build_constraint(ap[3], paths_sets[2])
    if P3 is None:
        omega3_basis = np.eye(len(ap.get(3, [])), dtype=np.int32)
        dim_omega3 = len(ap.get(3, []))
    else:
        rank_P3, nbasis3 = _gauss_nullbasis_i32(P3, na3, nc3, prime)
        dim_omega3 = nc3 - rank_P3
        omega3_basis = np.array(nbasis3, dtype=np.int32) if nbasis3 else None

    # Get Omega_4 basis
    P4, na4, nc4 = build_constraint(ap.get(4, []), paths_sets[3])
    if ap.get(4, []) and P4 is not None:
        rank_P4, nbasis4 = _gauss_nullbasis_i32(P4, na4, nc4, prime)
        dim_omega4 = nc4 - rank_P4
        omega4_basis = np.array(nbasis4, dtype=np.int32) if nbasis4 else None
    elif ap.get(4, []):
        omega4_basis = np.eye(len(ap[4]), dtype=np.int32)
        dim_omega4 = len(ap[4])
    else:
        omega4_basis = None
        dim_omega4 = 0

    # Build boundary d_3: A_3 -> A_2
    paths_2 = ap.get(2, [])
    paths_3 = ap.get(3, [])
    idx2 = {p: i for i, p in enumerate(paths_2)}
    bd3 = np.zeros((len(paths_2), len(paths_3)), dtype=np.int32)
    for j, path in enumerate(paths_3):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] = (bd3[idx2[face], j] + sign) % prime

    # Build boundary d_4: A_4 -> A_3
    paths_4 = ap.get(4, [])
    idx3 = {p: i for i, p in enumerate(paths_3)}
    bd4 = np.zeros((len(paths_3), len(paths_4)), dtype=np.int32)
    for j, path in enumerate(paths_4):
        for sign, face in boundary_faces(path):
            if face in idx3:
                bd4[idx3[face], j] = (bd4[idx3[face], j] + sign) % prime

    return {
        'ap': ap,
        'omega3_basis': omega3_basis,
        'omega4_basis': omega4_basis,
        'dim_omega3': dim_omega3,
        'dim_omega4': dim_omega4,
        'bd3': bd3,
        'bd4': bd4,
        'paths_3': paths_3,
        'paths_4': paths_4,
    }


def test_hyp398(A, n, v, data, prime=PRIME):
    """Test HYP-398 for vertex v in tournament A.

    Returns dict with:
      - old_ker_d3_dim: dim of ker(d_3) restricted to "old" (not involving v) Omega_3
      - old_imd4_dim: dim of im(d_4^old) restricted to old ker_d3
      - new_imd4_kills_old: dim of im(d_4^new) ∩ old_ker_d3 / im(d_4^old)
      - hyp398_holds: True if new boundaries don't kill old cycles
    """
    paths_3 = data['paths_3']
    paths_4 = data['paths_4']
    omega3_basis = data['omega3_basis']
    omega4_basis = data['omega4_basis']
    bd3 = data['bd3']
    bd4 = data['bd4']

    if omega3_basis is None or len(omega3_basis) == 0:
        return {'hyp398_holds': True, 'trivial': True}

    # Classify 3-paths as "old" (not involving v) or "new" (involving v)
    old_3_idx = [j for j, p in enumerate(paths_3) if v not in p]
    new_3_idx = [j for j, p in enumerate(paths_3) if v in p]

    # Classify 4-paths similarly
    old_4_idx = [j for j, p in enumerate(paths_4) if v not in p]
    new_4_idx = [j for j, p in enumerate(paths_4) if v in p]

    # Omega_3 basis vectors in A_3 coordinates
    # omega3_basis: (dim_omega3 x |A_3|) matrix

    # Compute d_3 restricted to Omega_3 (in old 3-path coordinates)
    d3_omega = (bd3.astype(np.int64) @ omega3_basis.T.astype(np.int64) % prime).astype(np.int32)

    # Compute ker(d_3) in Omega_3
    d3_list = [[int(d3_omega[i, j]) for j in range(d3_omega.shape[1])]
               for i in range(d3_omega.shape[0])]
    rk_d3, null_d3 = _gauss_nullbasis_i32(
        np.array(d3_list, dtype=np.int32) if d3_list else np.zeros((0, 0), dtype=np.int32),
        d3_omega.shape[0], d3_omega.shape[1], prime
    )
    ker_d3_dim = data['dim_omega3'] - rk_d3

    if ker_d3_dim == 0:
        return {'hyp398_holds': True, 'ker_d3_dim': 0}

    # Z_3 = ker(d_3) in Omega_3, in A_3 coordinates
    Z3_omega = np.array(null_d3, dtype=np.int32) if null_d3 else np.zeros((0, data['dim_omega3']), dtype=np.int32)
    Z3_A3 = (Z3_omega.astype(np.int64) @ omega3_basis.astype(np.int64) % prime).astype(np.int32)

    # Project Z_3 to "old" coordinates only
    Z3_old = Z3_A3[:, old_3_idx] if old_3_idx else np.zeros((Z3_A3.shape[0], 0), dtype=np.int32)

    # Compute im(d_4) restricted to Omega_4, in A_3 coordinates
    if omega4_basis is None or len(omega4_basis) == 0:
        im_d4_A3 = np.zeros((len(paths_3), 0), dtype=np.int32)
    else:
        im_d4_A3 = (bd4.astype(np.int64) @ omega4_basis.T.astype(np.int64) % prime).astype(np.int32)

    # Split im(d_4) by old/new 4-paths
    # old im(d_4): from 4-paths not involving v
    if omega4_basis is not None and len(omega4_basis) > 0:
        # We need to identify which Omega_4 basis vectors correspond to "old" 4-paths
        # Actually, the Omega_4 basis is a linear combination of ALL 4-paths,
        # so we can't simply split by old/new. Instead:

        # im(d_4^old): boundaries of 4-paths not involving v, restricted to Omega_4
        old_bd4 = bd4[:, old_4_idx] if old_4_idx else np.zeros((len(paths_3), 0), dtype=np.int32)
        new_bd4 = bd4[:, new_4_idx] if new_4_idx else np.zeros((len(paths_3), 0), dtype=np.int32)

        # Project old_bd4 to old 3-path coordinates
        old_bd4_old = old_bd4[old_3_idx, :] if old_3_idx else np.zeros((0, old_bd4.shape[1]), dtype=np.int32)
        new_bd4_old = new_bd4[old_3_idx, :] if old_3_idx else np.zeros((0, new_bd4.shape[1]), dtype=np.int32)

        # Now check: does new_bd4_old add anything to old_bd4_old in the Z3_old space?
        # Compute rank of old_bd4_old
        if old_bd4_old.shape[1] > 0:
            rk_old = _gauss_rank_i32(old_bd4_old.copy(), prime)
        else:
            rk_old = 0

        # Compute rank of [old_bd4_old | new_bd4_old]
        combined_old = np.concatenate([old_bd4_old, new_bd4_old], axis=1) if (
            old_bd4_old.shape[1] > 0 or new_bd4_old.shape[1] > 0
        ) else np.zeros((len(old_3_idx), 0), dtype=np.int32)
        if combined_old.shape[1] > 0:
            rk_combined = _gauss_rank_i32(combined_old.copy(), prime)
        else:
            rk_combined = 0

        new_kills = rk_combined - rk_old
    else:
        new_kills = 0
        rk_old = 0
        rk_combined = 0

    hyp398_holds = (new_kills == 0)

    return {
        'hyp398_holds': hyp398_holds,
        'ker_d3_dim': ker_d3_dim,
        'old_3_count': len(old_3_idx),
        'new_3_count': len(new_3_idx),
        'old_4_count': len(old_4_idx),
        'new_4_count': len(new_4_idx),
        'rk_old_bd4': rk_old,
        'rk_combined_bd4': rk_combined,
        'new_kills': new_kills,
    }


def main():
    print("=" * 70)
    print("HYP-398 TEST: NEW BOUNDARIES TARGET ONLY NEW CYCLES")
    print("=" * 70)

    # Phase 1: Test at n=8 for beta_3=2 tournaments
    print("\n--- Phase 1: n=8, finding beta_3=2 tournaments ---")
    n = 8
    rng = np.random.RandomState(12345)

    b3_2_tours = []
    b3_1_tours = []

    t0 = time.time()
    for trial in range(5000):
        A = random_tournament(n, rng)
        gc.collect()
        b3 = fast_beta3_lean(A, n)

        if b3 == 2 and len(b3_2_tours) < 3:
            b3_2_tours.append((trial, A.copy()))
            scores = sorted([int(sum(A[i])) for i in range(n)])
            print(f"  Found beta_3=2 at trial {trial}, scores={scores}")
        elif b3 == 1 and len(b3_1_tours) < 2:
            b3_1_tours.append((trial, A.copy()))

        if len(b3_2_tours) >= 3 and len(b3_1_tours) >= 2:
            break

    elapsed = time.time() - t0
    print(f"  Found {len(b3_2_tours)} b3=2 and {len(b3_1_tours)} b3=1 in {elapsed:.1f}s")

    # Phase 2: Test HYP-398 on beta_3=2 tournaments
    print(f"\n{'='*60}")
    print("TESTING HYP-398 ON BETA_3=2 TOURNAMENTS (n=8)")
    print(f"{'='*60}")

    for trial, A in b3_2_tours:
        scores = sorted([int(sum(A[i])) for i in range(n)])
        print(f"\n  Trial {trial}, scores={scores}")

        data = build_omega_and_boundary(A, n)
        print(f"    dim(Omega_3) = {data['dim_omega3']}, dim(Omega_4) = {data['dim_omega4']}")
        print(f"    |A_3| = {len(data['paths_3'])}, |A_4| = {len(data['paths_4'])}")

        for v in range(n):
            gc.collect()
            result = test_hyp398(A, n, v, data)

            out_deg = int(sum(A[v]))
            b3_v = fast_beta3_lean(
                np.array([[A[r][c] for c in range(n) if c != v]
                          for r in range(n) if r != v], dtype=np.int8), n-1
            )

            status = "HOLDS" if result.get('hyp398_holds', True) else "**FAILS**"
            if result.get('trivial'):
                status = "TRIVIAL"

            extra = ""
            if 'new_kills' in result:
                extra = (f" old_rk={result['rk_old_bd4']}, "
                        f"combined_rk={result['rk_combined_bd4']}, "
                        f"new_kills={result['new_kills']}")

            tag = "GOOD" if b3_v == 0 else f"BAD(b3={b3_v})"
            print(f"    v={v}(out={out_deg}): {status} [{tag}]{extra}")

    # Phase 3: Test HYP-398 on beta_3=1 tournaments for comparison
    print(f"\n{'='*60}")
    print("TESTING HYP-398 ON BETA_3=1 TOURNAMENTS (n=8, comparison)")
    print(f"{'='*60}")

    for trial, A in b3_1_tours[:1]:
        scores = sorted([int(sum(A[i])) for i in range(n)])
        print(f"\n  Trial {trial}, scores={scores}")

        data = build_omega_and_boundary(A, n)
        print(f"    dim(Omega_3) = {data['dim_omega3']}, dim(Omega_4) = {data['dim_omega4']}")

        for v in range(n):
            gc.collect()
            result = test_hyp398(A, n, v, data)

            out_deg = int(sum(A[v]))
            b3_v = fast_beta3_lean(
                np.array([[A[r][c] for c in range(n) if c != v]
                          for r in range(n) if r != v], dtype=np.int8), n-1
            )

            status = "HOLDS" if result.get('hyp398_holds', True) else "**FAILS**"
            if result.get('trivial'):
                status = "TRIVIAL"

            tag = "GOOD" if b3_v == 0 else f"BAD(b3={b3_v})"
            extra = ""
            if 'new_kills' in result:
                extra = f" new_kills={result['new_kills']}"
            print(f"    v={v}(out={out_deg}): {status} [{tag}]{extra}")

    print(f"\n{'='*70}")
    print("SUMMARY")
    print("  HYP-398 holds at n=7 (34/34 BAD vertices, opus-S55)")
    print("  Question: does HYP-398 FAIL at n=8 for beta_3=2 tournaments?")
    print("  If yes: new 5-paths through v kill old 3-cycles => beta_3 can exceed 1")
    print(f"{'='*70}")
    print("DONE.")


if __name__ == '__main__':
    main()
