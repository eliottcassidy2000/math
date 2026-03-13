#!/usr/bin/env python3
"""
eigenspace_uniformity_analysis.py

*** WARNING: THIS SCRIPT'S CONCLUSION IS INCORRECT ***
It uses interior-only boundary (TRH convention: removes vertices 1..m-1),
NOT the full GLMY boundary (removes ALL vertices 0..m including endpoints).
GLMY face_0 removes the first vertex, shifting the face to start at d_1 ≠ 0,
which introduces a k-dependent phase omega^{k*d_1}. Therefore B_k ≠ B_0.

Our per-eigenspace data CONFIRMS B_k ≠ B_0:
  P_11 k=0: rk(d) = [0, 0, 5, 15, 55, 150, 305, 390]
  P_11 k≠0: rk(d) = [0, 1, 4, 16, 54, 151, 309, 390]

See MISTAKE-021 for the TRH vs GLMY boundary confusion.

Original (incorrect) description:
Investigate WHY all eigenspaces k=0,...,p-1 of circulant tournaments on Z_p
have IDENTICAL Betti numbers (Eigenspace Betti Uniformity).

For k!=0: Galois symmetry explains it (gcd(j,p)=1 => eigenspaces k and jk isomorphic).
For k=0 matching k!=0: UNEXPLAINED.

Strategy:
1. For Interval n=7, explicitly compute projected boundary matrices B_k at each degree.
2. Compare B_0 vs B_1: singular values, ranks, structural relationships.
3. Look for a transformation relating B_0 to B_1.

Key notation:
- Orbit reps = paths starting from vertex 0
- B_k[beta, alpha] = sum over boundary faces of alpha, with phase omega^{k*shift}
  where shift rotates the face back to start from 0.
"""

import numpy as np
from numpy.linalg import svd, matrix_rank
from itertools import combinations

def circulant_tournament(n, S):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i][j] = 1
    return A

def get_regular_paths(A, m):
    """Get all regular m-paths (length-m chains of vertices with the regularity condition)."""
    n = A.shape[0]
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path:
                continue
            if not A[last][v]:
                continue
            if depth >= 1 and not A[prev][v]:
                continue
            path.append(v)
            dfs(path, depth + 1, last)
            path.pop()
    for start in range(n):
        dfs([start], 0, -1)
    return paths

def build_eigenspace_boundary(n, S, k, m, A=None):
    """Build the projected boundary matrix B_k for degree m.

    Returns B_k, paths_m_from0, paths_m1_from0.

    B_k maps m-path orbit reps to (m-1)-path orbit reps in eigenspace k.
    """
    if A is None:
        A = circulant_tournament(n, S)

    all_paths_m = get_regular_paths(A, m)
    all_paths_m1 = get_regular_paths(A, m - 1) if m >= 1 else []

    paths_m = [p for p in all_paths_m if p[0] == 0]
    paths_m1 = [p for p in all_paths_m1 if p[0] == 0]

    if not paths_m or not paths_m1:
        return None, paths_m, paths_m1

    omega = np.exp(2j * np.pi / n)
    path_to_idx = {p: i for i, p in enumerate(paths_m1)}

    B_k = np.zeros((len(paths_m1), len(paths_m)), dtype=complex)
    path_len = len(paths_m[0]) - 1  # = m

    for j, path in enumerate(paths_m):
        for i in range(1, path_len):  # interior vertices only
            face = path[:i] + path[i+1:]
            sign = (-1) ** i
            shift = face[0]
            canonical_face = tuple((v - shift) % n for v in face)
            phase = omega ** (k * shift)
            if canonical_face in path_to_idx:
                B_k[path_to_idx[canonical_face], j] += sign * phase

    return B_k, paths_m, paths_m1


def analyze_eigenspace_uniformity():
    n = 7
    S = {1, 2, 3}  # Interval tournament
    A = circulant_tournament(n, S)
    omega = np.exp(2j * np.pi / n)

    print("=" * 80)
    print("EIGENSPACE UNIFORMITY ANALYSIS: Interval n=7, S={1,2,3}")
    print("=" * 80)

    # ================================================================
    # PART 1: Compute and display B_k for all k at each degree
    # ================================================================
    print("\n" + "=" * 80)
    print("PART 1: Projected boundary matrices B_k")
    print("=" * 80)

    all_Bk = {}  # (k, m) -> matrix

    for m in range(1, n):
        print(f"\n--- Degree m={m} ---")

        # First compute B_0 and B_1 to compare
        B0, pm, pm1 = build_eigenspace_boundary(n, S, 0, m, A)
        B1, _, _ = build_eigenspace_boundary(n, S, 1, m, A)

        if B0 is None:
            print(f"  Empty (no paths)")
            continue

        nrows, ncols = B0.shape
        print(f"  Matrix size: {nrows} x {ncols}")
        print(f"  Paths at degree {m} (from 0): {[p for p in pm]}")
        print(f"  Paths at degree {m-1} (from 0): {[p for p in pm1]}")

        # Compute all B_k
        svd_data = {}
        for k in range(n):
            Bk, _, _ = build_eigenspace_boundary(n, S, k, m, A)
            all_Bk[(k, m)] = Bk
            if Bk is not None:
                U, s, Vh = svd(Bk, full_matrices=False)
                svd_data[k] = s
                rk = np.sum(s > 1e-8)
            else:
                svd_data[k] = np.array([])
                rk = 0

        # Print B_0 and B_1 entries explicitly
        print(f"\n  B_0 (k=0, real matrix since no phase):")
        B0_real = np.real(B0)
        for i in range(nrows):
            row_str = "    ["
            for j in range(ncols):
                val = B0_real[i, j]
                row_str += f" {val:6.2f}"
            row_str += " ]"
            print(row_str)

        print(f"\n  B_1 (k=1, complex):")
        for i in range(nrows):
            row_str = "    ["
            for j in range(ncols):
                val = B1[i, j]
                if abs(val.imag) < 1e-10:
                    row_str += f" {val.real:6.2f}"
                else:
                    row_str += f" {val:.2f}"
            row_str += " ]"
            print(row_str)

        # Compare singular values
        print(f"\n  Singular values comparison:")
        for k in range(n):
            s = svd_data[k]
            s_str = ", ".join(f"{v:.6f}" for v in s)
            print(f"    k={k}: [{s_str}]")

        # Check if all k have same singular values
        ref_svs = svd_data[0]
        all_same_sv = True
        for k in range(1, n):
            if len(svd_data[k]) != len(ref_svs):
                all_same_sv = False
                break
            if not np.allclose(sorted(svd_data[k], reverse=True),
                              sorted(ref_svs, reverse=True), atol=1e-6):
                all_same_sv = False
                break
        print(f"\n  All eigenspaces have same singular values? {all_same_sv}")

    # ================================================================
    # PART 2: Deeper structural comparison B_0 vs B_1
    # ================================================================
    print("\n" + "=" * 80)
    print("PART 2: Structural comparison B_0 vs B_k")
    print("=" * 80)

    for m in range(1, n):
        B0 = all_Bk.get((0, m))
        B1 = all_Bk.get((1, m))
        if B0 is None or B1 is None:
            continue

        print(f"\n--- Degree m={m} ---")
        nrows, ncols = B0.shape

        # Check 1: Are B_0^H B_0 and B_1^H B_1 the same? (Gram matrix)
        G0 = B0.conj().T @ B0
        G1 = B1.conj().T @ B1
        print(f"  B_0^H B_0 == B_1^H B_1 ? {np.allclose(G0, G1, atol=1e-8)}")

        # Check 2: Are B_0 B_0^H and B_1 B_1^H the same?
        G0r = B0 @ B0.conj().T
        G1r = B1 @ B1.conj().T
        print(f"  B_0 B_0^H == B_1 B_1^H ? {np.allclose(G0r, G1r, atol=1e-8)}")

        # Check 3: Is B_1 = D_row @ B_0 @ D_col for diagonal unitary matrices?
        # If B_1[i,j] = d_row[i] * B_0[i,j] * d_col[j], then
        # d_row[i] * d_col[j] = B_1[i,j] / B_0[i,j] wherever B_0[i,j] != 0
        print(f"\n  Looking for diagonal transformation B_1 = D_row * B_0 * D_col:")
        ratios = np.full_like(B0, np.nan, dtype=complex)
        for i in range(nrows):
            for j in range(ncols):
                if abs(B0[i, j]) > 1e-10:
                    ratios[i, j] = B1[i, j] / B0[i, j]

        # If diagonal transform exists, ratios[i,j] = d_row[i] * d_col[j]
        # So ratios[i,j] * ratios[i',j'] = ratios[i,j'] * ratios[i',j]
        # for all i,i',j,j' where defined
        found_diag = True
        for i1 in range(nrows):
            for i2 in range(nrows):
                for j1 in range(ncols):
                    for j2 in range(ncols):
                        if (not np.isnan(ratios[i1,j1]) and not np.isnan(ratios[i2,j2]) and
                            not np.isnan(ratios[i1,j2]) and not np.isnan(ratios[i2,j1])):
                            prod1 = ratios[i1,j1] * ratios[i2,j2]
                            prod2 = ratios[i1,j2] * ratios[i2,j1]
                            if abs(prod1 - prod2) > 1e-6:
                                found_diag = False

        if found_diag:
            print(f"    YES! B_1 = D_row * B_0 * D_col for diagonal unitaries")
            # Extract D_row and D_col
            # Fix d_col[0] = 1, then d_row[i] = ratios[i,0] (if defined)
            # and d_col[j] = ratios[i,j] / d_row[i] (for any valid i)
            d_col = np.ones(ncols, dtype=complex)
            d_row = np.ones(nrows, dtype=complex)
            # Find a reference column with known ratio
            ref_col = None
            for j in range(ncols):
                valid = [i for i in range(nrows) if not np.isnan(ratios[i, j])]
                if valid:
                    ref_col = j
                    break
            if ref_col is not None:
                for i in range(nrows):
                    if not np.isnan(ratios[i, ref_col]):
                        d_row[i] = ratios[i, ref_col]
                for j in range(ncols):
                    for i in range(nrows):
                        if not np.isnan(ratios[i, j]) and abs(d_row[i]) > 1e-10:
                            d_col[j] = ratios[i, j] / d_row[i]
                            break

                # Normalize so d_col[ref_col] = 1
                # already done since d_row was set from ref_col

                print(f"    D_row phases: {[f'{np.angle(d)/(np.pi):.4f}pi' for d in d_row]}")
                print(f"    D_col phases: {[f'{np.angle(d)/(np.pi):.4f}pi' for d in d_col]}")
                print(f"    D_row magnitudes: {[f'{abs(d):.6f}' for d in d_row]}")
                print(f"    D_col magnitudes: {[f'{abs(d):.6f}' for d in d_col]}")

                # Verify
                B1_reconstructed = np.diag(d_row) @ B0 @ np.diag(d_col)
                print(f"    Reconstruction error: {np.max(np.abs(B1_reconstructed - B1)):.2e}")
        else:
            print(f"    NO diagonal transformation exists")

        # Check 4: Is there a PERMUTATION relating them?
        # B_1 = P_row @ B_0 @ P_col ?
        # This would mean same sparsity pattern (up to permutation)

        # Check 5: Frobenius norm comparison
        print(f"\n  Frobenius norms: ||B_0||={np.linalg.norm(B0, 'fro'):.6f}, "
              f"||B_1||={np.linalg.norm(B1, 'fro'):.6f}")

    # ================================================================
    # PART 3: The algebraic structure of B_k
    # ================================================================
    print("\n" + "=" * 80)
    print("PART 3: Understanding B_k algebraically")
    print("=" * 80)

    print("\nThe boundary matrix B_k[beta, alpha] for eigenspace k has entries:")
    print("  B_k[beta, alpha] = sum over interior faces of alpha that canonicalize to beta")
    print("  Each term = (-1)^i * omega^{k * shift}")
    print("  where shift = starting vertex of the face before canonicalization")
    print("")
    print("For k=0: all phases = 1, so B_0 is a REAL INTEGER matrix")
    print("For k!=0: phases are powers of omega^k")
    print("")
    print("KEY QUESTION: Why does rank(B_0) = rank(B_k)?")

    # Let's look at the relationship more carefully
    # For each nonzero entry B_k[beta, alpha], the entry is sum of (-1)^i * omega^{k*s_i}
    # where s_i are the shifts for that (beta, alpha) pair.

    # Decompose each entry into its shift contributions
    for m in [2, 3, 4]:
        print(f"\n--- Degree m={m}: Entry decomposition ---")

        all_paths_m = get_regular_paths(A, m)
        all_paths_m1 = get_regular_paths(A, m - 1) if m >= 1 else []
        paths_m = [p for p in all_paths_m if p[0] == 0]
        paths_m1 = [p for p in all_paths_m1 if p[0] == 0]

        if not paths_m or not paths_m1:
            print(f"  Empty")
            continue

        path_to_idx = {p: i for i, p in enumerate(paths_m1)}
        path_len = len(paths_m[0]) - 1

        # For each (row, col) pair, collect ALL (sign, shift) contributions
        entry_data = {}  # (row_idx, col_idx) -> list of (sign, shift)

        for j, path in enumerate(paths_m):
            for i in range(1, path_len):
                face = path[:i] + path[i+1:]
                sign = (-1) ** i
                shift = face[0]
                canonical_face = tuple((v - shift) % n for v in face)
                if canonical_face in path_to_idx:
                    row = path_to_idx[canonical_face]
                    if (row, j) not in entry_data:
                        entry_data[(row, j)] = []
                    entry_data[(row, j)].append((sign, shift))

        # Show the shift structure
        print(f"  Nonzero entries and their shift decomposition:")
        for (row, col), contribs in sorted(entry_data.items()):
            beta_path = paths_m1[row]
            alpha_path = paths_m[col]

            # Compute entry for each k
            entries_by_k = []
            for k in range(n):
                val = sum(sign * omega**(k * shift) for sign, shift in contribs)
                entries_by_k.append(val)

            contribs_str = ", ".join(f"({s:+d}, shift={sh})" for s, sh in contribs)
            print(f"    ({row},{col}): beta={beta_path}, alpha={alpha_path}")
            print(f"      contributions: [{contribs_str}]")

            # Key observation: how many distinct shifts per entry?
            shifts = set(sh for _, sh in contribs)
            print(f"      distinct shifts: {sorted(shifts)} (count={len(shifts)})")

            # Show value at k=0 and k=1
            v0 = sum(sign for sign, shift in contribs)
            print(f"      B_0 entry = {v0}")
            v1 = entries_by_k[1]
            print(f"      B_1 entry = {v1:.4f}")

    # ================================================================
    # PART 4: Does B_k = D * B_0 * E for ALL k simultaneously?
    # ================================================================
    print("\n" + "=" * 80)
    print("PART 4: Universal diagonal gauge transformation")
    print("=" * 80)

    print("\nIf B_k = D_k * B_0 * E_k for diagonal D_k, E_k, then:")
    print("  rank(B_k) = rank(B_0) automatically (diagonal matrices are invertible)")
    print("  This would PROVE eigenspace Betti uniformity!")

    for m in range(1, n):
        Bks = [all_Bk.get((k, m)) for k in range(n)]
        if Bks[0] is None:
            continue

        print(f"\n--- Degree m={m} ---")
        B0 = Bks[0]

        for k in range(1, n):
            Bk = Bks[k]
            if Bk is None:
                continue

            # Check diagonal gauge: B_k[i,j] / B_0[i,j] = d_row[i] * d_col[j]?
            nrows, ncols = B0.shape
            valid_ratios = []
            for i in range(nrows):
                for j in range(ncols):
                    if abs(B0[i, j]) > 1e-10:
                        ratio = Bk[i, j] / B0[i, j]
                        valid_ratios.append((i, j, ratio))

            if not valid_ratios:
                print(f"  k={k}: No nonzero entries in B_0 to compare")
                continue

            # Check rank-1 structure of ratio matrix
            # If ratios[i,j] = d[i]*e[j], then the ratio matrix has rank 1
            ratio_matrix = np.zeros((nrows, ncols), dtype=complex)
            mask = np.zeros((nrows, ncols), dtype=bool)
            for i, j, r in valid_ratios:
                ratio_matrix[i, j] = r
                mask[i, j] = True

            # Extract submatrix of ratios where both B_0 and B_k are nonzero
            rows_used = sorted(set(i for i, j, r in valid_ratios))
            cols_used = sorted(set(j for i, j, r in valid_ratios))

            sub = ratio_matrix[np.ix_(rows_used, cols_used)]
            sub_mask = mask[np.ix_(rows_used, cols_used)]

            # Check if ratio matrix has rank 1 (ignoring zeros from masking)
            # A cleaner check: for all nonzero entries, r[i,j]*r[i',j'] = r[i,j']*r[i',j]
            is_rank1 = True
            for idx1, (i1, j1, r1) in enumerate(valid_ratios):
                for idx2, (i2, j2, r2) in enumerate(valid_ratios):
                    if (i1, j2) in {(i, j) for i, j, _ in valid_ratios} and \
                       (i2, j1) in {(i, j) for i, j, _ in valid_ratios}:
                        r12 = ratio_matrix[i1, j2]
                        r21 = ratio_matrix[i2, j1]
                        if abs(r1 * r2 - r12 * r21) > 1e-6:
                            is_rank1 = False
                            break
                if not is_rank1:
                    break

            if is_rank1:
                # Extract the diagonal factors
                # Fix e[cols_used[0]] = 1
                ref_j = cols_used[0]
                d_row = {}
                for i, j, r in valid_ratios:
                    if j == ref_j:
                        d_row[i] = r

                d_col = {}
                d_col[ref_j] = 1.0
                for i, j, r in valid_ratios:
                    if i in d_row and abs(d_row[i]) > 1e-10:
                        d_col[j] = r / d_row[i]

                # Show the phases
                d_row_phases = {i: np.angle(v) * n / (2 * np.pi) for i, v in d_row.items()}
                d_col_phases = {j: np.angle(v) * n / (2 * np.pi) for j, v in d_col.items()}

                print(f"  k={k}: DIAGONAL GAUGE EXISTS!")
                print(f"    row phases (units of 2pi/n): {d_row_phases}")
                print(f"    col phases (units of 2pi/n): {d_col_phases}")

                # Are these phases related to the path structure?
                # The phase on a path should be omega^{k * something_about_the_path}
                all_paths_m_from0 = [p for p in get_regular_paths(A, m) if p[0] == 0]
                all_paths_m1_from0 = [p for p in get_regular_paths(A, m - 1) if p[0] == 0]

                print(f"    Checking if col phases = k * (sum of vertices mod n):")
                for j, phase in d_col.items():
                    if j < len(all_paths_m_from0):
                        path = all_paths_m_from0[j]
                        vertex_sum = sum(path) % n
                        expected_phase = omega ** (k * vertex_sum)
                        match = abs(phase - expected_phase) < 1e-6
                        print(f"      col {j}: path={path}, sum_mod_n={vertex_sum}, "
                              f"phase={phase:.4f}, omega^(k*sum)={expected_phase:.4f}, match={match}")

                print(f"    Checking if row phases = k * (sum of vertices mod n):")
                for i, phase in d_row.items():
                    if i < len(all_paths_m1_from0):
                        path = all_paths_m1_from0[i]
                        vertex_sum = sum(path) % n
                        expected_phase = omega ** (k * vertex_sum)
                        match = abs(phase - expected_phase) < 1e-6
                        print(f"      row {i}: path={path}, sum_mod_n={vertex_sum}, "
                              f"phase={phase:.4f}, omega^(k*sum)={expected_phase:.4f}, match={match}")
            else:
                print(f"  k={k}: NO diagonal gauge (ratio matrix not rank 1)")
                # Show the problematic ratios
                print(f"    Sample ratios:")
                for i, j, r in valid_ratios[:10]:
                    print(f"      ({i},{j}): {r:.4f} (|r|={abs(r):.4f}, arg={np.angle(r)*n/(2*np.pi):.4f}*2pi/n)")

    # ================================================================
    # PART 5: Alternative — check if each entry has a SINGLE shift
    # ================================================================
    print("\n" + "=" * 80)
    print("PART 5: Single-shift entries analysis")
    print("=" * 80)
    print("\nIf every nonzero entry of B has exactly ONE contributing (sign, shift) pair,")
    print("then B_k[i,j] = sign * omega^{k*shift} = omega^{k*shift} * B_0[i,j]")
    print("which is exactly the diagonal gauge with d[i,j] depending only on shift.")

    for m in range(1, n):
        all_paths_m = get_regular_paths(A, m)
        all_paths_m1 = get_regular_paths(A, m - 1) if m >= 1 else []
        paths_m = [p for p in all_paths_m if p[0] == 0]
        paths_m1 = [p for p in all_paths_m1 if p[0] == 0]

        if not paths_m or not paths_m1:
            continue

        path_to_idx = {p: i for i, p in enumerate(paths_m1)}
        path_len = len(paths_m[0]) - 1

        entry_data = {}
        for j, path in enumerate(paths_m):
            for i in range(1, path_len):
                face = path[:i] + path[i+1:]
                sign = (-1) ** i
                shift = face[0]
                canonical_face = tuple((v - shift) % n for v in face)
                if canonical_face in path_to_idx:
                    row = path_to_idx[canonical_face]
                    if (row, j) not in entry_data:
                        entry_data[(row, j)] = []
                    entry_data[(row, j)].append((sign, shift))

        multi_shift = 0
        single_shift = 0
        for (row, col), contribs in entry_data.items():
            if len(contribs) == 1:
                single_shift += 1
            else:
                multi_shift += 1

        print(f"\n  Degree m={m}: {single_shift} single-shift entries, {multi_shift} multi-shift entries")

        if multi_shift > 0:
            print(f"    Multi-shift entries:")
            for (row, col), contribs in sorted(entry_data.items()):
                if len(contribs) > 1:
                    shifts = [sh for _, sh in contribs]
                    signs = [s for s, _ in contribs]
                    print(f"      ({row},{col}): shifts={shifts}, signs={signs}")
                    # For single-shift to work as diagonal gauge, we need
                    # sum(sign * omega^{k*shift}) / sum(sign) to factor as product of path-dependent phases
                    # This fails when there are multiple shifts... unless the shifts differ by a path-dependent amount

    # ================================================================
    # PART 6: The vertex-sum gauge conjecture
    # ================================================================
    print("\n" + "=" * 80)
    print("PART 6: Vertex-sum gauge conjecture")
    print("=" * 80)
    print("\nConjecture: B_k = D_k^row * B_0 * D_k^col")
    print("where D_k^row[i,i] = omega^{k * f(path_i)}")
    print("and D_k^col[j,j] = omega^{k * g(path_j)}")
    print("for some functions f, g of the path.")
    print("")
    print("Natural candidates for f, g:")
    print("  (a) sum of vertices mod n")
    print("  (b) sum of interior vertices mod n")
    print("  (c) last vertex mod n")
    print("  (d) first + last vertex mod n")

    for m in [2, 3, 4, 5]:
        B0 = all_Bk.get((0, m))
        B1 = all_Bk.get((1, m))
        if B0 is None or B1 is None:
            continue

        all_paths_m = get_regular_paths(A, m)
        all_paths_m1 = get_regular_paths(A, m - 1)
        paths_m = [p for p in all_paths_m if p[0] == 0]
        paths_m1 = [p for p in all_paths_m1 if p[0] == 0]

        nrows, ncols = B0.shape

        print(f"\n--- Degree m={m} ---")

        # Try each candidate for col gauge (m-paths)
        candidates_col = {
            "vertex_sum": lambda p: sum(p) % n,
            "interior_sum": lambda p: sum(p[1:-1]) % n if len(p) > 2 else 0,
            "last_vertex": lambda p: p[-1] % n,
            "sum_minus_first": lambda p: (sum(p) - p[0]) % n,
        }

        candidates_row = {
            "vertex_sum": lambda p: sum(p) % n,
            "interior_sum": lambda p: sum(p[1:-1]) % n if len(p) > 2 else 0,
            "last_vertex": lambda p: p[-1] % n,
            "sum_minus_first": lambda p: (sum(p) - p[0]) % n,
        }

        for cname_r, cfunc_r in candidates_row.items():
            for cname_c, cfunc_c in candidates_col.items():
                # Build gauge-transformed B_0
                d_row = np.array([omega ** (1 * cfunc_r(paths_m1[i])) for i in range(nrows)])
                d_col = np.array([omega ** (1 * cfunc_c(paths_m[j])) for j in range(ncols)])
                B1_test = np.diag(d_row) @ B0 @ np.diag(d_col)

                err = np.max(np.abs(B1_test - B1))
                if err < 1e-6:
                    print(f"  MATCH: row={cname_r}, col={cname_c}, error={err:.2e}")

    # ================================================================
    # PART 7: More general gauge — try omega^{k * last_vertex} on cols
    # and omega^{k * last_vertex} on rows
    # ================================================================
    print("\n" + "=" * 80)
    print("PART 7: Exhaustive gauge search (all linear combinations of vertex indices)")
    print("=" * 80)

    for m in [2, 3, 4, 5]:
        B0 = all_Bk.get((0, m))
        B1 = all_Bk.get((1, m))
        if B0 is None or B1 is None:
            continue

        all_paths_m = get_regular_paths(A, m)
        all_paths_m1 = get_regular_paths(A, m - 1)
        paths_m = [p for p in all_paths_m if p[0] == 0]
        paths_m1 = [p for p in all_paths_m1 if p[0] == 0]

        nrows, ncols = B0.shape

        print(f"\n--- Degree m={m} (ncols={ncols}, nrows={nrows}) ---")

        # For paths from vertex 0: path = (0, v_1, ..., v_m)
        # Try gauge: col_phase = sum of c_i * v_i for coefficients c_i in Z_n
        # row_phase = sum of r_i * v_i for coefficients r_i in Z_n
        # Since first vertex = 0, c_0 term vanishes.

        # For m-paths: try col_phase = a*v_last for a in Z_n
        # For (m-1)-paths: try row_phase = b*v_last for b in Z_n

        found = False
        for a in range(n):
            for b in range(n):
                d_col = np.array([omega ** (a * paths_m[j][-1]) for j in range(ncols)])
                d_row = np.array([omega ** (b * paths_m1[i][-1]) for i in range(nrows)])
                B1_test = np.diag(d_row) @ B0 @ np.diag(d_col)
                err = np.max(np.abs(B1_test - B1))
                if err < 1e-6:
                    print(f"  MATCH: row_phase=omega^({b}*v_last), col_phase=omega^({a}*v_last), err={err:.2e}")
                    found = True

        if not found:
            # Try with second-to-last vertex too
            for a1 in range(n):
                for a2 in range(n):
                    for b1 in range(n):
                        for b2 in range(n):
                            d_col = np.array([omega ** (a1 * paths_m[j][-1] + a2 * paths_m[j][-2] if len(paths_m[j]) >= 2 else 0)
                                             for j in range(ncols)])
                            d_row = np.array([omega ** (b1 * paths_m1[i][-1] + b2 * paths_m1[i][-2] if len(paths_m1[i]) >= 2 else 0)
                                             for i in range(nrows)])
                            B1_test = np.diag(d_row) @ B0 @ np.diag(d_col)
                            err = np.max(np.abs(B1_test - B1))
                            if err < 1e-6:
                                print(f"  MATCH: row_phase=omega^({b1}*v_last+{b2}*v_{-2}), "
                                      f"col_phase=omega^({a1}*v_last+{a2}*v_{-2}), err={err:.2e}")
                                found = True
                                break
                        if found:
                            break
                    if found:
                        break
                if found:
                    break

        if not found:
            print(f"  No simple linear gauge found")
            # Fall back: show the actual diagonal factors from Part 4 analysis
            # and try to identify them as a function of the path
            print(f"  Trying fully general per-path gauge extraction...")
            # Use Part 2 method
            nrows, ncols = B0.shape
            valid_ratios = []
            for i in range(nrows):
                for j in range(ncols):
                    if abs(B0[i, j]) > 1e-10:
                        valid_ratios.append((i, j, B1[i, j] / B0[i, j]))

            # Check rank-1
            is_r1 = True
            vr_set = {(i, j) for i, j, _ in valid_ratios}
            for i1, j1, r1 in valid_ratios:
                for i2, j2, r2 in valid_ratios:
                    if (i1, j2) in vr_set and (i2, j1) in vr_set:
                        r12 = [r for i, j, r in valid_ratios if i == i1 and j == j2][0]
                        r21 = [r for i, j, r in valid_ratios if i == i2 and j == j1][0]
                        if abs(r1 * r2 - r12 * r21) > 1e-6:
                            is_r1 = False
                            break
                if not is_r1:
                    break

            if is_r1:
                print(f"    Ratio matrix IS rank-1 => diagonal gauge exists")
            else:
                print(f"    Ratio matrix is NOT rank-1 => no diagonal gauge")

    # ================================================================
    # PART 8: Check n=5 too for comparison
    # ================================================================
    print("\n" + "=" * 80)
    print("PART 8: Same analysis for n=5 (Interval)")
    print("=" * 80)

    n5 = 5
    S5 = {1, 2}
    A5 = circulant_tournament(n5, S5)
    omega5 = np.exp(2j * np.pi / n5)

    for m in range(1, n5):
        B0, pm, pm1 = build_eigenspace_boundary(n5, S5, 0, m, A5)
        B1, _, _ = build_eigenspace_boundary(n5, S5, 1, m, A5)

        if B0 is None:
            continue

        nrows, ncols = B0.shape
        print(f"\n--- n=5, Degree m={m}, size {nrows}x{ncols} ---")

        # Show matrices
        print(f"  B_0:")
        for i in range(nrows):
            row = [f"{B0[i,j].real:6.2f}" for j in range(ncols)]
            print(f"    [{', '.join(row)}]")

        print(f"  B_1:")
        for i in range(nrows):
            row = [f"{B1[i,j]:.4f}" for j in range(ncols)]
            print(f"    [{', '.join(row)}]")

        # SVDs
        U0, s0, V0h = svd(B0, full_matrices=False)
        U1, s1, V1h = svd(B1, full_matrices=False)
        print(f"  SVD(B_0): {[f'{v:.6f}' for v in s0]}")
        print(f"  SVD(B_1): {[f'{v:.6f}' for v in s1]}")
        print(f"  Same SVDs? {np.allclose(sorted(s0, reverse=True), sorted(s1, reverse=True), atol=1e-6)}")

        # Diagonal gauge check
        paths_m = [p for p in get_regular_paths(A5, m) if p[0] == 0]
        paths_m1 = [p for p in get_regular_paths(A5, m - 1) if p[0] == 0]

        for a in range(n5):
            for b in range(n5):
                d_col = np.array([omega5 ** (a * paths_m[j][-1]) for j in range(ncols)])
                d_row = np.array([omega5 ** (b * paths_m1[i][-1]) for i in range(nrows)])
                B1_test = np.diag(d_row) @ B0 @ np.diag(d_col)
                err = np.max(np.abs(B1_test - B1))
                if err < 1e-6:
                    print(f"  GAUGE MATCH: row=omega^({b}*v_last), col=omega^({a}*v_last), err={err:.2e}")

    # ================================================================
    # PART 9: Check the difference-sequence representation
    # ================================================================
    print("\n" + "=" * 80)
    print("PART 9: Boundary in difference-sequence representation")
    print("=" * 80)
    print("\nPaths from 0 are determined by difference sequence (d_1,...,d_m), d_i in S.")
    print("Vertex sequence: (0, d_1, d_1+d_2, ..., d_1+...+d_m) mod n.")
    print("Removing interior vertex i (1<=i<=m-1) merges d_i and d_{i+1}.")
    print("The resulting face starts at vertex 0 if i>0, or at vertex d_1 if i=0.")
    print("Wait — we also remove i=0 (first vertex) and i=m (last vertex)?")
    print("NO: for GLMY path homology, boundary only removes INTERIOR vertices.")
    print("Interior = positions 1,...,m-1 (for an m-path = (m+1)-tuple).")
    print("")
    print("So: removing vertex at position i (1<=i<=m-1):")
    print("  - The face always starts at vertex 0 (shift=0)!")
    print("  - The diff sequence becomes (..., d_i + d_{i+1}, ...)")
    print("  - The face is STILL from vertex 0, so canonical = face itself!")
    print("")
    print("WAIT. This means shift is ALWAYS 0 for GLMY boundary.")
    print("If shift=0, then phase = omega^{k*0} = 1 for ALL k.")
    print("So B_k = B_0 for ALL k?!")
    print("")
    print("Let me verify this claim...")

    n = 7
    S = {1, 2, 3}
    A = circulant_tournament(n, S)

    for m in range(1, n):
        all_paths_m = get_regular_paths(A, m)
        paths_m = [p for p in all_paths_m if p[0] == 0]

        if not paths_m:
            continue

        path_len = len(paths_m[0]) - 1

        print(f"\n  Degree m={m}:")
        shifts_seen = set()
        for path in paths_m:
            for i in range(1, path_len):
                face = path[:i] + path[i+1:]
                shift = face[0]
                shifts_seen.add(shift)
                if shift != 0:
                    print(f"    PATH {path}, remove pos {i}: face={face}, shift={shift} != 0!")

        if shifts_seen == {0}:
            print(f"    All shifts are 0 => B_k = B_0 for all k at this degree")
        else:
            print(f"    Shifts seen: {sorted(shifts_seen)}")

    # ================================================================
    # PART 10: Direct comparison B_k matrices
    # ================================================================
    print("\n" + "=" * 80)
    print("PART 10: Direct equality check B_k == B_0")
    print("=" * 80)

    for m in range(1, n):
        B0 = all_Bk.get((0, m))
        if B0 is None:
            continue

        all_equal = True
        for k in range(1, n):
            Bk = all_Bk.get((k, m))
            if Bk is None:
                all_equal = False
                break
            if not np.allclose(B0, Bk, atol=1e-10):
                all_equal = False
                diff = np.max(np.abs(B0 - Bk))
                print(f"  m={m}, k={k}: NOT equal, max diff = {diff:.2e}")

        if all_equal:
            print(f"  m={m}: B_k == B_0 for ALL k=0,...,{n-1}")

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("""
The key insight about GLMY path homology boundary:

For an m-path (v_0, v_1, ..., v_m), the GLMY boundary removes INTERIOR vertices
only (positions 1 through m-1). The resulting face always starts at v_0 and ends
at v_m. When we use orbit representatives starting at vertex 0, the face ALSO
starts at vertex 0, so the shift needed to canonicalize is ALWAYS 0.

Therefore: phase = omega^{k * 0} = 1 for ALL k.

Consequence: B_k = B_0 for ALL k = 0, ..., p-1.

This means all eigenspace boundary matrices are IDENTICAL (not just isomorphic),
and hence all ranks, kernels, and Betti numbers are identical.

The uniformity is NOT a deep algebraic coincidence — it's a direct consequence of
the GLMY boundary only removing interior vertices (not endpoints), combined with
the convention of using orbit reps starting at a fixed vertex.

Algebraic proof:
  Let alpha = (0, a_1, ..., a_m) be an m-path orbit representative (starting at 0).
  The i-th face (1 <= i <= m-1) is:
    face_i = (0, a_1, ..., a_{i-1}, a_{i+1}, ..., a_m)
  which ALSO starts at 0 (no rotation needed to canonicalize).

  Therefore B_k[face_i_index, alpha_index] += (-1)^i * omega^{k*0} = (-1)^i

  This is INDEPENDENT of k.  QED.
""")


if __name__ == "__main__":
    analyze_eigenspace_uniformity()
