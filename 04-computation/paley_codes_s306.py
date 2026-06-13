#!/usr/bin/env python3
"""
paley_codes_s306.py — Paley tournaments as error-correcting codes

Claim: The adjacency matrix of the Paley tournament P_p, used as a
parity-check matrix over GF(2), yields classical codes:
  P_7  -> Hamming [7, 4, 3]
  P_23 -> binary Golay [23, 12, 7]

We test multiple constructions for p = 7, 11, 13, 23:
  (1) Raw adjacency matrix A as parity-check matrix
  (2) A + I  (include diagonal)
  (3) Bordered QR matrix (p+1 x p+1)
  (4) First (p-1)/2 rows of A
  (5) QR circulant Q (same as A for tournaments)
  (6) QR+I circulant

The key insight: for the standard QR code construction, one uses the
QR idempotent as a generator polynomial in GF(2)[x]/(x^p - 1).
The null space of the adjacency matrix gives a code, but it may be
the DUAL of the famous code.

opus-2026-03-24-S306
"""

import numpy as np
import time


def quadratic_residues(p):
    """Return set of quadratic residues mod p (excluding 0)."""
    qr = set()
    for a in range(1, p):
        qr.add((a * a) % p)
    return qr


def paley_adjacency(p):
    """
    Construct the Paley tournament on p vertices (p prime, p = 3 mod 4).
    Arc i->j iff (j-i) mod p is a quadratic residue.
    Returns adjacency matrix A (p x p) over integers.
    """
    qr = quadratic_residues(p)
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i, j] = 1
    return A


def gf2_rank(M):
    """Compute rank of matrix M over GF(2) using Gaussian elimination."""
    m, n = M.shape
    mat = M.copy() % 2
    rank = 0
    for col in range(n):
        pivot = -1
        for row in range(rank, m):
            if mat[row, col] == 1:
                pivot = row
                break
        if pivot == -1:
            continue
        mat[[rank, pivot]] = mat[[pivot, rank]]
        for row in range(m):
            if row != rank and mat[row, col] == 1:
                mat[row] = (mat[row] + mat[rank]) % 2
        rank += 1
    return rank


def gf2_nullspace(M):
    """
    Compute a basis for the null space of M over GF(2).
    Returns list of vectors (as numpy arrays) such that M @ v = 0 mod 2.
    """
    m, n = M.shape
    mat = M.copy() % 2
    pivot_cols = []
    row = 0
    for col in range(n):
        pivot = -1
        for r in range(row, m):
            if mat[r, col] == 1:
                pivot = r
                break
        if pivot == -1:
            continue
        mat[[row, pivot]] = mat[[pivot, row]]
        for r in range(m):
            if r != row and mat[r, col] == 1:
                mat[r] = (mat[r] + mat[row]) % 2
        pivot_cols.append(col)
        row += 1

    free_cols = [c for c in range(n) if c not in pivot_cols]

    basis = []
    for fc in free_cols:
        v = np.zeros(n, dtype=int)
        v[fc] = 1
        for i, pc in enumerate(pivot_cols):
            v[pc] = mat[i, fc]
        basis.append(v)
    return basis


def hamming_weight(v):
    return int(np.sum(v % 2))


def minimum_distance_bruteforce(null_basis, n, max_codewords=50_000_000):
    """Exact minimum distance via brute-force enumeration."""
    k = len(null_basis)
    if k == 0:
        return None
    total = (1 << k) - 1
    if total > max_codewords:
        return None

    min_d = n + 1
    for mask in range(1, (1 << k)):
        cw = np.zeros(n, dtype=int)
        for i in range(k):
            if mask & (1 << i):
                cw = (cw + null_basis[i]) % 2
        w = hamming_weight(cw)
        if w < min_d:
            min_d = w
    return min_d


def weight_distribution(null_basis, n, max_codewords=50_000_000):
    """Full weight distribution if feasible."""
    k = len(null_basis)
    total = 1 << k
    if total > max_codewords:
        return None
    dist = {}
    for mask in range(total):
        cw = np.zeros(n, dtype=int)
        for i in range(k):
            if mask & (1 << i):
                cw = (cw + null_basis[i]) % 2
        w = hamming_weight(cw)
        dist[w] = dist.get(w, 0) + 1
    return dist


def analyze_code(H_matrix, code_length, label, verbose=True):
    """
    Analyze a code with parity-check matrix H over GF(2).
    Returns (k, d, weight_dist_or_None).
    """
    H = H_matrix % 2
    r = gf2_rank(H)
    k = code_length - r

    basis = gf2_nullspace(H)
    assert len(basis) == k, f"Null space dim {len(basis)} != k={k}"
    for v in basis:
        assert np.all((H @ v) % 2 == 0), "Null space verification failed!"

    if k == 0:
        if verbose:
            print(f"    {label}: [{code_length}, 0, -] (trivial, rank={r})")
        return 0, None, {0: 1}

    t0 = time.time()
    d = minimum_distance_bruteforce(basis, code_length)
    elapsed = time.time() - t0

    dist = weight_distribution(basis, code_length)

    if d is not None:
        d_str = str(d)
        method = f"exact, {elapsed:.2f}s"
    else:
        d_str = "?"
        method = "too large for brute force"

    if verbose:
        print(f"    {label}: [{code_length}, {k}, {d_str}]  (rank={r}, {method})")
        if dist is not None and k <= 12:
            weights = sorted(dist.keys())
            nonzero = [(w, dist[w]) for w in weights if w > 0]
            if len(nonzero) <= 10:
                wt_str = ", ".join(f"w{w}:{c}" for w, c in nonzero)
            else:
                wt_str = ", ".join(f"w{w}:{c}" for w, c in nonzero[:5])
                wt_str += f", ... ({len(nonzero)} nonzero weights)"
            print(f"      Weight dist: {wt_str}")

    return k, d, dist


def bordered_qr_matrix(p):
    """
    Bordered QR matrix (p+1)x(p+1):
    H[0,0] = 0, H[0,j] = 1 for j>0, H[i,0] = 1 for i>0
    H[i,j] = Q[i-1,j-1] for i,j > 0
    where Q is the QR circulant.
    """
    qr = quadratic_residues(p)
    n = p + 1
    H = np.zeros((n, n), dtype=int)
    H[0, 1:] = 1
    H[1:, 0] = 1
    for i in range(p):
        for j in range(p):
            if ((j - i) % p) in qr:
                H[i + 1, j + 1] = 1
    return H


def main():
    print("=" * 70)
    print("  Paley Tournaments -> Error-Correcting Codes")
    print("  Testing p = 7, 11, 13, 23")
    print("=" * 70)

    all_results = {}

    for p in [7, 11, 13, 23]:
        print(f"\n{'='*70}")
        print(f"  p = {p}")
        print(f"{'='*70}")

        qr = quadratic_residues(p)
        print(f"  QR mod {p}: {sorted(qr)}")
        print(f"  p mod 4 = {p % 4}", end="")
        if p % 4 == 3:
            print("  =>  Paley TOURNAMENT (antisymmetric)")
        else:
            print("  =>  NOT a tournament (p=1 mod 4, -1 is QR, graph is undirected)")

        A = paley_adjacency(p)

        # Check tournament
        is_tourn = all(A[i, j] + A[j, i] == 1
                       for i in range(p) for j in range(i + 1, p))
        scores = sorted(A.sum(axis=1))
        print(f"  Tournament: {is_tourn}, scores: {list(scores)}")

        # Self-orthogonality of A over GF(2)
        prod = (A @ A.T) % 2
        aat_zero = np.all(prod == 0)
        print(f"  A*A^T = 0 mod 2: {aat_zero}")

        # Integer product A*A^T (to understand the structure)
        prod_int = A @ A.T
        diag = prod_int[0, 0]
        offdiag = prod_int[0, 1]
        print(f"  A*A^T over Z: diagonal={diag}, off-diagonal={offdiag}")
        # For a regular tournament with out-degree (p-1)/2:
        # (A A^T)_{ii} = sum_k A_{ik}^2 = (p-1)/2 (each row has (p-1)/2 ones)
        # (A A^T)_{ij} = sum_k A_{ik}*A_{jk} = # vertices k with i->k and j->k
        # For Paley: this equals (p-1)/4 - 1/4 if i->j, or (p-1)/4 + 1/4 if j->i
        # Actually for Paley: A*A^T = ((p-1)/4)*J + ((p+1)/4)*I - (1/2)*(A-A^T) mod something
        # The key: (A+A^T) is the Paley GRAPH adjacency (for p=3 mod 4, A+A^T = J-I)

        print()

        results = {}

        # Construction 1: Raw A
        k, d, dist = analyze_code(A, p, "C1: H=A (raw adjacency)")
        results['raw'] = (p, k, d)

        # Construction 2: A + I
        AI = (A + np.eye(p, dtype=int)) % 2
        k, d, dist = analyze_code(AI, p, "C2: H=A+I")
        results['A+I'] = (p, k, d)

        # Construction 3: Bordered QR matrix
        HB = bordered_qr_matrix(p)
        k, d, dist = analyze_code(HB, p + 1, f"C3: Bordered QR ({p+1}x{p+1})")
        results['bordered'] = (p + 1, k, d)

        # Construction 4: First (p-1)/2 rows of A
        half = (p - 1) // 2
        A_half = A[:half, :] % 2
        k, d, dist = analyze_code(A_half, p, f"C4: First {half} rows of A")
        results['half'] = (p, k, d)

        # Construction 5: Transpose of A (= complement tournament)
        AT = A.T % 2
        k, d, dist = analyze_code(AT, p, "C5: H=A^T (complement)")
        results['A^T'] = (p, k, d)

        # Construction 6: A + I bordered
        # (p+1)x(p+1) with A+I in the lower-right
        n2 = p + 1
        HBI = np.zeros((n2, n2), dtype=int)
        HBI[0, :] = 1  # top row all ones
        HBI[:, 0] = 1  # left col all ones
        HBI[0, 0] = 1  # corner
        HBI[1:, 1:] = (A + np.eye(p, dtype=int)) % 2
        k, d, dist = analyze_code(HBI % 2, p + 1, f"C6: Bordered (A+I) ({p+1}x{p+1})")
        results['bordered_AI'] = (p + 1, k, d)

        # Construction 7: Systematic form [I | A_half^T]
        # Use the rowspace of A as generator matrix instead
        print(f"\n    --- Dual code (use A as GENERATOR matrix) ---")
        # If A is generator, code = rowspace(A), parity check = left null space
        # Code C = rowspace of A over GF(2) has dimension = rank(A)
        r = gf2_rank(A % 2)
        print(f"    C7: G=A (generator): [{p}, {r}, ?]", end="")
        # For generator matrix, codewords = all linear combos of rows
        # This is equivalent to null space of null space... let's just enumerate
        # rows of A form a generating set
        if r <= 20:
            # Build the code from rows
            rows_gf2 = A % 2
            # Get a basis via row reduction
            mat = rows_gf2.copy()
            basis_rows = []
            used = 0
            for col in range(p):
                piv = -1
                for row in range(used, p):
                    if mat[row, col] == 1:
                        piv = row
                        break
                if piv == -1:
                    continue
                mat[[used, piv]] = mat[[piv, used]]
                for row in range(p):
                    if row != used and mat[row, col] == 1:
                        mat[row] = (mat[row] + mat[used]) % 2
                basis_rows.append(mat[used].copy())
                used += 1

            # Minimum distance of row space
            min_d = p + 1
            total = (1 << len(basis_rows)) - 1
            if total <= 50_000_000:
                for mask in range(1, total + 1):
                    cw = np.zeros(p, dtype=int)
                    for i in range(len(basis_rows)):
                        if mask & (1 << i):
                            cw = (cw + basis_rows[i]) % 2
                    w = hamming_weight(cw)
                    if w < min_d:
                        min_d = w
                print(f"  d={min_d}")
            else:
                print(f"  (too large)")
                min_d = None
            results['generator'] = (p, r, min_d)

        all_results[p] = results

    # =========================================================
    # GRAND SUMMARY
    # =========================================================
    print(f"\n\n{'='*70}")
    print(f"  GRAND SUMMARY")
    print(f"{'='*70}")

    target_codes = {
        7: "[7,4,3] Hamming",
        23: "[23,12,7] Golay",
    }

    for p in [7, 11, 13, 23]:
        print(f"\n  p = {p}:")
        results = all_results[p]
        for key, (n, k, d) in results.items():
            d_str = str(d) if d is not None else "?"
            code_str = f"[{n},{k},{d_str}]"
            tag = ""
            if p == 7 and n == 7 and k == 4 and d == 3:
                tag = " <-- HAMMING [7,4,3] !"
            if p == 7 and n == 7 and k == 3 and d == 4:
                tag = " <-- dual of Hamming = simplex code"
            if p == 23 and n == 23 and k == 12 and d == 7:
                tag = " <-- GOLAY [23,12,7] !"
            if p == 23 and n == 23 and k == 11 and d == 8:
                tag = " <-- dual of Golay (expurgated)"
            if p == 23 and n == 24 and k == 12 and d == 8:
                tag = " <-- EXTENDED GOLAY [24,12,8] !"
            print(f"    {key:20s}  {code_str:15s}{tag}")

    print(f"\n{'='*70}")
    print(f"  VERDICT ON ORIGINAL CLAIMS")
    print(f"{'='*70}")

    # P_7 analysis
    print(f"\n  P_7 -> Hamming [7,4,3]?")
    r7 = all_results[7]
    raw_n, raw_k, raw_d = r7['raw']
    print(f"    Raw A as parity-check: [{raw_n},{raw_k},{raw_d}]")
    if raw_k == 3 and raw_d == 4:
        print(f"    This is the [7,3,4] SIMPLEX code = dual of Hamming [7,4,3].")
        print(f"    The claim is HALF-RIGHT: A defines the Hamming code, but as")
        print(f"    GENERATOR matrix (rowspace), not parity-check (nullspace).")
    gen_n, gen_k, gen_d = r7.get('generator', (None, None, None))
    if gen_k == 4 and gen_d == 3:
        print(f"    A as GENERATOR matrix: [{gen_n},{gen_k},{gen_d}] = Hamming [7,4,3]. CONFIRMED.")
    elif gen_k is not None:
        print(f"    A as generator matrix: [{gen_n},{gen_k},{gen_d}]")

    # P_23 analysis
    print(f"\n  P_23 -> Golay [23,12,7]?")
    r23 = all_results[23]
    raw_n, raw_k, raw_d = r23['raw']
    print(f"    Raw A as parity-check: [{raw_n},{raw_k},{raw_d}]")
    gen_n, gen_k, gen_d = r23.get('generator', (None, None, None))
    if gen_k is not None:
        print(f"    A as generator matrix: [{gen_n},{gen_k},{gen_d}]")
    bord = r23.get('bordered', (None, None, None))
    if bord[0] is not None:
        print(f"    Bordered QR matrix:    [{bord[0]},{bord[1]},{bord[2]}]")
    bord_ai = r23.get('bordered_AI', (None, None, None))
    if bord_ai[0] is not None:
        print(f"    Bordered (A+I):        [{bord_ai[0]},{bord_ai[1]},{bord_ai[2]}]")

    # Check which construction actually gives Golay
    found_golay = False
    for key, (n, k, d) in r23.items():
        if n == 23 and k == 12 and d == 7:
            print(f"    GOLAY [23,12,7] found via: {key}")
            found_golay = True
        if n == 24 and k == 12 and d == 8:
            print(f"    EXTENDED GOLAY [24,12,8] found via: {key}")
            found_golay = True
    if not found_golay:
        if raw_k == 11 and raw_d == 8:
            print(f"    Null space of A gives [{raw_n},{raw_k},{raw_d}]:")
            print(f"    This is the EXPURGATED Golay code (dual perspective).")
            print(f"    The [23,12,7] Golay is the ROWSPACE of A (dual code).")

    # P_11 analysis
    print(f"\n  P_11 analysis:")
    r11 = all_results[11]
    for key, (n, k, d) in r11.items():
        d_str = str(d) if d is not None else "?"
        if k > 0:
            print(f"    {key:20s}  [{n},{k},{d_str}]")

    # P_13 analysis
    print(f"\n  P_13 analysis (p=1 mod 4, NOT a tournament):")
    r13 = all_results[13]
    for key, (n, k, d) in r13.items():
        d_str = str(d) if d is not None else "?"
        if k > 0:
            print(f"    {key:20s}  [{n},{k},{d_str}]")

    print()


if __name__ == "__main__":
    main()
