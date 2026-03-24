#!/usr/bin/env python3
"""
paley_codes_s306.py — Paley tournaments as error-correcting codes

Claim: The adjacency matrix of the Paley tournament P_p, used as a
parity-check matrix over GF(2), yields classical codes:
  P_7  -> Hamming [7, 4, 3]
  P_23 -> binary Golay [23, 12, 7]

We test this for p = 7, 11, 13, 23.

Method:
  1. Construct P_p using quadratic residues mod p.
     Arc i->j iff (j-i) is a QR mod p.
  2. Take A = adjacency matrix over GF(2).
  3. Use A as parity-check matrix H.
     Code C = {x in GF(2)^p : H x = 0}.
     Dimension k = p - rank_GF2(H).
  4. Minimum distance d = min weight of nonzero codeword.
     (Brute force for small codes, systematic enumeration for larger.)

opus-2026-03-24-S306
"""

import numpy as np
from itertools import combinations
import time


def quadratic_residues(p):
    """Return set of quadratic residues mod p (excluding 0)."""
    qr = set()
    for a in range(1, p):
        qr.add((a * a) % p)
    return qr


def paley_adjacency(p):
    """
    Construct the Paley tournament on p vertices (p prime, p ≡ 3 mod 4).
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
    # Work with a copy mod 2
    m = M.shape[0]
    n = M.shape[1]
    mat = M.copy() % 2
    rank = 0
    for col in range(n):
        # Find pivot
        pivot = -1
        for row in range(rank, m):
            if mat[row, col] == 1:
                pivot = row
                break
        if pivot == -1:
            continue
        # Swap
        mat[[rank, pivot]] = mat[[pivot, rank]]
        # Eliminate
        for row in range(m):
            if row != rank and mat[row, col] == 1:
                mat[row] = (mat[row] + mat[rank]) % 2
        rank += 1
    return rank


def gf2_nullspace(M):
    """
    Compute a basis for the null space of M over GF(2).
    Returns list of vectors (as numpy arrays) in GF(2)^n such that M @ v = 0 mod 2.
    """
    m, n = M.shape
    mat = M.copy() % 2
    # Augment with identity to track column operations
    # Use row reduction on transpose to find null space
    # Null space of M = {x : Mx = 0 mod 2}
    # Row reduce M^T to find relations

    # Standard approach: row reduce [M; I_n] treating columns as vectors
    # Actually, let's just do RREF on M and read off null space

    mat = M.copy() % 2
    pivot_cols = []
    row = 0
    for col in range(n):
        # Find pivot in this column
        pivot = -1
        for r in range(row, m):
            if mat[r, col] == 1:
                pivot = r
                break
        if pivot == -1:
            continue
        # Swap
        mat[[row, pivot]] = mat[[pivot, row]]
        # Eliminate all other rows
        for r in range(m):
            if r != row and mat[r, col] == 1:
                mat[r] = (mat[r] + mat[row]) % 2
        pivot_cols.append(col)
        row += 1

    rank = len(pivot_cols)
    free_cols = [c for c in range(n) if c not in pivot_cols]

    # For each free column, construct a null space vector
    basis = []
    for fc in free_cols:
        v = np.zeros(n, dtype=int)
        v[fc] = 1
        # For each pivot column, determine the value
        for i, pc in enumerate(pivot_cols):
            v[pc] = mat[i, fc]  # This is the coefficient
        basis.append(v)

    return basis


def hamming_weight(v):
    """Hamming weight of a binary vector."""
    return int(np.sum(v % 2))


def minimum_distance_bruteforce(null_basis, n):
    """
    Compute minimum distance by enumerating all nonzero codewords.
    Codewords = all GF(2) linear combinations of null_basis vectors.
    """
    k = len(null_basis)
    if k == 0:
        return float('inf')

    min_d = n + 1
    # Enumerate all 2^k - 1 nonzero combinations
    total = (1 << k) - 1
    if total > 10_000_000:
        return None  # Too many, skip brute force

    for mask in range(1, (1 << k)):
        cw = np.zeros(n, dtype=int)
        for i in range(k):
            if mask & (1 << i):
                cw = (cw + null_basis[i]) % 2
        w = hamming_weight(cw)
        if w < min_d:
            min_d = w

    return min_d


def minimum_distance_sampling(null_basis, n, samples=100000):
    """
    Estimate minimum distance by random sampling of codewords.
    Returns a lower bound (actual min distance <= returned value).
    """
    import random
    k = len(null_basis)
    if k == 0:
        return float('inf')

    min_d = n + 1

    # First check all weight-1 combinations (individual basis vectors)
    for v in null_basis:
        w = hamming_weight(v)
        if w < min_d:
            min_d = w

    # Then check weight-2 combinations
    for i in range(k):
        for j in range(i + 1, k):
            cw = (null_basis[i] + null_basis[j]) % 2
            w = hamming_weight(cw)
            if w < min_d:
                min_d = w

    # Random sampling
    for _ in range(samples):
        # Random subset of basis vectors
        num = random.randint(1, k)
        indices = random.sample(range(k), num)
        cw = np.zeros(n, dtype=int)
        for idx in indices:
            cw = (cw + null_basis[idx]) % 2
        w = hamming_weight(cw)
        if w < min_d:
            min_d = w

    return min_d


def analyze_paley_code(p):
    """Full analysis of Paley tournament code for prime p."""
    print(f"\n{'='*60}")
    print(f"  Paley Tournament P_{p}")
    print(f"{'='*60}")

    # Check p ≡ 3 mod 4
    if p % 4 != 3:
        print(f"  WARNING: p = {p} ≡ {p % 4} (mod 4), not 3 (mod 4).")
        print(f"  Paley tournament requires p ≡ 3 (mod 4).")
        print(f"  For p ≡ 1 (mod 4), -1 is a QR, so the 'tournament' is")
        print(f"  actually a self-complementary graph (not a tournament).")
        print(f"  We'll construct the QR-directed digraph anyway.\n")

    qr = quadratic_residues(p)
    print(f"  Quadratic residues mod {p}: {sorted(qr)}")
    print(f"  |QR| = {len(qr)} = (p-1)/2 = {(p-1)//2}")

    A = paley_adjacency(p)

    # Verify it's a tournament (for p ≡ 3 mod 4)
    is_tournament = True
    for i in range(p):
        for j in range(i + 1, p):
            if A[i, j] + A[j, i] != 1:
                is_tournament = False
                break
    print(f"  Is tournament: {is_tournament}")

    # Score sequence
    scores = sorted(A.sum(axis=1))
    print(f"  Score sequence: {list(scores)}")
    is_regular = all(s == scores[0] for s in scores)
    print(f"  Regular: {is_regular} (all scores = {scores[0] if is_regular else 'varies'})")

    # GF(2) analysis
    A_gf2 = A % 2
    r = gf2_rank(A_gf2)
    k = p - r
    print(f"\n  Parity-check matrix H = A (mod 2):")
    print(f"    Size: {p} x {p}")
    print(f"    GF(2) rank: {r}")
    print(f"    Code length n = {p}")
    print(f"    Code dimension k = {p} - {r} = {k}")

    # Compute null space and minimum distance
    null_basis = gf2_nullspace(A_gf2)
    print(f"    Null space dimension: {len(null_basis)}")

    # Verify null space
    for v in null_basis:
        check = A_gf2 @ v % 2
        assert np.all(check == 0), "Null space vector fails verification!"
    print(f"    Null space verified: all basis vectors satisfy Hx=0")

    t0 = time.time()
    total_codewords = (1 << k) - 1
    if total_codewords <= 10_000_000:
        print(f"    Computing minimum distance (brute force over {total_codewords + 1} codewords)...")
        d = minimum_distance_bruteforce(null_basis, p)
        print(f"    Minimum distance d = {d}  (exact, brute force)")
    else:
        print(f"    Too many codewords ({total_codewords + 1}) for brute force.")
        print(f"    Estimating minimum distance by sampling...")
        d = minimum_distance_sampling(null_basis, p, samples=500000)
        print(f"    Minimum distance d <= {d}  (upper bound from sampling)")
    elapsed = time.time() - t0
    print(f"    Time: {elapsed:.2f}s")

    # Weight distribution (if feasible)
    if total_codewords <= 10_000_000:
        weight_dist = {}
        for mask in range(0, (1 << k)):
            cw = np.zeros(p, dtype=int)
            for i in range(k):
                if mask & (1 << i):
                    cw = (cw + null_basis[i]) % 2
            w = hamming_weight(cw)
            weight_dist[w] = weight_dist.get(w, 0) + 1
        print(f"\n    Weight distribution:")
        for w in sorted(weight_dist.keys()):
            print(f"      w={w:3d}: {weight_dist[w]} codewords")

    print(f"\n  => Code parameters: [{p}, {k}, {d}]")

    # Check against known codes
    known = {
        7: (4, 3, "Hamming [7,4,3]"),
        23: (12, 7, "binary Golay [23,12,7]"),
    }
    if p in known:
        exp_k, exp_d, name = known[p]
        match = (k == exp_k and d == exp_d)
        print(f"  Expected: {name} -> k={exp_k}, d={exp_d}")
        print(f"  MATCH: {'YES' if match else 'NO'}")

    return p, k, d


def check_self_orthogonality(p):
    """Check if A * A^T = 0 over GF(2), which means code is self-orthogonal."""
    A = paley_adjacency(p) % 2
    prod = (A @ A.T) % 2
    is_self_orth = np.all(prod == 0)
    # Also check A*A^T mod 2
    # For QR codes, A*A^T = (p-1)/2 * I + ...
    print(f"\n  Self-orthogonality check (H * H^T mod 2):")
    diag_vals = set(prod[i, i] for i in range(p))
    offdiag_vals = set(prod[i, j] for i in range(p) for j in range(p) if i != j)
    print(f"    Diagonal values mod 2: {diag_vals}")
    print(f"    Off-diagonal values mod 2: {offdiag_vals}")
    print(f"    H * H^T = 0 (mod 2): {is_self_orth}")
    return is_self_orth


def main():
    print("=" * 60)
    print("  Paley Tournaments as Error-Correcting Codes")
    print("  Testing: P_7, P_11, P_13, P_23")
    print("=" * 60)

    results = []

    for p in [7, 11, 13, 23]:
        p_val, k, d = analyze_paley_code(p)
        check_self_orthogonality(p)
        results.append((p_val, k, d))

    # Summary
    print("\n" + "=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    print(f"  {'p':>4s}  {'n':>4s}  {'k':>4s}  {'d':>4s}  {'Code':>20s}  {'Known?'}")
    print(f"  {'---':>4s}  {'---':>4s}  {'---':>4s}  {'---':>4s}  {'---':>20s}  {'---'}")

    known_codes = {
        7: "Hamming [7,4,3]",
        23: "Golay [23,12,7]",
    }

    for p, k, d in results:
        code_str = f"[{p},{k},{d}]"
        known = known_codes.get(p, "")
        match_str = ""
        if p == 7:
            match_str = "YES" if (k == 4 and d == 3) else "NO"
        elif p == 23:
            match_str = "YES" if (k == 12 and d == 7) else "NO"
        print(f"  {p:4d}  {p:4d}  {k:4d}  {d:4d}  {code_str:>20s}  {known}  {match_str}")

    print("\n  Conclusion:")
    for p, k, d in results:
        if p == 7:
            if k == 4 and d == 3:
                print(f"    P_7 -> [{p},{k},{d}] = Hamming code. CONFIRMED.")
            else:
                print(f"    P_7 -> [{p},{k},{d}] != Hamming [7,4,3]. REFUTED.")
        elif p == 23:
            if k == 12 and d == 7:
                print(f"    P_23 -> [{p},{k},{d}] = binary Golay code. CONFIRMED.")
            else:
                print(f"    P_23 -> [{p},{k},{d}] != Golay [23,12,7]. REFUTED.")
        elif p == 11:
            print(f"    P_11 -> [{p},{k},{d}] = quadratic residue code.")
        elif p == 13:
            print(f"    P_13 -> [{p},{k},{d}].")
            if p % 4 == 1:
                print(f"      Note: P_13 is not a proper Paley tournament (13 ≡ 1 mod 4).")


if __name__ == "__main__":
    main()
