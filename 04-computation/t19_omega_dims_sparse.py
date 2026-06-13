"""
t19_omega_dims_sparse.py - Sparse computation of T_19 Omega dimensions.

ENGINEERING MOTIVATION:
  The dense constraint matrix at degree 6 is (166428 x 277236) x 4 bytes = 172 GB.
  But it has only ~7 nonzeros per column (degree m has m+1 faces).
  Sparse column reduction uses ~1-5 MB of storage instead.

  This script implements a general sparse modular rank computation using
  the column reduction algorithm (same as persistent homology software).

KNOWN RESULTS (from t19_omega_dims.py, which ran degrees 0-5):
  m=0: Omega=1
  m=1: Omega=9
  m=2: Omega=72   (|A_2|=81,   rank=9)
  m=3: Omega=540  (|A_3|=684,  rank=144)
  m=4: Omega=3753 (|A_4|=5436, rank=1683)
  m=5: Omega=23832 (|A_5|=40284, rank=16452)
  m=6: ??? (OOM with dense; this script solves it)

APPROACH:
  Sparse column reduction over F_191.
  Each degree-m column has at most m+1 nonzeros.
  Total nnz at degree 6: ~7 x 277236 = ~2M (vs 172 GB dense).

  Algorithm: process columns left to right. Maintain a dict pivot_of[row] -> column.
  For each column v, reduce against existing pivots. If v != 0 after reduction,
  it gives a new pivot. Rank = number of pivots found.

  Expected performance: fast if pivot columns stay sparse (low fill-in).
  Fill-in monitor: tracks max pivot density to detect performance issues.

Author: kind-pasteur-2026-03-10-S54
"""
import sys
import time
import os
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

PRIME = 191   # 191-1=190=10*19, so 19 | 191-1
P_PALEY = 19

KNOWN_OMEGA = {0: 1, 1: 9, 2: 72, 3: 540, 4: 3753, 5: 23832}
RESULTS_PATH = "05-knowledge/results/t19_omega_dims_sparse.out"


# ---------------------------------------------------------------------------
# Diff-seq enumeration (same as t19_omega_dims.py)
# ---------------------------------------------------------------------------

def qr_set(p):
    return set((a * a) % p for a in range(1, p))


def enumerate_diff_seqs(S, n, max_deg):
    """Enumerate all allowed diff-seqs up to length max_deg."""
    S_sorted = sorted(S)
    seqs = {0: [()]}
    partial_sums = {(): frozenset([0])}
    partial_last = {(): 0}
    for m in range(1, max_deg + 1):
        prev = seqs[m - 1]
        new, nps, nl = [], {}, {}
        for seq in prev:
            ps = partial_sums[seq]
            last = partial_last[seq]
            for s in S_sorted:
                nls = (last + s) % n
                if nls not in ps:
                    ns = seq + (s,)
                    new.append(ns)
                    nps[ns] = ps | frozenset([nls])
                    nl[ns] = nls
        seqs[m] = new
        partial_sums.update(nps)
        partial_last.update(nl)
    return seqs


def compute_face(D, face_idx, n):
    """Compute the face-idx face of diff-seq D (returns (m-1)-tuple)."""
    m = len(D)
    if face_idx == 0:
        return D[1:]
    elif face_idx == m:
        return D[:m - 1]
    else:
        merged = (D[face_idx - 1] + D[face_idx]) % n
        return D[:face_idx - 1] + (merged,) + D[face_idx + 1:]


# ---------------------------------------------------------------------------
# Sparse constraint matrix builder
# ---------------------------------------------------------------------------

def build_constraint_sparse(A_m_list, allowed_lower_set, n, prime):
    """
    Build constraint matrix as list of sparse column dicts over F_prime.

    For each D in A_m_list:
      column[D] = {row: sign mod prime} where row is a "junk" (m-1)-seq
      that appears as a face of D but is NOT in allowed_lower_set.

    Returns: (columns, n_junk, junk_idx)
      columns: list of dicts {row_int: val}
      n_junk: number of junk rows
      junk_idx: dict mapping junk_tuple -> row_int
    """
    junk_idx = {}   # junk_tuple -> row_index
    columns = []

    for D in A_m_list:
        col = {}  # row_idx -> accumulated value mod prime
        for fi in range(len(D) + 1):
            fd = compute_face(D, fi, n)
            if fd not in allowed_lower_set:
                if fd not in junk_idx:
                    junk_idx[fd] = len(junk_idx)
                row = junk_idx[fd]
                sign = 1 if fi % 2 == 0 else -1
                col[row] = (col.get(row, 0) + sign) % prime

        # Drop explicit zeros (can happen if two faces cancel mod prime)
        col = {r: v for r, v in col.items() if v != 0}
        columns.append(col)

    return columns, len(junk_idx), junk_idx


# ---------------------------------------------------------------------------
# Sparse column reduction rank computation over F_prime
# ---------------------------------------------------------------------------

def sparse_rank_mod(columns, prime, verbose=False, verbose_interval=50000):
    """
    Compute rank of a sparse matrix over F_prime using column reduction.

    columns: list of dicts {row_idx: value in [0, prime-1]}
    prime: a prime number

    Returns: (rank, max_pivot_density)

    Algorithm: process columns left to right. Reduce each column against
    existing pivot columns. If the reduced column is nonzero, it defines
    a new pivot (and rank increases by 1).

    Normalization: each pivot column is stored with its minimum-row entry = 1.

    Fill-in: each elimination step can add entries from the pivot column.
    For our boundary matrices (<=m+1 nonzeros per column), fill-in is bounded.
    """
    pivot_of = {}   # row_idx -> column dict (normalized: pivot_of[r][r] = 1)
    rank = 0
    n_cols = len(columns)
    t0 = time.time()
    max_pivot_density = 0
    total_elimination_steps = 0

    for j, col in enumerate(columns):
        if verbose and j % verbose_interval == 0 and j > 0:
            elapsed = time.time() - t0
            rate = j / elapsed if elapsed > 0 else 0
            eta = (n_cols - j) / rate if rate > 0 else float('inf')
            print(f"      col {j:7d}/{n_cols} ({100*j/n_cols:.1f}%), "
                  f"rank={rank:6d}, max_piv_den={max_pivot_density:4d}, "
                  f"{elapsed:6.0f}s elapsed, ETA={eta:.0f}s", flush=True)

        # Work on a mutable copy of the column
        v = dict(col)

        # Column reduction loop
        while v:
            r0 = min(v.keys())   # leftmost (smallest) nonzero row
            val0 = v[r0]

            if r0 not in pivot_of:
                # New pivot: normalize so v[r0] = 1
                inv = pow(int(val0), prime - 2, prime)
                v = {r: (val * inv) % prime for r, val in v.items() if (val * inv) % prime != 0}
                # r0 should still be in v (since val0 != 0 and inv != 0 mod prime)
                pivot_of[r0] = v
                rank += 1
                d = len(v)
                if d > max_pivot_density:
                    max_pivot_density = d
                break
            else:
                # Eliminate: v -= val0 * pivot_of[r0]  (since pivot_of[r0][r0] = 1)
                piv = pivot_of[r0]
                total_elimination_steps += 1
                for r, pv in piv.items():
                    new_val = (v.get(r, 0) - val0 * pv) % prime
                    if new_val == 0:
                        v.pop(r, None)
                    else:
                        v[r] = new_val
                # r0 is now eliminated (new_val for r0 was (val0 - val0*1) = 0)

        # If v is empty after reduction, column is in the span of pivots

    if verbose:
        elapsed = time.time() - t0
        print(f"    Rank={rank}, max_pivot_density={max_pivot_density}, "
              f"total_elim_steps={total_elimination_steps}, time={elapsed:.1f}s")

    return rank, max_pivot_density


# ---------------------------------------------------------------------------
# Dense rank (for verification on small cases)
# ---------------------------------------------------------------------------

def gauss_rank_dense_mod(C_dict_cols, n_rows, prime):
    """
    Dense Gaussian elimination for verification.
    Converts sparse column dicts to dense numpy array first.
    Only use for small matrices.
    """
    n_cols = len(C_dict_cols)
    if n_rows == 0 or n_cols == 0:
        return 0
    C = np.zeros((n_rows, n_cols), dtype=np.int32)
    for j, col in enumerate(C_dict_cols):
        for r, v in col.items():
            C[r, j] = v

    # Gaussian elimination mod prime
    rank = 0
    pivot_row = 0
    for col in range(n_cols):
        found = next((r for r in range(pivot_row, n_rows) if C[r, col] != 0), -1)
        if found < 0:
            continue
        C[[pivot_row, found]] = C[[found, pivot_row]]
        inv = pow(int(C[pivot_row, col]), prime - 2, prime)
        C[pivot_row] = (C[pivot_row].astype(np.int64) * inv % prime).astype(np.int32)
        fc = C[:, col].copy()
        fc[pivot_row] = 0
        for row in np.where(fc != 0)[0]:
            f = int(fc[row])
            C[row] = ((C[row].astype(np.int64) - f * C[pivot_row].astype(np.int64)) % prime
                      ).astype(np.int32)
        rank += 1
        pivot_row += 1
    return rank


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def compute_omega_dims_t19_sparse(max_deg, verify_dense_up_to=4):
    """
    Compute Omega dims for T_19 using sparse column reduction.

    verify_dense_up_to: for degrees <= this, also compute dense rank and compare.
    """
    p = P_PALEY
    prime = PRIME
    S = qr_set(p)

    print(f"T_19 Omega dimension computation (SPARSE METHOD)")
    print(f"  p={p}, prime={prime} ({prime-1} = {(prime-1)//p}*{p})")
    print(f"  QR_{p} = {sorted(S)}")
    print(f"  |S| = {len(S)} = (p-1)/2 = {(p-1)//2}")
    print()

    t0_total = time.time()
    print(f"Enumerating diff-seqs up to max_deg={max_deg}...")
    t_enum = time.time()
    diff_seqs = enumerate_diff_seqs(S, p, max_deg)
    print(f"  Done in {time.time()-t_enum:.1f}s")
    for m in range(max_deg + 1):
        print(f"  |A_{m}| = {len(diff_seqs.get(m, []))}")
    print()

    allowed_sets = {m: set(diff_seqs.get(m, [])) for m in range(max_deg + 2)}

    print(f"{'m':>3}  {'|A_m|':>8}  {'n_junk':>8}  {'nnz':>8}  "
          f"{'rank':>8}  {'Omega':>8}  {'time':>8}  {'max_piv_den':>12}  check")
    print("-" * 90)

    dims = []
    all_ok = True

    for m in range(max_deg + 1):
        A_m = diff_seqs.get(m, [])
        n_Am = len(A_m)

        if m == 0:
            dims.append(1)
            check = "OK" if KNOWN_OMEGA.get(0) == 1 else "FAIL"
            print(f"  0  {'1':>8}  {'0':>8}  {'0':>8}  {'0':>8}  {'1':>8}  "
                  f"{'0.0s':>8}  {'N/A':>12}  {check}")
            continue

        t1 = time.time()

        # Build sparse constraint matrix
        columns, n_junk, _ = build_constraint_sparse(
            A_m, allowed_sets[m - 1], p, prime
        )
        t_build = time.time() - t1

        total_nnz = sum(len(c) for c in columns)
        mem_cols_mb = total_nnz * 130 / 1e6   # rough bytes per dict entry

        # Estimate dense memory for comparison
        dense_gb = n_junk * n_Am * 4 / 1e9

        # Compute rank (sparse)
        verbose = (n_Am > 5000)
        t_rank = time.time()
        if n_junk == 0:
            rank = 0
            max_piv_den = 0
        else:
            rank, max_piv_den = sparse_rank_mod(
                columns, prime, verbose=verbose, verbose_interval=50000
            )
        t_rank_elapsed = time.time() - t_rank

        omega_dim = n_Am - rank
        dims.append(omega_dim)
        total_elapsed = time.time() - t1

        # Verification against known results
        known = KNOWN_OMEGA.get(m)
        if known is not None:
            check = "OK" if omega_dim == known else f"FAIL(exp {known})"
            if omega_dim != known:
                all_ok = False
        else:
            check = "NEW"

        # Dense verification for small degrees
        if m <= verify_dense_up_to and n_junk > 0:
            rank_dense = gauss_rank_dense_mod(columns, n_junk, prime)
            if rank_dense != rank:
                check += f" SPARSE/DENSE MISMATCH(dense={rank_dense})"
                all_ok = False

        print(f"{m:3d}  {n_Am:>8d}  {n_junk:>8d}  {total_nnz:>8d}  "
              f"{rank:>8d}  {omega_dim:>8d}  {total_elapsed:>7.1f}s  "
              f"{max_piv_den:>12d}  {check}")

        # Save intermediate result
        with open(RESULTS_PATH, 'a') as f:
            f.write(f"m={m:2d}: |A_m|={n_Am:7d}, n_junk={n_junk:7d}, nnz={total_nnz:7d}, "
                    f"rank={rank:7d}, Omega={omega_dim:7d}, "
                    f"max_piv_den={max_piv_den}, time={total_elapsed:.1f}s\n")

    print()
    print(f"Omega dims: {dims}")
    chi = sum((-1)**m * d for m, d in enumerate(dims))
    print(f"chi(T_19)  = {chi}  (per-eigenspace chi; total = 19*chi if all equal)")
    print(f"Known Paley chi: T_3=1, T_7=1, T_11=1 => expect chi=1 for T_19")
    print(f"Verification: {'ALL OK' if all_ok else 'SOME FAILURES'}")
    print(f"Total time: {time.time()-t0_total:.1f}s")

    return dims


def main():
    print("=" * 70)
    print("T_19 OMEGA DIMS -- SPARSE COLUMN REDUCTION")
    print("=" * 70)
    print()
    print("ENGINEERING NOTE:")
    print("  Dense method OOMed at degree 6: shape (166428, 277236), 172 GB int32")
    print("  Sparse method: ~7 nonzeros/col => ~2M total nonzeros => ~260 MB")
    print("  This is a 660x memory reduction.")
    print()
    print("KNOWN RESULTS (dense computation, kind-pasteur-S52):")
    print("  m=0: Omega=1")
    print("  m=1: Omega=9")
    print("  m=2: Omega=72  (|A_2|=81,   rank=9)")
    print("  m=3: Omega=540 (|A_3|=684,  rank=144)")
    print("  m=4: Omega=3753 (|A_4|=5436, rank=1683)")
    print("  m=5: Omega=23832 (|A_5|=40284, rank=16452)")
    print()

    # Initialize output file
    os.makedirs("05-knowledge/results", exist_ok=True)
    with open(RESULTS_PATH, 'w') as f:
        f.write("T_19 Omega dims (sparse computation, kind-pasteur-S54)\n")
        f.write("=" * 60 + "\n")
        f.write(f"PRIME={PRIME}, P={P_PALEY}\n")
        f.write("Known: m=0:1, m=1:9, m=2:72, m=3:540, m=4:3753, m=5:23832\n")
        f.write("=" * 60 + "\n\n")

    dims = compute_omega_dims_t19_sparse(max_deg=9, verify_dense_up_to=3)

    print()
    print("=" * 70)
    print("COMPARISON TABLE: Paley T_p Omega dims (per eigenspace, k=0)")
    print("=" * 70)
    print("  T_3  (chi=1, p=3):  [1, 1, 0]")
    print("  T_7  (chi=1, p=7):  [1, 3, 6, 9, 9, 6, 3]")
    print("  T_11 (chi=1, p=11): [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]")
    print(f"  T_19 (chi=?, p=19): {dims}")
    print()
    print(f"Results saved to: {RESULTS_PATH}")

    # Final save
    with open(RESULTS_PATH, 'a') as f:
        f.write(f"\nFinal Omega dims: {dims}\n")
        chi = sum((-1)**m * d for m, d in enumerate(dims))
        f.write(f"chi = {chi}\n")
        f.write("\nComparison:\n")
        f.write("  T_3:  [1, 1, 0]\n")
        f.write("  T_7:  [1, 3, 6, 9, 9, 6, 3]\n")
        f.write("  T_11: [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]\n")
        f.write(f"  T_19: {dims}\n")


if __name__ == '__main__':
    main()
    print("\nDONE.")
