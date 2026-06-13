#!/usr/bin/env python3
"""
integer_betti_n11.py — Integer Betti numbers for I_11 and P_11

Computes GLMY path homology Betti numbers using modular arithmetic
(Gaussian elimination mod small primes) to avoid floating-point errors.

KEY INSIGHT (stacking trick):
  rank(d_m restricted to Omega_m) = rank([C_m; D_m]) - rank(C_m)

  where C_m is the constraint matrix (junk faces) and D_m is the
  boundary matrix (allowed faces). This avoids computing null-space bases.

Uses SPARSE column reduction for memory efficiency on large matrices.

Author: opus-2026-03-13
"""
import sys
import time
import gc

sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from circulant_homology import (
    CirculantHomology, find_nth_root_of_unity,
    enumerate_diff_seqs, compute_face_with_offset,
    is_prime
)
import numpy as np


def find_primes_for_n(n, count=3):
    """Find several primes q with n | (q-1) for multi-prime verification."""
    primes = []
    q = n + 1
    while len(primes) < count:
        if is_prime(q) and (q - 1) % n == 0:
            primes.append(q)
        q += 1
    return primes


def sparse_rank_from_cols(columns, n_rows, prime):
    """
    Compute rank of a sparse matrix given as list of column dicts.
    columns: list of {row: value} dicts
    Returns: rank
    """
    pivot_of = {}  # row -> reduced column dict
    rank = 0

    for col in columns:
        v = dict(col)
        while v:
            r0 = min(v.keys())
            val0 = v[r0]
            if r0 not in pivot_of:
                inv = pow(int(val0), prime - 2, prime)
                v = {r: (val * inv) % prime for r, val in v.items()}
                v = {r: val for r, val in v.items() if val != 0}
                pivot_of[r0] = v
                rank += 1
                break
            else:
                piv = pivot_of[r0]
                # v = v - val0 * piv
                for r, pv in piv.items():
                    new_val = (v.get(r, 0) - val0 * pv) % prime
                    if new_val == 0:
                        v.pop(r, None)
                    else:
                        v[r] = new_val
    return rank


def build_stacked_sparse_cols(face_data_list, junk_idx, A_m1_idx, n_junk, n, m, prime, omega_k):
    """
    Build sparse column representation of stacked matrix [C_m; D_m].

    Top n_junk rows: junk constraint entries
    Bottom len(A_m1_idx) rows: allowed boundary entries

    Returns list of column dicts {row: value mod prime}.
    """
    columns = []
    for j, faces in enumerate(face_data_list):
        col = {}
        for fd, fi, offset, is_allowed in faces:
            sign = 1 if fi % 2 == 0 else -1
            w = pow(omega_k, offset, prime) if offset != 0 else 1
            entry = (sign * w) % prime
            if entry == 0:
                continue

            if not is_allowed:
                row = junk_idx[fd]
            elif fd in A_m1_idx:
                row = n_junk + A_m1_idx[fd]
            else:
                continue

            old = col.get(row, 0)
            new = (old + entry) % prime
            if new == 0:
                col.pop(row, None)
            else:
                col[row] = new
        columns.append(col)
    return columns


def build_constraint_sparse_cols(face_data_list, junk_idx, n, m, prime, omega_k):
    """Build sparse column representation of constraint matrix C_m only."""
    columns = []
    for j, faces in enumerate(face_data_list):
        col = {}
        for fd, fi, offset, is_allowed in faces:
            if is_allowed:
                continue
            sign = 1 if fi % 2 == 0 else -1
            w = pow(omega_k, offset, prime) if offset != 0 else 1
            entry = (sign * w) % prime
            if entry == 0:
                continue
            row = junk_idx[fd]
            old = col.get(row, 0)
            new = (old + entry) % prime
            if new == 0:
                col.pop(row, None)
            else:
                col[row] = new
        columns.append(col)
    return columns


def compute_full_betti(n, S, name, primes):
    """Full Betti computation with stacking trick and sparse rank."""
    max_deg = n - 1

    print(f"\n{'='*70}")
    print(f"BETTI COMPUTATION: {name}")
    print(f"n={n}, S={sorted(S)}")
    print(f"{'='*70}")

    # Enumerate diff-seqs
    print(f"\nEnumerating diff-seqs up to degree {max_deg+1}...")
    t0 = time.time()
    diff_seqs = enumerate_diff_seqs(S, n, max_deg + 1)
    allowed_sets = {m: set(diff_seqs.get(m, [])) for m in range(max_deg + 3)}
    print(f"  Time: {time.time()-t0:.1f}s")

    for m in range(max_deg + 2):
        cnt = len(diff_seqs.get(m, []))
        if cnt > 0:
            print(f"    |A_{m}| = {cnt}")

    # Compute Omega dims (per-eigenspace, same for all k by THM-125)
    print(f"\nComputing Omega dims...")
    h_tmp = CirculantHomology(n=n, S=S, prime=primes[0])
    omega_per_eig = h_tmp.omega_dims(max_degree=max_deg, verbose=True)
    print(f"  Per-eigenspace Omega = {omega_per_eig}")
    del h_tmp
    gc.collect()

    # Precompute face data for each degree (shared across eigenspaces and primes)
    print(f"\nPrecomputing face data...")
    all_face_data = {}  # m -> (face_data_list, junk_idx, A_m1_idx)

    for m in range(1, max_deg + 2):
        A_m = diff_seqs.get(m, [])
        A_m1 = diff_seqs.get(m - 1, [])
        if not A_m:
            all_face_data[m] = ([], {}, {})
            continue

        allowed_lower = allowed_sets.get(m - 1, set())
        A_m1_idx = {d: i for i, d in enumerate(A_m1)} if A_m1 else {}
        junk_idx = {}

        face_data_list = []
        for D in A_m:
            faces = []
            for fi in range(m + 1):
                fd, offset = compute_face_with_offset(D, fi, n)
                is_allowed = (fd in allowed_lower)
                if not is_allowed and fd not in junk_idx:
                    junk_idx[fd] = len(junk_idx)
                faces.append((fd, fi, offset, is_allowed))
            face_data_list.append(faces)

        all_face_data[m] = (face_data_list, junk_idx, A_m1_idx)
        n_junk = len(junk_idx)
        n_Am1 = len(A_m1_idx)
        print(f"    m={m}: |A_m|={len(A_m)}, |junk|={n_junk}, |A_{m-1}|={n_Am1}")

    # Compute boundary ranks for each prime
    results_by_prime = {}

    for prime in primes[:2]:
        omega_p = find_nth_root_of_unity(n, prime)
        print(f"\n--- Computing boundary ranks with prime={prime}, omega_p={omega_p} ---")

        boundary_ranks = {}

        for m in range(max_deg + 2):
            if m == 0 or m not in all_face_data or not all_face_data[m][0]:
                boundary_ranks[m] = [0] * n
                print(f"  m={m:2d}: rank(d|Omega)=    0 (trivial)")
                sys.stdout.flush()
                continue

            face_data_list, junk_idx, A_m1_idx = all_face_data[m]
            n_junk = len(junk_idx)
            n_Am1 = len(A_m1_idx)

            t0 = time.time()
            ranks_m = []

            for k in range(n):
                omega_k = pow(omega_p, k, prime)

                # Build sparse constraint columns (C_m)
                C_cols = build_constraint_sparse_cols(
                    face_data_list, junk_idx, n, m, prime, omega_k
                )
                rank_C = sparse_rank_from_cols(C_cols, n_junk, prime)

                # Build sparse stacked columns [C_m; D_m]
                CD_cols = build_stacked_sparse_cols(
                    face_data_list, junk_idx, A_m1_idx, n_junk,
                    n, m, prime, omega_k
                )
                rank_CD = sparse_rank_from_cols(CD_cols, n_junk + n_Am1, prime)

                rank_d = rank_CD - rank_C
                ranks_m.append(rank_d)

                del C_cols, CD_cols
                gc.collect()

            dt = time.time() - t0
            boundary_ranks[m] = ranks_m

            all_same = all(r == ranks_m[0] for r in ranks_m)
            if all_same:
                print(f"  m={m:2d}: rank(d|Omega)={ranks_m[0]:5d} (all k same)  ({dt:.1f}s)")
            else:
                rk0 = ranks_m[0]
                rest_same = all(r == ranks_m[1] for r in ranks_m[1:])
                if rest_same:
                    print(f"  m={m:2d}: rank(d|Omega) k=0:{rk0:5d}, k!=0:{ranks_m[1]:5d}  ({dt:.1f}s)")
                else:
                    print(f"  m={m:2d}: ranks={ranks_m}  ({dt:.1f}s)")
            sys.stdout.flush()
            gc.collect()

        results_by_prime[prime] = boundary_ranks

    # Verify agreement across primes
    if len(results_by_prime) >= 2:
        p0, p1 = primes[0], primes[1]
        print(f"\n  Cross-prime verification (prime {p0} vs {p1}):")
        all_match = True
        for m in range(max_deg + 2):
            r0 = results_by_prime[p0][m]
            r1 = results_by_prime[p1][m]
            if r0 != r1:
                print(f"    m={m}: MISMATCH! {r0} vs {r1}")
                all_match = False
        print(f"    All match: {all_match}")

    # Compute Betti numbers
    boundary_ranks = results_by_prime[primes[0]]
    print(f"\n  Betti numbers:")
    betti = []
    for m in range(max_deg + 1):
        total_beta_m = 0
        details = []
        for k in range(n):
            omega_mk = omega_per_eig[m]  # per-eigenspace Omega dim (same for all k)
            rank_dm = boundary_ranks[m][k]
            rank_dm1 = boundary_ranks.get(m + 1, [0] * n)[k]
            beta_mk = omega_mk - rank_dm - rank_dm1
            total_beta_m += beta_mk
            details.append(beta_mk)
        betti.append(total_beta_m)

        all_same = all(d == details[0] for d in details)
        if all_same:
            print(f"    beta_{m:2d} = {total_beta_m:5d}  (per-eig: {details[0]})")
        else:
            rest_same = all(d == details[1] for d in details[1:])
            if rest_same:
                print(f"    beta_{m:2d} = {total_beta_m:5d}  (k=0:{details[0]}, k!=0:{details[1]})")
            else:
                print(f"    beta_{m:2d} = {total_beta_m:5d}  per-eig={details}")

    chi = sum((-1)**i * b for i, b in enumerate(betti))
    print(f"\n  beta = {betti}")
    print(f"  chi  = {chi}")
    print(f"  Expected chi = {n}")

    all_div = all(b % n == 0 for b in betti)
    print(f"  All beta divisible by {n}: {all_div}")
    if all_div:
        print(f"  beta/{n} = {[b // n for b in betti]}")

    # Cross-check with known Paley
    known_paley = [1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0]
    if "Paley" in name and len(betti) == len(known_paley):
        print(f"\n  Known Paley: {known_paley}")
        print(f"  Match: {betti == known_paley}")

    return betti


def main():
    print("="*70)
    print("INTEGER BETTI COMPUTATION FOR n=11 TOURNAMENTS")
    print("Using modular arithmetic with stacking trick (sparse)")
    print("rank(d_m|_Omega) = rank([C_m; D_m]) - rank(C_m)")
    print("="*70)

    n = 11
    primes = find_primes_for_n(n, count=3)
    print(f"\nPrimes with {n} | (q-1): {primes}")

    # First validate with Paley (known result)
    S_paley = set((a*a) % n for a in range(1, n))
    print(f"\nPaley QR_11 = {sorted(S_paley)}")
    betti_paley = compute_full_betti(n, S_paley, "Paley P_11", primes)

    # Then Interval
    S_interval = set(range(1, 6))
    print(f"\nInterval S = {sorted(S_interval)}")
    betti_interval = compute_full_betti(n, S_interval, "Interval I_11", primes)

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"  Paley P_11:    beta = {betti_paley}")
    print(f"  Interval I_11: beta = {betti_interval}")
    print(f"\nDONE.")


if __name__ == '__main__':
    main()
