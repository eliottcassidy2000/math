#!/usr/bin/env python3
"""
paley19_streaming_omega.py — opus-2026-03-13-S71b

Stream-compute Omega dims for P_19 using depth-first diff-seq enumeration.
Instead of storing all diff-seqs, enumerate them on-the-fly and build
the constraint matrix column by column.

Goal: compute chi(P_19) = sum (-1)^d Omega_d / p

Known Omega/p through degree 8: [1, 9, 72, 540, 3753, 23832, 136260, 688266, 2987622]
Need degrees 9-18 to determine chi.
"""

import sys, time
sys.path.insert(0, '04-computation')
from circulant_homology import (
    find_prime_for_roots, find_nth_root_of_unity,
    compute_face, sparse_rank_mod
)

def streaming_omega(n, S, prime, target_degree, verbose=True):
    """
    Compute Omega_{target_degree}/n using streaming enumeration.

    Instead of storing all diff-seqs of degree target_degree,
    enumerate them via DFS and build constraint columns one at a time.
    Also needs the set of allowed (target_degree-1)-diff-seqs for face checking.
    """
    S_sorted = sorted(S)
    d = target_degree

    # First, enumerate ALL (d-1)-diff-seqs and store as a set
    # This is the "allowed lower set" for face checking
    if verbose:
        print(f"  Enumerating degree {d-1} diff-seqs...", flush=True)

    lower_set = set()
    lower_set.add(())  # degree 0

    # Build degree 1 through d-1 iteratively
    current = {(): (frozenset([0]), 0)}  # seq -> (visited, last_vertex)
    for deg in range(1, d):
        next_level = {}
        for seq, (visited, last) in current.items():
            for s in S_sorted:
                nv = (last + s) % n
                if nv not in visited:
                    ns = seq + (s,)
                    next_level[ns] = (visited | {nv}, nv)
        current = next_level
        if deg == d - 1:
            lower_set = set(current.keys())

    n_lower = len(lower_set)
    if verbose:
        print(f"    |A_{d-1}| = {n_lower}", flush=True)

    # Now stream through degree-d diff-seqs
    # Build constraint columns on-the-fly
    if verbose:
        print(f"  Streaming degree {d} diff-seqs and building constraints...", flush=True)

    junk_idx = {}
    columns = []
    n_Am = 0

    # DFS enumeration of degree-d diff-seqs
    # Stack: (partial_seq, visited_set, last_vertex)
    t0 = time.time()

    # Build degree d-1 first, then extend by one step
    for seq, (visited, last) in current.items():
        for s in S_sorted:
            nv = (last + s) % n
            if nv not in visited:
                D = seq + (s,)
                n_Am += 1

                # Build constraint column for this diff-seq
                col = {}
                for fi in range(len(D) + 1):
                    fd = compute_face(D, fi, n)
                    if fd not in lower_set:
                        if fd not in junk_idx:
                            junk_idx[fd] = len(junk_idx)
                        row = junk_idx[fd]
                        sign = 1 if fi % 2 == 0 else -1
                        col[row] = (col.get(row, 0) + sign) % prime
                col = {r: v for r, v in col.items() if v != 0}
                columns.append(col)

                if verbose and n_Am % 5000000 == 0:
                    elapsed = time.time() - t0
                    print(f"    {n_Am/1e6:.1f}M diff-seqs, "
                          f"{len(junk_idx)} junk, {elapsed:.0f}s", flush=True)

    elapsed = time.time() - t0
    n_junk = len(junk_idx)
    if verbose:
        print(f"    Done: |A_{d}| = {n_Am}, junk rows = {n_junk}, "
              f"time = {elapsed:.1f}s", flush=True)

    # Compute rank
    if verbose:
        print(f"  Computing rank of {n_junk} x {n_Am} constraint matrix...", flush=True)

    if n_junk == 0:
        rank = 0
    else:
        t1 = time.time()
        rank, mpd = sparse_rank_mod(columns, prime, verbose=verbose and n_Am > 1000000)
        elapsed2 = time.time() - t1
        if verbose:
            print(f"    Rank = {rank}, mpd = {mpd}, time = {elapsed2:.1f}s", flush=True)

    omega = n_Am - rank
    return omega, n_Am

def main():
    p = 19
    S = set()
    for x in range(1, p):
        S.add((x * x) % p)

    prime = 191  # Known working prime for P_19

    print("=" * 70)
    print(f"P_{p} STREAMING OMEGA COMPUTATION")
    print(f"S = {sorted(S)}, working mod {prime}")
    print("=" * 70)

    # Known values
    known = [1, 9, 72, 540, 3753, 23832, 136260, 688266, 2987622]
    print(f"Known Omega/p (d=0..{len(known)-1}): {known}")
    chi_partial = sum((-1)**d * o for d, o in enumerate(known))
    print(f"Partial chi/p = {chi_partial}")
    print(f"Need remaining = {1 - chi_partial} for chi/p=1")

    # Try degree 9
    for d in range(9, 19):
        print(f"\n{'='*50}")
        print(f"DEGREE {d}")
        print(f"{'='*50}")

        t0 = time.time()
        try:
            omega, n_Am = streaming_omega(p, S, prime, d, verbose=True)
            elapsed = time.time() - t0

            omega_per_p = omega
            A_per_p = n_Am
            known.append(omega_per_p)
            chi_partial = sum((-1)**dd * o for dd, o in enumerate(known))

            print(f"  Omega_{d}/p = {omega_per_p}")
            print(f"  A_{d}/p = {A_per_p}")
            print(f"  Omega/A = {omega_per_p/A_per_p:.6f}")
            print(f"  chi/p (d=0..{d}) = {chi_partial}")
            print(f"  Total time: {elapsed:.1f}s")

        except MemoryError:
            print(f"  OOM at degree {d}")
            break
        except Exception as e:
            print(f"  Error at degree {d}: {e}")
            break

    print(f"\nFinal Omega/p: {known}")
    chi = sum((-1)**d * o for d, o in enumerate(known))
    print(f"chi/p = {chi}")
    print(f"chi = {chi * p}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
