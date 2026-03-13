#!/usr/bin/env python3
"""
Compute Betti numbers for each eigenspace of P_11.
Uses same approach as eigenspace_betti_p7.py but with sparse rank for larger matrices.

opus-2026-03-13-S71c
"""

import sys, time
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

def find_mod_prime(p, min_val=100000):
    """Find a prime q > min_val with q ≡ 1 mod p."""
    from sympy import isprime
    q = min_val + (p - min_val % p) % p + 1  # first ≡ 1 mod p above min_val
    while not isprime(q):
        q += p
    return q

MODP = None  # set per p

def sparse_rank_complex(cols_data, n_rows, prime):
    """Sparse Gaussian elimination mod prime, with complex phase factors
    rounded to nearest integer mod prime."""
    pivots = {}
    rank = 0
    for col_idx, col in enumerate(cols_data):
        work = dict(col)
        changed = True
        while changed:
            changed = False
            for r in sorted(work.keys()):
                if r in pivots and r in work:
                    piv_data = pivots[r]
                    factor = work[r]
                    for pr, pv in piv_data.items():
                        work[pr] = (work.get(pr, 0) - factor * pv) % prime
                        if work[pr] == 0:
                            del work[pr]
                    changed = True
                    break
        if not work:
            continue
        pivot_row = min(work.keys())
        inv = pow(work[pivot_row], prime - 2, prime)
        normalized = {r: (v * inv) % prime for r, v in work.items()}
        pivots[pivot_row] = normalized
        rank += 1
    return rank

def compute_eigenspace_betti_p11():
    global MODP
    p = 11
    MODP = find_mod_prime(p)
    print(f"Using modular prime MODP = {MODP} (≡ 1 mod {p})")

    QR = sorted(set((a*a)%p for a in range(1,p)))
    m = len(QR)
    print(f"P_{p}: m={m}, QR={QR}")

    t0 = time.time()

    # Enumerate diff-sequences
    prev_seqs = [()]
    prev_ps = {(): frozenset([0])}
    prev_last = {(): 0}
    all_seqs = {0: [()]}

    for d in range(1, p):
        curr_seqs = []
        curr_ps = {}
        curr_last = {}
        for seq in prev_seqs:
            ps_set = prev_ps[seq]
            last = prev_last[seq]
            for s in QR:
                ns = (last + s) % p
                if ns not in ps_set:
                    nseq = seq + (s,)
                    curr_seqs.append(nseq)
                    curr_ps[nseq] = ps_set | {ns}
                    curr_last[nseq] = ns
        if not curr_seqs:
            break
        all_seqs[d] = curr_seqs
        prev_seqs = curr_seqs
        prev_ps = curr_ps
        prev_last = curr_last
        print(f"  d={d}: A={len(curr_seqs)} ({time.time()-t0:.1f}s)")

    max_d = max(all_seqs.keys())
    allowed = {d: set(all_seqs.get(d, [])) for d in range(max_d + 1)}

    # For each eigenspace k, compute using sparse modular rank
    # Phase factor for d_0 face: omega^{k*s1} mapped to integer mod MODP
    # omega = g^((MODP-1)/p) where g is a generator of (Z/MODP)*

    # Find generator
    g = 2
    while pow(g, (MODP-1)//2, MODP) == 1:
        g += 1
    omega_mod = pow(g, (MODP-1)//p, MODP)
    # Verify: omega_mod^p = 1 mod MODP
    assert pow(omega_mod, p, MODP) == 1, f"omega^p = {pow(omega_mod, p, MODP)} != 1"

    results = {}

    # Only need to compute k=0 (done), one k∈QR, one k∈NQR
    # By Galois symmetry, all k∈QR give same Betti, all k∈NQR give same
    # By complex conjugation, QR and NQR give same Betti
    # So actually all k≠0 are the same!
    # But let's verify by computing k=0, k=1 (NQR), k=3 (QR)
    test_ks = [0, 1, 3]

    for k in test_ks:
        omega_k = pow(omega_mod, k, MODP)
        print(f"\n--- Eigenspace k={k} (omega^k = {omega_k} mod {MODP}) ---")

        omega_dims = [1]
        boundary_ranks = {0: 0}

        for d in range(1, max_d + 1):
            A_d = all_seqs[d]
            A_dm1 = all_seqs.get(d-1, [])
            n_Ad = len(A_d)
            n_dm1 = len(A_dm1)

            if n_Ad == 0:
                omega_dims.append(0)
                boundary_ranks[d] = 0
                continue

            # Identify junk and allowed faces, with phase factors
            junk_set = set()
            face_data_junk = []
            face_data_allowed = []

            for D in A_d:
                junk_faces = []
                allowed_faces = []
                for fi in range(d + 1):
                    if fi == 0:
                        fd = D[1:]
                        phase = pow(omega_k, D[0], MODP)  # omega^{k*s1}
                    elif fi == d:
                        fd = D[:d-1]
                        phase = 1
                    else:
                        merged = (D[fi-1] + D[fi]) % p
                        fd = D[:fi-1] + (merged,) + D[fi+1:]
                        phase = 1

                    sign_int = (1 if fi % 2 == 0 else MODP - 1)
                    coeff = (sign_int * phase) % MODP

                    is_allowed = (fd in allowed[d-1])
                    if not is_allowed:
                        junk_set.add(fd)
                        junk_faces.append((fd, coeff))
                    else:
                        allowed_faces.append((fd, coeff))
                face_data_junk.append(junk_faces)
                face_data_allowed.append(allowed_faces)

            junk_list = sorted(junk_set)
            n_junk = len(junk_list)
            junk_idx = {j: i for i, j in enumerate(junk_list)}
            idx_dm1 = {s: i for i, s in enumerate(A_dm1)}

            # Build constraint columns
            constraint_cols = []
            for j, faces in enumerate(face_data_junk):
                col = {}
                for fd, coeff in faces:
                    row = junk_idx[fd]
                    col[row] = (col.get(row, 0) + coeff) % MODP
                    if col[row] == 0:
                        del col[row]
                constraint_cols.append(col)

            rank_constraint = sparse_rank_complex(constraint_cols, n_junk, MODP)
            omega_d = n_Ad - rank_constraint
            omega_dims.append(omega_d)

            # Build combined columns [constraint; boundary]
            combined_cols = []
            for j in range(n_Ad):
                col = {}
                for fd, coeff in face_data_junk[j]:
                    row = junk_idx[fd]
                    col[row] = (col.get(row, 0) + coeff) % MODP
                    if col[row] == 0:
                        del col[row]
                for fd, coeff in face_data_allowed[j]:
                    row = n_junk + idx_dm1[fd]
                    col[row] = (col.get(row, 0) + coeff) % MODP
                    if col[row] == 0:
                        del col[row]
                combined_cols.append(col)

            rank_combined = sparse_rank_complex(combined_cols, n_junk + n_dm1, MODP)
            bd_rank = rank_combined - rank_constraint
            boundary_ranks[d] = bd_rank

            dt = time.time() - t0
            print(f"  d={d}: Ω={omega_d}, rank(∂)={bd_rank} ({dt:.1f}s)")

        # Compute Betti
        betti = []
        for d in range(max_d + 1):
            od = omega_dims[d]
            rd = boundary_ranks.get(d, 0)
            rd1 = boundary_ranks.get(d+1, 0)
            bd = od - rd - rd1
            betti.append(bd)

        chi = sum((-1)**d * b for d, b in enumerate(betti))
        print(f"  β = {betti}")
        print(f"  χ = {chi}")
        results[k] = betti

    print("\n" + "=" * 60)
    print("SUMMARY:")
    for k in test_ks:
        print(f"  k={k}: β = {results[k]}")

    # Reconstruct total
    b0 = results[0]  # k=0
    b1 = results[1]  # k∈NQR representative
    b3 = results[3]  # k∈QR representative

    print(f"\n  k=0 (1 eigenspace): {b0}")
    print(f"  k∈QR (5 eigenspaces, rep k=3): {b3}")
    print(f"  k∈NQR (5 eigenspaces, rep k=1): {b1}")

    total = [b0[d] + 5*b3[d] + 5*b1[d] for d in range(len(b0))]
    print(f"\n  Total β = {total}")
    print(f"  Expected: [1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0]")
    print("=" * 60)

compute_eigenspace_betti_p11()
