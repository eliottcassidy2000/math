#!/usr/bin/env python3
"""
Compute k=0 eigenspace Betti numbers for Paley P_p using sparse modular arithmetic.

Strategy: Instead of computing null space bases, use the rank formula.
For the chain complex Omega_d -> Omega_{d-1} -> ...:
  beta_d = dim(Omega_d) - rank(∂_d) - rank(∂_{d+1})

The key insight is that we need rank(∂_d) for each d.
∂_d: Omega_d -> Omega_{d-1} is the boundary map on the Omega complex.

Since Omega_d = ker(constraint_d), we have:
  rank(∂_d) = rank of the composition: ker(constraint_d) -> A_d -> A_{d-1} -> ...

We can compute rank(∂_d) as:
  rank(∂_d) = dim(Omega_d) - dim(ker(∂_d|_{Omega_d}))
  = dim(Omega_d) - dim(ker(∂_d) ∩ Omega_d)

Alternatively, using the exact sequence structure:
  rank(∂_d) = dim(Omega_d) - dim(ker(∂_d) ∩ Omega_d)

The combined system [constraint_d; boundary_d] has null space = ker(∂_d) ∩ Omega_d.
So rank(∂_d) = dim(Omega_d) - nullity([constraint_d; boundary_d])
= dim(Omega_d) - (|A_d| - rank([constraint_d; boundary_d]))
= dim(Omega_d) - |A_d| + rank([constraint_d; boundary_d])

And dim(Omega_d) = |A_d| - rank(constraint_d), so:
rank(∂_d) = rank([constraint_d; boundary_d]) - rank(constraint_d)

This is MUCH easier — we just need ranks of two matrices!

opus-2026-03-13-S71c
"""

import sys, time
sys.stdout.reconfigure(line_buffering=True)

PRIME = 104729

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

def sparse_rank(cols_data, n_rows, prime):
    """Sparse Gaussian elimination mod prime. Returns rank."""
    pivots = {}
    rank = 0
    for col_idx, col in enumerate(cols_data):
        work = dict(col)
        for r in sorted(work.keys()):
            if r in pivots and r in work:
                piv_data = pivots[r]
                factor = work[r]
                for pr, pv in piv_data.items():
                    work[pr] = (work.get(pr, 0) - factor * pv) % prime
                    if work[pr] == 0:
                        del work[pr]
        if not work:
            continue
        pivot_row = min(work.keys())
        inv = pow(work[pivot_row], prime - 2, prime)
        normalized = {r: (v * inv) % prime for r, v in work.items()}
        pivots[pivot_row] = normalized
        rank += 1
    return rank

def compute_k0_betti(p, max_deg=None):
    S = qr(p)
    m = len(S)
    if max_deg is None:
        max_deg = p - 1

    print(f"P_{p}: m={m}, QR={S}")

    t0 = time.time()

    # Enumerate diff-sequences
    prev_seqs = [()]
    prev_ps = {(): frozenset([0])}
    prev_last = {(): 0}
    all_seqs = {0: [()]}

    for d in range(1, max_deg + 1):
        curr_seqs = []
        curr_ps = {}
        curr_last = {}
        for seq in prev_seqs:
            ps_set = prev_ps[seq]
            last = prev_last[seq]
            for s in S:
                ns = (last + s) % p
                if ns not in ps_set:
                    nseq = seq + (s,)
                    curr_seqs.append(nseq)
                    curr_ps[nseq] = ps_set | {ns}
                    curr_last[nseq] = ns
        all_seqs[d] = curr_seqs
        prev_seqs = curr_seqs
        prev_ps = curr_ps
        prev_last = curr_last
        print(f"  d={d}: A={len(curr_seqs)} ({time.time()-t0:.1f}s)")

    allowed = {d: set(all_seqs.get(d, [])) for d in range(max_deg + 1)}

    # For each degree d, compute:
    # rank(constraint_d) -> gives Omega_d = A_d - rank(constraint_d)
    # rank([constraint_d; boundary_d]) -> gives rank(∂_d) = rank_combined - rank_constraint

    omega_dims = [1]
    boundary_ranks = {0: 0}  # ∂_0 doesn't exist (or has rank 0)

    for d in range(1, max_deg + 1):
        A_d = all_seqs[d]
        A_dm1 = all_seqs.get(d-1, [])
        n_Ad = len(A_d)
        n_dm1 = len(A_dm1)

        if n_Ad == 0:
            omega_dims.append(0)
            boundary_ranks[d] = 0
            continue

        # Build constraint columns (junk faces)
        junk_set = set()
        face_data_junk = []  # For constraint matrix
        face_data_allowed = []  # For boundary matrix

        for D in A_d:
            junk_faces = []
            allowed_faces = []
            for fi in range(d + 1):
                if fi == 0: fd = D[1:]
                elif fi == d: fd = D[:d-1]
                else:
                    merged = (D[fi-1] + D[fi]) % p
                    fd = D[:fi-1] + (merged,) + D[fi+1:]
                sign = 1 if fi % 2 == 0 else -1
                is_allowed = (fd in allowed[d-1])
                if not is_allowed:
                    junk_set.add(fd)
                    junk_faces.append((fd, sign))
                else:
                    allowed_faces.append((fd, sign))
            face_data_junk.append(junk_faces)
            face_data_allowed.append(allowed_faces)

        junk_list = sorted(junk_set)
        n_junk = len(junk_list)
        junk_idx = {j: i for i, j in enumerate(junk_list)}

        # Build sparse columns for constraint matrix
        constraint_cols = []
        for j, faces in enumerate(face_data_junk):
            col = {}
            for fd, sign in faces:
                row = junk_idx[fd]
                col[row] = (col.get(row, 0) + sign) % PRIME
                if col[row] == 0:
                    del col[row]
            constraint_cols.append(col)

        rank_constraint = sparse_rank(constraint_cols, n_junk, PRIME)
        omega_d = n_Ad - rank_constraint
        omega_dims.append(omega_d)

        if n_junk == 0:
            # No constraint, so Omega_d = A_d
            # boundary_d just maps to A_{d-1}
            pass

        # Build combined columns: [constraint; boundary]
        # The boundary part maps to A_{d-1}
        idx_dm1 = {s: i for i, s in enumerate(A_dm1)}

        combined_cols = []
        for j in range(n_Ad):
            col = {}
            # Constraint part (rows 0 to n_junk-1)
            for fd, sign in face_data_junk[j]:
                row = junk_idx[fd]
                col[row] = (col.get(row, 0) + sign) % PRIME
                if col[row] == 0:
                    del col[row]
            # Boundary part (rows n_junk to n_junk + n_dm1 - 1)
            for fd, sign in face_data_allowed[j]:
                row = n_junk + idx_dm1[fd]
                col[row] = (col.get(row, 0) + sign) % PRIME
                if col[row] == 0:
                    del col[row]
            combined_cols.append(col)

        rank_combined = sparse_rank(combined_cols, n_junk + n_dm1, PRIME)
        bd_rank = rank_combined - rank_constraint
        boundary_ranks[d] = bd_rank

        dt = time.time() - t0
        print(f"  d={d}: Omega={omega_d}, rank(∂)={bd_rank} ({dt:.1f}s)")

    # Compute Betti numbers
    betti = []
    for d in range(max_deg + 1):
        od = omega_dims[d]
        rd = boundary_ranks.get(d, 0)  # rank of ∂_d (from Omega_d to Omega_{d-1})
        rd1 = boundary_ranks.get(d+1, 0)  # rank of ∂_{d+1} (image into Omega_d)
        beta_d = od - rd - rd1
        betti.append(beta_d)

    chi = sum((-1)**d * b for d, b in enumerate(betti))

    print(f"\n  Omega_k0 = {omega_dims}")
    print(f"  Boundary ranks = {[boundary_ranks.get(d, 0) for d in range(max_deg + 1)]}")
    print(f"  Betti_k0 = {betti}")
    print(f"  chi_k0 = {chi}")
    return omega_dims, betti

# P_7 (verification)
print("=" * 60)
o7, b7 = compute_k0_betti(7)

# P_11
print("\n" + "=" * 60)
o11, b11 = compute_k0_betti(11)
