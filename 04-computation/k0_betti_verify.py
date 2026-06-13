#!/usr/bin/env python3
"""
Compute k=0 eigenspace boundary ranks for Paley tournaments.
Verify β_d^(0) = 0 for intermediate degrees.

Key insight: rank(∂_d) + rank(∂_{d+1}) = Ω_d - β_d
Starting from rank(∂_1)=0 and known β_1=β_2=0,
we can determine rank(∂_d) up to the first unknown β.

This script computes boundary ranks via [constraint; boundary] combined matrix
to verify β_3=0, β_4=0, etc. for P_19.

opus-2026-03-13-S71c
"""

import sys, time
sys.stdout.reconfigure(line_buffering=True)

PRIME = 104729

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

def sparse_rank(cols_data, n_rows, prime):
    """Sparse Gaussian elimination mod prime."""
    pivots = {}
    rank = 0
    for col_idx, col in enumerate(cols_data):
        if col_idx % 20000 == 0 and col_idx > 0:
            print(f"    col {col_idx}/{len(cols_data)}, rank={rank}", flush=True)
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

def compute_k0_betti(p, max_d):
    """Compute k=0 eigenspace Betti numbers for P_p up to degree max_d."""
    S = qr(p)
    m = len(S)
    print(f"\nP_{p}: m={m}, QR={S}")
    t0 = time.time()

    # Enumerate diff-sequences
    prev_seqs = [()]
    prev_ps = {(): frozenset([0])}
    prev_last = {(): 0}
    all_seqs = {0: [()]}

    for d in range(1, max_d + 1):
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
        if not curr_seqs:
            break
        all_seqs[d] = curr_seqs
        prev_seqs = curr_seqs
        prev_ps = curr_ps
        prev_last = curr_last
        print(f"  d={d}: A={len(curr_seqs)} ({time.time()-t0:.1f}s)")

    actual_max_d = max(all_seqs.keys())
    max_d = min(max_d, actual_max_d)
    allowed = {d: set(all_seqs.get(d, [])) for d in range(max_d + 1)}

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

        junk_set = set()
        face_data_junk = []
        face_data_allowed = []

        for D in A_d:
            junk_faces = []
            allowed_faces = []
            for fi in range(d + 1):
                if fi == 0: fd = D[1:]
                elif fi == d: fd = D[:d-1]
                else:
                    merged = (D[fi-1] + D[fi]) % p
                    fd = D[:fi-1] + (merged,) + D[fi+1:]
                sign = 1 if fi % 2 == 0 else PRIME - 1
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
        idx_dm1 = {s: i for i, s in enumerate(A_dm1)}

        # Constraint columns
        constraint_cols = []
        for j, faces in enumerate(face_data_junk):
            col = {}
            for fd, coeff in faces:
                row = junk_idx[fd]
                col[row] = (col.get(row, 0) + coeff) % PRIME
                if col[row] == 0: del col[row]
            constraint_cols.append(col)

        print(f"  d={d}: computing constraint rank (A={n_Ad}, junk={n_junk})...", end=" ", flush=True)
        rank_constraint = sparse_rank(constraint_cols, n_junk, PRIME)
        omega_d = n_Ad - rank_constraint
        omega_dims.append(omega_d)
        print(f"Ω={omega_d}", flush=True)

        # Combined [constraint; boundary] columns
        combined_cols = []
        for j in range(n_Ad):
            col = {}
            for fd, coeff in face_data_junk[j]:
                row = junk_idx[fd]
                col[row] = (col.get(row, 0) + coeff) % PRIME
                if col[row] == 0: del col[row]
            for fd, coeff in face_data_allowed[j]:
                row = n_junk + idx_dm1[fd]
                col[row] = (col.get(row, 0) + coeff) % PRIME
                if col[row] == 0: del col[row]
            combined_cols.append(col)

        print(f"         computing combined rank ({n_junk}+{n_dm1} rows)...", end=" ", flush=True)
        rank_combined = sparse_rank(combined_cols, n_junk + n_dm1, PRIME)
        bd_rank = rank_combined - rank_constraint
        boundary_ranks[d] = bd_rank
        dt = time.time() - t0
        print(f"rank(∂)={bd_rank} ({dt:.1f}s)", flush=True)

    # Compute Betti
    betti = []
    for d in range(max_d + 1):
        od = omega_dims[d]
        rd = boundary_ranks.get(d, 0)
        rd1 = boundary_ranks.get(d+1, 0)
        bd = od - rd - rd1
        betti.append(bd)

    chi = sum((-1)**d * b for d, b in enumerate(betti))
    print(f"\n  Ω  = {omega_dims}")
    print(f"  ∂r = {[boundary_ranks.get(d,0) for d in range(max_d+1)]}")
    print(f"  β  = {betti}")
    print(f"  χ  = {chi}")

    # Verify rank formulas
    print(f"\n  Rank predictions (assuming β_d=0):")
    print(f"    rank(∂_1) = 0: actual {boundary_ranks.get(1,0)}")
    print(f"    rank(∂_2) = m = {m}: actual {boundary_ranks.get(2,0)}")
    print(f"    rank(∂_3) = m(m-2) = {m*(m-2)}: actual {boundary_ranks.get(3,0)}")
    if max_d >= 4:
        pred_r4 = omega_dims[3] - m*(m-2)
        print(f"    rank(∂_4) = Ω_3 - m(m-2) = {pred_r4}: actual {boundary_ranks.get(4,0)}")
    if max_d >= 5:
        pred_r5 = omega_dims[4] - boundary_ranks.get(4,0)
        print(f"    rank(∂_5) = Ω_4 - rank(∂_4) = {pred_r5}: actual {boundary_ranks.get(5,0)}")

    return omega_dims, boundary_ranks, betti

print("=" * 60)
print("k=0 EIGENSPACE BETTI NUMBERS — VERIFICATION")
print("=" * 60)

# P_7 (quick verification)
compute_k0_betti(7, 6)

# P_11 (verification)
compute_k0_betti(11, 10)

# P_19 — the key test! Compute up to d=8 if feasible (m=9)
print("\n" + "=" * 60)
print("P_19: TESTING β_3=β_4=...=β_8=0")
print("=" * 60)
compute_k0_betti(19, 8)

print("\n" + "=" * 60)
