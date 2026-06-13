#!/usr/bin/env python3
"""
Compute Omega_6 (k=0 eigenspace) for Paley P_p using sparse modular arithmetic.

Uses the rank formula: Omega_d = A_d - rank(constraint_d)
where constraint_d maps d-sequences to their non-allowed (d-1)-faces.

opus-2026-03-13-S71c
"""

import sys, time
sys.stdout.reconfigure(line_buffering=True)

PRIME = 104729

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

def sparse_rank(cols_data, n_rows, prime):
    """Sparse Gaussian elimination mod prime with correct fill-in handling."""
    pivots = {}
    rank = 0
    for col_idx, col in enumerate(cols_data):
        if col_idx % 50000 == 0 and col_idx > 0:
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

def compute_omega_up_to(p, max_d):
    """Compute Omega_d for d = 0, 1, ..., max_d."""
    S = qr(p)
    m = len(S)
    print(f"\nP_{p}: m={m}, QR={S[:5]}...")

    t0 = time.time()

    # Enumerate diff-sequences up to max_d
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
        all_seqs[d] = curr_seqs
        print(f"  d={d}: A={len(curr_seqs)} ({time.time()-t0:.1f}s)")
        prev_seqs = curr_seqs
        prev_ps = curr_ps
        prev_last = curr_last

    allowed = {d: set(all_seqs.get(d, [])) for d in range(max_d + 1)}

    # Compute Omega_d for each d
    omega_dims = [1]  # Omega_0 = 1

    for d in range(1, max_d + 1):
        A_d = all_seqs[d]
        n_Ad = len(A_d)

        if n_Ad == 0:
            omega_dims.append(0)
            continue

        # Build constraint columns (junk faces)
        junk_set = set()
        face_data_junk = []

        for D in A_d:
            junk_faces = []
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
            face_data_junk.append(junk_faces)

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

        dt = time.time() - t0
        print(f"  d={d}: A={n_Ad}, junk={n_junk}, rank={rank_constraint}, Omega={omega_d} ({dt:.1f}s)")

    return omega_dims

# Known values for verification
print("=" * 60)
print("Computing Omega dimensions for k=0 eigenspace")
print("=" * 60)

# P_7 (verification)
o7 = compute_omega_up_to(7, 6)
print(f"\n  P_7: Omega = {o7}")
print(f"  Expected: [1, 3, 6, 10, 10, 6, 1]")

# P_11
o11 = compute_omega_up_to(11, 6)
print(f"\n  P_11: Omega[0:7] = {o11}")
print(f"  Expected: [1, 5, 20, 70, 205, 460, 700]")

# P_19 — this is the new computation
o19 = compute_omega_up_to(19, 6)
print(f"\n  P_19: Omega[0:7] = {o19}")

# P_23
o23 = compute_omega_up_to(23, 6)
print(f"\n  P_23: Omega[0:7] = {o23}")

# P_31 — try if feasible
o31 = compute_omega_up_to(31, 6)
print(f"\n  P_31: Omega[0:7] = {o31}")

print("\n" + "=" * 60)
print("SUMMARY — Omega_6 values:")
print(f"  P_7  (m=3):  {o7[6] if len(o7)>6 else '?'}")
print(f"  P_11 (m=5):  {o11[6] if len(o11)>6 else '?'}")
print(f"  P_19 (m=9):  {o19[6] if len(o19)>6 else '?'}")
print(f"  P_23 (m=11): {o23[6] if len(o23)>6 else '?'}")
print(f"  P_31 (m=15): {o31[6] if len(o31)>6 else '?'}")
print("=" * 60)
