#!/usr/bin/env python3
"""
Compute Betti numbers for EACH eigenspace of the Paley tournament P_7.

Uses the rank formula: rank(∂_d^(k) on Omega) = rank([C; B^(k)]) - rank(C)
where C is the constraint matrix and B^(k) is the eigenspace-k boundary.

Key insight: For face d_0 in eigenspace k, the boundary includes a phase factor
omega^{k*s_1} (from shifting the path to start at vertex 0).

opus-2026-03-13-S71c
"""

import sys
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

def compute_eigenspace_betti(p, max_deg=None):
    """Compute Betti numbers for each eigenspace of P_p."""
    QR = sorted(set((a*a)%p for a in range(1,p)))
    m = len(QR)
    print(f"P_{p}: m={m}, QR={QR}")

    if max_deg is None:
        max_deg = p - 1

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
        print(f"  d={d}: A={len(curr_seqs)}")

    max_d = max(all_seqs.keys())
    allowed = {d: set(all_seqs.get(d, [])) for d in range(max_d + 1)}

    results = {}

    for k in range(p):
        omega_k = np.exp(2j * np.pi * k / p)
        print(f"\n--- Eigenspace k={k} ---")

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

            # Identify junk and allowed faces
            junk_set = set()
            face_data_junk = []
            face_data_allowed = []

            for D in A_d:
                junk_faces = []
                allowed_faces = []
                for fi in range(d + 1):
                    if fi == 0:
                        fd = D[1:]
                        # Phase for eigenspace k: omega^{k * s_1}
                        phase = omega_k ** D[0]
                    elif fi == d:
                        fd = D[:d-1]
                        phase = 1.0
                    else:
                        merged = (D[fi-1] + D[fi]) % p
                        fd = D[:fi-1] + (merged,) + D[fi+1:]
                        phase = 1.0

                    sign = (-1)**fi * phase
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

            # Build constraint matrix (junk part)
            C = np.zeros((n_junk, n_Ad), dtype=complex)
            for j, faces in enumerate(face_data_junk):
                for fd, sign in faces:
                    row = junk_idx[fd]
                    C[row, j] += sign

            rank_constraint = np.linalg.matrix_rank(C, tol=1e-8)
            omega_d = n_Ad - rank_constraint
            omega_dims.append(omega_d)

            # Build combined matrix [C; B^(k)]
            # B^(k) maps to A_{d-1} using allowed faces with phase
            combined = np.zeros((n_junk + n_dm1, n_Ad), dtype=complex)
            # Constraint part
            combined[:n_junk, :] = C
            # Boundary part
            for j in range(n_Ad):
                for fd, sign in face_data_allowed[j]:
                    row = n_junk + idx_dm1[fd]
                    combined[row, j] += sign

            rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)
            bd_rank = rank_combined - rank_constraint
            boundary_ranks[d] = bd_rank

        # Compute Betti
        betti = []
        for d in range(max_d + 1):
            od = omega_dims[d]
            rd = boundary_ranks.get(d, 0)
            rd1 = boundary_ranks.get(d+1, 0)
            bd = od - rd - rd1
            betti.append(bd)

        chi = sum((-1)**d * b for d, b in enumerate(betti))
        print(f"  Omega = {omega_dims}")
        print(f"  ∂ ranks = {[boundary_ranks.get(d, 0) for d in range(max_d + 1)]}")
        print(f"  β = {betti}")
        print(f"  χ = {chi}")

        results[k] = betti

    return results

# Main
p = 7
print("=" * 60)
print(f"EIGENSPACE BETTI NUMBERS FOR P_{p}")
print("=" * 60)

results = compute_eigenspace_betti(p)

print("\n" + "=" * 60)
print("SUMMARY:")
print("=" * 60)
for k in range(p):
    print(f"  k={k}: β = {results[k]}")

# Total Betti
max_len = max(len(v) for v in results.values())
total = [sum(results[k][d] if d < len(results[k]) else 0 for k in range(p)) for d in range(max_len)]
print(f"\n  Total β = {total}")
print(f"  Expected: [1, 0, 0, 0, 6, 0, 0]")

# Check: all k≠0 should have β_4^(k) = 1 only
for k in range(1, p):
    b = results[k]
    expected = [0] * len(b)
    expected[4] = 1  # m+1 = 3+1 = 4
    match = (b == expected)
    print(f"  k={k}: β={b}, expected {expected}, match={match}")
print("=" * 60)
