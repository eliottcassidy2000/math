#!/usr/bin/env python3
"""
Compute k=0 eigenspace Omega dims for Paley tournaments P_p.
Uses constraint matrix approach. Memory-efficient version.
opus-2026-03-13-S71c
"""
import numpy as np
import time, sys
from math import comb
sys.stdout.reconfigure(line_buffering=True)

PRIME = 104729  # smaller prime for faster modular arithmetic

def qr(p):
    return sorted(set((a*a) % p for a in range(1, p)))

def compute_k0_omega(p, max_deg):
    """Compute k=0 Omega dims using constraint matrix."""
    S = qr(p)
    m = len(S)
    S_set = set(S)
    print(f"P_{p}: m={m}, QR={S}")

    # Enumerate diff-seqs incrementally
    # Store as list of tuples + partial sum info
    prev_seqs = [()]
    prev_ps = {(): frozenset([0])}
    prev_last = {(): 0}

    omega_dims = [1]
    print(f"  d=0: Omega=1")

    # Also need to track allowed sets for junk checking
    allowed_prev = set([()])  # d=0 allowed

    for d in range(1, max_deg + 1):
        t0 = time.time()

        # Generate A_d
        curr_seqs = []
        curr_ps = {}
        curr_last = {}
        for seq in prev_seqs:
            ps = prev_ps[seq]
            last = prev_last[seq]
            for s in S:
                ns = (last + s) % p
                if ns not in ps:
                    nseq = seq + (s,)
                    curr_seqs.append(nseq)
                    curr_ps[nseq] = ps | {ns}
                    curr_last[nseq] = ns

        n_Ad = len(curr_seqs)
        if n_Ad == 0:
            omega_dims.append(0)
            print(f"  d={d}: A=0, Omega=0")
            break

        # Also need A_{d+1} for face checking? No - we need A_{d-1} = allowed_prev
        # Find junk faces
        junk_set = set()
        face_data = []

        for D in curr_seqs:
            faces = []
            for fi in range(d + 1):
                if fi == 0:
                    fd = D[1:]
                elif fi == d:
                    fd = D[:d-1]
                else:
                    merged = (D[fi-1] + D[fi]) % p
                    fd = D[:fi-1] + (merged,) + D[fi+1:]

                sign = 1 if fi % 2 == 0 else -1
                is_allowed = (fd in allowed_prev)
                faces.append((fd, sign, is_allowed))
                if not is_allowed:
                    junk_set.add(fd)
            face_data.append(faces)

        junk_list = sorted(junk_set)
        n_junk = len(junk_list)

        if n_junk == 0:
            omega_dims.append(n_Ad)
            dt = time.time() - t0
            print(f"  d={d}: A={n_Ad}, junk=0, Omega={n_Ad} ({dt:.1f}s)")
        else:
            # Build constraint matrix
            junk_idx = {j: i for i, j in enumerate(junk_list)}
            C = np.zeros((n_junk, n_Ad), dtype=np.int64)

            for j, faces in enumerate(face_data):
                for fd, sign, is_allowed in faces:
                    if not is_allowed:
                        row = junk_idx[fd]
                        entry = sign % PRIME
                        C[row, j] = (C[row, j] + entry) % PRIME

            # Gaussian elimination
            rk = 0
            pivot_row = 0
            for col in range(n_Ad):
                found = -1
                for row in range(pivot_row, n_junk):
                    if C[row, col] != 0:
                        found = row
                        break
                if found == -1:
                    continue
                C[[pivot_row, found]] = C[[found, pivot_row]]
                inv = pow(int(C[pivot_row, col]), PRIME - 2, PRIME)
                C[pivot_row] = (C[pivot_row] * inv) % PRIME
                for row in range(n_junk):
                    if row != pivot_row and C[row, col] != 0:
                        factor = C[row, col]
                        C[row] = (C[row] - factor * C[pivot_row]) % PRIME
                rk += 1
                pivot_row += 1

            omega_d = n_Ad - rk
            omega_dims.append(omega_d)
            dt = time.time() - t0
            print(f"  d={d}: A={n_Ad}, junk={n_junk}, rank={rk}, Omega={omega_d} ({dt:.1f}s)")

        # Update for next iteration
        allowed_prev = set(curr_seqs)
        prev_seqs = curr_seqs
        prev_ps = curr_ps
        prev_last = curr_last

        # Memory cleanup
        del face_data

    print(f"  Omega^(0) = {omega_dims}")
    return omega_dims

# P_7
print("=" * 60)
print("P_7 (expected [1, 3, 6, 9, 9, 6, 3])")
print("=" * 60)
o7 = compute_k0_omega(7, 6)

# P_11
print("\n" + "=" * 60)
print("P_11 (expected [1, 5, 20, 70, 205, 460, 700, ...])")
print("=" * 60)
o11 = compute_k0_omega(11, 7)

# P_19
print("\n" + "=" * 60)
print("P_19 (m=9) — NEW")
print("=" * 60)
o19 = compute_k0_omega(19, 5)

# P_23
print("\n" + "=" * 60)
print("P_23 (m=11) — NEW")
print("=" * 60)
o23 = compute_k0_omega(23, 4)

# Analysis
print("\n" + "=" * 60)
print("FORMULA SEARCH")
print("=" * 60)
for label, omega, p_val in [("P_7", o7, 7), ("P_11", o11, 11),
                             ("P_19", o19, 19), ("P_23", o23, 23)]:
    m = (p_val - 1) // 2
    print(f"\n{label} (m={m}): Omega = {omega}")
    for d, od in enumerate(omega):
        b = comb(m, d)
        ratio = od / b if b > 0 else None
        print(f"  d={d}: Omega={od}, C({m},{d})={b}, ratio={ratio}")

# Consecutive ratios
print("\n" + "=" * 60)
print("CONSECUTIVE RATIOS Omega_{d+1}/Omega_d")
print("=" * 60)
for label, omega, p_val in [("P_7", o7, 7), ("P_11", o11, 11),
                             ("P_19", o19, 19), ("P_23", o23, 23)]:
    m = (p_val - 1) // 2
    rs = [omega[d+1]/omega[d] if omega[d] > 0 else None for d in range(len(omega)-1)]
    print(f"{label} (m={m}): {[f'{r:.4f}' if r else 'N/A' for r in rs]}")
