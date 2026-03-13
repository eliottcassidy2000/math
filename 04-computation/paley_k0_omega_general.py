#!/usr/bin/env python3
"""
Compute k=0 eigenspace Omega dims for Paley tournaments P_p.

Uses the constraint matrix approach (correct GLMY definition):
Omega_d = ker(constraint_map) where constraint_map sends A_d -> non-A_{d-1}.

The k=0 eigenspace = diff-sequence complex (all Z_p translates identified).

opus-2026-03-13-S71c
"""

import numpy as np
import time
from math import comb

PRIME = 2147483647  # Mersenne prime for modular arithmetic

def quadratic_residues(p):
    qr = set()
    for a in range(1, p):
        qr.add((a * a) % p)
    return sorted(qr)

def enumerate_diff_seqs(S, p, max_deg):
    """Enumerate all valid diff sequences up to max_deg.
    A diff-seq (s1,...,sd) is valid iff partial sums 0, s1, s1+s2, ...
    are all distinct mod p.
    """
    S_sorted = sorted(S)
    seqs = {0: [()]}
    partial_sums = {(): frozenset([0])}
    partial_last = {(): 0}

    for d in range(1, max_deg + 1):
        prev = seqs[d - 1]
        new = []
        new_ps = {}
        new_last = {}
        for seq in prev:
            ps = partial_sums[seq]
            last = partial_last[seq]
            for s in S_sorted:
                new_last_sum = (last + s) % p
                if new_last_sum not in ps:
                    new_seq = seq + (s,)
                    new.append(new_seq)
                    new_ps[new_seq] = ps | frozenset([new_last_sum])
                    new_last[new_seq] = new_last_sum
        seqs[d] = new
        partial_sums.update(new_ps)
        partial_last.update(new_last)

    return seqs

def compute_face_diff(D, face_idx, p):
    """Compute face of diff-seq D at position face_idx.
    face_0: drop first vertex -> D[1:], offset = D[0]
    face_i (0 < i < d): merge D[i-1]+D[i] mod p, offset = 0
    face_d: drop last vertex -> D[:d-1], offset = 0
    """
    d = len(D)
    if face_idx == 0:
        return D[1:], D[0]
    elif face_idx == d:
        return D[:d - 1], 0
    else:
        merged = (D[face_idx - 1] + D[face_idx]) % p
        face = D[:face_idx - 1] + (merged,) + D[face_idx + 1:]
        return face, 0

def gauss_rank_mod(C, prime):
    """Gaussian elimination rank mod prime."""
    C = C.copy() % prime
    rows, cols = C.shape
    rank = 0
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if C[row, col] != 0:
                found = row
                break
        if found == -1:
            continue
        C[[pivot_row, found]] = C[[found, pivot_row]]
        inv = pow(int(C[pivot_row, col]), prime - 2, prime)
        C[pivot_row] = (C[pivot_row] * inv) % prime
        for row in range(rows):
            if row != pivot_row and C[row, col] != 0:
                factor = C[row, col]
                C[row] = (C[row] - factor * C[pivot_row]) % prime
        rank += 1
        pivot_row += 1
    return rank

def compute_k0_omega(p, max_deg=None):
    """Compute Omega dims for k=0 eigenspace of Paley P_p."""
    qr = quadratic_residues(p)
    m = len(qr)
    if max_deg is None:
        max_deg = min(m + 2, p - 1)

    print(f"P_{p}: m={m}, QR={qr}")

    # Enumerate diff-seqs up to max_deg + 1 (need d+1 for face computations at d)
    seqs = enumerate_diff_seqs(set(qr), p, max_deg + 1)

    # Build allowed sets
    allowed = {}
    for d in range(max_deg + 2):
        allowed[d] = set(seqs.get(d, []))

    omega_dims = []
    t0 = time.time()

    for d in range(max_deg + 1):
        A_d = seqs.get(d, [])
        n_Ad = len(A_d)

        if d == 0:
            omega_dims.append(1)
            print(f"  d=0: A=1, junk=0, Omega=1")
            continue

        if n_Ad == 0:
            omega_dims.append(0)
            print(f"  d={d}: A=0, Omega=0")
            continue

        # Find junk: non-allowed (d-1)-seqs that appear as faces of A_d elements
        junk_set = set()
        face_data = []

        for D in A_d:
            faces = []
            for fi in range(d + 1):
                fd, offset = compute_face_diff(D, fi, p)
                sign = 1 if fi % 2 == 0 else -1
                is_allowed = (fd in allowed[d - 1])
                faces.append((fd, sign, offset, is_allowed))
                if not is_allowed:
                    junk_set.add(fd)
            face_data.append(faces)

        junk_list = sorted(junk_set)
        n_junk = len(junk_list)

        if n_junk == 0:
            omega_dims.append(n_Ad)
            print(f"  d={d}: A={n_Ad}, junk=0, Omega={n_Ad}")
            continue

        # Build constraint matrix (n_junk x n_Ad)
        # For k=0 eigenspace: omega_k = 1, so phase = omega^offset = 1
        junk_idx = {j: i for i, j in enumerate(junk_list)}
        C = np.zeros((n_junk, n_Ad), dtype=np.int64)

        for j, faces in enumerate(face_data):
            for fd, sign, offset, is_allowed in faces:
                if not is_allowed:
                    row = junk_idx[fd]
                    # k=0: all phases are 1
                    entry = sign % PRIME
                    if entry < 0:
                        entry += PRIME
                    C[row, j] = (C[row, j] + entry) % PRIME

        rk = gauss_rank_mod(C, PRIME)
        omega_d = n_Ad - rk
        omega_dims.append(omega_d)
        elapsed = time.time() - t0
        print(f"  d={d}: A={n_Ad}, junk_rows={n_junk}, rank={rk}, Omega={omega_d} ({elapsed:.1f}s)")

    print(f"\n  Omega^(0) = {omega_dims}")
    return omega_dims

# Verification
print("=" * 70)
print("VERIFICATION: P_7 (expected [1, 3, 6, 9, 9, 6, 3])")
print("=" * 70)
o7 = compute_k0_omega(7)

print("\n" + "=" * 70)
print("VERIFICATION: P_11 (expected [1, 5, 20, 70, 205, 460])")
print("=" * 70)
o11 = compute_k0_omega(11, max_deg=6)

# New computations
print("\n" + "=" * 70)
print("NEW: P_19 (m=9)")
print("=" * 70)
o19 = compute_k0_omega(19, max_deg=6)

print("\n" + "=" * 70)
print("NEW: P_23 (m=11)")
print("=" * 70)
o23 = compute_k0_omega(23, max_deg=5)

# P_43 at low degrees
print("\n" + "=" * 70)
print("NEW: P_43 (m=21)")
print("=" * 70)
o43 = compute_k0_omega(43, max_deg=4)

# Analysis
print("\n" + "=" * 70)
print("RATIO ANALYSIS: Omega_d / C(m,d)")
print("=" * 70)
for label, omega, p_val in [("P_7", o7, 7), ("P_11", o11, 11),
                              ("P_19", o19, 19), ("P_23", o23, 23),
                              ("P_43", o43, 43)]:
    m = (p_val - 1) // 2
    ratios = []
    for d, od in enumerate(omega):
        b = comb(m, d)
        if b > 0:
            ratios.append(f"{od/b:.4f}")
        else:
            ratios.append("N/A")
    print(f"{label} (m={m}): {ratios}")

# Check consecutive ratios
print("\n" + "=" * 70)
print("CONSECUTIVE RATIOS: Omega_{d+1} / Omega_d")
print("=" * 70)
for label, omega, p_val in [("P_7", o7, 7), ("P_11", o11, 11),
                              ("P_19", o19, 19), ("P_23", o23, 23),
                              ("P_43", o43, 43)]:
    m = (p_val - 1) // 2
    ratios = []
    for d in range(len(omega) - 1):
        if omega[d] > 0:
            ratios.append(f"{omega[d+1]/omega[d]:.4f}")
        else:
            ratios.append("N/A")
    print(f"{label} (m={m}): {ratios}")

# Check if Omega_d matches C(m,d) * m^d / d! or similar
print("\n" + "=" * 70)
print("FORMULA CANDIDATES")
print("=" * 70)
for label, omega, p_val in [("P_7", o7, 7), ("P_11", o11, 11),
                              ("P_19", o19, 19), ("P_23", o23, 23),
                              ("P_43", o43, 43)]:
    m = (p_val - 1) // 2
    print(f"\n{label} (m={m}):")
    for d, od in enumerate(omega):
        candidates = []
        # C(m,d)
        candidates.append(f"C({m},{d})={comb(m,d)}")
        # C(2m,d)/C(m,0) = C(2m,d)
        candidates.append(f"C({2*m},{d})={comb(2*m,d)}")
        # d! * C(m,d)^2 / C(2m,d)?
        # Catalan
        if d > 0:
            candidates.append(f"C({m},{d})*{m}^{d}/{m**d}={comb(m,d)}")

        print(f"  d={d}: Omega={od}  | {', '.join(candidates)}")
