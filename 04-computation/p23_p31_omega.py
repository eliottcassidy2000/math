#!/usr/bin/env python3
"""Compute P_23 and P_31 k=0 Omega dims. opus-2026-03-13-S71c"""
import numpy as np
import sys, time
sys.stdout.reconfigure(line_buffering=True)

PRIME = 104729

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

def gauss_rank(C, prime):
    C = C.copy() % prime
    rows, cols = C.shape
    rank = 0
    pr = 0
    for col in range(cols):
        found = -1
        for row in range(pr, rows):
            if C[row, col] != 0:
                found = row
                break
        if found == -1:
            continue
        C[[pr, found]] = C[[found, pr]]
        inv = pow(int(C[pr, col]), prime - 2, prime)
        C[pr] = (C[pr] * inv) % prime
        for row in range(rows):
            if row != pr and C[row, col] != 0:
                factor = C[row, col]
                C[row] = (C[row] - factor * C[pr]) % prime
        rank += 1
        pr += 1
    return rank

def compute_k0_omega(p, max_deg):
    S = qr(p)
    m = len(S)
    print(f"P_{p}: m={m}")

    prev_seqs = [()]
    prev_ps = {(): frozenset([0])}
    prev_last = {(): 0}
    allowed_prev = set([()])
    omega_dims = [1]
    print(f"  d=0: Omega=1")

    for d in range(1, max_deg + 1):
        t0 = time.time()
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
            break

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
                is_allowed = (fd in allowed_prev)
                faces.append((fd, 1 if fi % 2 == 0 else -1, is_allowed))
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
            junk_idx = {j: i for i, j in enumerate(junk_list)}
            C = np.zeros((n_junk, n_Ad), dtype=np.int64)
            for j, faces in enumerate(face_data):
                for fd, sign, is_allowed in faces:
                    if not is_allowed:
                        row = junk_idx[fd]
                        C[row, j] = (C[row, j] + sign) % PRIME

            rk = gauss_rank(C, PRIME)
            omega_d = n_Ad - rk
            omega_dims.append(omega_d)
            dt = time.time() - t0
            print(f"  d={d}: A={n_Ad}, junk={n_junk}, rank={rk}, Omega={omega_d} ({dt:.1f}s)")

        allowed_prev = set(curr_seqs)
        prev_seqs = curr_seqs
        prev_ps = curr_ps
        prev_last = curr_last
        del face_data

    print(f"  Omega = {omega_dims}")
    return omega_dims

from math import comb

print("=" * 60)
o23 = compute_k0_omega(23, 4)
m23 = 11
print(f"  Formula check d=3: m(m-1)(2m-3)/2 = {m23*(m23-1)*(2*m23-3)//2}")
print(f"  A_4 formula: m(m-1)(2m^2-m-2)/2 = {m23*(m23-1)*(2*m23**2-m23-2)//2}")

print()
print("=" * 60)
o31 = compute_k0_omega(31, 4)
m31 = 15
print(f"  Formula check d=3: m(m-1)(2m-3)/2 = {m31*(m31-1)*(2*m31-3)//2}")
print(f"  A_4 formula: m(m-1)(2m^2-m-2)/2 = {m31*(m31-1)*(2*m31**2-m31-2)//2}")

# Junk rank analysis
print("\n" + "=" * 60)
print("JUNK RANK AT d=4")
print("=" * 60)
known = [
    (3, 7, 39, 9),
    (5, 11, 430, 205),
    (9, 19, 5436, 3753),
]
# Add P_23 and P_31 if computed
if len(o23) > 4:
    known.append((11, 23, 12595, o23[4]))
if len(o31) > 4:
    known.append((15, 31, 45465, o31[4]))

for m, p, A4, O4 in known:
    jr = A4 - O4
    print(f"  m={m:2d}: A_4={A4:6d}, Omega_4={O4:6d}, junk_rank={jr:5d}, jr/m={jr/m:.2f}, jr/(m(m-1))={jr/(m*(m-1)):.4f}")

# Try polynomial fit for junk_rank
print("\nTrying polynomial fit for junk_rank_4:")
ms = [row[0] for row in known]
jrs = [row[2] - row[3] for row in known]
for deg in range(1, min(len(known), 5)):
    coeffs = np.polyfit(ms, jrs, deg)
    print(f"  Degree {deg}: coeffs = {[f'{c:.2f}' for c in coeffs]}")
    fitted = [int(round(np.polyval(coeffs, m))) for m in ms]
    print(f"    Fitted: {fitted}")
    print(f"    Actual: {jrs}")
    print(f"    Match: {fitted == jrs}")
