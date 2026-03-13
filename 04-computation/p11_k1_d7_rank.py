#!/usr/bin/env python3
"""
P_11 eigenspace k=1, degree 7: compute rank of boundary map ∂_7^{(1)}|_{Ω_7}
using the stacking trick over GF(q).

Uses scipy sparse matrices + sparse GF rank via dense subblocks.
Key insight: use column-pivoting Gaussian elimination on the TRANSPOSE
(work with the smaller dimension).

Actually, the simplest reliable approach: build the matrices as scipy sparse,
convert to dense column-by-column blocks, and do GF(q) Gaussian elimination
using vectorized numpy operations on int16 arrays (q=23 fits in int16).
"""

import numpy as np
from scipy import sparse
from time import time
import sys

p = 11
q = 23  # working prime

# QR set for p=11
QR = set()
for x in range(1, p):
    QR.add(pow(x, 2, p))
QR_list = sorted(QR)
print(f"QR = {QR_list}")
sys.stdout.flush()

# Primitive root and omega
def primitive_root(q):
    for g in range(2, q):
        if all(pow(g, (q-1)//pf, q) != 1 for pf in [2, 11]):  # factors of q-1=22
            return g
    raise ValueError

g = primitive_root(q)
omega = pow(g, (q - 1) // p, q)
print(f"omega = {omega}, omega^11 mod 23 = {pow(omega, 11, q)}")
assert pow(omega, p, q) == 1 and omega != 1
sys.stdout.flush()

# Build diff-seqs
def build_diffseqs(d):
    results = []
    if d == 0:
        return [()]
    def backtrack(seq, ps_list, ps_set):
        if len(seq) == d:
            results.append(tuple(seq))
            return
        for s in QR_list:
            new_ps = (ps_list[-1] + s) % p
            if new_ps in ps_set:
                continue
            seq.append(s)
            ps_list.append(new_ps)
            ps_set.add(new_ps)
            backtrack(seq, ps_list, ps_set)
            seq.pop()
            ps_list.pop()
            ps_set.remove(new_ps)
    backtrack([], [0], {0})
    return results

t0 = time()
diffseqs_7 = build_diffseqs(7)
diffseqs_6 = build_diffseqs(6)
print(f"|A_7| = {len(diffseqs_7)}, |A_6| = {len(diffseqs_6)}  ({time()-t0:.1f}s)")
sys.stdout.flush()

idx7 = {s: i for i, s in enumerate(diffseqs_7)}
idx6 = {s: i for i, s in enumerate(diffseqs_6)}
n7 = len(diffseqs_7)
n6 = len(diffseqs_6)

# Build constraint matrix C_7 and boundary matrix B_7^{(1)} as sparse COO
# then stack and compute rank

print("Building matrices as sparse...")
sys.stdout.flush()
t0 = time()

# C_7: junk interior faces
# Collect junk 6-seqs and their indices
c_rows, c_cols, c_vals = [], [], []
junk_6seqs = {}
junk_count = 0

for col, seq in enumerate(diffseqs_7):
    for i in range(1, 7):  # interior face i (merges seq[i-1] and seq[i])
        merged = (seq[i - 1] + seq[i]) % p
        if merged not in QR:
            result = seq[:i-1] + (merged,) + seq[i+1:]
            if result not in junk_6seqs:
                junk_6seqs[result] = junk_count
                junk_count += 1
            row = junk_6seqs[result]
            sign = (1 if i % 2 == 0 else -1) % q
            c_rows.append(row)
            c_cols.append(col)
            c_vals.append(sign)

n_junk = junk_count
print(f"  Junk faces: {len(c_rows)} entries, {n_junk} distinct junk 6-seqs")

# B_7^{(1)}: boundary with phase on face 0
b_rows, b_cols, b_vals = [], [], []

for col, seq in enumerate(diffseqs_7):
    # Face 0: (s_2,...,s_7), coeff = omega^{s_1}
    face0 = seq[1:]
    if face0 in idx6:
        row = idx6[face0]
        coeff = pow(omega, seq[0], q)
        b_rows.append(row)
        b_cols.append(col)
        b_vals.append(coeff)

    # Interior faces 1..6 (valid only)
    for i in range(1, 7):
        merged = (seq[i - 1] + seq[i]) % p
        if merged in QR:
            result = seq[:i-1] + (merged,) + seq[i+1:]
            if result in idx6:
                row = idx6[result]
                sign = pow(-1, i, q)
                b_rows.append(row)
                b_cols.append(col)
                b_vals.append(sign)

    # Face 7: (s_1,...,s_6), coeff = (-1)^7 = -1
    face7 = seq[:6]
    if face7 in idx6:
        row = idx6[face7]
        b_rows.append(row)
        b_cols.append(col)
        b_vals.append((-1) % q)

print(f"  Boundary entries: {len(b_rows)}")
print(f"  Built in {time()-t0:.1f}s")
sys.stdout.flush()

# For the stacked matrix, shift B rows by n_junk
# Stacked = [C_7 (n_junk rows); B_7 (n6 rows)] with n7 columns
total_rows = n_junk + n6

all_rows = c_rows + [r + n_junk for r in b_rows]
all_cols = c_cols + b_cols
all_vals = c_vals + b_vals

print(f"\nStacked matrix: {total_rows} x {n7}")
print(f"Total nonzero entries: {len(all_rows)}")
sys.stdout.flush()

# Convert to dense int16 matrix — size = total_rows * n7 * 2 bytes
# = ~12000 * 8735 * 2 ≈ 200 MB — should be OK
print(f"Allocating dense matrix ({total_rows} x {n7}, int16, ~{total_rows * n7 * 2 // 1024 // 1024} MB)...")
sys.stdout.flush()

# Build stacked matrix
stacked = np.zeros((total_rows, n7), dtype=np.int16)
for r, c, v in zip(all_rows, all_cols, all_vals):
    stacked[r, c] = (stacked[r, c] + v) % q

# Also need C_7 alone for its rank
C7_dense = np.zeros((n_junk, n7), dtype=np.int16)
for r, c, v in zip(c_rows, c_cols, c_vals):
    C7_dense[r, c] = (C7_dense[r, c] + v) % q

print("Matrices built.")
sys.stdout.flush()

def gf_rank_fast(M, q):
    """GF(q) rank via vectorized Gaussian elimination. Works on int16/int32."""
    M = M.astype(np.int32) % q  # int32 for intermediate products
    nrows, ncols = M.shape
    pivot_row = 0

    for col in range(ncols):
        if pivot_row >= nrows:
            break

        # Find nonzero in column at or below pivot_row
        col_below = M[pivot_row:, col]
        nonzero = np.flatnonzero(col_below)
        if len(nonzero) == 0:
            continue

        found = pivot_row + nonzero[0]
        if found != pivot_row:
            M[[pivot_row, found]] = M[[found, pivot_row]]

        inv = pow(int(M[pivot_row, col]), q - 2, q)
        M[pivot_row] = (M[pivot_row] * inv) % q

        # Eliminate all other rows with nonzero in this column
        factors = M[:, col].copy()
        factors[pivot_row] = 0
        nz = np.flatnonzero(factors)
        if len(nz) > 0:
            # M[nz] = (M[nz] - outer(factors[nz], M[pivot_row])) % q
            # Do in chunks to avoid huge temporaries
            chunk = 2000
            for start in range(0, len(nz), chunk):
                end = min(start + chunk, len(nz))
                idx = nz[start:end]
                M[idx] = (M[idx] - np.outer(factors[idx], M[pivot_row])) % q

        pivot_row += 1

        if pivot_row % 200 == 0:
            print(f"  rank: {pivot_row} pivots found, col {col}/{ncols}")
            sys.stdout.flush()

    return pivot_row

print(f"\nComputing rank(C_7) over GF({q})...")
sys.stdout.flush()
t0 = time()
rank_C7 = gf_rank_fast(C7_dense, q)
print(f"rank(C_7) = {rank_C7}  ({time()-t0:.1f}s)")
print(f"Omega_7 = {n7} - {rank_C7} = {n7 - rank_C7}")
sys.stdout.flush()

print(f"\nComputing rank([C_7; B_7^{{(1)}}]) over GF({q})...")
sys.stdout.flush()
t0 = time()
rank_stacked = gf_rank_fast(stacked, q)
print(f"rank(stacked) = {rank_stacked}  ({time()-t0:.1f}s)")
sys.stdout.flush()

R7_k1 = rank_stacked - rank_C7
print(f"\nR_7^{{(1)}} = {rank_stacked} - {rank_C7} = {R7_k1}")
print(f"\nPrediction: R_7^{{(1)}} = 390")
print(f"If R_6^{{(1)}} = 309 and Omega_6 = 700:")
print(f"  beta_6^{{(1)}} = 700 - 309 - {R7_k1} = {700 - 309 - R7_k1}")
sys.stdout.flush()

# Verify with q2=67
print(f"\n--- Verification with q=67 ---")
sys.stdout.flush()

q2 = 67
g2 = primitive_root(q2)
omega2 = pow(g2, (q2 - 1) // p, q2)
assert pow(omega2, p, q2) == 1 and omega2 != 1
print(f"omega2 = {omega2}")

# Rebuild stacked mod q2
stacked2 = np.zeros((total_rows, n7), dtype=np.int16)
# C_7 part (signs only, no omega)
for r, c, v in zip(c_rows, c_cols, c_vals):
    # v was computed mod 23, recompute sign
    pass

# Actually let's recompute from scratch for q2
C7_v2 = np.zeros((n_junk, n7), dtype=np.int16)
for r, c, v in zip(c_rows, c_cols, c_vals):
    # signs are ±1, recompute mod q2
    # v was (1 or q-1) mod q. Map back to ±1 then mod q2
    if v == q - 1:
        sv = q2 - 1
    else:
        sv = v  # 1
    C7_v2[r, c] = (C7_v2[r, c] + sv) % q2

# B_7 for q2: need to recompute omega values
b_rows2, b_cols2, b_vals2 = [], [], []
for col, seq in enumerate(diffseqs_7):
    face0 = seq[1:]
    if face0 in idx6:
        row = idx6[face0]
        coeff = pow(omega2, seq[0], q2)
        b_rows2.append(row)
        b_cols2.append(col)
        b_vals2.append(coeff)

    for i in range(1, 7):
        merged = (seq[i - 1] + seq[i]) % p
        if merged in QR:
            result = seq[:i-1] + (merged,) + seq[i+1:]
            if result in idx6:
                row = idx6[result]
                sign = pow(-1, i, q2)
                b_rows2.append(row)
                b_cols2.append(col)
                b_vals2.append(sign)

    face7 = seq[:6]
    if face7 in idx6:
        row = idx6[face7]
        b_rows2.append(row)
        b_cols2.append(col)
        b_vals2.append((-1) % q2)

B7_v2 = np.zeros((n6, n7), dtype=np.int16)
for r, c, v in zip(b_rows2, b_cols2, b_vals2):
    B7_v2[r, c] = (B7_v2[r, c] + v) % q2

stacked2 = np.vstack([C7_v2, B7_v2])

t0 = time()
rank_C7_v2 = gf_rank_fast(C7_v2.astype(np.int32), q2)
print(f"rank(C_7) over GF({q2}) = {rank_C7_v2}  ({time()-t0:.1f}s)")
sys.stdout.flush()

t0 = time()
rank_stacked_v2 = gf_rank_fast(stacked2.astype(np.int32), q2)
print(f"rank(stacked) over GF({q2}) = {rank_stacked_v2}  ({time()-t0:.1f}s)")

R7_v2 = rank_stacked_v2 - rank_C7_v2
print(f"R_7^{{(1)}} over GF(67) = {R7_v2}")

if R7_k1 == R7_v2:
    print(f"\n*** CONFIRMED: R_7^{{(1)}} = {R7_k1} ***")
else:
    print(f"\n*** MISMATCH: GF(23)={R7_k1}, GF(67)={R7_v2} ***")

print("\nDone.")
sys.stdout.flush()
