"""
P_11 k=0 Betti via boundary rank stacking trick.
Only compute boundary ranks for the degrees that matter for Betti.

The Omega values are KNOWN:
  [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]

For Betti: β_d = Omega_d - R_d - R_{d+1}
By HYP-672: β_d = 0 for d ∉ {0, 5, 6} (m=5).
So we expect R_d + R_{d+1} = Omega_d for d ∉ {0, 5, 6}.

Known: R_0 = 0.
From β_0 = 1: R_0 + R_1 = Omega_0 - β_0 = 1 - 1 = 0. So R_1 = 0. ✓

From β_1 = 0: R_1 + R_2 = Omega_1 = 5. So R_2 = 5.
From β_2 = 0: R_2 + R_3 = Omega_2 = 20. So R_3 = 15.
From β_3 = 0: R_3 + R_4 = Omega_3 = 70. So R_4 = 55.
From β_4 = 0: R_4 + R_5 = Omega_4 = 205. So R_5 = 150.
β_5 = Omega_5 - R_5 - R_6 = 460 - 150 - R_6 = 5.  So R_6 = 305.
β_6 = Omega_6 - R_6 - R_7 = 700 - 305 - R_7 = 5.  So R_7 = 390.
From β_7 = 0: R_7 + R_8 = Omega_7 = 690. So R_8 = 300.
From β_8 = 0: R_8 + R_9 = Omega_8 = 450. So R_9 = 150.
From β_9 = 0: R_9 + R_10 = Omega_9 = 180. So R_10 = 30.
From β_10 = 0: R_10 = Omega_10 - β_10 = 30 - 0 = 30. ✓ (and R_11 = 0).

So the PREDICTED boundary ranks are:
R = [0, 0, 5, 15, 55, 150, 305, 390, 300, 150, 30, 0]

These are consistent if β = [1, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0].

But wait: β_5 = 5 = m(m-3)/2 = 5*2/2 = 5. ✓
And β_6 = 5. But HYP-672 says β_6 = C(6,2) = 15. ✗!!

HOLD ON. HYP-672 says β_{m+1} = C(m+1,2) = 15 for total Betti.
But β_6^{(0)} = 5, and β_6 = β_6^{(0)} + sum_{k≠0} β_6^{(k)} = 5 + 10*1 = 15. ✓

OK so I was confusing per-eigenspace (k=0) with total Betti.
The per-eigenspace (k=0) Betti should be:
  β_0^{(0)} = 1, β_5^{(0)} = 5 = m(m-3)/2, β_6^{(0)} = 5 = m(m-3)/2
  All others 0.

This is ALREADY KNOWN from eigenspace_contractibility.out. The predicted
boundary ranks are:
"""
import sys

Omega = [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]
m = 5
p = 11

# Predicted boundary ranks from β^{(0)} = [1,0,0,0,0,5,5,0,0,0,0]:
R_pred = [0, 0, 5, 15, 55, 150, 305, 390, 300, 150, 30, 0]
beta_pred = [1, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0]

print("PREDICTED k=0 eigenspace for P_11 (m=5):")
print(f"{'d':>3} {'Omega':>8} {'R_d':>8} {'R_{d+1}':>8} {'beta_d':>8}")
for d in range(len(Omega)):
    print(f"{d:>3} {Omega[d]:>8} {R_pred[d]:>8} {R_pred[d+1]:>8} {beta_pred[d]:>8}")

# Verify: R_d + R_{d+1} + beta_d = Omega_d
for d in range(len(Omega)):
    check = R_pred[d] + R_pred[d+1] + beta_pred[d]
    ok = "✓" if check == Omega[d] else f"✗ (got {check})"
    print(f"  d={d}: {R_pred[d]} + {R_pred[d+1]} + {beta_pred[d]} = {check} = Omega_{d} {ok}")

# Now let me verify a few of these boundary ranks computationally
print("\n\nVerifying boundary ranks computationally (small degrees)...")

import numpy as np

def get_QR(p):
    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))
    return QR

def build_diffseqs(p, d):
    QR = get_QR(p)
    QR_list = sorted(QR)
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

def mod_rank(matrix, prime=3):
    M = matrix.copy().astype(np.int64) % prime
    rows, cols = M.shape
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if M[row, col] % prime != 0:
                pivot = row
                break
        if pivot is None:
            continue
        M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = (M[rank] * inv) % prime
        for row in range(rows):
            if row != rank and M[row, col] % prime != 0:
                M[row] = (M[row] - M[row, col] * M[rank]) % prime
        rank += 1
    return rank

QR = get_QR(p)

# Precompute diff-seqs for small degrees
diffseqs = {}
for d in range(8):
    diffseqs[d] = build_diffseqs(p, d)
    print(f"  d={d}: |A_d| = {len(diffseqs[d])}", flush=True)

def build_constraint(d, Ad):
    junk_faces = {}
    for seq in Ad:
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:
                face = list(seq)
                face[i-1] = merged
                del face[i]
                ft = tuple(face)
                if ft not in junk_faces:
                    junk_faces[ft] = len(junk_faces)
    n_rows = len(junk_faces)
    if n_rows == 0:
        return np.zeros((0, len(Ad)), dtype=np.int8)
    C = np.zeros((n_rows, len(Ad)), dtype=np.int8)
    for j, seq in enumerate(Ad):
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:
                face = list(seq)
                face[i-1] = merged
                del face[i]
                ft = tuple(face)
                C[junk_faces[ft], j] += (-1) ** i
    return C

def build_boundary(d, Ad_from, Ad_to):
    to_idx = {seq: i for i, seq in enumerate(Ad_to)}
    B = np.zeros((len(Ad_to), len(Ad_from)), dtype=np.int8)
    for j, seq in enumerate(Ad_from):
        for i in range(d + 1):
            if i == 0:
                face = seq[1:]
            elif i == d:
                face = seq[:-1]
            else:
                face = seq[:i-1] + ((seq[i-1] + seq[i]) % p,) + seq[i+1:]
            if face in to_idx:
                B[to_idx[face], j] += (-1) ** i
    return B

# Compute R_d for d = 1..7
for d in range(1, 8):
    print(f"\n  R_{d}:", flush=True)
    C_d = build_constraint(d, diffseqs[d])
    B_d = build_boundary(d, diffseqs[d], diffseqs[d-1])

    stacked = np.vstack([C_d, B_d]) if C_d.shape[0] > 0 else B_d
    print(f"    stacked matrix: {stacked.shape}", flush=True)

    rank_C = mod_rank(C_d) if C_d.shape[0] > 0 else 0
    rank_stacked = mod_rank(stacked)
    R_d = rank_stacked - rank_C
    print(f"    rank([C;B]) = {rank_stacked}, rank(C) = {rank_C}, R_{d} = {R_d}")
    print(f"    predicted: {R_pred[d]}  {'✓' if R_d == R_pred[d] else '✗'}")
