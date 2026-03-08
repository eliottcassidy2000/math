#!/usr/bin/env python3
"""
beta2_surplus_formula.py - Search for a formula for surplus in terms of tournament invariants

surplus = dim(Omega_3) - dim(Z_2)
       = dim(Omega_3) - [dim(Omega_2) - rk(d_2|Omega_2)]

Known: surplus >= 0 for all tournaments (equivalent to beta_2 = 0).

Can surplus be expressed in terms of:
- Score sequence (d_0, ..., d_{n-1})
- 3-cycle count t_3
- TT triple count (= dim of TT subspace of Omega_2)
- Number of DT 4-paths
- Edge counts between vertex pairs

Let's tabulate surplus against all available invariants.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_3cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                c = A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i]
                count += c
    return count

def count_TT(A, n):
    count = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                if A[a][c]:
                    count += 1
    return count

def count_DT(A, n):
    count = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                if not A[a][c]: continue
                for d in range(n):
                    if d == a or d == b or d == c or not A[c][d]: continue
                    if A[b][d]:
                        count += 1
    return count

def compute_surplus(A, n):
    allowed = {}
    for p in range(5):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break

    omega_basis = {}
    for p in range(4):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p >= 1 and p-1 in allowed else [])
        omega_basis[p] = basis

    dim_O2 = omega_basis[2].shape[1] if omega_basis[2].ndim == 2 else 0
    dim_O3 = omega_basis[3].shape[1] if omega_basis[3].ndim == 2 else 0

    if dim_O2 == 0:
        Z2 = 0
        rk2 = 0
    else:
        bd2 = build_full_boundary_matrix(allowed[2], allowed.get(1, []))
        bd2_om = bd2 @ omega_basis[2]
        S_v = np.linalg.svd(bd2_om, compute_uv=False)
        rk2 = int(np.sum(np.abs(S_v) > 1e-8))
        Z2 = dim_O2 - rk2

    surplus = dim_O3 - Z2
    A2 = len(allowed.get(2, []))
    A3 = len(allowed.get(3, []))
    return surplus, dim_O2, dim_O3, Z2, rk2, A2, A3


n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"SURPLUS FORMULA SEARCH AT n={n}")
print("=" * 70)

data = []
for bits in range(total):
    A = build_adj(n, bits)
    surplus, O2, O3, Z2, rk2, A2, A3 = compute_surplus(A, n)
    t3 = count_3cycles(A, n)
    tt = count_TT(A, n)
    dt = count_DT(A, n)
    scores = tuple(sorted(sum(row) for row in A))
    score_var = sum((s - (n-1)/2)**2 for s in scores)

    data.append({
        'bits': bits, 'surplus': surplus, 'O2': O2, 'O3': O3,
        'Z2': Z2, 'rk2': rk2, 'A2': A2, 'A3': A3,
        't3': t3, 'TT': tt, 'DT': dt, 'scores': scores,
        'score_var': score_var
    })

# Check: is surplus determined by (t3)?
print(f"\n  Surplus by t3:")
by_t3 = defaultdict(set)
for d in data:
    by_t3[d['t3']].add(d['surplus'])
for t3 in sorted(by_t3.keys()):
    print(f"    t3={t3}: surplus in {sorted(by_t3[t3])}")

determined_by_t3 = all(len(v) == 1 for v in by_t3.values())
print(f"  Determined by t3: {determined_by_t3}")

# Check: is surplus determined by (TT)?
print(f"\n  Surplus by TT:")
by_TT = defaultdict(set)
for d in data:
    by_TT[d['TT']].add(d['surplus'])
for tt in sorted(by_TT.keys()):
    print(f"    TT={tt}: surplus in {sorted(by_TT[tt])}")

# Check: is surplus determined by (DT)?
print(f"\n  Surplus by DT:")
by_DT = defaultdict(set)
for d in data:
    by_DT[d['DT']].add(d['surplus'])
for dt in sorted(by_DT.keys()):
    print(f"    DT={dt}: surplus in {sorted(by_DT[dt])}")

# Check: is surplus = DT - TT + something?
print(f"\n  Surplus vs DT - TT:")
residual = Counter(d['surplus'] - (d['DT'] - d['TT']) for d in data)
print(f"    surplus - (DT - TT): {dict(sorted(residual.items()))}")

# Check: is surplus determined by score sequence?
print(f"\n  Surplus by score sequence:")
by_score = defaultdict(set)
for d in data:
    by_score[d['scores']].add(d['surplus'])
for scores in sorted(by_score.keys()):
    vals = sorted(by_score[scores])
    if len(vals) > 1:
        print(f"    {scores}: surplus in {vals}")

determined_by_score = all(len(v) == 1 for v in by_score.values())
print(f"  Determined by score sequence: {determined_by_score}")

# Try formulas
print(f"\n  Testing formulas:")

# Formula 1: surplus = DT - Z_2 = DT - (O2 - rk2)
# No, surplus = O3 - Z2, and O3 != DT always.

# Formula 2: try surplus = A3 - A2 + something
for d in data[:5]:
    print(f"    A2={d['A2']}, A3={d['A3']}, O2={d['O2']}, O3={d['O3']}, Z2={d['Z2']}, surplus={d['surplus']}, t3={d['t3']}, TT={d['TT']}, DT={d['DT']}")

# Check: is O2 = TT + f(t3)?
print(f"\n  O2 - TT by t3:")
by_t3_O2 = defaultdict(set)
for d in data:
    by_t3_O2[d['t3']].add(d['O2'] - d['TT'])
for t3 in sorted(by_t3_O2.keys()):
    print(f"    t3={t3}: O2-TT in {sorted(by_t3_O2[t3])}")

# Check: is O3 = DT + f(t3)?
print(f"\n  O3 - DT by t3:")
by_t3_O3 = defaultdict(set)
for d in data:
    by_t3_O3[d['t3']].add(d['O3'] - d['DT'])
for t3 in sorted(by_t3_O3.keys()):
    print(f"    t3={t3}: O3-DT in {sorted(by_t3_O3[t3])}")

# Key: what about O2 - TT (= number of "NT cancellation" elements in Omega_2)?
print(f"\n  O2 - TT distribution:")
nt_dist = Counter(d['O2'] - d['TT'] for d in data)
for k in sorted(nt_dist.keys()):
    print(f"    O2-TT={k}: {nt_dist[k]}")

# And O3 - DT?
print(f"\n  O3 - DT distribution:")
ndt_dist = Counter(d['O3'] - d['DT'] for d in data)
for k in sorted(ndt_dist.keys()):
    print(f"    O3-DT={k}: {ndt_dist[k]}")

print("\nDone.")
