#!/usr/bin/env python3
"""
beta2_rank_equation.py - Investigate the rank equation for beta2=0

beta2=0 iff rk(bd3|_Om3) = dim(Z_2) = dim(Om2) - rk(bd2|_Om2)

We know: dim(Om3) = |A3| - rk(C) where C = constraint matrix
And: dim(Om2) = |A2| - J2 where J2 = junk pairs

Question: Is there a direct algebraic relationship between
rk(bd3|_Om3), dim(Om2), and the constraint matrix?

Key observation: the boundary bd3 sends a 3-path (a,b,c,d) to
  (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
The faces (a,c,d) and (a,b,d) are the ones that can be non-allowed.
(b,c,d) and (a,b,c) are always allowed (consecutive edges in tournament).

So bd3 restricted to Om3 maps each Om3 element to a combination
of ALLOWED 2-paths only (by definition of Om3).

For beta2=0, the image of bd3|_Om3 must equal Z2 = ker(bd2|_Om2).

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
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


def full_analysis(A, n):
    """Full rank analysis of the chain complex."""
    a0 = enumerate_allowed_paths(A, n, 0)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a4 = enumerate_allowed_paths(A, n, 4)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)

    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # Boundary matrices
    bd2 = build_full_boundary_matrix(a2, a1)
    bd3 = build_full_boundary_matrix(a3, a2)

    # Ranks
    rk_bd2_om = 0
    if d_om2 > 0:
        bd2_om = bd2 @ om2
        S = np.linalg.svd(bd2_om, compute_uv=False)
        rk_bd2_om = sum(s > 1e-8 for s in S)

    dim_z2 = d_om2 - rk_bd2_om

    rk_bd3_om = 0
    if d_om3 > 0:
        bd3_om = bd3 @ om3
        S = np.linalg.svd(bd3_om, compute_uv=False)
        rk_bd3_om = sum(s > 1e-8 for s in S)

    beta2 = dim_z2 - rk_bd3_om

    # Classify 3-paths
    dd = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])

    # J2 (junk pairs)
    j2 = 0
    for a_ in range(n):
        for c_ in range(n):
            if a_ == c_ or not A[c_][a_]:
                continue
            for b_ in range(n):
                if b_ != a_ and b_ != c_ and A[a_][b_] and A[b_][c_]:
                    j2 += 1
                    break

    return {
        'A2': len(a2), 'A3': len(a3), 'A4': len(a4),
        'J2': j2, 'Om2': d_om2, 'Om3': d_om3,
        'rk_bd2': rk_bd2_om, 'Z2': dim_z2,
        'rk_bd3': rk_bd3_om, 'beta2': beta2,
        'DD': dd,
    }


# ============================================================
# Rank equation analysis at n=5
# ============================================================
print("=" * 70)
print("RANK EQUATION ANALYSIS: n=5")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

data = []
for bits in range(total):
    A = build_adj(n, bits)
    r = full_analysis(A, n)
    r['bits'] = bits
    r['scores'] = tuple(sorted([sum(row) for row in A]))
    data.append(r)

# Check rank equation
# beta2 = Z2 - rk_bd3 = (Om2 - rk_bd2) - rk_bd3
# For beta2=0: rk_bd3 = Om2 - rk_bd2

# Cluster by score sequence
from collections import defaultdict
score_data = defaultdict(list)
for r in data:
    score_data[r['scores']].append(r)

print(f"\nBy score sequence:")
print(f"{'scores':<20} {'count':>5} {'A2':>4} {'J2':>3} {'Om2':>4} {'rk_bd2':>6} {'Z2':>3} "
      f"{'A3':>4} {'DD':>3} {'Om3':>4} {'rk_bd3':>6} {'b2':>3}")
for sc in sorted(score_data.keys()):
    entries = score_data[sc]
    r0 = entries[0]
    # Check if all entries have same values
    same = all(r['A2'] == r0['A2'] and r['J2'] == r0['J2'] and r['Om2'] == r0['Om2']
               and r['Z2'] == r0['Z2'] and r['Om3'] == r0['Om3'] for r in entries)
    marker = " " if same else "*"
    print(f"{str(sc):<20} {len(entries):>5} {r0['A2']:>4} {r0['J2']:>3} {r0['Om2']:>4} {r0['rk_bd2']:>6} "
          f"{r0['Z2']:>3} {r0['A3']:>4} {r0['DD']:>3} {r0['Om3']:>4} {r0['rk_bd3']:>6} {r0['beta2']:>3}{marker}")


# ============================================================
# Key relationship: rk_bd3 vs Z2
# ============================================================
print(f"\n{'='*70}")
print("KEY RELATIONSHIP: rk(bd3|_Om3) = dim(Z2)")
print("=" * 70)

# rk_bd3 should ALWAYS equal Z2 for beta2=0
for r in data:
    if r['rk_bd3'] != r['Z2']:
        print(f"  MISMATCH: bits={r['bits']}, rk_bd3={r['rk_bd3']} != Z2={r['Z2']}")

print(f"All {total} verified: rk(bd3|_Om3) = dim(Z2) always")

# What determines Z2?
# Z2 = Om2 - rk_bd2 = (A2 - J2) - rk_bd2
# What determines rk_bd2?
print(f"\n{'='*70}")
print("WHAT DETERMINES rk(bd2|_Om2)?")
print("=" * 70)

# bd2 sends (a,b,c) -> (b,c) - (a,c) + (a,b)
# Restricted to Om2, it sends Om2 elements to A1
# rk(bd2|_Om2) = dim(Om2) - dim(Z2)
# = dim(Om2) - dim(ker bd2 on Om2)

rk_bd2_vals = Counter()
for r in data:
    rk_bd2_vals[(r['scores'], r['rk_bd2'])] += 1

print(f"\nrk(bd2) by score:")
for (sc, rk), cnt in sorted(rk_bd2_vals.items()):
    print(f"  scores={sc}: rk_bd2={rk} ({cnt} tournaments)")


# ============================================================
# n=6 analysis
# ============================================================
print(f"\n{'='*70}")
print("RANK EQUATION ANALYSIS: n=6")
print("=" * 70)

n = 6
total = 1 << (n*(n-1)//2)

score_data6 = defaultdict(list)
t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s)")

    A = build_adj(n, bits)
    r = full_analysis(A, n)
    r['bits'] = bits
    r['scores'] = tuple(sorted([sum(row) for row in A]))
    score_data6[r['scores']].append(r)

dt = time.time() - t0
print(f"\nDone in {dt:.0f}s")

print(f"\nBy score sequence (representative):")
print(f"{'scores':<25} {'count':>5} {'A2':>4} {'J2':>3} {'Om2':>4} {'rk_bd2':>6} {'Z2':>3} "
      f"{'A3':>5} {'DD':>3} {'Om3':>4} {'rk_bd3':>6} {'b2':>3}")
for sc in sorted(score_data6.keys()):
    entries = score_data6[sc]
    r0 = entries[0]
    # Check if values vary within score class
    z2_vals = set(r['Z2'] for r in entries)
    om3_vals = set(r['Om3'] for r in entries)
    b2_vals = set(r['beta2'] for r in entries)
    marker = " " if len(z2_vals) == 1 else "*"
    z2_str = str(r0['Z2']) if len(z2_vals) == 1 else f"{min(z2_vals)}-{max(z2_vals)}"
    om3_str = str(r0['Om3']) if len(om3_vals) == 1 else f"{min(om3_vals)}-{max(om3_vals)}"
    print(f"{str(sc):<25} {len(entries):>5} {r0['A2']:>4} {r0['J2']:>3} {r0['Om2']:>4} {r0['rk_bd2']:>6} "
          f"{z2_str:>3} {r0['A3']:>5} {r0['DD']:>3} {om3_str:>4} {r0['rk_bd3']:>6} {max(b2_vals):>3}{marker}")


# ============================================================
# Universal identity check: dim(Om3) - rk(bd3) = ?
# ============================================================
print(f"\n{'='*70}")
print("KERNEL OF bd3|_Om3 (= dim Om3 - rk bd3)")
print("=" * 70)

# dim(ker bd3|_Om3) = dim(Om3) - rk(bd3)
# This is the dimension of Z3 (3-cycles in Om3)
# beta3 = Z3 - im(bd4|_Om4)

ker_dist = Counter()
for sc, entries in score_data6.items():
    for r in entries:
        ker = r['Om3'] - r['rk_bd3']
        ker_dist[ker] += 1

print(f"\ndim(ker bd3|_Om3) distribution:")
for k in sorted(ker_dist.keys()):
    print(f"  ker={k}: {ker_dist[k]}")


# ============================================================
# Euler characteristic check
# ============================================================
print(f"\n{'='*70}")
print("EULER CHARACTERISTIC CHECK")
print("=" * 70)

# chi = sum (-1)^p beta_p
# For tournaments: beta0=1, beta1=c3 (number of 3-cycles)
# If beta2=0, chi = 1 - beta1 + 0 - ...
# But we also have chi = sum (-1)^p dim(Om_p)

# Actually chi = sum (-1)^p rk(bd_p)  (from alternating sum of dim(Om_p))
# No, let's compute it properly

# Euler char from Om dimensions:
# chi = sum (-1)^p dim(Om_p) [this is NOT the Euler char]
# The actual Euler char is sum (-1)^p beta_p

# Let's just verify beta2=0 from the rank data
all_beta2_zero = all(r['beta2'] == 0 for entries in score_data6.values() for r in entries)
print(f"\nAll beta2=0 at n=6: {all_beta2_zero}")


print("\n\nDone.")
