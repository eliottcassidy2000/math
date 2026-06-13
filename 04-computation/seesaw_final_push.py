#!/usr/bin/env python3
"""
SEESAW: FINAL PUSH for the proof.

Strategy: understand WHAT structural property separates β₁=1 from β₃>0.

The key: β₁ = 1 iff rank(∂₂|Ω₂) = C(n,2) - n (one below maximum).
This means exactly ONE linear dependency among the transitive triple constraints.

β₃ > 0 iff some 3-chains in Ω₃ have zero boundary but don't come from Ω₄.

HYPOTHESIS: The "one missing constraint" that creates β₁=1 is INCOMPATIBLE
with the structure needed for β₃>0. Specifically:

When rank(∂₂) drops by 1 (β₁=1), the chain complex "tightens" at dimension 2,
which FORCES ∂₃ to have higher rank (fills more of Ω₂), which leaves LESS
room for ker(∂₃) (which contains β₃), pushing β₃ toward 0.

This is the SEESAW MECHANISM at the algebraic level.

Let me verify: when β₁=1, does rank(∂₃) increase compared to β₁=0 tournaments?

kind-pasteur-2026-03-25-S28
"""
import sys, time, random
import numpy as np
from itertools import permutations
from collections import Counter, defaultdict

sys.stdout.reconfigure(line_buffering=True)
random.seed(20260325)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def enum_allowed(A, n, p):
    paths = []
    for perm in permutations(range(n), p+1):
        ok = True
        for i in range(p):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False; break
        if ok: paths.append(perm)
    return paths

def compute_omega_rank(A, n, p):
    """Return (dim_Omega_p, rank_d_p_on_Omega_p, Omega_basis)."""
    ap = enum_allowed(A, n, p)
    if not ap: return 0, 0, None
    if p == 0: return len(ap), 0, np.eye(len(ap))

    ap_m1 = enum_allowed(A, n, p-1)
    ap_m1_set = set(ap_m1)

    na = {}; ni = 0
    for path in ap:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in ap_m1_set and face not in na:
                na[face] = ni; ni += 1

    if not na:
        om_dim = len(ap); om_basis = np.eye(len(ap))
    else:
        P = np.zeros((len(na), len(ap)))
        for j, path in enumerate(ap):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in na: P[na[face], j] += sign
        _, S, Vt = np.linalg.svd(P)
        rk = sum(s > 1e-8 for s in S)
        om_dim = len(ap) - rk
        om_basis = Vt[rk:].T if om_dim > 0 else np.zeros((len(ap), 0))

    # Boundary rank
    if ap_m1 and om_dim > 0:
        idx = {path: i for i, path in enumerate(ap_m1)}
        D = np.zeros((len(ap_m1), len(ap)))
        for j, path in enumerate(ap):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in idx: D[idx[face], j] += sign
        Do = D @ om_basis
        _, S2, _ = np.linalg.svd(Do)
        rk_d = sum(s > 1e-8 for s in S2)
    else:
        rk_d = 0

    return om_dim, rk_d, om_basis

# ============================================================
# COLLECT STATISTICS: β₁ vs rank(∂₃) at n=7
# ============================================================

n = 7
print("=" * 80)
print(f"  SEESAW FINAL PUSH: Chain complex statistics at n={n}")
print("=" * 80)

data = defaultdict(list)
n_samples = 300

t0 = time.time()
for trial in range(n_samples):
    A = random_tournament(n)

    o1, r1, _ = compute_omega_rank(A, n, 1)  # Ω₁, rank(∂₁)
    o2, r2, _ = compute_omega_rank(A, n, 2)  # Ω₂, rank(∂₂)
    o3, r3, _ = compute_omega_rank(A, n, 3)  # Ω₃, rank(∂₃)
    o4, r4, _ = compute_omega_rank(A, n, 4)  # Ω₄, rank(∂₄)

    ker_d1 = n*(n-1)//2 - (n-1)  # always 15
    b1 = max(0, ker_d1 - r2)

    # β₂ = dim(Ω₂) - r2 - r3 ... wait, need to be more careful
    # β₂ = ker(∂₂|Ω₂) - im(∂₃|Ω₃→Ω₂)
    # ker(∂₂) = Ω₂ - rank(∂₂) = o2 - r2
    # im(∂₃ into Ω₂) needs projection... approximate with r3
    b2 = max(0, (o2 - r2) - r3)

    # β₃ = ker(∂₃) - im(∂₄)
    # ker(∂₃) = o3 - r3
    # im(∂₄ into Ω₃) ≈ r4
    b3 = max(0, (o3 - r3) - r4)

    profile = (b1, b2, b3)
    data[profile].append({
        'o2': o2, 'r2': r2, 'o3': o3, 'r3': r3, 'o4': o4, 'r4': r4,
        'ker_d3': o3 - r3, 'im_d4': r4,
    })

    if (trial+1) % 100 == 0:
        print(f"  [{trial+1}/{n_samples}] {time.time()-t0:.1f}s")

print(f"\n  Profile distribution:")
for prof in sorted(data.keys()):
    count = len(data[prof])
    print(f"  (β₁,β₂,β₃) = {prof}: {count} ({count/n_samples*100:.1f}%)")

# ============================================================
# COMPARE CHAIN COMPLEX STRUCTURE: β₁=1 vs β₃>0
# ============================================================

print(f"\n  === CHAIN COMPLEX COMPARISON ===")
print(f"  {'Statistic':<15} {'β₁=0,β₃=0':>12} {'β₁=1,β₃=0':>12} {'β₁=0,β₃=1':>12}")
print("  " + "-" * 55)

groups = {
    '(0,0,0)': data.get((0,0,0), []),
    '(1,0,0)': data.get((1,0,0), []),
    '(0,0,1)': data.get((0,0,1), []),
}

for stat in ['o2', 'r2', 'o3', 'r3', 'o4', 'r4', 'ker_d3', 'im_d4']:
    vals = {}
    for gname, gdata in groups.items():
        if gdata:
            vals[gname] = np.mean([d[stat] for d in gdata])
        else:
            vals[gname] = 0
    print(f"  {stat:<15} {vals.get('(0,0,0)', 0):>12.1f} {vals.get('(1,0,0)', 0):>12.1f} {vals.get('(0,0,1)', 0):>12.1f}")

# THE KEY COMPARISON
print(f"\n  === THE SEESAW MECHANISM ===")
for gname, gdata in groups.items():
    if not gdata: continue
    mean_o2 = np.mean([d['o2'] for d in gdata])
    mean_r2 = np.mean([d['r2'] for d in gdata])
    mean_r3 = np.mean([d['r3'] for d in gdata])
    mean_ker3 = np.mean([d['ker_d3'] for d in gdata])
    mean_r4 = np.mean([d['im_d4'] for d in gdata])

    # β₂ = 0 means: ker(∂₂) = im(∂₃)
    # So: o2 - r2 = r3 (exactness at Ω₂)
    exactness_gap = mean_o2 - mean_r2 - mean_r3

    print(f"\n  {gname}:")
    print(f"    dim(Ω₂) = {mean_o2:.1f}, rank(∂₂) = {mean_r2:.1f}")
    print(f"    ker(∂₂) = dim(Ω₂) - rank(∂₂) = {mean_o2-mean_r2:.1f}")
    print(f"    rank(∂₃) = {mean_r3:.1f} (should equal ker(∂₂) by β₂=0)")
    print(f"    Exactness gap: {exactness_gap:.1f} (should be ~0)")
    print(f"    ker(∂₃) = {mean_ker3:.1f}, im(∂₄) = {mean_r4:.1f}")
    print(f"    β₃ ≈ ker(∂₃) - im(∂₄) = {mean_ker3-mean_r4:.1f}")

print(f"""
  THE SEESAW IN NUMBERS:
  When β₁ goes from 0 to 1:
    rank(∂₂) drops by 1 (one fewer constraint)
    → ker(∂₂) increases by 1
    → by β₂=0 exactness: rank(∂₃) must INCREASE by 1
    → ker(∂₃) = Ω₃ - rank(∂₃) DECREASES by 1
    → β₃ = ker(∂₃) - im(∂₄) DECREASES by 1

  IF ker(∂₃) starts at exactly im(∂₄) + 1 when β₃=1,
  THEN decreasing ker(∂₃) by 1 gives β₃=0. SEESAW!

  This requires: β₃ = 1 is the MAXIMUM possible (no β₃≥2 at n=7).
  And when β₁=1, the rank increase in ∂₃ is EXACTLY 1.

  CHECK: Is rank(∂₃) exactly 1 higher in β₁=1 tournaments vs β₁=0?
""")

# Verify: rank(∂₃) difference between β₁=0 and β₁=1
if groups['(0,0,0)'] and groups['(1,0,0)']:
    r3_b1_0 = [d['r3'] for d in groups['(0,0,0)']]
    r3_b1_1 = [d['r3'] for d in groups['(1,0,0)']]
    print(f"  rank(∂₃) when β₁=0: mean={np.mean(r3_b1_0):.2f}, values={Counter(r3_b1_0).most_common(5)}")
    print(f"  rank(∂₃) when β₁=1: mean={np.mean(r3_b1_1):.2f}, values={Counter(r3_b1_1).most_common(5)}")
    print(f"  Difference: {np.mean(r3_b1_1) - np.mean(r3_b1_0):.2f}")

# Also check: does ker(∂₃) = im(∂₄) exactly when β₃=0?
print(f"\n  ker(∂₃) - im(∂₄) when β₃=0:")
for gname in ['(0,0,0)', '(1,0,0)']:
    if groups[gname]:
        gaps = [d['ker_d3'] - d['im_d4'] for d in groups[gname]]
        print(f"    {gname}: {Counter(gaps).most_common(5)}")

print(f"\n  ker(∂₃) - im(∂₄) when β₃=1:")
if groups.get('(0,0,1)'):
    gaps = [d['ker_d3'] - d['im_d4'] for d in groups['(0,0,1)']]
    print(f"    (0,0,1): {Counter(gaps).most_common(5)}")
