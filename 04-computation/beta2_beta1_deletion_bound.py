#!/usr/bin/env python3
"""
beta2_beta1_deletion_bound.py — Why β₁(T)=0 implies Σ_v β₁(T\v) ≤ 3

KEY DISCOVERY: The proof of β₂=0 reduces to showing:
  β₁(T) = 0 ⟹ at most 3 vertex-deletions have β₁(T\v) > 0

At n=4: β₁ = 1 iff t₃ = 2. So β₁(T\v) = 1 iff t₃(T\v) = 2.
t₃(T\v) = t₃(T) - c_v where c_v = #{3-cycles through v}.

Study:
1. Distribution of (t₃, c_v, β₁(T), β₁(T\v)) at n=5,6
2. Why does β₁(T)=0 constrain the count of v with c_v = t₃-2?
3. Is there a purely combinatorial proof?

Also verify: β₁(T)=0 ⟺ Σβ₁(T\v) ≤ 3 at n=7 (sampling)

Author: opus-2026-03-08-S49
"""
import sys, time, random
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def build_random_adj(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def count_3cycles(A, n):
    """Count 3-cycles in tournament."""
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check if {i,j,k} forms a 3-cycle
                edges = (A[i][j], A[j][k], A[i][k],
                         A[j][i], A[k][j], A[k][i])
                # 3-cycle: i→j→k→i or i→k→j→i
                if A[i][j] and A[j][k] and A[k][i]:
                    t3 += 1
                elif A[i][k] and A[k][j] and A[j][i]:
                    t3 += 1
    return t3

def compute_beta1(A, n):
    """Compute β₁ of tournament."""
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

    d1 = dim_om(om1)
    if d1 == 0:
        return 0

    bd1 = build_full_boundary_matrix(ap1, ap0)
    rk1 = np.linalg.matrix_rank(bd1 @ om1, tol=1e-8)
    z1 = d1 - rk1

    if dim_om(om2) > 0:
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2om = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
        b1 = np.linalg.matrix_rank(bd2om, tol=1e-8)
    else:
        b1 = 0

    return z1 - b1


# ============================================================
# n=5: Full analysis
# ============================================================
print("=" * 70)
print("n=5: β₁(T) vs Σβ₁(T\\v) — 3-CYCLE MECHANISM")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

joint = Counter()  # (beta1_T, sigma_beta1_del, t3) -> count
cycle_details = defaultdict(list)

for bits in range(total):
    A = build_adj(n, bits)
    t3 = count_3cycles(A, n)
    beta1_T = compute_beta1(A, n)

    beta1_del = []
    cv_list = []
    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
        b1v = compute_beta1(A_sub, n-1)
        beta1_del.append(b1v)
        # Count 3-cycles through v
        cv = sum(1 for i in range(n) for j in range(n)
                 if i != v and j != v and i != j
                 and A[v][i] and A[i][j] and A[j][v])
        cv_list.append(cv)

    sigma = sum(beta1_del)
    joint[(beta1_T, sigma, t3)] += 1

    # Detailed analysis for sigma=3, beta1=0
    if sigma == 3 and beta1_T == 0 and len(cycle_details[(3,0)]) < 3:
        scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
        critical = [v for v in range(n) if beta1_del[v] == 1]
        cycle_details[(3,0)].append({
            'bits': bits, 't3': t3, 'scores': scores,
            'cv': cv_list, 'critical': critical,
            'beta1_del': beta1_del
        })

print(f"\nJoint (β₁(T), Σβ₁(T\\v), t₃):")
for (b1, sb, t3), cnt in sorted(joint.items()):
    print(f"  β₁={b1}, Σβ₁_del={sb}, t₃={t3}: {cnt}")

print(f"\nExamples with Σβ₁_del=3, β₁(T)=0:")
for info in cycle_details[(3,0)]:
    print(f"  bits={info['bits']}: t₃={info['t3']}, scores={info['scores']}")
    print(f"    c_v={info['cv']}, critical={info['critical']}, β₁_del={info['beta1_del']}")
    print(f"    t₃(T\\v) = {[info['t3']-info['cv'][v] for v in range(n)]}")

# ============================================================
# KEY TEST: At n=4, β₁=1 iff t₃=2
# So β₁(T\v) = 1 iff c_v = t₃-2
# Count of critical v = #{v: c_v = t₃-2}
# ============================================================
print(f"\n{'='*70}")
print("KEY: #{v: c_v = t₃-2} when β₁(T)=0")
print("=" * 70)

# For each (t3, beta1_T): what's the distribution of #{v: cv = t3-2}?
for t3_target in range(8):
    for b1 in [0, 1]:
        count = 0
        ncrit_dist = Counter()
        for bits in range(total):
            A = build_adj(n, bits)
            t3 = count_3cycles(A, n)
            if t3 != t3_target:
                continue
            beta1_T = compute_beta1(A, n)
            if beta1_T != b1:
                continue
            count += 1

            # Count critical vertices
            ncrit = 0
            for v in range(n):
                cv = sum(1 for i in range(n) for j in range(n)
                         if i != v and j != v and i != j
                         and A[v][i] and A[i][j] and A[j][v])
                if cv == t3 - 2:
                    ncrit += 1
            ncrit_dist[ncrit] += 1

        if count > 0:
            print(f"  t₃={t3_target}, β₁={b1}: {count} tournaments, #{'{v:cv=t₃-2}'} dist = {dict(sorted(ncrit_dist.items()))}")


# ============================================================
# n=6,7: Verify β₁(T)=0 ⟺ Σβ₁(T\v) ≤ 3
# ============================================================
print(f"\n{'='*70}")
print("n=6: EXHAUSTIVE CHECK")
print("=" * 70)

n = 6
m = n*(n-1)//2
total = 1 << m
t0 = time.time()

iff_holds = True
joint6 = Counter()  # (beta1_T, sigma_beta1_del) -> count
max_sigma_when_b1_0 = 0

for bits in range(total):
    A = build_adj(n, bits)
    beta1_T = compute_beta1(A, n)

    sigma = 0
    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
        b1v = compute_beta1(A_sub, n-1)
        sigma += b1v

    joint6[(beta1_T, sigma)] += 1

    if beta1_T == 0:
        max_sigma_when_b1_0 = max(max_sigma_when_b1_0, sigma)
        if sigma > 3:
            iff_holds = False
            scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
            print(f"  FAIL: bits={bits}, β₁=0, Σβ₁_del={sigma}, scores={scores}")
    elif beta1_T > 0 and sigma <= 3:
        iff_holds = False
        scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
        print(f"  FAIL (reverse): bits={bits}, β₁={beta1_T}, Σβ₁_del={sigma}, scores={scores}")

    if (bits+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total} ({elapsed:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=6: {total} tournaments in {elapsed:.0f}s")
print(f"  Max Σβ₁_del when β₁=0: {max_sigma_when_b1_0}")
print(f"  β₁=0 ⟺ Σβ₁_del≤3: {'✓ HOLDS' if iff_holds else '✗ FAILS'}")

print(f"\n  Joint (β₁(T), Σβ₁(T\\v)):")
for (b1, sb), cnt in sorted(joint6.items()):
    print(f"    β₁={b1}, Σβ₁_del={sb}: {cnt}")


# ============================================================
# n=7: Sampled check
# ============================================================
print(f"\n{'='*70}")
print("n=7: SAMPLED CHECK (2000 samples)")
print("=" * 70)

n = 7
samples = 2000
t0 = time.time()
joint7 = Counter()
max_sigma_b1_0 = 0
counterexample = False

for trial in range(samples):
    A = build_random_adj(n)
    beta1_T = compute_beta1(A, n)

    sigma = 0
    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
        b1v = compute_beta1(A_sub, n-1)
        sigma += b1v

    joint7[(beta1_T, sigma)] += 1

    if beta1_T == 0 and sigma > 3:
        counterexample = True
        scores = sorted([sum(A[i][j] for j in range(n) if j!=i) for i in range(n)])
        print(f"  COUNTEREXAMPLE: β₁=0, Σβ₁_del={sigma}, scores={scores}")
    if beta1_T > 0 and sigma <= 3:
        counterexample = True
        scores = sorted([sum(A[i][j] for j in range(n) if j!=i) for i in range(n)])
        print(f"  COUNTEREXAMPLE (rev): β₁={beta1_T}, Σβ₁_del={sigma}, scores={scores}")

    if beta1_T == 0:
        max_sigma_b1_0 = max(max_sigma_b1_0, sigma)

    if (trial+1) % 500 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{samples} ({elapsed:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=7: {samples} tournaments in {elapsed:.0f}s")
print(f"  Max Σβ₁_del when β₁=0: {max_sigma_b1_0}")
print(f"  Counterexample found: {counterexample}")
print(f"\n  Joint (β₁(T), Σβ₁(T\\v)):")
for (b1, sb), cnt in sorted(joint7.items()):
    print(f"    β₁={b1}, Σβ₁_del={sb}: {cnt}")

print("\nDone.")
