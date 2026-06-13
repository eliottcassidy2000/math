"""
Test beta_1 * beta_3 = 0 seesaw at n=8.
THM-095 proves this for n<=7. Does it hold at n=8?

The seesaw mechanism: beta_2=0 forces im(d_2) to take only 2 values,
which mediates exclusion between beta_1 and beta_3.

At n=8, beta_2=0 still holds (HYP confirmed). So the seesaw mechanism
should still apply. But we need to verify computationally.

kind-pasteur-2026-03-25-S1
"""
import sys
import time
import random
import numpy as np
from collections import Counter
from itertools import permutations

sys.stdout.reconfigure(line_buffering=True)
random.seed(2026)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def enum_allowed(A, n, p):
    paths = []
    for perm in permutations(range(n), p+1):
        ok = True
        for i in range(p):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            paths.append(perm)
    return paths

def compute_omega_and_boundary(A, n, p, ap, ap_m1):
    """Compute Omega_p basis and boundary map d_p: Omega_p -> A_{p-1}."""
    if not ap:
        return np.zeros((0, 0)), np.zeros((len(ap_m1) if ap_m1 else 0, 0))

    ap_m1_set = set(ap_m1) if ap_m1 else set()

    # Find non-allowed faces
    non_allowed = {}
    na_idx = 0
    for path in ap:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in ap_m1_set and face not in non_allowed:
                non_allowed[face] = na_idx
                na_idx += 1

    if not non_allowed:
        omega_basis = np.eye(len(ap))
    else:
        P = np.zeros((len(non_allowed), len(ap)))
        for j, path in enumerate(ap):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in non_allowed:
                    P[non_allowed[face], j] += sign

        _, S, Vt = np.linalg.svd(P)
        rank_P = sum(s > 1e-8 for s in S)
        if rank_P >= len(ap):
            return np.zeros((len(ap), 0)), np.zeros((len(ap_m1) if ap_m1 else 0, 0))
        omega_basis = Vt[rank_P:].T

    # Build boundary d_p: A_p -> A_{p-1}
    if p >= 1 and ap_m1:
        idx = {path: i for i, path in enumerate(ap_m1)}
        D = np.zeros((len(ap_m1), len(ap)))
        for j, path in enumerate(ap):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in idx:
                    D[idx[face], j] += sign
        D_omega = D @ omega_basis
    else:
        D_omega = np.zeros((0, omega_basis.shape[1]))

    return omega_basis, D_omega

def compute_betti_full(A, n, max_dim=None):
    """Compute all Betti numbers up to max_dim."""
    if max_dim is None:
        max_dim = n - 1

    # Enumerate all allowed paths
    allowed = {}
    for p in range(max_dim + 2):
        allowed[p] = enum_allowed(A, n, p)

    betti = []
    for p in range(max_dim + 1):
        omega_p, D_p = compute_omega_and_boundary(A, n, p, allowed[p], allowed.get(p-1, []))

        if omega_p.shape[1] == 0:
            betti.append(0)
            continue

        # ker(d_p)
        if D_p.shape[0] > 0 and D_p.shape[1] > 0:
            _, S_dp, _ = np.linalg.svd(D_p)
            rank_dp = sum(s > 1e-8 for s in S_dp)
        else:
            rank_dp = 0
        ker_dp = omega_p.shape[1] - rank_dp

        # im(d_{p+1})
        omega_p1, D_p1 = compute_omega_and_boundary(A, n, p+1, allowed.get(p+1, []), allowed[p])
        if omega_p1.shape[1] == 0 or D_p1.shape[1] == 0:
            im_dp1 = 0
        else:
            # Project D_p1 into Omega_p coordinates
            im_matrix = omega_p.T @ D_p1
            _, S_im, _ = np.linalg.svd(im_matrix)
            im_dp1 = sum(s > 1e-8 for s in S_im)

        betti.append(max(0, ker_dp - im_dp1))

    return betti

# ===== MAIN =====
print("=" * 60)
print("SEESAW TEST: beta_1 * beta_3 = 0 at n=8")
print("kind-pasteur-2026-03-25-S1")
print("=" * 60)

# n=8: compute beta_1, beta_2, beta_3 for many random tournaments
n = 8
print(f"\n--- n={n}: seesaw test ---")

results = []
seesaw_violations = 0
beta2_violations = 0
t0 = time.time()
target = 500

for trial in range(target):
    A = random_tournament(n)
    betti = compute_betti_full(A, n, max_dim=4)

    b0, b1, b2, b3 = betti[0], betti[1], betti[2], betti[3]
    b4 = betti[4] if len(betti) > 4 else 0

    results.append((b0, b1, b2, b3, b4))

    if b1 > 0 and b3 > 0:
        seesaw_violations += 1
    if b2 > 0:
        beta2_violations += 1

    if (trial+1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  [{trial+1}/{target}] {elapsed:.1f}s, seesaw_violations={seesaw_violations}, b2>0={beta2_violations}")

elapsed = time.time() - t0
print(f"\n  Total: {len(results)} tournaments, {elapsed:.1f}s")

# Analyze
b1_hist = Counter(r[1] for r in results)
b2_hist = Counter(r[2] for r in results)
b3_hist = Counter(r[3] for r in results)
b4_hist = Counter(r[4] for r in results)

print(f"\n  beta_1 distribution: {dict(sorted(b1_hist.items()))}")
print(f"  beta_2 distribution: {dict(sorted(b2_hist.items()))}")
print(f"  beta_3 distribution: {dict(sorted(b3_hist.items()))}")
print(f"  beta_4 distribution: {dict(sorted(b4_hist.items()))}")

print(f"\n  beta_1 * beta_3 = 0 seesaw: {'HOLDS' if seesaw_violations == 0 else 'VIOLATED'}")
print(f"  Seesaw violations: {seesaw_violations}/{len(results)}")
print(f"  beta_2 = 0: {'ALWAYS' if beta2_violations == 0 else f'VIOLATED ({beta2_violations} times)'}")

# Cross-tabulate beta_1 vs beta_3
print(f"\n  Cross-tabulation beta_1 vs beta_3:")
cross = Counter((r[1], r[3]) for r in results)
b1_vals = sorted(set(r[1] for r in results))
b3_vals = sorted(set(r[3] for r in results))
print(f"  {'':>8}", end="")
for b3 in b3_vals:
    print(f"  b3={b3:>2}", end="")
print()
for b1 in b1_vals:
    print(f"  b1={b1:>2}", end="")
    for b3 in b3_vals:
        print(f"  {cross.get((b1,b3), 0):>5}", end="")
    print()

# Also cross-tabulate beta_3 vs beta_4
print(f"\n  Cross-tabulation beta_3 vs beta_4:")
cross34 = Counter((r[3], r[4]) for r in results)
b3v = sorted(set(r[3] for r in results))
b4v = sorted(set(r[4] for r in results))
print(f"  {'':>8}", end="")
for b4 in b4v:
    print(f"  b4={b4:>2}", end="")
print()
for b3 in b3v:
    print(f"  b3={b3:>2}", end="")
    for b4 in b4v:
        print(f"  {cross34.get((b3,b4), 0):>5}", end="")
    print()

# Betti profile histogram
print(f"\n  Full Betti profile histogram:")
profiles = Counter(tuple(r) for r in results)
for prof, count in sorted(profiles.items(), key=lambda x: -x[1])[:15]:
    print(f"    {list(prof)}: {count} ({count/len(results)*100:.1f}%)")
