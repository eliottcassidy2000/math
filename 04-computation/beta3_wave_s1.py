"""
Beta_3 wave computation: extend the defect rate pattern to n=8,9.
Verify HYP-402 (beta_3 <= 2) with larger sample.

The "defect wave" (HYP-338):
  n=5: beta_1=29.7%, beta_3=0%
  n=6: beta_1=14.6%, beta_3=1.2%
  n=7: beta_1=5.8%, beta_3=8-11%
  n=8: beta_1=~1%, beta_3=~21%

Goal: compute beta_1 and beta_3 rates at n=8 (verify) and n=9 (new data).

kind-pasteur-2026-03-25-S1
"""
import sys
import time
import random
import numpy as np
from collections import defaultdict, Counter
from itertools import permutations
from math import comb

sys.stdout.reconfigure(line_buffering=True)
random.seed(2025)

def random_tournament(n):
    """Generate random tournament on n vertices."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def is_sc_score(s, n):
    return all(s[i] + s[n-1-i] == n-1 for i in range(n))

def enum_allowed(A, n, p):
    """Enumerate allowed p-paths (length p, p+1 vertices)."""
    if p == 0:
        return list(range(n))  # Just vertices
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

def compute_beta_fast(A, n, target_dim):
    """Compute beta_{target_dim} efficiently.

    Uses the formula: beta_p = dim(ker d_p on Omega_p) - dim(im d_{p+1} from Omega_{p+1})

    Needs: Omega_p, Omega_{p+1}, and boundary maps between them.
    """
    p = target_dim

    # Enumerate allowed paths at dimensions p-1, p, p+1, p+2
    ap_m1 = enum_allowed(A, n, p-1) if p >= 1 else [()]
    ap = enum_allowed(A, n, p)
    ap_p1 = enum_allowed(A, n, p+1)

    if not ap:
        return 0

    # Compute Omega_p: subspace of A_p whose boundary is in A_{p-1}
    # For each allowed p-path, check which faces are NOT allowed
    ap_m1_set = set(ap_m1) if p >= 1 else set()

    # Find non-allowed faces
    non_allowed_faces = {}
    na_idx = 0

    for path in ap:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in ap_m1_set and face not in non_allowed_faces:
                non_allowed_faces[face] = na_idx
                na_idx += 1

    if not non_allowed_faces:
        # Omega_p = A_p (all boundaries land in A_{p-1})
        omega_p_dim = len(ap)
        omega_p_to_ap = np.eye(len(ap))
    else:
        # Build non-allowed face projection matrix
        P = np.zeros((len(non_allowed_faces), len(ap)), dtype=np.float64)
        for j, path in enumerate(ap):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in non_allowed_faces:
                    P[non_allowed_faces[face], j] += sign

        # Omega_p = null space of P
        _, S, Vt = np.linalg.svd(P)
        rank_P = sum(s > 1e-8 for s in S)
        omega_p_dim = len(ap) - rank_P
        if omega_p_dim == 0:
            return 0
        omega_p_to_ap = Vt[rank_P:].T  # columns are basis of null space

    # Build boundary d_p: A_p -> A_{p-1}
    if p >= 1 and ap_m1:
        ap_m1_idx = {path: i for i, path in enumerate(ap_m1)}
        D_p = np.zeros((len(ap_m1), len(ap)), dtype=np.float64)
        for j, path in enumerate(ap):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in ap_m1_idx:
                    D_p[ap_m1_idx[face], j] += sign

        # Restrict d_p to Omega_p
        D_p_omega = D_p @ omega_p_to_ap
        _, S_dp, _ = np.linalg.svd(D_p_omega)
        ker_dp_dim = omega_p_dim - sum(s > 1e-8 for s in S_dp)
    else:
        ker_dp_dim = omega_p_dim

    # Now compute Omega_{p+1} and im(d_{p+1}) restricted to Omega_{p+1}
    if not ap_p1:
        im_dp1_dim = 0
    else:
        # Compute Omega_{p+1}
        ap_set = set(ap)
        non_allowed_p = {}
        na_idx2 = 0
        for path in ap_p1:
            for i in range(p+2):
                face = path[:i] + path[i+1:]
                if face not in ap_set and face not in non_allowed_p:
                    non_allowed_p[face] = na_idx2
                    na_idx2 += 1

        if not non_allowed_p:
            omega_p1_dim = len(ap_p1)
            omega_p1_to_ap1 = np.eye(len(ap_p1))
        else:
            Q = np.zeros((len(non_allowed_p), len(ap_p1)), dtype=np.float64)
            for j, path in enumerate(ap_p1):
                for i in range(p+2):
                    face = path[:i] + path[i+1:]
                    sign = (-1)**i
                    if face in non_allowed_p:
                        Q[non_allowed_p[face], j] += sign

            _, Sq, Vqt = np.linalg.svd(Q)
            rank_Q = sum(s > 1e-8 for s in Sq)
            omega_p1_dim = len(ap_p1) - rank_Q
            if omega_p1_dim == 0:
                return ker_dp_dim  # beta = ker_dp_dim - 0
            omega_p1_to_ap1 = Vqt[rank_Q:].T

        # Build d_{p+1}: A_{p+1} -> A_p
        ap_idx = {path: i for i, path in enumerate(ap)}
        D_p1 = np.zeros((len(ap), len(ap_p1)), dtype=np.float64)
        for j, path in enumerate(ap_p1):
            for i in range(p+2):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in ap_idx:
                    D_p1[ap_idx[face], j] += sign

        # Restrict d_{p+1} to Omega_{p+1}, landing in A_p
        D_p1_omega = D_p1 @ omega_p1_to_ap1

        # Project into Omega_p basis
        # The image of d_{p+1} should land in Omega_p (d^2 = 0)
        # So we compute rank of D_p1_omega projected into Omega_p
        # = rank of omega_p_to_ap.T @ D_p1_omega
        im_matrix = omega_p_to_ap.T @ D_p1_omega
        _, S_im, _ = np.linalg.svd(im_matrix)
        im_dp1_dim = sum(s > 1e-8 for s in S_im)

    beta = ker_dp_dim - im_dp1_dim
    return max(0, beta)  # Should be non-negative

def compute_beta1(A, n):
    """Compute beta_1 for a tournament."""
    return compute_beta_fast(A, n, 1)

def compute_beta3(A, n):
    """Compute beta_3 for a tournament."""
    return compute_beta_fast(A, n, 3)

# ===== MAIN =====

print("=" * 60)
print("BETA WAVE COMPUTATION")
print("kind-pasteur-2026-03-25-S1")
print("=" * 60)

# Test at n=5 first (known: beta_1 rate ~29.7%)
n = 5
print(f"\n--- n={n}: calibration ---")
b1_count = 0
total = 0
t0 = time.time()
for _ in range(200):
    A = random_tournament(n)
    b1 = compute_beta1(A, n)
    if b1 > 0:
        b1_count += 1
    total += 1
rate = b1_count / total * 100
print(f"  beta_1 > 0: {b1_count}/{total} = {rate:.1f}% (expected ~29.7%)")
print(f"  Time: {time.time()-t0:.1f}s")

# n=7 calibration
n = 7
print(f"\n--- n={n}: calibration ---")
b1_count = 0
b3_count = 0
total = 0
t0 = time.time()
for _ in range(100):
    A = random_tournament(n)
    b1 = compute_beta1(A, n)
    b3 = compute_beta3(A, n)
    if b1 > 0:
        b1_count += 1
    if b3 > 0:
        b3_count += 1
    total += 1
print(f"  beta_1 > 0: {b1_count}/{total} = {b1_count/total*100:.1f}% (expected ~5.8%)")
print(f"  beta_3 > 0: {b3_count}/{total} = {b3_count/total*100:.1f}% (expected ~8-11%)")
print(f"  Time: {time.time()-t0:.1f}s ({(time.time()-t0)/total:.2f}s per tournament)")

# n=8
n = 8
print(f"\n--- n={n}: main computation ---")
b1_hist = Counter()
b3_hist = Counter()
total = 0
t0 = time.time()
target = 200

for trial in range(target):
    A = random_tournament(n)
    b1 = compute_beta1(A, n)
    b3 = compute_beta3(A, n)
    b1_hist[b1] += 1
    b3_hist[b3] += 1
    total += 1
    if total % 50 == 0:
        elapsed = time.time() - t0
        print(f"  [{total}/{target}] {elapsed:.1f}s, b1>0={sum(v for k,v in b1_hist.items() if k>0)}, b3>0={sum(v for k,v in b3_hist.items() if k>0)}")

elapsed = time.time() - t0
print(f"\n  n={n} Results ({total} tournaments, {elapsed:.1f}s):")
print(f"  beta_1 distribution: {dict(sorted(b1_hist.items()))}")
print(f"  beta_3 distribution: {dict(sorted(b3_hist.items()))}")
b1_rate = sum(v for k,v in b1_hist.items() if k>0) / total * 100
b3_rate = sum(v for k,v in b3_hist.items() if k>0) / total * 100
print(f"  beta_1 > 0: {b1_rate:.1f}%")
print(f"  beta_3 > 0: {b3_rate:.1f}%")
max_b3 = max(b3_hist.keys())
print(f"  max beta_3 = {max_b3}")

# n=9 (slower, fewer samples)
n = 9
print(f"\n--- n={n}: new territory ---")
b1_hist9 = Counter()
b3_hist9 = Counter()
total = 0
t0 = time.time()
target = 50  # Fewer due to size

for trial in range(target):
    A = random_tournament(n)
    try:
        b1 = compute_beta1(A, n)
        b3 = compute_beta3(A, n)
    except Exception as e:
        print(f"  Trial {trial}: ERROR {e}")
        continue
    b1_hist9[b1] += 1
    b3_hist9[b3] += 1
    total += 1
    if total % 10 == 0:
        elapsed = time.time() - t0
        print(f"  [{total}/{target}] {elapsed:.1f}s")

elapsed = time.time() - t0
print(f"\n  n={n} Results ({total} tournaments, {elapsed:.1f}s):")
print(f"  beta_1 distribution: {dict(sorted(b1_hist9.items()))}")
print(f"  beta_3 distribution: {dict(sorted(b3_hist9.items()))}")
b1_rate9 = sum(v for k,v in b1_hist9.items() if k>0) / max(total,1) * 100
b3_rate9 = sum(v for k,v in b3_hist9.items() if k>0) / max(total,1) * 100
print(f"  beta_1 > 0: {b1_rate9:.1f}%")
print(f"  beta_3 > 0: {b3_rate9:.1f}%")

# Summary: the defect wave
print(f"\n{'='*60}")
print("DEFECT WAVE SUMMARY")
print(f"{'='*60}")
print(f"  n  | beta_1>0  | beta_3>0  | max(beta_3)")
print(f"  5  |   ~30%    |     0%    |     0")
print(f"  6  |  ~14.6%   |   ~1.2%   |     1")
print(f"  7  |   ~5.8%   |   ~8-11%  |     1")
print(f"  8  |  {b1_rate:5.1f}%   |  {b3_rate:5.1f}%   |     {max_b3}")
if total > 0:
    max_b3_9 = max(b3_hist9.keys()) if b3_hist9 else 0
    print(f"  9  |  {b1_rate9:5.1f}%   |  {b3_rate9:5.1f}%   |     {max_b3_9}")
