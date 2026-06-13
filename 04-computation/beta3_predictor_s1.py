"""
Beta_3 Predictor: which tournament invariants determine beta_3?

Cross-correlate beta_3 with: H, score entropy, alpha_1, alpha_2,
regularity, SC status, c_3, c_5, Pfaffian, etc.

kind-pasteur-2026-03-25-S1
"""
import sys
import time
import random
import numpy as np
from collections import Counter, defaultdict
from itertools import permutations, combinations
from math import comb, factorial

sys.stdout.reconfigure(line_buffering=True)
random.seed(314159)

def random_tournament(n):
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

def score_variance(A, n):
    scores = [sum(A[i]) for i in range(n)]
    mean = sum(scores) / n
    return sum((s - mean)**2 for s in scores) / n

def is_sc_score(s, n):
    return all(s[i] + s[n-1-i] == n-1 for i in range(n))

def count_3cycles(A, n):
    """Count directed 3-cycle vertex sets."""
    count = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        # Check if cyclic: i->j->k->i or i->k->j->i
        if (A[i][j] and A[j][k] and A[k][i]):
            count += 1
        if (A[i][k] and A[k][j] and A[j][i]):
            count += 1
    return count  # This counts directed 3-cycles

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1<<v, v)] = 1
    for sz in range(2, n+1):
        for mask in range(1<<n):
            if bin(mask).count('1') != sz:
                continue
            for v in range(n):
                if not (mask & (1<<v)):
                    continue
                prev = mask ^ (1<<v)
                for u in range(n):
                    if u == v or not (prev & (1<<u)):
                        continue
                    if A[u][v]:
                        dp[(mask,v)] = dp.get((mask,v), 0) + dp.get((prev,u), 0)
    full = (1<<n)-1
    return sum(dp.get((full,v), 0) for v in range(n))

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

def compute_beta_p(A, n, p):
    """Compute beta_p."""
    ap_m1 = enum_allowed(A, n, p-1) if p >= 1 else [()]
    ap = enum_allowed(A, n, p)
    ap_p1 = enum_allowed(A, n, p+1)
    if not ap:
        return 0

    ap_m1_set = set(ap_m1)
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
            return 0
        omega_basis = Vt[rank_P:].T

    # ker(d_p)
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
        _, S_dp, _ = np.linalg.svd(D_omega)
        rank_dp = sum(s > 1e-8 for s in S_dp)
    else:
        rank_dp = 0
    ker_dp = omega_basis.shape[1] - rank_dp

    # im(d_{p+1})
    if not ap_p1:
        return ker_dp

    ap_set = set(ap)
    non_allowed2 = {}
    na2 = 0
    for path in ap_p1:
        for i in range(p+2):
            face = path[:i] + path[i+1:]
            if face not in ap_set and face not in non_allowed2:
                non_allowed2[face] = na2
                na2 += 1

    if not non_allowed2:
        omega_p1 = np.eye(len(ap_p1))
    else:
        Q = np.zeros((len(non_allowed2), len(ap_p1)))
        for j, path in enumerate(ap_p1):
            for i in range(p+2):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in non_allowed2:
                    Q[non_allowed2[face], j] += sign
        _, Sq, Vqt = np.linalg.svd(Q)
        rQ = sum(s > 1e-8 for s in Sq)
        if rQ >= len(ap_p1):
            return ker_dp
        omega_p1 = Vqt[rQ:].T

    ap_idx = {path: i for i, path in enumerate(ap)}
    D1 = np.zeros((len(ap), len(ap_p1)))
    for j, path in enumerate(ap_p1):
        for i in range(p+2):
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            if face in ap_idx:
                D1[ap_idx[face], j] += sign

    D1_omega = D1 @ omega_p1
    im_matrix = omega_basis.T @ D1_omega
    _, S_im, _ = np.linalg.svd(im_matrix)
    im_dp1 = sum(s > 1e-8 for s in S_im)

    return max(0, ker_dp - im_dp1)

# ===== MAIN =====
print("=" * 60)
print("BETA_3 PREDICTOR ANALYSIS")
print("kind-pasteur-2026-03-25-S1")
print("=" * 60)

n = 7
target = 300
print(f"\n--- n={n}: {target} samples ---")

data = []
t0 = time.time()

for trial in range(target):
    A = random_tournament(n)
    s = score_seq(A, n)
    sv = score_variance(A, n)
    sc = is_sc_score(s, n)
    c3 = count_3cycles(A, n)
    H = count_ham_paths(A, n)
    b1 = compute_beta_p(A, n, 1)
    b3 = compute_beta_p(A, n, 3)

    data.append({
        'score': s,
        'score_var': sv,
        'sc_score': sc,
        'c3': c3,
        'H': H,
        'b1': b1,
        'b3': b3,
    })

    if (trial+1) % 100 == 0:
        print(f"  [{trial+1}/{target}] {time.time()-t0:.1f}s")

print(f"  Done: {time.time()-t0:.1f}s")

# Analyze: what predicts beta_3 > 0?
b3_pos = [d for d in data if d['b3'] > 0]
b3_zero = [d for d in data if d['b3'] == 0]

print(f"\n  beta_3 > 0: {len(b3_pos)}/{len(data)} = {len(b3_pos)/len(data)*100:.1f}%")

# Compare means
for key in ['H', 'c3', 'score_var']:
    mean_pos = np.mean([d[key] for d in b3_pos]) if b3_pos else 0
    mean_zero = np.mean([d[key] for d in b3_zero]) if b3_zero else 0
    print(f"\n  {key}:")
    print(f"    beta_3=0: mean={mean_zero:.2f}")
    print(f"    beta_3>0: mean={mean_pos:.2f}")
    if b3_pos and b3_zero:
        # t-test
        from scipy import stats
        t, p = stats.ttest_ind([d[key] for d in b3_pos], [d[key] for d in b3_zero])
        print(f"    t-test: t={t:.3f}, p={p:.4f} {'***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'}")

# SC score vs beta_3
sc_pos = sum(1 for d in b3_pos if d['sc_score'])
sc_zero = sum(1 for d in b3_zero if d['sc_score'])
sc_total_pos = sum(1 for d in data if d['sc_score'])
sc_total = len(data)
print(f"\n  SC score:")
print(f"    beta_3=0: SC rate = {sc_zero/(len(b3_zero) if b3_zero else 1)*100:.1f}%")
print(f"    beta_3>0: SC rate = {sc_pos/(len(b3_pos) if b3_pos else 1)*100:.1f}%")
print(f"    Overall SC rate: {sc_total_pos/sc_total*100:.1f}%")

# Score sequence distribution
print(f"\n  Score sequences with beta_3 > 0:")
score_counts = Counter(d['score'] for d in b3_pos)
for s, c in score_counts.most_common(10):
    total_with_score = sum(1 for d in data if d['score'] == s)
    print(f"    {s}: {c}/{total_with_score} = {c/total_with_score*100:.0f}%")

# H histogram by beta_3
print(f"\n  H distribution:")
H_b3_0 = sorted([d['H'] for d in b3_zero])
H_b3_1 = sorted([d['H'] for d in b3_pos])
if H_b3_0:
    print(f"    beta_3=0: min={H_b3_0[0]}, median={H_b3_0[len(H_b3_0)//2]}, max={H_b3_0[-1]}")
if H_b3_1:
    print(f"    beta_3>0: min={H_b3_1[0]}, median={H_b3_1[len(H_b3_1)//2]}, max={H_b3_1[-1]}")

# c3 histogram by beta_3
print(f"\n  c3 (directed 3-cycles) distribution:")
c3_b3_0 = sorted([d['c3'] for d in b3_zero])
c3_b3_1 = sorted([d['c3'] for d in b3_pos])
if c3_b3_0:
    print(f"    beta_3=0: min={c3_b3_0[0]}, median={c3_b3_0[len(c3_b3_0)//2]}, max={c3_b3_0[-1]}")
if c3_b3_1:
    print(f"    beta_3>0: min={c3_b3_1[0]}, median={c3_b3_1[len(c3_b3_1)//2]}, max={c3_b3_1[-1]}")

# beta_1 vs beta_3 (seesaw)
print(f"\n  Seesaw at n={n}:")
both = sum(1 for d in data if d['b1'] > 0 and d['b3'] > 0)
print(f"    beta_1>0 AND beta_3>0: {both}/{len(data)}")
