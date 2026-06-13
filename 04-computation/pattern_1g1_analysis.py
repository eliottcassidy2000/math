#!/usr/bin/env python3
"""
Analyze the (1,g,1) pattern sigmas that contribute to w_{n-5}.

At n=7:
  (1,1,1): sigma_1111 — 4 consecutive positions
  (1,2,1): sigma_121 — positions {i,i+1,i+3,i+4}
  (1,3,1): sigma_131 — positions {i,i+1,i+4,i+5}

w_2 = 3*sigma_1111 + 2*sigma_121 + 1*sigma_131

Known: w_2 = -60*t_3 + 12*t_5 + 24*alpha_2 + 231

So: 3*sigma_1111 + 2*sigma_121 + sigma_131 = -60*t3 + 12*t5 + 24*a2 + 231

What does each sigma depend on?

sigma_1111 = sigma({0,1,2,3}) — this is the "4-path" sum
sigma_121 = sigma({0,1,3,4}) — "2-path, gap, 2-path"
sigma_131 = sigma({0,1,4,5}) — "2-path, big gap, 2-path"

Hypothesis: sigma_121 and sigma_131 factor into sigma_adj^2 somehow?
Or: sigma_121 depends on t5 and alpha_2?

kind-pasteur-2026-03-06-S25g
"""
from itertools import permutations, combinations
from math import factorial, comb
import random
import numpy as np

def compute_sigma(A, n, positions):
    total = 0.0
    for p in permutations(range(n)):
        prod = 1.0
        for i in positions:
            prod *= (A[p[i]][p[i+1]] - 0.5)
        total += prod
    return total

def count_directed_k_cycles(A, n, k):
    count = 0
    for verts in combinations(range(n), k):
        seen = set()
        for p in permutations(verts):
            if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                start = min(range(k), key=lambda i: p[i])
                canon = tuple(p[start:] + p[:start])
                if canon not in seen:
                    seen.add(canon)
                    count += 1
    return count

def count_alpha2(A, n):
    all_cycles = []
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            seen = set()
            for p in permutations(verts):
                if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                    start = min(range(k), key=lambda i: p[i])
                    canon = tuple(p[start:] + p[:start])
                    if canon not in seen:
                        seen.add(canon)
                        all_cycles.append(frozenset(verts))
    count = 0
    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            if all_cycles[i].isdisjoint(all_cycles[j]):
                count += 1
    return count

n = 7
random.seed(42)

data = []
for trial in range(15):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_directed_k_cycles(A, n, 3)
    t5 = count_directed_k_cycles(A, n, 5)
    a2 = count_alpha2(A, n)

    s1111 = compute_sigma(A, n, (0, 1, 2, 3))
    s121 = compute_sigma(A, n, (0, 1, 3, 4))
    s131 = compute_sigma(A, n, (0, 1, 4, 5))

    data.append({
        't3': t3, 't5': t5, 'a2': a2,
        's1111': s1111, 's121': s121, 's131': s131
    })

# Fit each sigma as linear function of (t3, t5, a2, 1)
for name, key in [('sigma_1111', 's1111'), ('sigma_121', 's121'), ('sigma_131', 's131')]:
    y = np.array([d[key] for d in data])
    X = np.array([[d['t3'], d['t5'], d['a2'], 1] for d in data])
    c = np.linalg.lstsq(X, y, rcond=None)[0]
    err = max(abs(y - X @ c))
    print(f"{name} = {c[0]:.4f}*t3 + {c[1]:.4f}*t5 + {c[2]:.4f}*a2 + {c[3]:.4f}  (max_err={err:.4f})")

# Check: do 3*s1111 + 2*s121 + s131 match w2?
print(f"\nVerify: 3*s1111 + 2*s121 + s131 = w2")
for d in data[:5]:
    w2 = 3*d['s1111'] + 2*d['s121'] + d['s131']
    w2_pred = -60*d['t3'] + 12*d['t5'] + 24*d['a2'] + 231
    print(f"  t3={d['t3']:2d} t5={d['t5']:2d} a2={d['a2']:2d}: w2={w2:.1f}, pred={w2_pred:.1f}")

print("\nDONE")
