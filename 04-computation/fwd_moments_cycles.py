#!/usr/bin/env python3
"""
fwd_moments_cycles.py — What cycle invariants determine E[fwd^r]?

We know:
  E[fwd^0] = 1 (trivial)
  E[fwd^1] = (n-1)/2 (universal, by palindrome)
  E[fwd^2] = (n-1)(3n-2)/12 + 4*t3/(n(n-1)) (from THM-089)

QUESTION: What determines E[fwd^3]?

fwd = X_0 + X_1 + ... + X_{n-2} where X_i = indicator(P[i]->P[i+1] in T)
fwd^3 = (sum X_i)^3 = sum X_i^3 + 3 sum_{i!=j} X_i^2 X_j + 6 sum_{i<j<k} X_i X_j X_k

Since X_i is binary: X_i^r = X_i for all r >= 1.
So: fwd^3 = sum X_i + 3 sum_{i!=j} X_i X_j + 6 sum_{i<j<k} X_i X_j X_k
     = fwd + 3(fwd^2 - fwd) + 6 sum_{i<j<k} X_i X_j X_k
     = -2*fwd + 3*fwd^2 + 6 * sum_{i<j<k} X_i X_j X_k

So E[fwd^3] = -2*(n-1)/2 + 3*E[fwd^2] + 6 * sum_{i<j<k} E[X_i X_j X_k]

The triple correlation E[X_i X_j X_k] is the probability that ALL three
consecutive pairs (P[i],P[i+1]), (P[j],P[j+1]), (P[k],P[k+1]) are
forward edges in T.

By the non-adjacent uncorrelatedness result (THM-089):
  E[X_i X_j] = E[X_i]*E[X_j] = 1/4 when |i-j| >= 2
  E[X_i X_{i+1}] = 1/4 + Cov = ...

For triples, we need to consider adjacency patterns:
  - All three non-adjacent: E[X_i X_j X_k] = E[X_i]E[X_j]E[X_k] = 1/8 (if all pairwise non-adjacent)
  - One adjacent pair: more complex
  - Two adjacent pairs (e.g., i,i+1,i+2): involves 4-vertex correlations

Let me compute E[fwd^3] directly and see what invariants determine it.

Author: opus-2026-03-07-S46c
"""
from itertools import permutations, combinations
from math import comb, factorial
from fractions import Fraction
from collections import defaultdict
import numpy as np

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def compute_F(adj, n):
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def count_3cycles(adj, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t3 += 1
    return t3

def count_kcycles(adj, n, k):
    count = 0
    for combo in combinations(range(n), k):
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%k]] for i in range(k)):
                count += 1
    return count // k

print("=" * 60)
print("E[fwd^r] AND CYCLE INVARIANTS")
print("=" * 60)

for n in [4, 5, 6]:
    m_vals = n*(n-1)//2
    seen = set()
    data = []

    for bits in range(1 << m_vals):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        t3 = count_3cycles(adj, n)
        t5 = count_kcycles(adj, n, 5) if n >= 5 else 0

        total = sum(F)
        moments = {}
        for r in range(6):
            moments[r] = Fraction(sum(k**r * F[k] for k in range(n)), total)

        data.append({
            't3': t3, 't5': t5, 'F': F,
            'moments': moments
        })

    print(f"\nn={n}: {len(data)} distinct F-classes")

    # Check which moments are determined by t3
    for r in range(5):
        # Group by t3
        t3_to_mom = defaultdict(set)
        for d in data:
            t3_to_mom[d['t3']].add(d['moments'][r])
        determined = all(len(v) == 1 for v in t3_to_mom.values())
        if determined:
            # Find linear formula
            items = [(t3, list(vals)[0]) for t3, vals in sorted(t3_to_mom.items())]
            if len(items) >= 2:
                slope = (items[1][1] - items[0][1]) / (items[1][0] - items[0][0])
                intercept = items[0][1] - slope * items[0][0]
                ok = all(abs(v - (intercept + slope * t)) < Fraction(1, 10000) for t, v in items)
                if ok:
                    print(f"  E[fwd^{r}] = {float(intercept):.6f} + {float(slope):.6f} * t3  [EXACT linear in t3]")
                else:
                    print(f"  E[fwd^{r}]: determined by t3 but NOT linear")
            else:
                val = list(t3_to_mom[items[0][0]])[0]
                print(f"  E[fwd^{r}] = {float(val):.6f}  [universal, only 1 t3 value]")
        else:
            n_amb = sum(1 for v in t3_to_mom.values() if len(v) > 1)
            print(f"  E[fwd^{r}]: NOT determined by t3 alone ({n_amb} ambiguous t3 values)")

            # Check if (t3, t5) determines it
            if n >= 5:
                t3t5_to_mom = defaultdict(set)
                for d in data:
                    t3t5_to_mom[(d['t3'], d['t5'])].add(d['moments'][r])
                det2 = all(len(v) == 1 for v in t3t5_to_mom.values())
                if det2:
                    print(f"    BUT determined by (t3, t5)")
                else:
                    n_amb2 = sum(1 for v in t3t5_to_mom.values() if len(v) > 1)
                    print(f"    NOT determined by (t3, t5) either ({n_amb2} ambiguous)")

    # Print moments table at n=6
    if n == 6:
        print(f"\n  Full moment table:")
        print(f"  {'t3':>3} {'t5':>3} {'E[fwd]':>10} {'E[fwd^2]':>12} {'E[fwd^3]':>14} {'E[fwd^4]':>16}")
        for d in sorted(data, key=lambda x: (x['t3'], x['t5'])):
            m1 = float(d['moments'][1])
            m2 = float(d['moments'][2])
            m3 = float(d['moments'][3])
            m4 = float(d['moments'][4])
            print(f"  {d['t3']:>3} {d['t5']:>3} {m1:>10.4f} {m2:>12.6f} {m3:>14.6f} {m4:>16.6f}")

# ============================================================
# KEY QUESTION: E[fwd^3] formula
# ============================================================
print("\n" + "=" * 60)
print("E[fwd^3] ANALYSIS AT n=6")
print("=" * 60)

n = 6
seen = set()
data6 = []
for bits in range(1 << (n*(n-1)//2)):
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    t3 = count_3cycles(adj, n)
    t5 = count_kcycles(adj, n, 5)

    # Also compute alpha_2
    cycles3 = []
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            cycles3.append((i, j, k))
    alpha_2 = 0
    for a in range(len(cycles3)):
        for b in range(a+1, len(cycles3)):
            if set(cycles3[a]).isdisjoint(set(cycles3[b])):
                alpha_2 += 1

    total = sum(F)
    m3 = Fraction(sum(k**3 * F[k] for k in range(n)), total)

    data6.append({
        't3': t3, 't5': t5, 'alpha_2': alpha_2,
        'E_fwd3': m3
    })

# Linear regression: E[fwd^3] = a + b*t3 + c*t5 + d*alpha_2
X = np.array([[1, d['t3'], d['t5'], d['alpha_2']] for d in data6], dtype=float)
y = np.array([float(d['E_fwd3']) for d in data6])
c, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
err = max(abs(X @ c - y))
print(f"E[fwd^3] = {c[0]:.6f} + {c[1]:.6f}*t3 + {c[2]:.6f}*t5 + {c[3]:.6f}*alpha_2")
print(f"max error: {err:.10f}")

# Try exact rational coefficients
from fractions import Fraction
# E[fwd^3] for transitive (t3=t5=alpha_2=0)
trans = [d for d in data6 if d['t3'] == 0][0]
print(f"\nTransitive: E[fwd^3] = {trans['E_fwd3']} = {float(trans['E_fwd3']):.6f}")

# Check if t3 alone works
X_t3 = np.array([[1, d['t3']] for d in data6], dtype=float)
c_t3, _, _, _ = np.linalg.lstsq(X_t3, y, rcond=None)
err_t3 = max(abs(X_t3 @ c_t3 - y))
print(f"\nt3 only: E[fwd^3] = {c_t3[0]:.6f} + {c_t3[1]:.6f}*t3, max_err={err_t3:.10f}")

if err_t3 < 1e-8:
    print("  --> E[fwd^3] is EXACTLY determined by t3!")
else:
    print(f"  --> E[fwd^3] NOT determined by t3 alone")
    if err < 1e-8:
        print(f"  --> BUT determined by (t3, t5, alpha_2)")
    else:
        print(f"  --> NOT determined by (t3, t5, alpha_2) either")
