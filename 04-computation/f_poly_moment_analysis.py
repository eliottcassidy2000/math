#!/usr/bin/env python3
"""
f_poly_moment_analysis.py — Moment generating function interpretation.

F(T,x) = sum_P x^{fwd(P)} = E_P[x^X] where X = fwd(P) is a random
variable on {0, 1, ..., n-1}.

The MGF is M_X(t) = E[e^{tX}] = F(T, e^t) / n!.
The cumulant generating function is K_X(t) = log M_X(t).

Moments: mu_k = E[X^k] = (1/n!) sum_P fwd(P)^k
Cumulants: kappa_1 = mean, kappa_2 = variance, kappa_3 = skewness*sigma^3

We know:
  mu_1 = (n-1)/2 (universal!)
  mu_2 = 2*D_2/n! + mu_1 (from D_2 = sum C(fwd,2))
  var = mu_2 - mu_1^2 = 2/3 + 2*t3/21 at n=7

QUESTION: Are higher cumulants also simple functions of t3, t5, etc.?

ALSO: The centered generating function
  G(T, x) = F(T, x) / (1+x)^{(n-1)/2}
or
  G(T, x) = F(T, x) / x^{(n-1)/2}
might have simpler structure.

Author: opus-2026-03-07-S44
"""
from itertools import permutations, combinations
import math
import random
import numpy as np

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F_dp(A, n):
    full = (1 << n) - 1
    dp = [[None]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = {0: 1}
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] is None:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                new_mask = mask | (1 << u)
                fwd_edge = A[v][u]
                if dp[new_mask][u] is None:
                    dp[new_mask][u] = {}
                for fwd_count, cnt in dp[mask][v].items():
                    new_fwd = fwd_count + fwd_edge
                    dp[new_mask][u][new_fwd] = dp[new_mask][u].get(new_fwd, 0) + cnt
    F = [0] * n
    for v in range(n):
        if dp[full][v] is not None:
            for fwd_count, cnt in dp[full][v].items():
                F[fwd_count] += cnt
    return F

def count_t3(A, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            t3 += 1
    return t3

def count_t5(A, n):
    t5 = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t5 += 1
                break
    return t5

# ============================================================
# MOMENTS OF fwd(P)
# ============================================================
print("=" * 60)
print("MOMENTS OF fwd(P) AS FUNCTION OF TOURNAMENT INVARIANTS")
print("=" * 60)

n = 7
m = n*(n-1)//2
nfact = math.factorial(n)
random.seed(42)

seen = set()
data = []

for trial in range(300):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    F = compute_F_dp(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    t3 = count_t3(A, n)

    # Moments
    mu = [(1/nfact) * sum(F[k] * k**p for k in range(n)) for p in range(6)]

    # Centered moments
    center = mu[1]
    cm = [(1/nfact) * sum(F[k] * (k - center)**p for k in range(n)) for p in range(6)]

    data.append({
        't3': t3, 'mu': mu, 'cm': cm, 'F': F, 'H': F[n-1]
    })

# Check which moments depend only on t3
print(f"\nn=7: {len(data)} distinct tournaments")
print(f"\nMoment analysis (checking if mu_k depends only on t3):")

for p in range(1, 6):
    t3_to_mu = {}
    for d in data:
        t3 = d['t3']
        mu_p = round(d['mu'][p], 8)
        if t3 not in t3_to_mu:
            t3_to_mu[t3] = set()
        t3_to_mu[t3].add(mu_p)

    # Is mu_p a function of t3 alone?
    unique_per_t3 = max(len(v) for v in t3_to_mu.values())
    if unique_per_t3 == 1:
        # Linear regression
        t3_vals = sorted(t3_to_mu.keys())
        mu_vals = [list(t3_to_mu[t])[0] for t in t3_vals if len(t3_to_mu[t]) == 1]
        if len(mu_vals) >= 2:
            slope = (mu_vals[-1] - mu_vals[0]) / (t3_vals[-1] - t3_vals[0])
            intercept = mu_vals[0] - slope * t3_vals[0]
            print(f"  mu_{p}: DEPENDS ONLY ON t3! ≈ {intercept:.6f} + {slope:.6f}*t3")
    else:
        print(f"  mu_{p}: depends on MORE than t3 (max {unique_per_t3} values per t3)")

# Centered moments
print(f"\nCentered moment analysis:")
for p in range(2, 6):
    t3_to_cm = {}
    for d in data:
        t3 = d['t3']
        cm_p = round(d['cm'][p], 8)
        if t3 not in t3_to_cm:
            t3_to_cm[t3] = set()
        t3_to_cm[t3].add(cm_p)

    unique_per_t3 = max(len(v) for v in t3_to_cm.values())
    if unique_per_t3 == 1:
        t3_vals = sorted(t3_to_cm.keys())
        cm_vals = [list(t3_to_cm[t])[0] for t in t3_vals]
        if len(cm_vals) >= 2:
            slope = (cm_vals[-1] - cm_vals[0]) / (t3_vals[-1] - t3_vals[0])
            intercept = cm_vals[0] - slope * t3_vals[0]
            print(f"  cm_{p}: DEPENDS ONLY ON t3! ≈ {intercept:.6f} + {slope:.6f}*t3")
    else:
        print(f"  cm_{p}: depends on MORE than t3 (max {unique_per_t3} values per t3)")

# ============================================================
# CUMULANTS
# ============================================================
print("\n" + "=" * 60)
print("CUMULANTS OF fwd(P)")
print("=" * 60)
# kappa_1 = mu_1
# kappa_2 = cm_2 = variance
# kappa_3 = cm_3
# kappa_4 = cm_4 - 3*cm_2^2
# kappa_5 = cm_5 - 10*cm_3*cm_2

for trial_idx, d in enumerate(data[:5]):
    cm = d['cm']
    kappa = [0] * 6
    kappa[1] = d['mu'][1]
    kappa[2] = cm[2]
    kappa[3] = cm[3]
    kappa[4] = cm[4] - 3 * cm[2]**2
    kappa[5] = cm[5] - 10 * cm[3] * cm[2]

    print(f"  t3={d['t3']}: kappa = {[f'{k:.6f}' for k in kappa[1:]]}")

# ============================================================
# NOVEL: F POLYNOMIAL DIVIDED BY (1+x)^{(n-1)/2}
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) / (1+x)^{(n-1)/2} — REDUCED POLYNOMIAL")
print("=" * 60)
# Since F has palindromic coefficients, we can write
# F(x) = (1+x)^a * Q(x) where Q might have nice properties.
# Let's compute the quotient when dividing by (1+x).

# At n=5: F is degree 4 (even). (1+x)^2 = 1+2x+x^2.
# If F(-1) = 0, then (1+x) | F. S(T)=0 means F(-1)=0.

n = 5
m = n*(n-1)//2
print(f"\nn=5:")

seen = set()
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F_dp(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]
    # Check if (1+x) divides F
    f_at_neg1 = sum(F[k] * (-1)**k for k in range(n))

    # How many times does (1+x) divide F?
    poly = list(F)  # coefficients from x^0 to x^{n-1}
    div_count = 0
    while len(poly) > 1:
        val = sum(poly[k] * (-1)**k for k in range(len(poly)))
        if val == 0:
            # Divide by (1+x): polynomial long division
            new_poly = [0] * (len(poly) - 1)
            new_poly[0] = poly[0]
            for i in range(1, len(poly) - 1):
                new_poly[i] = poly[i] - (-1) * new_poly[i-1]
                # Actually: divide [a0, a1, ..., an] by [1, 1]
                # Quotient: q0 = a0, q1 = a1 - q0, q2 = a2 - q1, ...
            # Synthetic division by (x + 1) = (x - (-1))
            quot = [0] * (len(poly) - 1)
            quot[-1] = poly[-1]
            for i in range(len(poly) - 2, 0, -1):
                quot[i-1] = poly[i] - (-1) * quot[i]  # wait, wrong direction

            # Let me just use numpy
            quotient, remainder = np.polydiv(
                [F[k] for k in range(n-1, -1, -1)] if div_count == 0 else list(reversed(poly)),
                [1, 1]  # x + 1
            )
            if abs(remainder[-1]) < 1e-6:
                div_count += 1
                poly = [int(round(c)) for c in reversed(quotient)]
            else:
                break
        else:
            break

    reduced = poly
    print(f"  H={H:3d}: F(-1)={f_at_neg1:4d}, (1+x) divides {div_count} times, "
          f"reduced={reduced}")

# ============================================================
# NOVEL: ENTROPY OF F(T,x)
# ============================================================
print("\n" + "=" * 60)
print("ENTROPY OF THE fwd(P) DISTRIBUTION")
print("=" * 60)

n = 7
random.seed(42)
seen = set()

for trial in range(200):
    bits = random.getrandbits(21)
    A = tournament_from_bits(bits, n)
    F = compute_F_dp(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H_val = F[n-1]
    t3 = count_t3(A, n)

    # Shannon entropy of F_k/n!
    probs = [F[k]/nfact for k in range(n)]
    entropy = -sum(p * math.log2(p) for p in probs if p > 0)

    # Compare with binomial B(n-1, 1/2) entropy
    binom_probs = [math.comb(n-1, k) / 2**(n-1) for k in range(n)]
    binom_entropy = -sum(p * math.log2(p) for p in binom_probs if p > 0)

    if trial < 15 or t3 in [0, 1, 13, 14, 15]:
        print(f"  H={H_val:3d} t3={t3:2d}: entropy={entropy:.6f}, "
              f"binom_entropy={binom_entropy:.6f}, ratio={entropy/binom_entropy:.6f}")
