#!/usr/bin/env python3
"""
fwd_cumulants.py — Cumulants of the forward-edge distribution.

The cumulant generating function is log(E[exp(t*fwd)]).
Cumulants kappa_r are often more natural than moments:
  kappa_1 = mean = (n-1)/2
  kappa_2 = variance = (n+1)/12 + 4t3/(n(n-1))
  kappa_3 = third central moment (= 0 by symmetry for transitive, but ?)
  kappa_4 = fourth cumulant (excess kurtosis * sigma^4)

For symmetric distributions: kappa_3 = 0.

QUESTION: Is kappa_3 = 0 for ALL tournaments (not just transitive)?
If so, that's equivalent to saying E[fwd^3] is determined by
E[fwd] and E[fwd^2] (which are both t3-determined), explaining
why E[fwd^3] depends only on t3.

Author: opus-2026-03-07-S46c
"""
from itertools import permutations, combinations
from math import comb, factorial
from fractions import Fraction

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
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t3 += 1
    return t3

def count_5cycles(adj, n):
    count = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                count += 1
    return count // 5

print("=" * 70)
print("CUMULANTS OF FORWARD-EDGE DISTRIBUTION")
print("=" * 70)

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
        t5 = count_5cycles(adj, n) if n >= 5 else 0

        total = sum(F)
        mu = Fraction(sum(k * F[k] for k in range(n)), total)
        # Central moments
        mu2 = Fraction(sum((k - mu)**2 * F[k] for k in range(n)), total)
        mu3 = Fraction(sum((k - mu)**3 * F[k] for k in range(n)), total)
        mu4 = Fraction(sum((k - mu)**4 * F[k] for k in range(n)), total)

        # Cumulants
        kappa1 = mu
        kappa2 = mu2
        kappa3 = mu3  # third cumulant = third central moment
        kappa4 = mu4 - 3 * mu2**2  # excess kurtosis * sigma^4

        data.append({
            't3': t3, 't5': t5,
            'kappa1': kappa1, 'kappa2': kappa2,
            'kappa3': kappa3, 'kappa4': kappa4,
            'F': F
        })

    print(f"\nn={n}: {len(data)} F-classes")
    print(f"  {'t3':>3} {'t5':>3} {'kappa1':>8} {'kappa2':>12} {'kappa3':>12} {'kappa4':>14}")
    for d in sorted(data, key=lambda x: (x['t3'], x['t5'])):
        print(f"  {d['t3']:>3} {d['t5']:>3} {str(d['kappa1']):>8} {str(d['kappa2']):>12} "
              f"{str(d['kappa3']):>12} {str(d['kappa4']):>14}")

    # Check: is kappa3 always 0?
    all_zero_k3 = all(d['kappa3'] == 0 for d in data)
    print(f"\n  kappa_3 = 0 for ALL tournaments: {all_zero_k3}")

    if all_zero_k3:
        print(f"  ==> F(T,x) always has symmetric (central) distribution!")
        print(f"      This means the palindrome symmetry forces zero skewness")
        print(f"      even though F(T,x) is NOT palindromic in general.")

    # Check kappa4
    # kappa4 = mu4 - 3*sigma^4
    # Is it determined by t3?
    t3_to_k4 = {}
    k4_by_t3 = True
    for d in data:
        t3 = d['t3']
        if t3 in t3_to_k4:
            if t3_to_k4[t3] != d['kappa4']:
                k4_by_t3 = False
        else:
            t3_to_k4[t3] = d['kappa4']

    print(f"  kappa_4 determined by t3: {k4_by_t3}")

    if not k4_by_t3:
        # Check (t3, t5)
        t3t5_to_k4 = {}
        k4_by_t3t5 = True
        for d in data:
            key = (d['t3'], d['t5'])
            if key in t3t5_to_k4:
                if t3t5_to_k4[key] != d['kappa4']:
                    k4_by_t3t5 = False
            else:
                t3t5_to_k4[key] = d['kappa4']
        print(f"  kappa_4 determined by (t3, t5): {k4_by_t3t5}")

# ============================================================
# WHY IS kappa_3 = 0?
# ============================================================
print("\n" + "=" * 70)
print("WHY kappa_3 = 0 (PROOF SKETCH)")
print("=" * 70)
print("""
kappa_3 = E[(fwd - mu)^3] = third central moment.

The palindrome symmetry says F_k = F_{n-1-k} for PALINDROMIC tournaments.
But F(T,x) is NOT palindromic in general (only for self-complementary T).

However, kappa_3 = 0 means the distribution is SYMMETRIC about its mean,
even though the distribution is not palindromic!

The resolution: fwd and (n-1)-fwd have the SAME distribution when we
average over all permutations. This is because:
  fwd(sigma) = number of edges i->j in T where i appears before j in sigma
  (n-1) - fwd(sigma) = number of edges j->i in T where i appears before j
                      = fwd(sigma) for T^op (the reverse tournament)

Wait, that gives fwd(sigma, T) + fwd(sigma, T^op) = n-1 for each sigma.
But this doesn't immediately give symmetry unless T = T^op.

Actually, let's think again. For the SAME sigma:
  fwd(sigma) = sum_{i<j in sigma} T[sigma_i][sigma_j]
  (n-1) - fwd(sigma) = sum_{i<j in sigma} T[sigma_j][sigma_i]
                      = fwd(sigma^{rev})

where sigma^{rev}(i) = sigma(n-1-i) is the reverse permutation.

So fwd(sigma) + fwd(sigma^{rev}) = n-1.
And sigma and sigma^{rev} have the same distribution (reversal is a bijection on S_n).
So fwd has the same distribution as (n-1) - fwd.
Therefore the distribution IS symmetric about (n-1)/2.
QED: kappa_3 = 0 (and all odd central moments are zero).
""")

print("WAIT — but this means E[fwd^3] is EXACTLY determined by")
print("E[fwd^2] and E[fwd] via the relation:")
print("  E[fwd^3] = 3*mu*E[fwd^2] - 2*mu^3")
print("(the third moment of a symmetric distribution)")
print()
print("And E[fwd^2] = Var + mu^2 = (n+1)/12 + 4t3/(n(n-1)) + ((n-1)/2)^2")
print("So E[fwd^3] = 3 * (n-1)/2 * [(n+1)/12 + 4t3/(n(n-1)) + ((n-1)/2)^2] - 2*((n-1)/2)^3")
print("             = A(n) + 3*(n-1)/2 * 4t3/(n(n-1))")
print("             = A(n) + 6t3/n")
print()
print("This proves E[fwd^3] = A(n) + 6t3/n ALGEBRAICALLY!")
print("The 6/n slope comes from 3 * (n-1)/2 * 4/(n(n-1)) = 6/n.")
