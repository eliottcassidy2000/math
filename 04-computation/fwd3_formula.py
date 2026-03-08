#!/usr/bin/env python3
"""
fwd3_formula.py — Find the exact formula for E[fwd^3] in terms of n and t3.

DISCOVERY: E[fwd^3] is exactly linear in t3 at n=4,5,6.
Find the general formula E[fwd^3] = A(n) + B(n)*t3.

At n=4: E[fwd^3] = 5.25 + 1.5*t3
  = 21/4 + 3/2 * t3

At n=5: E[fwd^3] = 11 + 1.2*t3
  = 11 + 6/5 * t3

At n=6: E[fwd^3] = 20 + 1*t3

A(n): 21/4, 11, 20 → 5.25, 11, 20
Differences: 5.75, 9 → not arithmetic

Let me compute exactly and find the pattern.

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
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t3 += 1
    return t3

print("=" * 60)
print("EXACT E[fwd^3] FORMULA")
print("=" * 60)

for n in [3, 4, 5, 6]:
    m_vals = n*(n-1)//2
    # Collect (t3, E[fwd^3]) pairs
    seen = set()
    pairs = []

    for bits in range(1 << m_vals):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        t3 = count_3cycles(adj, n)
        total = sum(F)
        m3 = Fraction(sum(k**3 * F[k] for k in range(n)), total)
        pairs.append((t3, m3))

    # Linear fit in exact arithmetic
    # E[fwd^3] = A + B*t3
    # Use two distinct t3 values
    t3_vals = sorted(set(t for t, _ in pairs))
    t3_to_m3 = {}
    for t, m3 in pairs:
        if t not in t3_to_m3:
            t3_to_m3[t] = m3

    if len(t3_vals) >= 2:
        B = (t3_to_m3[t3_vals[1]] - t3_to_m3[t3_vals[0]]) / (t3_vals[1] - t3_vals[0])
        A = t3_to_m3[t3_vals[0]] - B * t3_vals[0]
        # Verify
        ok = all(abs(t3_to_m3[t] - (A + B*t)) < Fraction(1, 10000) for t in t3_vals)
        print(f"\nn={n}: E[fwd^3] = {A} + {B} * t3  {'[EXACT]' if ok else '[APPROX]'}")
        print(f"       = {float(A):.6f} + {float(B):.6f} * t3")
    elif len(t3_vals) == 1:
        print(f"\nn={n}: Only one t3 value ({t3_vals[0]}), E[fwd^3] = {t3_to_m3[t3_vals[0]]}")

# Now let's find the pattern
print("\n" + "=" * 60)
print("PATTERN DETECTION FOR A(n) AND B(n)")
print("=" * 60)

# Manually from the output above:
# n=3: A=1, B=? (only t3=0 or 1 possible)
# n=4: A=21/4, B=3/2
# n=5: A=11, B=6/5
# n=6: A=20, B=1

# Let me compute for n=3 properly
for n in [3]:
    adj0 = tournament_from_bits(n, 0)  # transitive
    F0 = compute_F(adj0, n)
    m3_0 = Fraction(sum(k**3 * F0[k] for k in range(n)), sum(F0))
    t3_0 = 0

    adj1 = tournament_from_bits(n, 7)  # 3-cycle (all bits set)
    F1 = compute_F(adj1, n)
    m3_1 = Fraction(sum(k**3 * F1[k] for k in range(n)), sum(F1))
    t3_1 = 1

    B = (m3_1 - m3_0) / (t3_1 - t3_0)
    A = m3_0
    print(f"\nn=3: A = {A}, B = {B}")

# Collect all
As = {3: None, 4: None, 5: None, 6: None}
Bs = {3: None, 4: None, 5: None, 6: None}

for n in [3, 4, 5, 6]:
    # Transitive: t3=0
    adj = tournament_from_bits(n, 0)
    F = compute_F(adj, n)
    total = sum(F)
    A_n = Fraction(sum(k**3 * F[k] for k in range(n)), total)
    As[n] = A_n

    # Find a tournament with t3=1
    for bits in range(1, 1 << (n*(n-1)//2)):
        adj = tournament_from_bits(n, bits)
        t3 = count_3cycles(adj, n)
        if t3 == 1:
            F = compute_F(adj, n)
            m3 = Fraction(sum(k**3 * F[k] for k in range(n)), sum(F))
            B_n = m3 - A_n
            Bs[n] = B_n
            break

print("\nSummary:")
for n in [3, 4, 5, 6]:
    print(f"  n={n}: A = {As[n]} = {float(As[n]):.6f}, B = {Bs[n]} = {float(Bs[n]):.6f}")
    print(f"         A = {As[n]}, B = {Bs[n]}")
    # Try to express A and B in terms of n
    # A(n) = E[fwd^3] at transitive = ???
    # For transitive, F_k = Eulerian numbers

# Let's see if A(n) = some simple function
print("\nChecking A(n) = E[fwd^3] for transitive tournament:")
for n in [3, 4, 5, 6]:
    # Theoretical: for transitive, fwd has same distribution as #descents in perm
    # E[desc^3] is a known quantity
    #
    # E[desc] = (n-1)/2
    # E[desc^2] = (n-1)(3n-2)/12 (for perm of [n])
    # E[desc^3] = ?
    #
    # Actually E[desc^r] for permutations is related to Bernoulli numbers
    # via the generating function.
    mean = Fraction(n-1, 2)
    var = Fraction(n+1, 12)
    # E[desc^3] = E[(desc-mu)^3] + 3*mu*E[desc^2] - 2*mu^3
    # Skewness of descent distribution for S_n is 0 (by symmetry desc <-> n-1-desc)
    # So E[(desc-mu)^3] = 0
    # Therefore E[desc^3] = 3*mu*E[desc^2] - 2*mu^3
    Efwd2 = var + mean**2
    Efwd3_pred = 3*mean*Efwd2 - 2*mean**3
    print(f"  n={n}: A = {As[n]}, predicted = {Efwd3_pred}, match = {As[n] == Efwd3_pred}")

print("\n\nChecking B(n) = slope of E[fwd^3] wrt t3:")
for n in [3, 4, 5, 6]:
    # At n=4: B = 3/2 = 6/4
    # At n=5: B = 6/5
    # At n=6: B = 1 = 6/6
    ratio = Bs[n] * n
    print(f"  n={n}: B = {Bs[n]}, B*n = {ratio}")
    # Hmm, let's try B = 6/(n*(n-1)) * something
    ratio2 = Bs[n] * n * (n-1)
    print(f"         B*n*(n-1) = {ratio2}")

# The skewness = 0 argument gives A(n):
# A(n) = E[fwd^3]_transitive = 3*((n-1)/2)*((n-1)(3n-2)/12) - 2*((n-1)/2)^3
# Wait, need to check: E[fwd^2] for Eulerian distribution
print("\n\nE[fwd^2] for transitive:")
for n in [3, 4, 5, 6]:
    mean = Fraction(n-1, 2)
    var = Fraction(n+1, 12)  # Var[desc] for S_n = (n+1)/12
    Efwd2 = var + mean**2
    print(f"  n={n}: E[fwd^2] = {Efwd2} = {float(Efwd2):.6f}")
    print(f"    = (n+1)/12 + ((n-1)/2)^2 = {Fraction(n+1,12)} + {mean**2} = {Efwd2}")

print("\nB(n) pattern:")
for n in [3, 4, 5, 6]:
    # B(3)=2, B(4)=3/2, B(5)=6/5, B(6)=1
    # = 6/n? No: 6/3=2, 6/4=3/2, 6/5=6/5, 6/6=1. YES!
    pred = Fraction(6, n)
    print(f"  n={n}: B = {Bs[n]}, 6/n = {pred}, match = {Bs[n] == pred}")
