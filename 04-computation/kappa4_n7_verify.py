#!/usr/bin/env python3
"""
kappa4_n7_verify.py — Verify kappa_4 structure at n=7 via sampling.

Predicted:
  E[fwd^4] t3 coefficient = 116/21
  kappa_4 constant (t3=t5=a2=0) = -1/15

Author: opus-2026-03-07-S46d
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict
import random
import sys

n = 7

def tournament_from_bits(bits):
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

def count_3cycles(adj):
    t = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t += 1
    return t

def count_5cycles(adj):
    """Count directed 5-cycles (each counted once per direction)."""
    t = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t += 1
    return t // 5

def count_alpha2(adj):
    cycles_3 = []
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            cycles_3.append(frozenset(triple))
    count = 0
    for a in range(len(cycles_3)):
        for b in range(a+1, len(cycles_3)):
            if cycles_3[a].isdisjoint(cycles_3[b]):
                count += 1
    return count

def count_7cycles(adj):
    """Count directed 7-cycles (each counted once per direction)."""
    # All 7 vertices used (n=7), so just check all permutations
    t = 0
    for perm in permutations(range(n)):
        if all(adj[perm[i]][perm[(i+1)%n]] for i in range(n)):
            t += 1
    return t // n

# Use sampling approach with small perm sample for speed
ALL_PERMS = list(permutations(range(n)))
TOTAL = len(ALL_PERMS)  # 5040

m_edges = n*(n-1)//2
random.seed(42)
seen_F = {}
data = []

print(f"n={n}: Sampling tournaments...", flush=True)

# Include transitive tournament first
adj_trans = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
F_trans = [0]*n
for P in ALL_PERMS:
    fwd = sum(1 for i in range(n-1) if adj_trans[P[i]][P[i+1]])
    F_trans[fwd] += 1

t3 = 0; t5 = 0; a2 = 0; t7 = count_7cycles(adj_trans)
mu = Fraction(n-1, 2)
M4 = sum(Fraction(k**4 * F_trans[k], TOTAL) for k in range(n))
M3 = sum(Fraction(k**3 * F_trans[k], TOTAL) for k in range(n))
M2 = sum(Fraction(k**2 * F_trans[k], TOTAL) for k in range(n))
Var = M2 - mu**2
mu4_val = M4 - 4*mu*M3 + 6*mu**2*M2 - 3*mu**4
k4_val = mu4_val - 3*Var**2

print(f"Transitive: t3=0, t5=0, a2=0, t7={t7}")
print(f"  M4 = {M4}, kappa_4 = {k4_val}")
print(f"  Expected: -(n+1)/120 = {Fraction(-(n+1),120)}")
print(f"  Match: {k4_val == Fraction(-(n+1),120)}")

data.append({'t3': 0, 't5': 0, 'a2': 0, 't7': t7, 'M4': M4, 'kappa4': k4_val, 'Var': Var})
seen_F[tuple(F_trans)] = True

for trial in range(3000):
    bits = random.getrandbits(m_edges)
    adj = tournament_from_bits(bits)

    F = [0]*n
    for P in ALL_PERMS:
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1

    key = tuple(F)
    if key in seen_F:
        continue
    seen_F[key] = True

    t3 = count_3cycles(adj)
    t5 = count_5cycles(adj)
    a2 = count_alpha2(adj)
    t7 = count_7cycles(adj)

    mu = Fraction(n-1, 2)
    M4 = sum(Fraction(k**4 * F[k], TOTAL) for k in range(n))
    M3 = sum(Fraction(k**3 * F[k], TOTAL) for k in range(n))
    M2 = sum(Fraction(k**2 * F[k], TOTAL) for k in range(n))
    Var = M2 - mu**2
    mu4_val = M4 - 4*mu*M3 + 6*mu**2*M2 - 3*mu**4
    k4_val = mu4_val - 3*Var**2

    data.append({'t3': t3, 't5': t5, 'a2': a2, 't7': t7, 'M4': M4, 'kappa4': k4_val, 'Var': Var})

    if len(data) % 20 == 0:
        print(f"  {len(data)} F-classes found (trial {trial})", file=sys.stderr, flush=True)

print(f"\nTotal: {len(data)} distinct F-classes", flush=True)

# Check: is kappa_4 determined by (t3, t5, alpha_2)?
inv_to_k4 = defaultdict(set)
for d in data:
    key = (d['t3'], d['t5'], d['a2'])
    inv_to_k4[key].add(d['kappa4'])

ambiguous_3 = sum(1 for v in inv_to_k4.values() if len(v) > 1)
print(f"kappa_4 determined by (t3, t5, a2): {'YES' if ambiguous_3 == 0 else f'NO ({ambiguous_3} ambiguous)'}")

# Check: is kappa_4 determined by (t3, t5, alpha_2, t7)?
inv_to_k4_7 = defaultdict(set)
for d in data:
    key = (d['t3'], d['t5'], d['a2'], d['t7'])
    inv_to_k4_7[key].add(d['kappa4'])

ambiguous_7 = sum(1 for v in inv_to_k4_7.values() if len(v) > 1)
print(f"kappa_4 determined by (t3, t5, a2, t7): {'YES' if ambiguous_7 == 0 else f'NO ({ambiguous_7} ambiguous)'}")

# Check: is M4 determined by (t3, t5, alpha_2)?
inv_to_M4 = defaultdict(set)
for d in data:
    key = (d['t3'], d['t5'], d['a2'])
    inv_to_M4[key].add(d['M4'])

ambiguous_M4 = sum(1 for v in inv_to_M4.values() if len(v) > 1)
print(f"M4 determined by (t3, t5, a2): {'YES' if ambiguous_M4 == 0 else f'NO ({ambiguous_M4} ambiguous)'}")

# Check: is M4 determined by (t3, t5, alpha_2, t7)?
inv_to_M4_7 = defaultdict(set)
for d in data:
    key = (d['t3'], d['t5'], d['a2'], d['t7'])
    inv_to_M4_7[key].add(d['M4'])

ambiguous_M4_7 = sum(1 for v in inv_to_M4_7.values() if len(v) > 1)
print(f"M4 determined by (t3, t5, a2, t7): {'YES' if ambiguous_M4_7 == 0 else f'NO ({ambiguous_M4_7} ambiguous)'}")

# If M4 is determined by (t3, t5, a2), fit it
if ambiguous_M4 == 0:
    # Try linear fit: M4 = a + b*t3 + c*t5 + d*a2
    from sympy import Matrix, Rational
    rows = [(1, d['t3'], d['t5'], d['a2'], d['M4']) for d in data]
    A_mat = Matrix([[r[0], r[1], r[2], r[3]] for r in rows])
    b_vec = Matrix([r[4] for r in rows])
    ATA = A_mat.T * A_mat
    ATb = A_mat.T * b_vec
    x = ATA.solve(ATb)

    # Check residuals
    residuals = 0
    for d in data:
        pred = x[0] + x[1]*d['t3'] + x[2]*d['t5'] + x[3]*d['a2']
        if pred != d['M4']:
            residuals += 1

    print(f"\nLinear fit M4 = {x[0]} + ({x[1]})*t3 + ({x[2]})*t5 + ({x[3]})*a2")
    print(f"  Residuals: {residuals}/{len(data)}")

    if residuals == 0:
        print(f"  Predicted t3 coeff = 116/21 = {float(Fraction(116,21)):.6f}, actual = {float(x[1]):.6f}")
        print(f"  Match 116/21: {x[1] == Fraction(116,21)}")

# Also fit kappa_4
if ambiguous_3 == 0:
    # kappa_4 = c0 + c1*t5 + c2*a2 + c3*t3^2  (no linear t3)
    from sympy import Matrix
    rows = [(1, d['t5'], d['a2'], d['t3']**2, d['kappa4']) for d in data]
    A_mat = Matrix([[r[0], r[1], r[2], r[3]] for r in rows])
    b_vec = Matrix([r[4] for r in rows])
    ATA = A_mat.T * A_mat
    ATb = A_mat.T * b_vec
    x = ATA.solve(ATb)

    residuals = 0
    for d in data:
        pred = x[0] + x[1]*d['t5'] + x[2]*d['a2'] + x[3]*d['t3']**2
        if pred != d['kappa4']:
            residuals += 1

    print(f"\nFitted kappa_4 = {x[0]} + ({x[1]})*t5 + ({x[2]})*a2 + ({x[3]})*t3^2")
    print(f"  Residuals: {residuals}/{len(data)}")
    print(f"  Expected constant = -1/15, actual = {x[0]}, match = {x[0] == Fraction(-1,15)}")
    print(f"  Expected t3^2 coeff = -48/(7*6)^2 = {Fraction(-48, 42**2)}, actual = {x[3]}")
    print(f"  Match: {x[3] == Fraction(-48, 42**2)}")
