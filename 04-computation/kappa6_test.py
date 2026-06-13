#!/usr/bin/env python3
"""
kappa6_test.py — Test whether kappa_6 introduces t7 at n=7.

At n=7, we have all 7-cycles. If the hierarchy prediction holds,
kappa_6 should depend on (t3, t5, alpha_2, t7, ...) where t7 is new.

kappa_6 = mu_6 - 15*mu_4*mu_2 - 10*mu_3^2 + 30*mu_2^3
(since mu_3 = 0 by symmetry, simplifies to)
kappa_6 = mu_6 - 15*mu_4*mu_2 + 30*mu_2^3

Author: opus-2026-03-07-S46d
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict
import random
import sys

n = 7
ALL_PERMS = list(permutations(range(n)))
TOTAL = len(ALL_PERMS)

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
    t = 0
    for perm in ALL_PERMS:
        if all(adj[perm[i]][perm[(i+1)%n]] for i in range(n)):
            t += 1
    return t // n

m_edges = n*(n-1)//2
random.seed(42)
seen_F = {}
data = []

# Include transitive
adj_trans = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
F_trans = [0]*n
for P in ALL_PERMS:
    fwd = sum(1 for i in range(n-1) if adj_trans[P[i]][P[i+1]])
    F_trans[fwd] += 1

seen_F[tuple(F_trans)] = True
mu = Fraction(n-1, 2)
moments = {}
for r in range(9):
    moments[r] = sum(Fraction((k - mu)**r * F_trans[k], TOTAL) for k in range(n))

k6 = moments[6] - 15*moments[4]*moments[2] + 30*moments[2]**3
data.append({
    't3': 0, 't5': 0, 'a2': 0, 't7': count_7cycles(adj_trans),
    'kappa6': k6, 'kappa4': moments[4] - 3*moments[2]**2,
    'kappa2': moments[2]
})
print(f"Transitive: kappa_6 = {k6} = {float(k6):.8f}")

# Sample random tournaments
for trial in range(2000):
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

    moments = {}
    for r in range(9):
        moments[r] = sum(Fraction((k - mu)**r * F[k], TOTAL) for k in range(n))

    k6 = moments[6] - 15*moments[4]*moments[2] + 30*moments[2]**3
    k4 = moments[4] - 3*moments[2]**2
    k2 = moments[2]

    data.append({
        't3': t3, 't5': t5, 'a2': a2, 't7': t7,
        'kappa6': k6, 'kappa4': k4, 'kappa2': k2
    })

    if len(data) % 20 == 0:
        print(f"  {len(data)} F-classes (trial {trial})", file=sys.stderr, flush=True)

print(f"\nTotal: {len(data)} F-classes")

# Is kappa_6 determined by (t3, t5, a2)?
inv_to_k6 = defaultdict(set)
for d in data:
    key = (d['t3'], d['t5'], d['a2'])
    inv_to_k6[key].add(d['kappa6'])
amb3 = sum(1 for v in inv_to_k6.values() if len(v) > 1)
print(f"kappa_6 determined by (t3, t5, a2): {'YES' if amb3 == 0 else f'NO ({amb3} ambiguous)'}")

# Is kappa_6 determined by (t3, t5, a2, t7)?
inv_to_k6_7 = defaultdict(set)
for d in data:
    key = (d['t3'], d['t5'], d['a2'], d['t7'])
    inv_to_k6_7[key].add(d['kappa6'])
amb7 = sum(1 for v in inv_to_k6_7.values() if len(v) > 1)
print(f"kappa_6 determined by (t3, t5, a2, t7): {'YES' if amb7 == 0 else f'NO ({amb7} ambiguous)'}")

if amb7 == 0 and amb3 > 0:
    print("\nkappa_6 REQUIRES t7 (7-cycle count) — confirms hierarchy prediction!")

    # Fit: kappa_6 = c0 + c1*t5 + c2*a2 + c3*t7 + c4*t3^2 + c5*t3*t5 + c6*t3^3
    # Try simpler: c0 + c1*t5 + c2*a2 + c3*t7 + nonlinear t3 terms
    # Actually, by analogy with kappa_4 structure:
    # kappa_6 = kappa_6(trans) + new_7vertex_part + nonlinear_lower
    # where nonlinear_lower involves products of lower corrections

    # Let's first check what invariants are needed beyond (t3,t5,a2)
    # Show an ambiguous case
    for key, vals in inv_to_k6.items():
        if len(vals) > 1:
            t3, t5, a2 = key
            # Find corresponding t7 values
            matching = [d for d in data if d['t3']==t3 and d['t5']==t5 and d['a2']==a2]
            print(f"\n  Ambiguous: (t3,t5,a2) = ({t3},{t5},{a2})")
            for d in matching:
                print(f"    t7={d['t7']}: kappa_6 = {float(d['kappa6']):.8f}")
            break

if amb7 == 0:
    # Fit kappa_6
    # Try: kappa_6 = c0 + c1*t7 + c2*t5*t3 + c3*a2*t3 + c4*t3^3 + c5*t3*t5 + c6*t5^2
    # Actually first try the simplest: c0 + c1*t7 + terms involving t3,t5,a2
    from sympy import Matrix
    # Basis: 1, t3, t5, a2, t7, t3^2, t3*t5, t3*a2, t3^3
    rows = []
    for d in data:
        t3, t5, a2, t7 = d['t3'], d['t5'], d['a2'], d['t7']
        row = [1, t7, t3**2*t5, t3**2*a2, t3**3, t3*t5, t3*a2, t5*a2]
        rows.append(row + [d['kappa6']])

    A_mat = Matrix([[r[i] for i in range(len(rows[0])-1)] for r in rows])
    b_vec = Matrix([r[-1] for r in rows])
    ATA = A_mat.T * A_mat
    ATb = A_mat.T * b_vec
    try:
        x = ATA.solve(ATb)
        residuals = 0
        for r in rows:
            pred = sum(x[i]*r[i] for i in range(len(x)))
            if pred != r[-1]:
                residuals += 1
        print(f"\nFit with basis [1, t7, t3^2*t5, t3^2*a2, t3^3, t3*t5, t3*a2, t5*a2]:")
        print(f"  Coefficients: {[str(c) for c in x]}")
        print(f"  Residuals: {residuals}/{len(data)}")
    except Exception as e:
        print(f"Fit failed: {e}")
