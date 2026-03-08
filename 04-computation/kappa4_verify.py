#!/usr/bin/env python3
"""
kappa4_verify.py — Verify E[fwd^4] formulas and compute kappa_4 from scratch.

Author: opus-2026-03-07-S46d
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict

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

def count_3cycles(adj, n):
    t = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t += 1
    return t

def count_5cycles(adj, n):
    """Count directed 5-cycles (each cycle counted once per direction)."""
    t = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t += 1
    return t // 5  # each directed cycle counted 5x (start vertex)

def count_alpha2(adj, n):
    """Count unordered pairs of vertex-disjoint 3-cycles."""
    cycles_3 = []
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            cycles_3.append(set(triple))
    count = 0
    for a in range(len(cycles_3)):
        for b in range(a+1, len(cycles_3)):
            if cycles_3[a].isdisjoint(cycles_3[b]):
                count += 1
    return count

for n in [5, 6]:
    print(f"\n{'='*60}")
    print(f"n = {n}")
    print(f"{'='*60}")

    m_edges = n*(n-1)//2
    seen = {}
    data = []

    for bits in range(1 << m_edges):
        adj = tournament_from_bits(n, bits)
        total = Fraction(1, n.__class__(1))  # will use factorial

        from math import factorial
        total = factorial(n)

        fwd_dist = [0]*n
        for P in permutations(range(n)):
            fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
            fwd_dist[fwd] += 1

        key = tuple(fwd_dist)
        if key in seen:
            continue
        seen[key] = True

        t3_val = count_3cycles(adj, n)
        t5_val = count_5cycles(adj, n)
        a2_val = count_alpha2(adj, n) if n >= 6 else 0

        mu = Fraction(n-1, 2)
        M1 = sum(Fraction(k * fwd_dist[k], total) for k in range(n))
        M2 = sum(Fraction(k**2 * fwd_dist[k], total) for k in range(n))
        M3 = sum(Fraction(k**3 * fwd_dist[k], total) for k in range(n))
        M4 = sum(Fraction(k**4 * fwd_dist[k], total) for k in range(n))

        Var = M2 - mu**2
        mu4 = M4 - 4*mu*M3 + 6*mu**2*M2 - 3*mu**4
        k4 = mu4 - 3*Var**2

        data.append({
            't3': t3_val, 't5': t5_val, 'a2': a2_val,
            'M4': M4, 'Var': Var, 'mu4': mu4, 'kappa4': k4
        })

    # Check E[fwd^4] formula from THM-092
    print(f"\n{len(data)} F-classes")
    print(f"\nVerifying E[fwd^4] formula:")
    for d in sorted(data, key=lambda x: (x['t3'], x['t5'], x['a2'])):
        if n == 5:
            predicted_M4 = Fraction(287, 10) + Fraction(27, 5)*d['t3'] + Fraction(2, 5)*d['t5']
        elif n == 6:
            predicted_M4 = Fraction(619, 10) + Fraction(82, 15)*d['t3'] + Fraction(2, 15)*d['t5'] + Fraction(4, 15)*d['a2']

        match = (predicted_M4 == d['M4'])
        if not match:
            print(f"  t3={d['t3']}, t5={d['t5']}, a2={d['a2']}: predicted={predicted_M4}, actual={d['M4']} {'✓' if match else 'MISMATCH'}")

    # If all match, say so
    all_match = all(
        (Fraction(287, 10) + Fraction(27, 5)*d['t3'] + Fraction(2, 5)*d['t5'] == d['M4'] if n == 5
         else Fraction(619, 10) + Fraction(82, 15)*d['t3'] + Fraction(2, 15)*d['t5'] + Fraction(4, 15)*d['a2'] == d['M4'])
        for d in data
    )
    if all_match:
        print(f"  ALL {len(data)} classes match E[fwd^4] formula ✓")

    # Now fit kappa_4 as function of invariants
    print(f"\nkappa_4 data:")
    print(f"{'t3':>3} {'t5':>3} {'a2':>3} {'kappa_4':>18} {'Var':>12}")
    for d in sorted(data, key=lambda x: (x['t3'], x['t5'], x['a2'])):
        print(f"{d['t3']:>3} {d['t5']:>3} {d['a2']:>3} {str(d['kappa4']):>18} {str(d['Var']):>12}")

    # Try fitting: kappa_4 = c0 + c1*t3 + c2*t5 + c3*a2 + c4*t3^2 + c5*t3*t5 + ...
    # First check: is kappa_4 determined by (t3, t5, a2)?
    inv_to_k4 = defaultdict(set)
    for d in data:
        key = (d['t3'], d['t5'], d['a2'])
        inv_to_k4[key].add(d['kappa4'])

    ambiguous = sum(1 for v in inv_to_k4.values() if len(v) > 1)
    print(f"\nkappa_4 determined by (t3, t5, alpha_2): {'YES' if ambiguous == 0 else f'NO ({ambiguous} ambiguous)'}")

    if ambiguous == 0 and n == 5:
        # Fit kappa_4 = a + b*t3 + c*t5 + d*t3^2
        from fractions import Fraction
        # Use least squares with Fraction
        # Build system: for each data point, kappa_4 = a + b*t3 + c*t5 + d*t3^2
        rows = [(Fraction(1), Fraction(d['t3']), Fraction(d['t5']), Fraction(d['t3']**2), d['kappa4']) for d in data]

        # Use sympy for exact solve
        from sympy import Matrix, symbols, Rational
        A_mat = Matrix([[r[0], r[1], r[2], r[3]] for r in rows])
        b_vec = Matrix([r[4] for r in rows])

        # Solve A^T A x = A^T b
        ATA = A_mat.T * A_mat
        ATb = A_mat.T * b_vec
        x = ATA.solve(ATb)
        print(f"\nFitted kappa_4 = {x[0]} + ({x[1]})*t3 + ({x[2]})*t5 + ({x[3]})*t3^2")

        # Verify
        for d in data:
            pred = x[0] + x[1]*d['t3'] + x[2]*d['t5'] + x[3]*d['t3']**2
            if pred != d['kappa4']:
                print(f"  RESIDUAL: t3={d['t3']}, t5={d['t5']}: pred={pred}, actual={d['kappa4']}")

    if ambiguous == 0 and n == 6:
        # Fit kappa_4 = a + b*t3 + c*t5 + d*a2 + e*t3^2
        rows = [(Fraction(1), Fraction(d['t3']), Fraction(d['t5']), Fraction(d['a2']),
                 Fraction(d['t3']**2), d['kappa4']) for d in data]

        from sympy import Matrix
        A_mat = Matrix([[r[0], r[1], r[2], r[3], r[4]] for r in rows])
        b_vec = Matrix([r[5] for r in rows])

        ATA = A_mat.T * A_mat
        ATb = A_mat.T * b_vec
        x = ATA.solve(ATb)
        print(f"\nFitted kappa_4 = {x[0]} + ({x[1]})*t3 + ({x[2]})*t5 + ({x[3]})*a2 + ({x[4]})*t3^2")

        # Verify
        residuals = 0
        for d in data:
            pred = x[0] + x[1]*d['t3'] + x[2]*d['t5'] + x[3]*d['a2'] + x[4]*d['t3']**2
            if pred != d['kappa4']:
                residuals += 1
        if residuals == 0:
            print(f"  ALL {len(data)} classes match ✓")
        else:
            print(f"  {residuals} residuals — need more terms")
