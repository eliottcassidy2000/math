#!/usr/bin/env python3
"""
kappa4_general.py — Derive kappa_4 formula at general n.

Strategy:
kappa_4 = mu_4 - 3*mu_2^2  where mu_k = E[(fwd - mean)^k] = central moment

We know:
  mean = (n-1)/2
  Var = sigma^2 = (n+1)/12 + 4*t3/(n*(n-1))  [THM-089]
  E[fwd^3] = A(n) + 6*t3/n  [THM-090]
  E[fwd^4] = known at n=5,6

By symmetry (THM-091), all odd central moments vanish.
So kappa_4 = mu_4 - 3*sigma^4

mu_4 = E[(fwd - mu)^4] = E[fwd^4] - 4*mu*E[fwd^3] + 6*mu^2*E[fwd^2] - 4*mu^3*E[fwd] + mu^4
     = M4 - 4*mu*M3 + 6*mu^2*M2 - 3*mu^4  (since M1=mu, so last two terms combine)

Actually let's be careful:
mu_4 = E[fwd^4] - 4*mu*E[fwd^3] + 6*mu^2*E[fwd^2] - 4*mu^3*E[fwd] + mu^4
     = M4 - 4*mu*M3 + 6*mu^2*M2 - 4*mu^3*mu + mu^4
     = M4 - 4*mu*M3 + 6*mu^2*M2 - 3*mu^4

kappa_4 = mu_4 - 3*sigma^4 = mu_4 - 3*(M2 - mu^2)^2

Let's compute symbolically for n=5 and n=6 to find the pattern.

Author: opus-2026-03-07-S46d
"""
from fractions import Fraction
from sympy import symbols, simplify, collect, Rational, expand, factor

# Symbolic computation
t3, t5, a2, n_sym = symbols('t3 t5 alpha2 n', positive=True)

def compute_kappa4(n_val):
    """Compute kappa_4 symbolically for given n."""
    mu = Rational(n_val - 1, 2)

    # Variance formula: Var = (n+1)/12 + 4*t3/(n*(n-1))  [THM-089]
    Var = Rational(n_val + 1, 12) + 4*t3 / (n_val * (n_val - 1))

    # M2 = Var + mu^2
    M2 = Var + mu**2

    # M3 = A(n) + 6*t3/n, where A(n) = 3*mu*M2_transitive - 2*mu^3
    # For transitive tournament: t3=0, so Var_trans = (n+1)/12
    M2_trans = Rational(n_val + 1, 12) + mu**2
    A_n = 3*mu*M2_trans - 2*mu**3
    M3 = A_n + 6*t3 / n_val

    # M4: use exact formulas
    if n_val == 5:
        M4 = Rational(287, 10) + Rational(27, 5)*t3 + Rational(2, 5)*t5
    elif n_val == 6:
        M4 = Rational(619, 10) + Rational(82, 15)*t3 + Rational(2, 15)*t5 + Rational(4, 15)*a2
    else:
        return None

    # Central fourth moment
    mu4 = expand(M4 - 4*mu*M3 + 6*mu**2*M2 - 3*mu**4)

    # kappa_4 = mu4 - 3*Var^2
    k4 = expand(mu4 - 3*Var**2)

    return {
        'mu': mu, 'Var': Var, 'M2': M2, 'M3': M3, 'M4': M4,
        'mu4': mu4, 'kappa4': k4,
        'kappa4_simplified': simplify(k4),
        'kappa4_collected': collect(expand(k4), [t3, t5, a2])
    }

print("=" * 60)
print("KAPPA_4 FORMULA DERIVATION")
print("=" * 60)

for n_val in [5, 6]:
    print(f"\n{'='*40}")
    print(f"n = {n_val}")
    print(f"{'='*40}")

    result = compute_kappa4(n_val)

    print(f"  mu = {result['mu']}")
    print(f"  Var = {result['Var']}")
    print(f"  M2 = {result['M2']}")
    print(f"  M3 = {result['M3']}")
    print(f"  M4 = {result['M4']}")
    print(f"  mu4 = {result['mu4']}")
    print(f"  kappa_4 = {result['kappa4_collected']}")
    print(f"  kappa_4 (simplified) = {result['kappa4_simplified']}")

# Now let's try to find the PATTERN
print("\n" + "=" * 60)
print("PATTERN ANALYSIS")
print("=" * 60)

# For n=5: kappa_4 should be f(t3, t5)
# For n=6: kappa_4 should be f(t3, t5, alpha_2)
# Let's extract coefficients

for n_val in [5, 6]:
    result = compute_kappa4(n_val)
    k4 = expand(result['kappa4'])

    # Extract coefficient of each variable
    from sympy import Poly
    if n_val == 5:
        # k4 is polynomial in t3, t5
        coeff_const = k4.subs([(t3, 0), (t5, 0)])
        coeff_t3 = k4.coeff(t3, 1).subs([(t5, 0)])
        coeff_t5 = k4.coeff(t5, 1).subs([(t3, 0)])
        coeff_t3_sq = k4.coeff(t3, 2)
        print(f"\nn={n_val}:")
        print(f"  kappa_4 = {coeff_const} + ({coeff_t3})*t3 + ({coeff_t5})*t5 + ({coeff_t3_sq})*t3^2")
    elif n_val == 6:
        coeff_const = k4.subs([(t3, 0), (t5, 0), (a2, 0)])
        coeff_t3 = k4.coeff(t3, 1).subs([(t5, 0), (a2, 0)])
        coeff_t5 = k4.coeff(t5, 1).subs([(t3, 0), (a2, 0)])
        coeff_a2 = k4.coeff(a2, 1).subs([(t3, 0), (t5, 0)])
        coeff_t3_sq = k4.coeff(t3, 2).subs([(t5, 0), (a2, 0)])
        print(f"\nn={n_val}:")
        print(f"  kappa_4 = {coeff_const} + ({coeff_t3})*t3 + ({coeff_t5})*t5 + ({coeff_a2})*alpha_2 + ({coeff_t3_sq})*t3^2")

# Now let me also verify numerically against known E[fwd^4] data
print("\n" + "=" * 60)
print("NUMERICAL VERIFICATION AT n=5")
print("=" * 60)

from itertools import permutations, combinations
from math import comb, factorial

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
    t = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t += 1
    # Each 5-cycle counted 5 times (cyclic rotations) * 2 (two directions)
    # Actually each directed 5-cycle counted 5 times (start vertex)
    return t // 5  # = number of directed 5-cycles (both orientations)

def count_alpha2(adj, n):
    """Count vertex-disjoint 3-cycle pairs."""
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

n = 5
m_edges = n*(n-1)//2
seen = {}
all_data = []

for bits in range(1 << m_edges):
    adj = tournament_from_bits(n, bits)

    # Compute fwd distribution
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

    mu_val = Fraction(n-1, 2)
    M1 = sum(Fraction(k * fwd_dist[k], total) for k in range(n))
    M2_val = sum(Fraction(k**2 * fwd_dist[k], total) for k in range(n))
    M3_val = sum(Fraction(k**3 * fwd_dist[k], total) for k in range(n))
    M4_val = sum(Fraction(k**4 * fwd_dist[k], total) for k in range(n))

    Var_val = M2_val - mu_val**2
    mu4_val = M4_val - 4*mu_val*M3_val + 6*mu_val**2*M2_val - 3*mu_val**4
    k4_val = mu4_val - 3*Var_val**2

    all_data.append({
        't3': t3_val, 't5': t5_val,
        'M4': M4_val, 'Var': Var_val,
        'mu4': mu4_val, 'kappa4': k4_val
    })

print(f"\nn=5: {len(all_data)} F-classes")
print(f"{'t3':>3} {'t5':>3} {'Var':>12} {'kappa_4':>15}")
for d in sorted(all_data, key=lambda x: (x['t3'], x['t5'])):
    print(f"{d['t3']:>3} {d['t5']:>3} {str(d['Var']):>12} {str(d['kappa4']):>15}")

# Verify formula
print("\nVerification of symbolic formula:")
result5 = compute_kappa4(5)
k4_formula = result5['kappa4']
for d in sorted(all_data, key=lambda x: (x['t3'], x['t5'])):
    predicted = k4_formula.subs([(t3, d['t3']), (t5, d['t5'])])
    match = (Fraction(predicted) == d['kappa4'])
    if not match:
        print(f"  t3={d['t3']}, t5={d['t5']}: formula={predicted}, actual={d['kappa4']} MISMATCH")
    else:
        print(f"  t3={d['t3']}, t5={d['t5']}: {d['kappa4']} ✓")

# n=6 verification
print("\n" + "=" * 60)
print("NUMERICAL VERIFICATION AT n=6")
print("=" * 60)

n = 6
m_edges = n*(n-1)//2
seen = {}
all_data_6 = []

for bits in range(1 << m_edges):
    adj = tournament_from_bits(n, bits)

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
    a2_val = count_alpha2(adj, n)

    mu_val = Fraction(n-1, 2)
    M2_val = sum(Fraction(k**2 * fwd_dist[k], total) for k in range(n))
    M3_val = sum(Fraction(k**3 * fwd_dist[k], total) for k in range(n))
    M4_val = sum(Fraction(k**4 * fwd_dist[k], total) for k in range(n))

    Var_val = M2_val - mu_val**2
    mu4_val = M4_val - 4*mu_val*M3_val + 6*mu_val**2*M2_val - 3*mu_val**4
    k4_val = mu4_val - 3*Var_val**2

    all_data_6.append({
        't3': t3_val, 't5': t5_val, 'a2': a2_val,
        'M4': M4_val, 'Var': Var_val,
        'mu4': mu4_val, 'kappa4': k4_val
    })

print(f"\nn=6: {len(all_data_6)} F-classes")
print(f"{'t3':>3} {'t5':>3} {'a2':>3} {'kappa_4':>15}")
for d in sorted(all_data_6, key=lambda x: (x['t3'], x['t5'], x['a2'])):
    print(f"{d['t3']:>3} {d['t5']:>3} {d['a2']:>3} {str(d['kappa4']):>15}")

# Verify formula
print("\nVerification of symbolic formula:")
result6 = compute_kappa4(6)
k4_formula_6 = result6['kappa4']
mismatches = 0
for d in sorted(all_data_6, key=lambda x: (x['t3'], x['t5'], x['a2'])):
    predicted = k4_formula_6.subs([(t3, d['t3']), (t5, d['t5']), (a2, d['a2'])])
    match = (Fraction(predicted) == d['kappa4'])
    if not match:
        mismatches += 1
        if mismatches <= 5:
            print(f"  t3={d['t3']}, t5={d['t5']}, a2={d['a2']}: formula={predicted}, actual={d['kappa4']} MISMATCH")

if mismatches == 0:
    print(f"  ALL {len(all_data_6)} classes MATCH ✓")
else:
    print(f"  {mismatches}/{len(all_data_6)} MISMATCHES")
