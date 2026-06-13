#!/usr/bin/env python3
"""
free_cumulant_connection.py - Explore free probability connections.

In free probability, the FREE cumulants of a random variable are defined via
the moment-cumulant formula using NON-CROSSING partitions (Speicher).

Classical cumulants use ALL set partitions; free cumulants use non-crossing ones.

For our forward-edge distribution:
  Classical cumulants: kappa_2, kappa_4, kappa_6 (computed above)
  Free cumulants: k_2, k_4, k_6 (to be computed)

The moment-free cumulant relation:
  mu_r = sum_{pi in NC(r)} prod_{B in pi} k_{|B|}

For the transitive tournament (Eulerian distribution), the FREE cumulants
might have a cleaner structure than classical ones.

Also: the Catalan numbers appear in free probability (counting NC partitions).
The Eulerian numbers have known connections to Catalan numbers.
Is there a free-probability interpretation of the forward-edge distribution?

Author: opus-2026-03-07-S46d
"""
from fractions import Fraction
from math import factorial, comb
from itertools import permutations

def compute_free_cumulants(moments, max_r=8):
    """
    Compute free cumulants from moments using Mobius inversion on NC lattice.

    The moment-free cumulant relation for centered moments:
    mu_1 = k_1
    mu_2 = k_2 + k_1^2
    mu_3 = k_3 + 3*k_2*k_1 + k_1^3
    mu_4 = k_4 + 2*k_2^2 + 4*k_3*k_1 + 6*k_2*k_1^2 + k_1^4

    For centered variables (k_1 = 0):
    mu_2 = k_2
    mu_3 = k_3
    mu_4 = k_4 + 2*k_2^2
    mu_5 = k_5 + 5*k_3*k_2
    mu_6 = k_6 + 6*k_4*k_2 + 3*k_3^2 + 5*k_2^3

    (coefficients are Catalan-related)
    """
    k = {}
    # Centered: moments[1] = 0
    k[1] = Fraction(0)
    k[2] = moments[2]
    k[3] = moments[3]  # = 0 by symmetry
    k[4] = moments[4] - 2*k[2]**2
    k[5] = moments[5] - 5*k[3]*k[2]  # = moments[5] since k[3]=0
    k[6] = moments[6] - 6*k[4]*k[2] - 3*k[3]**2 - 5*k[2]**3
    # = moments[6] - 6*k[4]*k[2] - 5*k[2]^3  (since k[3]=0)
    return k

# Compute for Eulerian distribution at various n
print("FREE CUMULANTS OF EULERIAN (TRANSITIVE) DISTRIBUTION")
print("=" * 60)

for n in range(4, 10):
    # Get Eulerian numbers
    A = [0]*n
    A[0] = 1
    for step in range(2, n+1):
        new_A = [0]*n
        for k_val in range(n):
            new_A[k_val] = (k_val+1) * A[k_val] + (step - k_val) * (A[k_val-1] if k_val > 0 else 0)
        A = new_A

    total = factorial(n)
    mu = Fraction(n-1, 2)

    # Central moments
    moments = {}
    for r in range(1, 9):
        moments[r] = sum(Fraction((k_val - mu)**r * A[k_val], total) for k_val in range(n))

    free_k = compute_free_cumulants(moments)

    print(f"\nn={n}: mu = {mu}")
    print(f"  k_2 = {free_k[2]} = {float(free_k[2]):.8f}")
    print(f"  k_4 = {free_k[4]} = {float(free_k[4]):.8f}")
    print(f"  k_6 = {free_k[6]} = {float(free_k[6]):.8f}")

    # Classical cumulants for comparison
    kappa_2 = moments[2]
    kappa_4 = moments[4] - 3*moments[2]**2
    kappa_6 = moments[6] - 15*moments[4]*moments[2] + 30*moments[2]**3

    print(f"  kappa_2 = {kappa_2}")
    print(f"  kappa_4 = {kappa_4}")
    print(f"  kappa_6 = {kappa_6}")

    # Ratios
    if free_k[2] != 0:
        r42 = free_k[4] / free_k[2]
        r62 = free_k[6] / free_k[2]
        print(f"  k_4/k_2 = {r42}")
        print(f"  k_6/k_2 = {r62}")

# Now compute free cumulants for non-transitive tournaments at n=5
print("\n" + "=" * 60)
print("FREE CUMULANTS FOR n=5 TOURNAMENTS")
print("=" * 60)

n = 5

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
    for i, j, k in __import__('itertools').combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t += 1
    return t

def count_5cycles(adj, n):
    t = 0
    for combo in __import__('itertools').combinations(range(n), 5):
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t += 1
    return t // 5

m_edges = n*(n-1)//2
seen = {}
data = []

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

    t3 = count_3cycles(adj, n)
    t5 = count_5cycles(adj, n)

    mu_val = Fraction(n-1, 2)
    moments = {}
    for r in range(1, 9):
        moments[r] = sum(Fraction((k_val - mu_val)**r * fwd_dist[k_val], total) for k_val in range(n))

    free_k = compute_free_cumulants(moments)

    # Classical
    kappa_2 = moments[2]
    kappa_4 = moments[4] - 3*moments[2]**2

    data.append({
        't3': t3, 't5': t5,
        'k2': free_k[2], 'k4': free_k[4], 'k6': free_k[6],
        'kappa2': kappa_2, 'kappa4': kappa_4
    })

print(f"\nn=5: {len(data)} F-classes")
print(f"{'t3':>3} {'t5':>3} {'k_2':>12} {'k_4':>12} {'k_6':>12} {'kappa_4':>12}")
for d in sorted(data, key=lambda x: (x['t3'], x['t5'])):
    print(f"{d['t3']:>3} {d['t5']:>3} {str(d['k2']):>12} {str(d['k4']):>12} {str(d['k6']):>12} {str(d['kappa4']):>12}")

# Check: is free k_4 simpler than classical kappa_4?
# k_4 = mu_4 - 2*k_2^2 = mu_4 - 2*Var^2
# kappa_4 = mu_4 - 3*Var^2
# So k_4 = kappa_4 + Var^2!

print("\nVerification: k_4 = kappa_4 + Var^2?")
for d in data:
    check = d['kappa4'] + d['k2']**2
    match = (check == d['k4'])
    if not match:
        print(f"  FAIL: t3={d['t3']}, t5={d['t5']}")
        break
else:
    print("  YES, k_4 = kappa_4 + Var^2 for all F-classes ✓")

# So: k_4 = -(n+1)/120 + Var^2 + 2*t5/C(n,4) + 4*a2/C(n,4) - 48*t3^2/(n(n-1))^2
# = -(n+1)/120 + ((n+1)/12 + 4*t3/(n(n-1)))^2 + 2*t5/C(n,4) - 48*t3^2/(n(n-1))^2
# The Var^2 expansion: Var^2 = ((n+1)/12)^2 + 2*(n+1)/12*4*t3/(n(n-1)) + (4*t3/(n(n-1)))^2
# = (n+1)^2/144 + 8(n+1)*t3/(12*n(n-1)) + 16*t3^2/(n(n-1))^2

# So k_4 = -(n+1)/120 + (n+1)^2/144 + 2(n+1)/(3n(n-1))*t3 + 16t3^2/(n(n-1))^2 + 2t5/C(n,4) - 48t3^2/(n(n-1))^2
# = -(n+1)/120 + (n+1)^2/144 + 2(n+1)*t3/(3n(n-1)) + (16-48)*t3^2/(n(n-1))^2 + 2t5/C(n,4)
# = -(n+1)/120 + (n+1)^2/144 + 2(n+1)*t3/(3n(n-1)) - 32*t3^2/(n(n-1))^2 + 2t5/C(n,4)

# Hmm, the free cumulant k_4 has a LINEAR t3 term, unlike the classical kappa_4!
# This means kappa_4 is SIMPLER than k_4 for our purposes.
# The classical cumulants are the natural framework for the tournament hierarchy.

print("\nFree k_4 has LINEAR t3 term (non-zero), unlike classical kappa_4!")
print("Classical cumulants are the natural framework for tournaments.")
