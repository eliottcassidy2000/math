#!/usr/bin/env python3
"""
signed_F_analysis.py — The signed forward-edge polynomial SF(T,x).

DISCOVERY (S46): SF(T,x) = sum_{sigma in S_n} sgn(sigma) * x^{fwd_T(sigma)}
is PALINDROMIC and SF(T,1) = 0.

QUESTIONS:
1. SF(T,-1) = ? (always zero for odd n? Always a specific value?)
2. What determines SF? Is it a finer invariant than F?
3. Does (x-1) | SF(T,x)? What is SF(T,x)/(x-1)?
4. SF via path reversal: sgn(P^rev) = sgn(P) * (-1)^{C(n,2)} and fwd(P^rev)=n-1-fwd(P)
   So SF(T,x) = (-1)^{C(n,2)} * x^{n-1} * SF(T, 1/x)
   At n=5: C(5,2)=10, so SF(T,x) = x^4 SF(T,1/x) — palindromic. ✓
   At n=4: C(4,2)=6, so SF(T,x) = x^3 SF(T,1/x) — palindromic. ✓
   At n=3: C(3,2)=3, so SF(T,x) = -x^2 SF(T,1/x) — ANTI-palindromic!

5. Does SF factor nicely?

Author: opus-2026-03-07-S46
"""
from itertools import permutations, combinations
from math import comb, factorial, gcd
from functools import reduce

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

def compute_F_and_SF(adj, n):
    F = [0]*n
    SF = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        inv = sum(1 for i in range(n) for j in range(i+1, n) if P[i] > P[j])
        sgn = (-1)**inv
        F[fwd] += 1
        SF[fwd] += sgn
    return F, SF

def count_3cycles(adj, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t3 += 1
    return t3

def poly_eval(F, x):
    return sum(F[k] * x**k for k in range(len(F)))

# ============================================================
# SIGNED F POLYNOMIAL: STRUCTURE
# ============================================================
print("=" * 60)
print("SIGNED F(T,x) STRUCTURE")
print("=" * 60)

for n in [3, 4, 5, 6]:
    m = n*(n-1)//2
    seen = set()
    sf_data = []

    import random
    random.seed(42)
    num = min(1 << m, 100000)

    for trial in range(num):
        if n <= 5:
            bits = trial
        else:
            bits = random.getrandbits(m)

        adj = tournament_from_bits(n, bits)
        F, SF = compute_F_and_SF(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        t3 = count_3cycles(adj, n)

        # Check palindrome type
        cn2 = comb(n, 2)
        parity = (-1)**cn2
        is_palindromic = all(SF[k] == parity * SF[n-1-k] for k in range(n))

        sf_data.append({
            'F': F, 'SF': SF, 't3': t3,
            'sf_at_1': poly_eval(SF, 1),
            'sf_at_m1': poly_eval(SF, -1),
            'palindromic': is_palindromic,
        })

    print(f"\nn={n}: C(n,2)={comb(n,2)}, parity=(-1)^C(n,2)={(-1)**comb(n,2)}")
    print(f"  {len(sf_data)} distinct F-vectors")

    for d in sf_data[:8]:
        pal_type = "palindromic" if d['palindromic'] else "NOT palindromic"
        print(f"  F={d['F']}, SF={d['SF']}, SF(1)={d['sf_at_1']}, SF(-1)={d['sf_at_m1']}, {pal_type}")

    # Check: does SF determine F or vice versa?
    f_to_sf = {}
    sf_to_f = {}
    for d in sf_data:
        fk = tuple(d['F'])
        sfk = tuple(d['SF'])
        f_to_sf[fk] = sfk
        if sfk not in sf_to_f:
            sf_to_f[sfk] = set()
        sf_to_f[sfk].add(fk)

    ambiguous = sum(1 for v in sf_to_f.values() if len(v) > 1)
    print(f"  SF→F ambiguity: {ambiguous} SF values map to >1 F-vector")

# ============================================================
# SF DIVISIBILITY BY (x-1)
# ============================================================
print("\n" + "=" * 60)
print("SF(T,x) / (x-1): quotient")
print("=" * 60)

for n in [4, 5]:
    m = n*(n-1)//2
    seen = set()
    print(f"\nn={n}:")

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F, SF = compute_F_and_SF(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        # Divide SF by (x-1) using synthetic division
        sf_at_1 = poly_eval(SF, 1)
        if sf_at_1 != 0:
            print(f"  WARNING: SF(1)={sf_at_1} != 0 for F={F}")
            continue

        # Synthetic division by (x-1)
        coeffs = list(SF)
        quot = [0] * (n-1)
        quot[-1] = coeffs[-1]
        for i in range(n-3, -1, -1):
            quot[i] = coeffs[i+1] + quot[i+1]
        # Verify: quot * (x-1) should give SF
        check = [0]*n
        for i in range(n-1):
            check[i+1] += quot[i]
            check[i] -= quot[i]
        assert check == SF, f"Division check failed: {check} vs {SF}"

        t3 = count_3cycles(adj, n)
        print(f"  F={F}, SF={SF}, SF/(x-1)={quot}, t3={t3}")

        # Check if quotient is palindromic or anti-palindromic
        is_pal = all(quot[k] == quot[len(quot)-1-k] for k in range(len(quot)//2 + 1))
        is_apal = all(quot[k] == -quot[len(quot)-1-k] for k in range(len(quot)//2 + 1))
        if is_pal:
            print(f"    quotient is PALINDROMIC")
        elif is_apal:
            print(f"    quotient is ANTI-PALINDROMIC")

# ============================================================
# SF AND DETERMINANT CONNECTION
# ============================================================
print("\n" + "=" * 60)
print("SF(T,x) AND DETERMINANT")
print("=" * 60)

# SF(T,x) = sum sgn(σ) x^{fwd_T(σ)}
# = sum sgn(σ) prod_{i=0}^{n-2} (x if σ(i)→σ(i+1) else 1)
# This looks like a DETERMINANT of something...
# Actually the permanent of W(x) = J-I+(x-1)A is:
# per(W) = sum_σ prod_i W[i][σ(i)] — but this is a matching, not a path.
#
# For a path, we need: prod_{i=0}^{n-2} W[σ(i)][σ(i+1)]
# This is the "path permanent" not the matrix permanent.
# With sign: sum sgn(σ) prod_i W[σ(i)][σ(i+1)] — "path determinant"
#
# Is there a matrix whose determinant equals SF(T,x)?
# Not obviously. But maybe via the transfer matrix...

# Transfer matrix approach: SF(T,x) = sum over Ham paths of sgn(path) * x^fwd
# where sgn(path) = sgn of the permutation (i_0, i_1, ..., i_{n-1}).
#
# This is related to the IMMANANT of the weighted adjacency matrix.
# The permanent gives F(T,x), the determinant gives something else,
# and SF is like the "path immanant" for the sign character.

import numpy as np

n = 5
print(f"\nn={n}:")
seen = set()
for bits in range(1 << (n*(n-1)//2)):
    adj = tournament_from_bits(n, bits)
    F, SF = compute_F_and_SF(adj, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    # Try: is SF related to det(W(x))?
    # W(x) = J-I + (x-1)A = J-I + (x-1)A
    # det(W(x)) is a polynomial in x of degree n (since A has degree 1 entries)

    A = np.array(adj, dtype=float)
    J = np.ones((n,n))
    I = np.eye(n)

    # Compute det(W(x)) as polynomial: evaluate at n+1 points, interpolate
    x_pts = list(range(n+2))
    det_pts = [np.linalg.det((J-I) + (x-1)*A) for x in x_pts]
    det_coeffs = np.round(np.polyfit(x_pts, det_pts, n), 2).tolist()

    sf_poly = list(SF)  # degree n-1

    print(f"  F={F}, SF={SF}")
    print(f"  det(W(x)) coeffs (high to low): {[f'{c:.0f}' for c in det_coeffs]}")

    # Check: is det(W(x)) / (something) = SF?
    # det(W(x)) has degree n, SF has degree n-1.
    # Maybe det(W(x)) = (x-1) * SF(T,x) + something?
    # Or det(W(x)) = (something) * (x-1) since det(W(1)) = det(J-I) = (-1)^{n-1}(n-1)

    # At x=1: det = (-1)^{n-1}(n-1) = 4 for n=5
    # SF(1) = 0
    # So (x-1) divides both... but det(W(x))/(x-1) has degree n-1 like SF

    # Compute det(W(x))/(x-1) numerically
    # det at x=1 should be 4, and we need (det - 4) ... no.
    # Actually, (x-1) may or may not divide det(W(x)).
    # det(W(1)) = 4 ≠ 0, so (x-1) does NOT divide det(W(x)).

    if len(seen) >= 4:
        break

# ============================================================
# SF AND THE REGULAR REPRESENTATION
# ============================================================
print("\n" + "=" * 60)
print("SF(T,x) GCD AND STRUCTURE")
print("=" * 60)

for n in [4, 5]:
    m = n*(n-1)//2
    seen = set()
    all_SF = []

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F, SF = compute_F_and_SF(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        all_SF.append(SF)

    # GCD of all SF coefficients
    all_gcd = reduce(gcd, [abs(c) for sf in all_SF for c in sf if c != 0])
    print(f"\nn={n}: {len(all_SF)} distinct SF, coefficient GCD = {all_gcd}")

    # Coefficient-wise GCD
    for k in range(n):
        vals = [sf[k] for sf in all_SF]
        g = reduce(gcd, [abs(v) for v in vals if v != 0], 0)
        print(f"  SF[{k}] values: {sorted(set(vals))}, GCD={g}")

    # Check: SF(T,x) / GCD = ?
    print(f"  Normalized (SF/GCD):")
    for sf in all_SF:
        norm = [c // all_gcd for c in sf]
        print(f"    {norm}")
