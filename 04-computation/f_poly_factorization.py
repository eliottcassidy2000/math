#!/usr/bin/env python3
"""
f_poly_factorization.py — Factoring F(T,x) over Z[x] and studying divisibility.

QUESTIONS:
1. Does (1+x) always divide F(T,x) at even n? (YES by palindrome: F(-1)=0)
2. What is the max power of (1+x) dividing F(T,x)?
3. Does F(T,x)/(1+x)^k have any universal structure?
4. Over Q, what are the irreducible factors of F(T,x)?
5. What is gcd of all F(T,x) over all tournaments at fixed n?

Author: opus-2026-03-07-S46
"""
from itertools import permutations
from math import factorial, gcd
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

def compute_F(adj, n):
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def poly_eval(F, x):
    return sum(F[k] * x**k for k in range(len(F)))

def poly_div_by_xp1(F):
    """Divide polynomial F by (x+1). Returns (quotient, remainder)."""
    n = len(F)
    Q = [0]*(n-1)
    r = 0
    for i in range(n-1, -1, -1):
        r = F[i] - r
        if i < n-1:
            Q[i] = r
    # Remainder is -r (since we divide by x+1)
    remainder = F[0] + sum((-1)**k * F[k] for k in range(len(F)))
    # Actually let's just use synthetic division
    coeffs = list(F)  # F[0] + F[1]*x + ... + F[n-1]*x^{n-1}
    # Dividing by (x+1) = (x - (-1))
    n = len(coeffs)
    quot = [0]*(n-1)
    quot[n-2] = coeffs[n-1]  # leading coeff
    for i in range(n-3, -1, -1):
        quot[i] = coeffs[i+1] - (-1)*quot[i+1]  # wait, this isn't right
    # Let me just do it properly
    # f(x) = q(x)*(x+1) + r
    # Use Horner at x=-1 for remainder
    rem = 0
    for i in range(n-1, -1, -1):
        rem = coeffs[i] + (-1)*rem  # wrong direction
    # Actually remainder = F(-1)
    rem = poly_eval(coeffs, -1)
    if rem != 0:
        return None, rem
    # Compute quotient by deflation
    # F(x) = (x+1)*Q(x), so Q(x) = F(x)/(x+1)
    # Synthetic division by root -1
    q = [0]*(n-1)
    q[n-2] = coeffs[n-1]
    for i in range(n-3, -1, -1):
        q[i] = coeffs[i+1] + (-1)*q[i+1]
    # Verify: wait, synthetic division from high to low
    # coeffs = [c0, c1, ..., c_{n-1}]
    # Dividing by (x - r) where r = -1
    # Start from highest degree
    qq = [0]*(n-1)
    qq[-1] = coeffs[-1]
    for i in range(n-3, -1, -1):
        qq[i] = coeffs[i+1] + (-1)*qq[i+1]
    return qq, rem

def max_power_xp1(F):
    """Find max k such that (x+1)^k | F(x)."""
    k = 0
    cur = list(F)
    while True:
        val = poly_eval(cur, -1)
        if val != 0:
            return k, cur
        # Divide by (x+1)
        n = len(cur)
        if n <= 1:
            return k, cur
        qq = [0]*(n-1)
        qq[-1] = cur[-1]
        for i in range(n-3, -1, -1):
            qq[i] = cur[i+1] + (-1)*qq[i+1]
        # Verify
        check = poly_eval(qq, -1)  # This should be the next derivative-like thing
        k += 1
        cur = qq

# ============================================================
# (1+x) DIVISIBILITY
# ============================================================
print("=" * 60)
print("(1+x) DIVISIBILITY OF F(T,x)")
print("=" * 60)

for n in [3, 4, 5, 6]:
    m = n*(n-1)//2
    powers = {}
    seen = set()
    all_quotients = []

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        k, quot = max_power_xp1(F)
        if k not in powers:
            powers[k] = 0
        powers[k] += 1
        all_quotients.append((k, quot, F))

    min_pow = min(powers.keys())
    max_pow = max(powers.keys())
    print(f"\nn={n}: {len(seen)} distinct F-vectors")
    print(f"  (1+x) power distribution: {dict(sorted(powers.items()))}")
    print(f"  min power: {min_pow}, max power: {max_pow}")

    # Show some examples
    for k, quot, F in all_quotients[:3]:
        print(f"  F={F}, (1+x)^{k} divides, quotient={quot}")

# ============================================================
# GCD OF ALL F(T,x) AT FIXED n
# ============================================================
print("\n" + "=" * 60)
print("GCD OF ALL F(T,x) COEFFICIENTS")
print("=" * 60)

for n in [3, 4, 5, 6]:
    m = n*(n-1)//2
    seen = set()
    all_F = []

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        all_F.append(F)

    # GCD of corresponding coefficients
    coeff_gcds = []
    for k in range(n):
        vals = [F[k] for F in all_F]
        g = reduce(gcd, vals)
        coeff_gcds.append(g)
    print(f"\nn={n}: coeff GCDs = {coeff_gcds}")

    # Polynomial GCD (coefficients of F/gcd)
    overall_gcd = reduce(gcd, [c for F in all_F for c in F])
    print(f"  Overall coefficient GCD = {overall_gcd}")

    # After dividing by (1+x)^min_pow, what's the GCD?
    min_pow = min(max_power_xp1(F)[0] for F in all_F)
    if min_pow > 0:
        reduced = []
        for F in all_F:
            _, quot = max_power_xp1(F)
            # Get quotient at exactly min_pow
            cur = list(F)
            for _ in range(min_pow):
                nn = len(cur)
                qq = [0]*(nn-1)
                qq[-1] = cur[-1]
                for i in range(nn-3, -1, -1):
                    qq[i] = cur[i+1] + (-1)*qq[i+1]
                cur = qq
            reduced.append(cur)
        red_gcds = []
        for k in range(len(reduced[0])):
            vals = [r[k] for r in reduced]
            g = reduce(gcd, vals)
            red_gcds.append(g)
        print(f"  After dividing by (1+x)^{min_pow}: coeff GCDs = {red_gcds}")

# ============================================================
# F(-1) VALUES AT ODD n
# ============================================================
print("\n" + "=" * 60)
print("F(-1) AT ODD n")
print("=" * 60)

for n in [3, 5, 7]:
    m = n*(n-1)//2
    seen = set()
    fm1_values = {}

    num = min(1 << m, 200000)
    import random
    random.seed(42)

    for trial in range(num):
        if n <= 5:
            bits = trial
        else:
            bits = random.getrandbits(m)

        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        fm1 = poly_eval(F, -1)
        if fm1 not in fm1_values:
            fm1_values[fm1] = 0
        fm1_values[fm1] += 1

    print(f"\nn={n}: {len(seen)} distinct F-vectors")
    print(f"  F(-1) values: {dict(sorted(fm1_values.items()))}")
    # Is F(-1) always divisible by something?
    all_fm1 = list(fm1_values.keys())
    if all_fm1:
        g = reduce(gcd, [abs(v) for v in all_fm1 if v != 0], 0)
        print(f"  GCD of |F(-1)| values: {g}")

# ============================================================
# F(omega) = F(e^{2pi*i/3}) — mod 3 structure
# ============================================================
print("\n" + "=" * 60)
print("F(x) mod 3 — residue classes")
print("=" * 60)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    seen = set()
    mod3_classes = {}

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        F_mod3 = tuple(f % 3 for f in F)
        if F_mod3 not in mod3_classes:
            mod3_classes[F_mod3] = 0
        mod3_classes[F_mod3] += 1

    print(f"\nn={n}: {len(mod3_classes)} distinct F mod 3 classes")
    for cls, cnt in sorted(mod3_classes.items(), key=lambda x: -x[1])[:10]:
        print(f"  F mod 3 = {list(cls)}: {cnt} F-vectors")

# ============================================================
# RESULTANT: does F(T,x) share roots with F(T',x)?
# ============================================================
print("\n" + "=" * 60)
print("SHARED ROOTS BETWEEN F(T,x) AND F(T',x) UNDER FLIP")
print("=" * 60)

# If F(T) - F(T') = (x-1)*D(x), then any common root of F(T) and F(T')
# must be a root of (x-1)*D(x). So common roots are either x=1 (impossible,
# F(1)=n!>0) or roots of D(x).

# What if D(x) = 0 (identically)? Then F(T) = F(T'): flipping doesn't change F.
# When does this happen?

n = 5
m = n*(n-1)//2
zero_D = 0
total_arcs = 0

for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    F_T = compute_F(adj, n)

    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue
            total_arcs += 1

            flip_adj = [row[:] for row in adj]
            flip_adj[u][v] = 0
            flip_adj[v][u] = 1
            F_flip = compute_F(flip_adj, n)

            if F_T == F_flip:
                zero_D += 1

print(f"\nn={n}: D=0 (F unchanged by flip) for {zero_D}/{total_arcs} arcs = {100*zero_D/total_arcs:.1f}%")
print("These are 'F-invariant arcs' — flipping them preserves the full F polynomial.")

# Which tournaments have the most F-invariant arcs?
best = {}
for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    F_T = compute_F(adj, n)
    count = 0
    for u in range(n):
        for v in range(n):
            if u == v or not adj[u][v]:
                continue
            flip_adj = [row[:] for row in adj]
            flip_adj[u][v] = 0
            flip_adj[v][u] = 1
            if compute_F(flip_adj, n) == F_T:
                count += 1
    key = tuple(F_T)
    if key not in best or count > best[key][0]:
        best[key] = (count, bits)

print("\nTournaments with most F-invariant arcs:")
for F_key, (count, bits) in sorted(best.items(), key=lambda x: -x[1][0])[:5]:
    print(f"  F={list(F_key)}, invariant arcs={count}/{n*(n-1)//2*2}, bits={bits}")
