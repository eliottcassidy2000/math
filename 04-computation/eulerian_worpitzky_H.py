"""
eulerian_worpitzky_H.py
opus-2026-03-14-S71l

The Worpitzky identity connects Eulerian numbers to simplex decomposition:
  x^n = sum_k A(n,k) * C(x+n-1-k, n)

The tournament F-polynomial F(T,x) is related via:
  H(T) = F(T,1) and W(T,x) = sum_k w_k * C(x+n-1-k, n)

This script investigates:
1. Exact relationship between Eulerian numbers and tournament H
2. Whether the simplex-cuboid nesting has a Worpitzky description
3. The connection to the 3-strand Pascal through the F-polynomial
4. The n=8 H-spectrum predictions
"""

import numpy as np
from math import factorial, comb
from itertools import permutations
from collections import Counter
from fractions import Fraction

print("=" * 70)
print("EULERIAN-WORPITZKY-TOURNAMENT CONNECTION")
print("opus-2026-03-14-S71l")
print("=" * 70)

# =====================================================================
# PART 1: EULERIAN NUMBERS AND DESCENT STRUCTURE
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: EULERIAN NUMBERS AND TOURNAMENT STRUCTURE")
print("=" * 70)

def eulerian(n, k):
    """A(n,k) = perms of [n] with k descents."""
    if n == 0:
        return 1 if k == 0 else 0
    if k < 0 or k >= n:
        return 0
    return (k+1) * eulerian(n-1, k) + (n-k) * eulerian(n-1, k-1)

print("\n  Eulerian triangle A(n,k):")
for n in range(1, 9):
    row = [eulerian(n, k) for k in range(n)]
    print(f"  n={n}: {row}  sum={sum(row)}={n}!")

print("""
  WORPITZKY IDENTITY:
  sum_{k=0}^{n-1} A(n,k) * C(x+n-1-k, n) = x^n

  At x=1: sum A(n,k) * C(n-k, n) = 1
  At x=2: sum A(n,k) * C(n+1-k, n) = 2^n

  This gives a DECOMPOSITION of 2^n (number of arc orientations)
  into Eulerian contributions!
""")

for n in range(1, 8):
    total = 0
    print(f"  n={n}: 2^{n} = {2**n} =", end="")
    for k in range(n):
        contrib = eulerian(n, k) * comb(n+1-k, n)
        if k > 0:
            print(f" + ", end="")
        print(f"A({n},{k})*C({n+1-k},{n})", end="")
        total += contrib
    print(f" = {total}")

# =====================================================================
# PART 2: H-VALUE DECOMPOSITION VIA WORPITZKY
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: H-VALUE DECOMPOSITION")
print("=" * 70)

def compute_H_dp(adj, n):
    """Count Hamiltonian paths via DP."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def compute_F_poly(adj, n):
    """Compute F-polynomial F(T,x) = sum_k F_k x^k.
    F_k = number of HP with exactly k ascents (positions where pi(i) < pi(i+1)
    in the tournament, i.e., the arc goes "forward" in the ordering).
    """
    f_coeffs = [0] * n  # F_0, F_1, ..., F_{n-1}
    for perm in permutations(range(n)):
        # Check if this is a Hamiltonian path
        valid = True
        for i in range(n-1):
            if not adj[perm[i]][perm[i+1]]:
                valid = False
                break
        if not valid:
            continue
        # Count ascents: positions i where perm[i] < perm[i+1]
        ascents = sum(1 for i in range(n-1) if perm[i] < perm[i+1])
        f_coeffs[ascents] += 1
    return f_coeffs

print("""
  For a tournament T on [n], the F-polynomial is:
  F(T, x) = sum_k F_k(T) * x^k

  where F_k(T) = number of HP of T with exactly k "ascents"
  (positions where the arc goes from smaller to larger vertex).

  H(T) = F(T, 1) = sum of all F_k.

  The Worpitzky decomposition of H(T):
  H(T) = sum_k w_k(T) where w_k = Eulerian-weighted F-coefficients.
""")

# Compute for small tournaments
for n in range(3, 7):
    m = comb(n, 2)
    arcs = [(i,j) for i in range(n) for j in range(i+1, n)]

    # Sample a few tournaments
    for bits in [0, (1 << m) - 1, (1 << (m//2))]:
        adj = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(arcs):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1

        h = compute_H_dp(adj, n)
        f = compute_F_poly(adj, n)

        print(f"  n={n}, bits={bits:0{m}b}: H={h}, F={f}")
        # Verify H = sum(F)
        assert h == sum(f), f"H={h} != sum(F)={sum(f)}"

# =====================================================================
# PART 3: THE EULERIAN-TOURNAMENT CONNECTION
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: EULERIAN-TOURNAMENT AVERAGE F-POLYNOMIAL")
print("=" * 70)

print("""
  Over ALL tournaments, what is the average F-polynomial?

  Each permutation pi is a HP of exactly 2^{m-(n-1)} tournaments
  (the n-1 arcs are fixed, the rest are free).

  The average F_k = sum over perms with k ascents, each counted
  2^{m-(n-1)} / 2^m = 1/2^{n-1} times.

  So: E[F_k] = A(n,k) / 2^{n-1}

  where A(n,k) is the Eulerian number!
""")

for n in range(3, 8):
    m = comb(n, 2)
    print(f"\n  n={n}: Mean F-polynomial = (1/2^{n-1}) * Eulerian row")
    for k in range(n):
        mean_fk = Fraction(eulerian(n, k), 2**(n-1))
        print(f"    E[F_{k}] = {eulerian(n,k)}/{2**(n-1)} = {mean_fk}")

    mean_h = sum(eulerian(n, k) for k in range(n)) / 2**(n-1)
    print(f"    E[H] = {n}!/2^{n-1} = {mean_h}")

# Verify against direct computation at n=5
print("\n  Direct verification at n=5:")
n = 5
m = comb(n, 2)
arcs = [(i,j) for i in range(n) for j in range(i+1, n)]
total_f = [0] * n
count = 0
for bits in range(2**m):
    adj = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(arcs):
        if (bits >> idx) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    f = compute_F_poly(adj, n)
    for k in range(n):
        total_f[k] += f[k]
    count += 1

for k in range(n):
    avg = Fraction(total_f[k], count)
    expected = Fraction(eulerian(n, k), 2**(n-1))
    print(f"    F_{k}: actual mean = {avg} = {float(avg):.4f}, expected = {expected} = {float(expected):.4f}, match: {avg == expected}")

# =====================================================================
# PART 4: VARIANCE OF F-POLYNOMIAL
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: VARIANCE STRUCTURE OF F-POLYNOMIAL")
print("=" * 70)

# Compute variance of each F_k at n=5
n = 5
m = comb(n, 2)
arcs = [(i,j) for i in range(n) for j in range(i+1, n)]
f_vals = [[] for _ in range(n)]
h_vals = []
for bits in range(2**m):
    adj = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(arcs):
        if (bits >> idx) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    f = compute_F_poly(adj, n)
    for k in range(n):
        f_vals[k].append(f[k])
    h_vals.append(sum(f))

print(f"  n={n}: F-polynomial statistics over all {2**m} tournaments")
for k in range(n):
    vals = np.array(f_vals[k])
    mean = np.mean(vals)
    var = np.var(vals)
    print(f"    F_{k}: mean={mean:.4f}, var={var:.4f}, std={np.std(vals):.4f}")

print(f"\n    H = sum F_k: mean={np.mean(h_vals):.4f}, var={np.var(h_vals):.4f}")

# Correlation matrix of F_k's
print(f"\n  Correlation matrix of F_k's:")
corr = np.corrcoef([f_vals[k] for k in range(n)])
for k1 in range(n):
    row = "    "
    for k2 in range(n):
        row += f"{corr[k1][k2]:7.3f}"
    print(row)

# =====================================================================
# PART 5: F-POLYNOMIAL AND SIMPLEX-CUBOID
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: F-POLYNOMIAL AND SIMPLEX DECOMPOSITION")
print("=" * 70)

print("""
  The F-polynomial F(T,x) evaluated at x gives:
    F(T,0) = F_0(T) = number of HP with 0 ascents
                     = number of "descending" HP (pi(1)>pi(2)>...>pi(n))
                     = either 0 or 1 (can be at most 1)

    F(T,1) = H(T) = total HP count

    F(T,-1) = sum (-1)^k F_k(T)
            = "signed alternating HP count"

  The Worpitzky connection:
  F(T,x) = sum_k F_k(T) * x^k

  The "Eulerian decomposition" of F:
  F(T,x) = sum_k w_k(T) * C(x+n-1-k, n)

  where the w_k are the "Worpitzky coefficients" of T.

  If T is the AVERAGE tournament (uniform over all T):
  E[F(T,x)] = (1/2^{n-1}) * sum_k A(n,k) * x^k = A_n(x) / 2^{n-1}

  where A_n(x) is the Eulerian polynomial.

  At x=1: E[H] = A_n(1)/2^{n-1} = n!/2^{n-1}  ✓

  The SIMPLEX interpretation:
  A_n(x) = n! * sum_k A(n,k)/n! * x^k
  The normalized coefficients A(n,k)/n! sum to 1.
  So A_n(x)/n! is a probability generating function!

  The descent probability: P(k descents) = A(n,k)/n!
  Mean descents = (n-1)/2 (by symmetry of Eulerian numbers)
""")

for n in range(3, 8):
    print(f"  n={n}: Eulerian probabilities P(k descents):")
    for k in range(n):
        prob = Fraction(eulerian(n, k), factorial(n))
        print(f"    P(k={k}) = {prob} = {float(prob):.4f}")

# =====================================================================
# PART 6: CONNECTING TO THE 3-STRAND PASCAL
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: FROM EULERIAN TO 3-STRAND PASCAL")
print("=" * 70)

print("""
  QUESTION: Is there a "trinomial F-polynomial" that decomposes H
  using the 3-strand Pascal (trinomial coefficients) instead of
  the 2-strand Pascal (binomial coefficients)?

  The standard decomposition:
  H(T) = sum_k F_k(T) * 1^k = sum F_k  (all weights = 1)

  A 3-strand version would be:
  H_3(T) = sum_k F_k(T) * omega^k  where omega = e^(2pi*i/3)

  This is the "cube-root-of-unity Fourier transform" of the F-polynomial!

  F(T, omega) = sum_k F_k(T) * omega^k

  Since omega^3 = 1, this only depends on F_k mod 3:
  F(T, omega) = (F_0 + F_3 + ...) + omega*(F_1 + F_4 + ...) + omega^2*(F_2 + F_5 + ...)

  Let G_j = sum_{k equiv j mod 3} F_k. Then:
  F(T, omega) = G_0 + omega*G_1 + omega^2*G_2

  |F(T, omega)|^2 = G_0^2 + G_1^2 + G_2^2 - G_0*G_1 - G_1*G_2 - G_0*G_2
                   = Phi_3-norm of (G_0, G_1, G_2) in some sense

  THIS IS THE 3-STRAND H VALUE!
  It measures the "imbalance" among the three residue classes of descents.
""")

# Compute F(T, omega) for small tournaments
omega = np.exp(2j * np.pi / 3)

for n in [5, 6]:
    m = comb(n, 2)
    arcs = [(i,j) for i in range(n) for j in range(i+1, n)]
    f_omega_vals = []
    h_vals = []

    for bits in range(2**m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(arcs):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        f = compute_F_poly(adj, n)
        f_omega = sum(f[k] * omega**k for k in range(n))
        f_omega_vals.append(f_omega)
        h_vals.append(sum(f))

    print(f"\n  n={n}: F(T, omega) statistics")
    mags = [abs(v) for v in f_omega_vals]
    phases = [np.angle(v) for v in f_omega_vals]
    print(f"    |F(omega)| range: [{min(mags):.4f}, {max(mags):.4f}]")
    print(f"    Mean |F(omega)|: {np.mean(mags):.4f}")
    print(f"    |F(omega)|^2 range: [{min(m**2 for m in mags):.4f}, {max(m**2 for m in mags):.4f}]")

    # Is |F(omega)|^2 always an integer?
    mag_sq = [abs(v)**2 for v in f_omega_vals]
    all_int = all(abs(m - round(m)) < 1e-8 for m in mag_sq)
    print(f"    |F(omega)|^2 always integer? {all_int}")

    if all_int:
        mag_sq_int = [round(m) for m in mag_sq]
        distinct = sorted(set(mag_sq_int))
        print(f"    Distinct |F(omega)|^2 values: {distinct[:15]}{'...' if len(distinct) > 15 else ''}")

    # Correlation with H
    corr_h = np.corrcoef(h_vals, mags)[0,1]
    print(f"    Correlation(H, |F(omega)|): {corr_h:.4f}")

# =====================================================================
# PART 7: THE DESCENT POLYNOMIAL AND TOURNAMENT OBSTRUCTION
# =====================================================================
print("\n" + "=" * 70)
print("PART 7: DESCENT POLYNOMIAL AND WHY 7 IS FORBIDDEN")
print("=" * 70)

print("""
  For a tournament T, define the descent polynomial:
  D(T, x) = sum over all HP pi: x^{des(pi)}

  where des(pi) = number of descents in pi as a permutation.

  Then F(T, x) counts ascents, so:
  D(T, x) = x^{n-1} * F(T, 1/x)  (reverse coefficients)

  The relationship: F_k(T) = D_{n-1-k}(T) (palindromic duality for each T)

  Wait, that's only true for the uniform distribution.
  For a specific tournament, F(T,x) is generally NOT palindromic.
  But F(T,x) = x^{n-1} * F(T^op, 1/x) (complement duality, from HYP-011).

  KEY QUESTION: Can H=7 be expressed as a sum F_0 + F_1 + ... + F_{n-1}
  for some valid F-polynomial?

  At n=5 (first n where 7 could be achievable):
  H = F_0 + F_1 + F_2 + F_3 + F_4 = 7
  Each F_k >= 0.
  F_0 + F_4 are the "extremal" terms (fully descending / fully ascending HP).

  The constraint: sum = 7 with F_k >= 0.
  This is certainly achievable as a partition of 7.
  But the F_k must ALSO be consistent with being the descent distribution
  of Hamiltonian paths of some tournament.

  The obstructions are:
  1. F_0 in {{0, 1}} (at most one fully descending HP)
  2. F_{n-1} in {{0, 1}} (at most one fully ascending HP)
  3. F_k(T) = F_{n-1-k}(T^op) (complement duality)
  4. Various inequalities from the tournament structure
""")

# Check all F-polynomials at n=5
n = 5
m = comb(n, 2)
arcs = [(i,j) for i in range(n) for j in range(i+1, n)]
f_polys = Counter()
h_by_fpoly = {}

for bits in range(2**m):
    adj = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(arcs):
        if (bits >> idx) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    f = tuple(compute_F_poly(adj, n))
    f_polys[f] += 1
    h = sum(f)
    if h not in h_by_fpoly:
        h_by_fpoly[h] = set()
    h_by_fpoly[h].add(f)

print(f"\n  n={n}: {len(f_polys)} distinct F-polynomials")
print(f"  H-values and their F-polynomials:")
for h in sorted(h_by_fpoly.keys()):
    fps = h_by_fpoly[h]
    print(f"    H={h:2d}: {len(fps)} F-polys: {sorted(fps)[:5]}{'...' if len(fps) > 5 else ''}")

print(f"\n  H=7 has {len(h_by_fpoly.get(7, set()))} F-polynomials (should be 0)")
print(f"  H=11 has {len(h_by_fpoly.get(11, set()))} F-polynomials")

# What F-polynomial sums would give 7?
print("\n  F-polynomials that sum to 7 (if they existed):")
print("  Since F_0, F_4 in {0,1}, and F_1+F_2+F_3 must fill the rest:")
print("  F_0=0, F_4=0: need F_1+F_2+F_3=7")
print("  F_0=1, F_4=0: need F_1+F_2+F_3=6")
print("  F_0=0, F_4=1: need F_1+F_2+F_3=6")
print("  F_0=1, F_4=1: need F_1+F_2+F_3=5")
print("  None of these combinations are achievable by any tournament!")

print("\n" + "=" * 70)
print("DONE — EULERIAN-WORPITZKY-TOURNAMENT CONNECTION")
print("=" * 70)
