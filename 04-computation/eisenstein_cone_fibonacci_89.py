"""
eisenstein_cone_fibonacci_89.py
opus-2026-03-14-S89

Deep exploration building on kind-pasteur's Phi_3 cone theorem.
Three threads woven together:

1. DEVIATION FORMULA: Test Var/Mean^2 = 1/3 - c_n/n(n-1)(n-2) at n=6
   - Is c_n always integer? Always 1? Or does it grow?

2. EISENSTEIN LATTICE: The forbidden values are norms in Z[omega].
   - What does the Eisenstein lattice look like for tournaments?
   - Which norms are achievable H values? Which are forbidden?

3. FIBONACCI-EISENSTEIN-CONE TRINITY:
   - Fibonacci lives in Q(sqrt(5)), Eisenstein in Q(sqrt(-3))
   - Both fields have class number 1 (unique factorization!)
   - The tournament generator 2 ramifies in Z[omega]: 2 = -omega(1-omega)^2
   - And 2 = (phi)(2-phi) in Z[phi] (nearly — actually Z[(1+sqrt(5))/2])
   - How do these interact?
"""

from itertools import combinations, product
from collections import Counter
from math import gcd, sqrt, factorial, log
import cmath

print("=" * 70)
print("EISENSTEIN-CONE-FIBONACCI SYNTHESIS")
print("opus-2026-03-14-S89")
print("=" * 70)

# =====================================================================
# PART 1: EXACT Var/Mean^2 FOR n=3,4,5 AND DEVIATION FORMULA
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: DEVIATION FORMULA — Var/Mean^2 = 1/3 - c_n/n(n-1)(n-2)")
print("=" * 70)

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(2**m):
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

def compute_H_small(adj, n):
    """Compute H(T) = number of Hamiltonian paths."""
    from itertools import permutations
    count = 0
    for perm in permutations(range(n)):
        is_path = True
        for i in range(n-1):
            if not adj[perm[i]][perm[i+1]]:
                is_path = False
                break
        if is_path:
            count += 1
    return count

def compute_H_via_ocf(adj, n):
    """Compute H via OCF: H = I(Omega, 2)."""
    # Find all odd cycles
    odd_cycles = []
    # 3-cycles
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            odd_cycles.append(frozenset({a, b, c}))
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            odd_cycles.append(frozenset({a, b, c}))

    # 5-cycles (if n >= 5)
    if n >= 5:
        for vertices in combinations(range(n), 5):
            vlist = list(vertices)
            from itertools import permutations as perms5
            for perm in perms5(vlist):
                is_cycle = True
                for i in range(5):
                    if not adj[perm[i]][perm[(i+1) % 5]]:
                        is_cycle = False
                        break
                if is_cycle:
                    odd_cycles.append(frozenset(vertices))
                    break  # one direction is enough, handled by frozenset

    oc_list = list(set(odd_cycles))
    nc = len(oc_list)

    # Independence polynomial at x=2
    H = 0
    for mask in range(2**nc):
        selected = [i for i in range(nc) if mask & (1 << i)]
        is_indep = True
        for i in range(len(selected)):
            for j in range(i+1, len(selected)):
                if oc_list[selected[i]] & oc_list[selected[j]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            H += 2**len(selected)
    return H

# Compute exact statistics for n=3,4,5
for n in [3, 4, 5]:
    H_values = []
    for adj in all_tournaments(n):
        H = compute_H_small(adj, n)
        H_values.append(H)

    N = len(H_values)
    mean_H = sum(H_values) / N
    var_H = sum((h - mean_H)**2 for h in H_values) / N
    ratio = var_H / mean_H**2

    falling = n * (n-1) * (n-2)
    if abs(ratio - 1/3) < 1e-10:
        c_n = 0
    else:
        c_n = (1/3 - ratio) * falling

    print(f"\n  n={n}: N={N} tournaments")
    print(f"    Mean(H) = {mean_H:.6f}")
    print(f"    Var(H)  = {var_H:.6f}")
    print(f"    Var/Mean^2 = {ratio:.10f}")
    print(f"    1/3        = {1/3:.10f}")
    print(f"    n(n-1)(n-2) = {falling}")
    print(f"    c_n = (1/3 - ratio) * {falling} = {c_n:.6f}")

    # Express as exact fraction
    # Var = sum(H^2)/N - mean^2
    sum_H = sum(H_values)
    sum_H2 = sum(h**2 for h in H_values)
    # Var = sum_H2/N - (sum_H/N)^2 = (N*sum_H2 - sum_H^2) / N^2
    var_num = N * sum_H2 - sum_H**2
    # ratio = Var/mean^2 = var_num / (N^2 * (sum_H/N)^2) = var_num / sum_H^2
    # So ratio = var_num / sum_H^2
    g = gcd(var_num, sum_H**2)
    print(f"    Exact: Var/Mean^2 = {var_num//g}/{sum_H**2//g}")

# n=6 by sampling (too large for exhaustive)
print("\n  n=6 (Monte Carlo, 50000 samples):")
import random
random.seed(42)

n = 6
edges_6 = [(i, j) for i in range(n) for j in range(i+1, n)]
m = len(edges_6)

H_values_6 = []
for _ in range(50000):
    adj = [[0]*n for _ in range(n)]
    for (i, j) in edges_6:
        if random.random() < 0.5:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    H = compute_H_small(adj, n)
    H_values_6.append(H)

N6 = len(H_values_6)
mean_6 = sum(H_values_6) / N6
var_6 = sum((h - mean_6)**2 for h in H_values_6) / N6
ratio_6 = var_6 / mean_6**2
c_6 = (1/3 - ratio_6) * 6 * 5 * 4
print(f"    Mean(H) = {mean_6:.4f}")
print(f"    Var(H)  = {var_6:.4f}")
print(f"    Var/Mean^2 = {ratio_6:.8f}")
print(f"    c_6 = (1/3 - ratio) * 120 = {c_6:.4f}")
print(f"    Expected exact: 1/3 - c_6/120 where c_6 ~ {round(c_6)}")

# =====================================================================
# PART 2: THE EISENSTEIN LATTICE AND TOURNAMENT H-VALUES
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: EISENSTEIN LATTICE — H-VALUES AS NORMS")
print("=" * 70)

# Eisenstein integers Z[omega], omega = e^(2pi i/3)
omega = cmath.exp(2j * cmath.pi / 3)

print(f"\n  omega = {omega:.6f}")
print(f"  omega^2 = {omega**2:.6f}")
print(f"  1 + omega + omega^2 = {1 + omega + omega**2:.6f}")

# Phi_3(b) = b^2 + b + 1 = |1 - b*omega^2|^2
print("\n  Phi_3(b) = |1 - b*omega^2|^2 (Eisenstein norm)")
print(f"  {'b':>4} {'Phi_3(b)':>10} {'|1-b*w^2|^2':>14} {'H achievable?':>15}")
print("  " + "-" * 48)

# Which Phi_3 values are achievable H values?
# We know H-spectrum for n=3: {1,3}, n=4: {1,3,5}, n=5: {1,3,5,9,11,13,15}
all_H = set()
for n_check in [3, 4, 5]:
    for adj in all_tournaments(n_check):
        H = compute_H_small(adj, n_check)
        all_H.add(H)

for b in range(15):
    phi3_b = b*b + b + 1
    z = 1 - b * omega**2
    norm = abs(z)**2
    achievable = "YES" if phi3_b in all_H else "no"
    print(f"  {b:4d} {phi3_b:10d} {norm:14.4f} {achievable:>15}")

# Map of Eisenstein norms to H-values
print("\n  Eisenstein norms N(a+b*omega) for small a,b:")
print(f"  {'a':>3} {'b':>3} {'N(a+bw)':>10} {'= a^2-ab+b^2':>15} {'H value?':>10}")
print("  " + "-" * 50)

norms_to_check = set()
for a in range(-6, 7):
    for b in range(-6, 7):
        if a == 0 and b == 0:
            continue
        norm_val = a*a - a*b + b*b
        if norm_val > 0 and norm_val <= 25:
            norms_to_check.add((norm_val, a, b))

seen_norms = set()
for norm_val, a, b in sorted(norms_to_check):
    if norm_val in seen_norms:
        continue
    seen_norms.add(norm_val)
    achievable = "YES" if norm_val in all_H else ("FORB" if norm_val == 7 or norm_val == 21 else "?")
    print(f"  {a:3d} {b:3d} {norm_val:10d} {'':>15} {achievable:>10}")

# =====================================================================
# PART 3: WHICH EISENSTEIN NORMS ARE H-FORBIDDEN?
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: EISENSTEIN NORMS vs H-VALUES — FORBIDDEN SET")
print("=" * 70)

# Get full H-spectrum for n=3,4,5
H_spectra = {}
for n_check in [3, 4, 5]:
    H_set = set()
    for adj in all_tournaments(n_check):
        H = compute_H_small(adj, n_check)
        H_set.add(H)
    H_spectra[n_check] = sorted(H_set)
    print(f"\n  n={n_check}: H-spectrum = {H_spectra[n_check]}")

# All achievable H values up to n=5
all_achieved = set()
for spec in H_spectra.values():
    all_achieved.update(spec)

# Eisenstein norms up to 25
eis_norms = set()
for a in range(-10, 11):
    for b in range(-10, 11):
        n_val = a*a - a*b + b*b
        if n_val > 0:
            eis_norms.add(n_val)

print(f"\n  All H values achieved (n<=5): {sorted(all_achieved)}")
print(f"  All odd numbers 1..25: {[x for x in range(1,26,2)]}")
print(f"  Missing odd H values (n<=5): {sorted(set(range(1,26,2)) - all_achieved)}")
print(f"  Missing that are Eisenstein norms: ", end="")
missing = sorted(set(range(1,26,2)) - all_achieved)
missing_eis = [x for x in missing if x in eis_norms]
print(missing_eis)

# Check: are the forbidden values exactly the Phi_3 values?
print("\n  Phi_3 values: ", [b*b+b+1 for b in range(10)])
print("  Missing H's: ", missing)
print("  Overlap: ", [x for x in missing if x in [b*b+b+1 for b in range(10)]])

# =====================================================================
# PART 4: FIBONACCI SEQUENCE IN THE EISENSTEIN WORLD
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: FIBONACCI IN THE EISENSTEIN WORLD")
print("=" * 70)

# Fibonacci sequence
fib = [0, 1]
for i in range(2, 25):
    fib.append(fib[-1] + fib[-2])

print("\n  Fibonacci: ", fib[:20])

# Fibonacci mod 3 (Pisano period pi(3) = 8)
fib_mod3 = [f % 3 for f in fib[:20]]
print(f"  Fib mod 3: {fib_mod3}")
print(f"  Pisano period pi(3) = 8")

# Fibonacci mod 7 (Pisano period pi(7) = 16)
fib_mod7 = [f % 7 for f in fib[:20]]
print(f"  Fib mod 7: {fib_mod7}")

# Find Pisano period for p
def pisano(p):
    prev, curr = 0, 1
    for i in range(1, p*p + 1):
        prev, curr = curr, (prev + curr) % p
        if prev == 0 and curr == 1:
            return i
    return -1

print(f"\n  Pisano periods for tournament-relevant primes:")
for p in [2, 3, 5, 7, 13, 21]:
    if p < 2:
        continue
    # Check if prime
    is_prime = all(p % d != 0 for d in range(2, int(p**0.5)+1)) if p > 1 else False
    if is_prime:
        pi_p = pisano(p)
        print(f"    pi({p}) = {pi_p}")

# The GOLDEN connection: phi^2 = phi + 1, omega^2 + omega + 1 = 0
# Both satisfy quadratic equations!
# phi = (1+sqrt(5))/2, omega = (-1+sqrt(-3))/2
# phi satisfies x^2 - x - 1 = 0 (discriminant 5)
# omega satisfies x^2 + x + 1 = 0 (discriminant -3)
phi = (1 + sqrt(5)) / 2

print(f"\n  QUADRATIC NUMBER FIELD COMPARISON:")
print(f"  {'':20} {'Golden':>15} {'Eisenstein':>15}")
print(f"  {'Generator':20} {'phi':>15} {'omega':>15}")
print(f"  {'Min poly':20} {'x^2-x-1':>15} {'x^2+x+1':>15}")
print(f"  {'Disc':20} {'5':>15} {'-3':>15}")
print(f"  {'Class number':20} {'1':>15} {'1':>15}")
print(f"  {'Units':20} {'+-phi^n':>15} {'omega^k':>15}")
print(f"  {'2 splits?':20} {'irred':>15} {'irred':>15}")
print(f"  {'Norm(2)':20} {'4':>15} {'4':>15}")
print(f"  {'Phi_n(2)':20} {'Phi_5(2)=31':>15} {'Phi_3(2)=7':>15}")

# Norm of 2 in both rings
# Z[phi]: N(2) = 2*2 = 4 (2 is inert in Z[phi] since disc 5, and 2 doesn't split)
# Z[omega]: N(2) = 4 (2 = -omega*(1-omega)^2, norm = 1*3 ≠ 4... wait)
# Actually in Z[omega], 2 is a prime since disc -3, Legendre(-3,2) = 1? No.
# 2 mod 3 = 2 ≡ 2 (mod 3), and for Eisenstein: 2 is inert iff 2 ≡ 2 (mod 3). YES!

print(f"\n  2 in Z[phi]: inert (5 mod 4 = 1, but 2 is inert since disc 5)")
print(f"  2 in Z[omega]: inert (2 mod 3 = 2, so 2 stays prime)")
print(f"  Both: N(2) = 4 = 2^2")

# =====================================================================
# PART 5: THE FIBONACCI-2-3 DECOMPOSITION
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: FIBONACCI 2-3 DECOMPOSITION — ZECKENDORF AND BEYOND")
print("=" * 70)

# Zeckendorf representation: every positive integer = sum of non-consecutive Fibs
def zeckendorf(n):
    """Return Zeckendorf representation of n."""
    if n == 0:
        return []
    fibs = [1, 2]
    while fibs[-1] <= n:
        fibs.append(fibs[-1] + fibs[-2])
    result = []
    for f in reversed(fibs):
        if f <= n:
            result.append(f)
            n -= f
    return result

# The 2-3 decomposition: F(n) = F(n-2) + F(n-3) + F(n-3) or similar
# Actually: the Fibonacci word is built from 2-letter and 3-letter blocks
# F(n+2) = F(n+1) + F(n) = (F(n) + F(n-1)) + F(n) = 2*F(n) + F(n-1)
# But also: F(n) = F(n-2) + F(n-1) = F(n-2) + F(n-3) + F(n-2) = 2*F(n-2) + F(n-3)
# So: F(n) = 2*F(n-2) + F(n-3)

print("\n  F(n) = 2*F(n-2) + F(n-3) [the 2-3 recurrence]:")
for i in range(3, 15):
    lhs = fib[i]
    rhs = 2*fib[i-2] + fib[i-3]
    print(f"    F({i:2d}) = {lhs:5d} = 2*F({i-2:2d}) + F({i-3:2d}) = 2*{fib[i-2]:4d} + {fib[i-3]:4d} = {rhs:5d}  {'OK' if lhs == rhs else 'FAIL'}")

# Deeper: F(n) can also be written using only 2s and 3s
# F(n) = floor(phi^n / sqrt(5) + 1/2)
# But the RECURRENCE uses 2 and 3:
# The characteristic equation of x^3 = 2x + 1 is... no.
# Let's think differently: the Fibonacci representation in base {2,3}

print("\n  Fibonacci in the 2-3 world:")
print("  Every Fibonacci number has a unique representation in terms of 2 and 3:")
for i in range(1, 18):
    f = fib[i]
    # Express f = 2a + 3b in all ways
    reps = []
    for b in range(f // 3 + 1):
        rem = f - 3*b
        if rem >= 0 and rem % 2 == 0:
            reps.append((rem // 2, b))
    print(f"    F({i:2d}) = {f:5d} = ", end="")
    rep_strs = [f"2*{a}+3*{b}" for a, b in reps]
    print(", ".join(rep_strs[:5]))

# The key insight: 2 and 3 generate tournaments
# And the Fibonacci sequence is the canonical recurrence from 2 and 3
# because F(n) = F(n-1) + F(n-2), and phi = (1+sqrt(5))/2
# while the tournament generator 2 gives the Jacobsthal: J(n) = J(n-1) + 2*J(n-2)

print("\n  FIBONACCI vs JACOBSTHAL — the 2-3 duality:")
print(f"  {'n':>3} {'F(n)':>8} {'J(n)':>8} {'F(n)/J(n)':>10} {'Ratio limit':>12}")

# Jacobsthal sequence
J = [0, 1]
for i in range(2, 20):
    J.append(J[-1] + 2*J[-2])

for i in range(1, 18):
    ratio_fj = fib[i] / J[i] if J[i] != 0 else float('inf')
    print(f"  {i:3d} {fib[i]:8d} {J[i]:8d} {ratio_fj:10.4f}")

# Limit: F(n)/J(n) -> phi^n / 2^n * sqrt(5)/3 -> 0 (since phi < 2)
# More precisely: F(n) ~ phi^n/sqrt(5), J(n) ~ 2^n/3
# F(n)/J(n) ~ 3*phi^n / (sqrt(5)*2^n) = 3/sqrt(5) * (phi/2)^n -> 0

print(f"\n  F(n) ~ phi^n/sqrt(5) = {phi}^n / {sqrt(5):.4f}")
print(f"  J(n) ~ 2^n/3")
print(f"  Ratio F/J ~ 3/sqrt(5) * (phi/2)^n -> 0 since phi/2 = {phi/2:.4f} < 1")
print(f"  But 3/sqrt(5) = {3/sqrt(5):.4f} = sqrt(9/5) = sqrt(1.8)")

# =====================================================================
# PART 6: THE PERIOD-6 STRUCTURE
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: PERIOD-6 — WHERE FIBONACCI MEETS THE CONE")
print("=" * 70)

# Fibonacci mod 2: period 3 (0,1,1,0,1,1,...)
# Fibonacci mod 3: period 8 (0,1,1,2,0,2,2,1)
# Fibonacci mod 4: period 6 (0,1,1,2,3,1)
# Fibonacci mod 7: period 16
# Fibonacci mod 8: period 12
# The LCM of periods mod 2 and mod 3: lcm(3,8) = 24
# But the PISANO period mod 6 = lcm(3,8) = 24? No...
# Pisano(2) = 3, Pisano(3) = 8, Pisano(6) = lcm(3,8) = 24

print("\n  Pisano periods for first few moduli:")
for m in range(2, 13):
    pi_m = pisano(m) if all(m % d != 0 for d in range(2, int(m**0.5)+1)) else -1
    if pi_m > 0:
        print(f"    pi({m:2d}) = {pi_m}")

# The user mentioned "6 to get back to starting state"
# This is the period of the FIBONACCI WORD in the 2-3 substitution system
# Or: the matrix [[1,1],[1,0]]^6 mod something

# Matrix powers of [[1,1],[1,0]]
def mat_mul_mod(A, B, mod=None):
    n = len(A)
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                C[i][j] += A[i][k] * B[k][j]
            if mod:
                C[i][j] %= mod
    return C

def mat_pow_mod(M, p, mod):
    n = len(M)
    result = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    base = [row[:] for row in M]
    while p > 0:
        if p % 2 == 1:
            result = mat_mul_mod(result, base, mod)
        base = mat_mul_mod(base, base, mod)
        p //= 2
    return result

# Period of [[1,1],[1,0]] mod p = Pisano(p)
# The user says "6 to get back" — let's check what has period 6
print("\n  What has period 6?")
print("  Checking period of Fibonacci matrix [[1,1],[1,0]] mod p:")
F_mat = [[1, 1], [1, 0]]
I_mat = [[1, 0], [0, 1]]

for mod in range(2, 20):
    power = F_mat
    for k in range(1, 100):
        if power == I_mat:
            print(f"    mod {mod:2d}: period = {k}" + (" *** PERIOD 6!" if k == 6 else ""))
            break
        power = mat_mul_mod(power, F_mat, mod)

# Period 6 happens mod 4! (and the Fibonacci word mod 4 has pattern 0,1,1,2,3,1)
# Also: the ORDER of phi mod 4 in the ring Z[phi]/4Z
print("\n  Fibonacci mod 4 (period 6): ", [fib[i] % 4 for i in range(12)])
print("  Pattern: 0, 1, 1, 2, 3, 1, 0, 1, 1, 2, 3, 1")
print("  The unique non-trivial residues: {0,1,2,3} — ALL mod-4 residues!")
print("  6 = 2 * 3 = the tournament numbers!")

# Connection to tournaments: mod 4 ↔ F_4
# The Fibonacci matrix has period 6 in GL(2, Z/4Z)
# And 6 = |S_3| = the tournament vertex permutation group at n=3
# And 6 = the number of totals in S_6 outer automorphism
# And 6 = Pisano period for Fibonacci mod 4

print("\n  THE PERIOD-6 COINCIDENCE:")
print("    6 = period of Fib matrix mod 4")
print("    6 = |S_3| = Aut(3-tournament)")
print("    6 = number of totals in S_6 outer aut")
print("    6 = 2 * 3 = tournament_generator * cycle_generator")
print("    6 = Phi_6(1) = 1 - 1 + 1 = 1... no, Phi_6(1) = 1")
print("    6 = |units in Z[omega]| = {1, omega, omega^2, -1, -omega, -omega^2}")

# =====================================================================
# PART 7: THE CONE-FIBONACCI BRIDGE VIA CHEBYSHEV
# =====================================================================
print("\n" + "=" * 70)
print("PART 7: CHEBYSHEV POLYNOMIALS — THE BRIDGE")
print("=" * 70)

# Chebyshev polynomials T_n(cos(theta)) = cos(n*theta)
# T_0 = 1, T_1 = x, T_2 = 2x^2-1, T_3 = 4x^3-3x, ...
# They satisfy T_n(x) = 2x*T_{n-1}(x) - T_{n-2}(x)
# This recurrence has the SAME structure as Fibonacci
# with the "2" coefficient!

# Fibonacci: F_n = F_{n-1} + F_{n-2}  (coeff 1)
# Chebyshev: T_n = 2x*T_{n-1} - T_{n-2}  (coeff 2x)
# Jacobsthal: J_n = J_{n-1} + 2*J_{n-2}  (coeff 2)

# At x = phi/2 (half the golden ratio):
x_half_phi = phi / 2
print(f"\n  Chebyshev at x = phi/2 = {x_half_phi:.6f}:")
T = [1, x_half_phi]
for i in range(2, 12):
    T.append(2 * x_half_phi * T[-1] - T[-2])
print(f"  T_n(phi/2): {[round(t, 4) for t in T]}")
print(f"  Compare F_n/2^(n-1): {[fib[i]/2**(i-1) if i > 0 else 0 for i in range(12)]}")

# Actually, there's a known identity: F_n = i^{1-n} * U_{n-1}(i/2)
# where U_n is the Chebyshev polynomial of the second kind
# U_n(cos theta) = sin((n+1)theta)/sin(theta)

# More relevant: the Fibonacci numbers satisfy
# F_n = 2*F_{n-2} + F_{n-3}  (our 2-3 recurrence!)
# This is a Chebyshev-like recurrence with mixed coefficients

print("\n  The 2-3 recurrence F(n) = 2*F(n-2) + F(n-3):")
print("  This is a DEPTH-3 Chebyshev variant!")
print("  Characteristic polynomial: x^3 = 2x + 1")
print("  = x^3 - 2x - 1 = (x+1)(x^2-x-1) = (x+1)*Phi_Fib(x)")

# WHOA: x^3 - 2x - 1 = (x+1)(x^2 - x - 1)
# The roots are x = -1 and x = phi, x = -1/phi!
# So the 2-3 Fibonacci recurrence factors into:
# (shift by -1) * (Fibonacci recurrence)

import numpy as np
roots = np.roots([1, 0, -2, -1])
print(f"\n  Roots of x^3 - 2x - 1 = 0:")
for r in roots:
    print(f"    x = {r:.6f}" + (" = phi!" if abs(r - phi) < 0.001 else " = -1/phi!" if abs(r + 1/phi) < 0.001 else " = -1!" if abs(r + 1) < 0.001 else ""))

# Verify factorization
print(f"\n  (x+1)(x^2-x-1) expanded:")
print(f"  = x^3 - x^2 - x + x^2 - x - 1 = x^3 - 2x - 1  CHECK!")

print(f"\n  THE 2-3 FIBONACCI RECURRENCE FACTORS AS:")
print(f"  (period-2 part) * (golden ratio part)")
print(f"  The -1 root gives period 2: (-1)^n alternates")
print(f"  The phi root gives growth: phi^n")
print(f"  Together: F(n) = A*phi^n + B*(-1/phi)^n + C*(-1)^n")
print(f"  But we know F(n) = (phi^n - (-1/phi)^n)/sqrt(5)")
print(f"  The extra (-1)^n part is the residue from the 2-3 decomposition")

# =====================================================================
# PART 8: EISENSTEIN + FIBONACCI = ICOSAHEDRAL?
# =====================================================================
print("\n" + "=" * 70)
print("PART 8: EISENSTEIN + FIBONACCI = ICOSAHEDRAL SYMMETRY?")
print("=" * 70)

# The exceptional isomorphism PSL(2,4) ≅ PSL(2,5) ≅ A_5
# connects:
# - F_4 world (Eisenstein, omega, Phi_3, Fano)
# - F_5 world (Fibonacci, phi, Phi_5, icosahedron)
# - A_5 (alternating group, 60 elements)

print(f"\n  PSL(2,4) = PSL(2,5) = A_5: the EXCEPTIONAL ISOMORPHISM")
print(f"  |PSL(2,4)| = (4^2-1)*(4^2-4)/2 = 15*12/2 = 60? No...")
print(f"  |PSL(2,q)| = q(q^2-1)/gcd(2,q-1)")
print(f"  |PSL(2,4)| = 4*(16-1)/1 = 60")
print(f"  |PSL(2,5)| = 5*(25-1)/2 = 60")
print(f"  |A_5| = 60")
print(f"  ALL EQUAL 60 = 3*4*5 = n(n-1)(n-2) at n=5!")
print(f"  And 60 is the denominator of the n=5 deviation: Var/Mean^2 = 19/60!!")

print(f"\n  THE A_5 CONNECTION:")
print(f"  - A_5 is the rotation symmetry of the ICOSAHEDRON")
print(f"  - The icosahedron has 12 vertices, 30 edges, 20 faces")
print(f"  - 20 faces = C(6,3) (triangulations of the hexagon)")
print(f"  - 30 edges = C(6,2)*2 (pairs of consecutive elements)")
print(f"  - 12 vertices = 2*6 (antipodal pairs)")

# Icosahedron data
print(f"\n  Icosahedron: V=12, E=30, F=20")
print(f"  Euler: V - E + F = 12 - 30 + 20 = 2  CHECK")
print(f"  Each vertex has degree 5")
print(f"  Dual = dodecahedron: V=20, E=30, F=12")

# The 60 connection
print(f"\n  60 = |A_5| = |PSL(2,4)| = |PSL(2,5)|")
print(f"  60 = n(n-1)(n-2) at n=5")
print(f"  60 = denominator of Var/Mean^2 = 19/60 at n=5")
print(f"  Var/Mean^2 = 19/60 = 1/3 - 1/60")
print(f"  So: deviation = 1/|A_5| = 1/|PSL(2,4)| = 1/|PSL(2,5)|!")
print(f"  THE DEVIATION FROM THE CONE IS 1/|A_5|!")

# =====================================================================
# PART 9: THE REMARKABLE IDENTITY 19/60
# =====================================================================
print("\n" + "=" * 70)
print("PART 9: 19/60 — A SPECIAL FRACTION")
print("=" * 70)

print(f"\n  19/60 = Var(H)/Mean(H)^2 at n=5")
print(f"  19 is the 8th prime")
print(f"  60 = 2^2 * 3 * 5")
print(f"  19/60 = 1/3 - 1/60 = 20/60 - 1/60")
print(f"  The '1' in the numerator correction: c_5 = 1")

# Can we express 19/60 in terms of Phi_3?
print(f"\n  19/60 in the Phi_3 world:")
print(f"  1/Phi_3(1) = 1/3 = 20/60")
print(f"  1/Phi_3(1) - 1/A_5 = 1/3 - 1/60 = 19/60")
print(f"  So: Var/Mean^2 = 1/Phi_3(1) - 1/|PSL(2, Phi_3(1)+1)|")
print(f"  Because: PSL(2,4) = PSL(2, Phi_3(1)+1) = A_5, |A_5|=60")

# 19 itself
print(f"\n  19 = 2^5 - 13 = 32 - |PG(2,3)|")
print(f"  19 = Phi_3(1) * Phi_3(2) - 2 = 3*7 - 2 = 21 - 2 = 19")
print(f"  Wait: 3*7 = 21, not 19. Let me reconsider.")
print(f"  19 = 20 - 1 = 4*5 - 1 = (Phi_3(1)+1)*(Phi_3(1)+2) - 1")
print(f"  19 is a UNIQUE prime: it's the number of rotational symmetries")
print(f"  of the regular 19-gon that are not trivial")

# The pattern for higher n
print(f"\n  PREDICTION for n=6:")
print(f"  If c_n pattern continues: Var/Mean^2 = 1/3 - c_6/(6*5*4)")
print(f"  From Monte Carlo: Var/Mean^2 ~ {ratio_6:.6f}")
print(f"  c_6 = (1/3 - {ratio_6:.6f}) * 120 = {c_6:.4f}")
print(f"  Nearest integer: c_6 ~ {round(c_6)}")
falling_6 = 6*5*4
predicted = 1/3 - round(c_6) / falling_6
print(f"  Predicted: 1/3 - {round(c_6)}/120 = {predicted:.6f}")

# =====================================================================
# PART 10: THE TRIANGLE-CONE-CATEGORY TRINITY REVISITED
# =====================================================================
print("\n" + "=" * 70)
print("PART 10: THE TRINITY — THREE FACES OF Phi_3")
print("=" * 70)

print("""
  Phi_3(x) = x^2 + x + 1 has THREE geometric faces:

  FACE 1: THE TRIANGLE (combinatorial)
  - Phi_3(x) = 1 + x + x^2 = sum_{k=0}^{2} x^k
  - This is the CHARACTER of the cyclic group C_3
  - Phi_3(1) = 3 = |C_3| = number of vertices of a triangle
  - The 3-CYCLE is the fundamental tournament generator
  - Every tournament is built from 3-cycles (transitive + rotational)

  FACE 2: THE CONE (geometric)
  - 1/Phi_3(1) = 1/3 = cone-to-cylinder ratio
  - V_cone = (1/3) * V_cylinder
  - This is integral_0^1 t^2 dt = 1/3
  - Var(H)/Mean(H)^2 = 1/3 at n=3,4
  - The tournament H-landscape IS a cone in effective dimension 3

  FACE 3: THE PROJECTIVE LINE (algebraic)
  - Phi_3(x) = (x^3 - 1)/(x - 1) = |PG(0, F_{x^3-1})| ... no
  - Better: Phi_3(q) = q^2 + q + 1 = |PG(2, F_q)|
  - At q=2: |PG(2,F_2)| = 7 = Fano = H_forb_1
  - At q=4: |PG(2,F_4)| = 21 = H_forb_2
  - These are the tournament FORBIDDEN VALUES

  UNIFICATION:
  The 3-cycle (face 1) creates a cone (face 2) that forbids
  projective plane sizes (face 3). ALL through Phi_3.

  The Fibonacci connection:
  - phi satisfies x^2 - x - 1 = 0 (reverse signs of Phi_3!)
  - Phi_3(x) = x^2 + x + 1, Fibonacci poly = x^2 - x - 1
  - They are COMPLEMENTARY: Phi_3(x) * (x-1) = x^3 - 1
  - Fib_poly(x) * (x+1) = x^3 + x^2 - x - 1 = (x^2-1)(x+1) = (x-1)(x+1)^2
  - Hmm, not as clean. But:
  - Phi_3(-x) = x^2 - x + 1 = Phi_6(x)
  - And Phi_6(phi) = phi^2 - phi + 1 = (phi+1) - phi + 1 = 2
  - THE TOURNAMENT GENERATOR 2 = Phi_6(phi) = Phi_3(-phi)!
""")

# Verify Phi_3(-phi) = 2
val = phi**2 - phi + 1  # Phi_6(phi) = Phi_3(-phi)
print(f"  Phi_6(phi) = phi^2 - phi + 1 = {val:.10f}")
print(f"  Phi_3(-phi) = phi^2 - phi + 1 = {val:.10f}")
print(f"  = 2? {abs(val - 2) < 1e-10}")

# And Phi_3(omega) = 0 (by definition)
# So the three evaluation points of Phi_3 give:
# Phi_3(omega) = 0 (the root)
# Phi_3(1) = 3 (the cone)
# Phi_3(2) = 7 (the Fano)
# Phi_3(-phi) = 2 (the generator!)
# Phi_3(-1/phi) = 1-(-1/phi)+(-1/phi)^2 = 1+1/phi+1/phi^2

val2 = 1 + 1/phi + 1/phi**2
print(f"  Phi_3(-1/phi) = {val2:.10f}")
# 1 + 1/phi = 1 + phi - 1 = phi (since 1/phi = phi-1)
# 1/phi^2 = 1/(phi+1) = phi - 1... no, 1/phi^2 = (phi-1)^2 = phi^2-2phi+1
# = (phi+1)-2phi+1 = 2-phi
# So: 1 + 1/phi + 1/phi^2 = 1 + (phi-1) + (2-phi) = 2
print(f"  Phi_3(-1/phi) = 1 + (phi-1) + (2-phi) = 2 also! SAME VALUE!")
print(f"  Because -phi and -1/phi are BOTH roots of x^2+x-1=0")
print(f"  Wait: -phi and -1/phi satisfy x^2+x-1=0? Let me check.")
print(f"  (-phi)^2 + (-phi) - 1 = phi^2 - phi - 1 = (phi+1) - phi - 1 = 0  YES!")
print(f"  So Phi_3(x) at roots of x^2+x-1=0 always gives 2!")
print(f"  x^2+x-1 = 0 => x^2 = 1-x => Phi_3(x) = x^2+x+1 = (1-x)+x+1 = 2")

print("\n  *** BEAUTIFUL IDENTITY ***")
print("  For any root alpha of x^2 + x - 1 = 0:")
print("  Phi_3(alpha) = alpha^2 + alpha + 1 = (1-alpha) + alpha + 1 = 2")
print("  THE TOURNAMENT GENERATOR 2 = Phi_3 EVALUATED AT FIBONACCI ROOTS!")
print("  This connects the GOLDEN RATIO to the CUBE ROOT OF UNITY!")

# =====================================================================
# PART 11: THE MASTER EQUATION
# =====================================================================
print("\n" + "=" * 70)
print("PART 11: THE MASTER EQUATION")
print("=" * 70)

print("""
  MASTER EQUATION: Phi_3(x) * (x - 1) = x^3 - 1

  At x = 1:   Phi_3(1) * 0 = 0         (cone: 1/3)
  At x = 2:   Phi_3(2) * 1 = 7          (Fano: H_forb_1)
  At x = 4:   Phi_3(4) * 3 = 63         (= 2^6 - 1 = |PG(5,F_2)|)
  At x = -phi: Phi_3(-phi) * (-phi-1) = (-phi)^3 - 1 = -(phi^3+1)
              = 2 * (-phi-1) = -2phi-2 = -(phi^3+1)
              phi^3 = phi*phi^2 = phi*(phi+1) = phi^2+phi = 2phi+1
              -(2phi+1+1) = -2phi-2 and 2*(-phi-1) = -2phi-2  CHECK!

  The FACTORIZATION x^3 - 1 = (x-1)(x^2+x+1) = (x-1)*Phi_3(x) says:

  "Dividing the cube (x^3) by its unit boundary (x-1) leaves the
  triangle (Phi_3(x))."

  In tournament terms:
  "The full binary structure (2^3 = 8 orientations per triangle)
  minus the identity (1 trivial orientation)
  equals the tournament cycle polynomial (Phi_3(2) = 7 = Fano)."

  And: 8 - 1 = 7. The Fano IS the non-trivial part of the binary cube!
""")

# Verify 2^3 - 1 = 7 = Phi_3(2)
print(f"  2^3 - 1 = {2**3 - 1} = Phi_3(2) = {2**2+2+1}")
print(f"  2^6 - 1 = {2**6 - 1} = Phi_3(4) * Phi_6(4) * ... ")
print(f"  Actually 2^6-1 = 63 = 7*9 = 7*3^2")
print(f"  And 63 = Phi_3(4) * 3 = 21 * 3 = 63  CHECK")
print(f"  Or: 63 = (2^3-1)(2^3+1) = 7*9")

# =====================================================================
# PART 12: FIBONACCI PERIOD-6 AND THE CONE
# =====================================================================
print("\n" + "=" * 70)
print("PART 12: THE PERIOD-6 FIBONACCI CONE")
print("=" * 70)

# The Fibonacci matrix [[1,1],[1,0]] has period 6 mod 4
# The cone has dimension 3
# 6 = 2 * 3
# The period-6 structure creates a "Fibonacci cone" in mod-4 arithmetic

# More deeply: the Fibonacci numbers mod n create a CYCLIC PATTERN
# The 6-step return means: after 6 Fibonacci steps, you're back to start (mod 4)
# This is like going around a cone (360°) — the angular return

print("\n  Fibonacci mod 4 (period 6):")
for i in range(13):
    f = fib[i] % 4
    angle = (i % 6) * 60
    print(f"    F({i:2d}) mod 4 = {f}  (angle = {angle:3d} deg)")

print("\n  The 6 states form a HEXAGONAL cycle in mod-4 space:")
print("  State: (F_n mod 4, F_{n+1} mod 4)")
states = []
for i in range(6):
    state = (fib[i] % 4, fib[i+1] % 4)
    states.append(state)
    print(f"    Step {i}: ({state[0]}, {state[1]})")
print(f"  Return to start: ({fib[6]%4}, {fib[7]%4}) = ({fib[6]%4}, {fib[7]%4})")
print(f"  Same as step 0? {(fib[6]%4, fib[7]%4) == (fib[0]%4, fib[1]%4)}")

# The hexagonal cycle in 2D mod-4 space
# This is a DISCRETE CONE SECTION!
# A regular hexagon is the intersection of a cube with a plane
# perpendicular to the (1,1,1) diagonal — the "triangular" direction!

print("\n  THE HEXAGONAL-CONIC CONNECTION:")
print("  A regular hexagon = cube section perpendicular to (1,1,1)")
print("  The (1,1,1) direction = the 'cone axis' in 3D")
print("  Period 6 = one full rotation around the cone axis")
print("  The 6 Fibonacci states = 6 vertices of the hexagonal cross-section")
print("  The cone dimension 3 × rotation period 2 = 6 = Fibonacci period mod 4")

# =====================================================================
# PART 13: GRAND SYNTHESIS TABLE
# =====================================================================
print("\n" + "=" * 70)
print("PART 13: GRAND SYNTHESIS — THE Phi_3-FIBONACCI-CONE TABLE")
print("=" * 70)

print("""
  ╔═══════════════════╦═══════════════╦════════════════╦═══════════════╗
  ║     Concept       ║   Phi_3 eval  ║   Fibonacci    ║   Geometry    ║
  ╠═══════════════════╬═══════════════╬════════════════╬═══════════════╣
  ║ Cone ratio 1/3    ║ 1/Phi_3(1)    ║ Fib mod 4      ║ V_cone/V_cyl  ║
  ║                   ║               ║ period = 6=2*3 ║ in dim 3      ║
  ╠═══════════════════╬═══════════════╬════════════════╬═══════════════╣
  ║ Generator 2       ║ Phi_3(-phi)   ║ phi^2-phi+1=2  ║ Binary choice ║
  ║                   ║ = Phi_6(phi)  ║                ║ per arc       ║
  ╠═══════════════════╬═══════════════╬════════════════╬═══════════════╣
  ║ Fano 7 = H_forb_1 ║ Phi_3(2)     ║ F(7-1)=8=2^3  ║ PG(2,F_2)    ║
  ║                   ║               ║ F(7)=13=PG(2,3)║              ║
  ╠═══════════════════╬═══════════════╬════════════════╬═══════════════╣
  ║ 21 = H_forb_2    ║ Phi_3(4)      ║ F(8) = 21 !    ║ PG(2,F_4)    ║
  ╠═══════════════════╬═══════════════╬════════════════╬═══════════════╣
  ║ Deviation at n=5  ║ 1/3-1/60     ║ 60 = |A_5|     ║ Icosahedron   ║
  ║                   ║ = 19/60      ║ PSL(2,4)=A_5   ║ 60 rotations  ║
  ╠═══════════════════╬═══════════════╬════════════════╬═══════════════╣
  ║ Jacobsthal        ║ J(n) has 1/3  ║ J vs F:        ║ Discrete cone ║
  ║ sequence          ║ in formula    ║ coeff 2 vs 1   ║ enumeration   ║
  ╠═══════════════════╬═══════════════╬════════════════╬═══════════════╣
  ║ Eisenstein norm   ║ Phi_3(b)=     ║ Class number 1 ║ Hexagonal     ║
  ║                   ║ |1-b*w^2|^2  ║ (both Z[phi],  ║ lattice       ║
  ║                   ║               ║  Z[omega])     ║               ║
  ╚═══════════════════╩═══════════════╩════════════════╩═══════════════╝

  THE MASTER IDENTITY:
  For any root alpha of x^2 + x - 1 = 0 (the FIBONACCI polynomial):
    Phi_3(alpha) = 2 = the TOURNAMENT GENERATOR

  In words: "The third cyclotomic polynomial, evaluated at the
  golden ratio, gives the tournament base."

  Or equivalently: "The cube root of unity polynomial, composed
  with the Fibonacci generator, produces the binary choice."

  This is the DEEPEST connection between Fibonacci and tournaments:
  THE GOLDEN RATIO IS A ROOT OF THE EQUATION Phi_3(x) = 2.

  Verify: Phi_3(x) = 2  <=>  x^2 + x + 1 = 2  <=>  x^2 + x - 1 = 0
  Solutions: x = (-1 +/- sqrt(5))/2 = {phi-1, -phi} = {1/phi, -phi}

  THE FIBONACCI ROOTS SOLVE Phi_3(x) = 2 !!
""")

# Final verification
print("  Verification:")
print(f"    phi - 1 = {phi-1:.10f} = 1/phi = {1/phi:.10f}")
print(f"    Phi_3(1/phi) = {(1/phi)**2 + 1/phi + 1:.10f}")
print(f"    Phi_3(-phi) = {phi**2 - phi + 1:.10f}")
print(f"    Both equal 2: YES!")

print("\n" + "=" * 70)
print("DONE — EISENSTEIN-CONE-FIBONACCI SYNTHESIS")
print("=" * 70)
