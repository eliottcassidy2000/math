#!/usr/bin/env python3
"""flip_markov_s116n.py — Spectral theory of the tournament flip Markov chain.

The flip chain: random walk on {0,1}^m where m = C(n-1,2) tiling bits.
Each step flips one random bit. Transition eigenvalues: lambda_k = 1-2k/m.

KEY OBSERVATION: The denominators of these eigenvalues (in lowest terms)
are the divisors of m (odd m) or m/2 (even m). These denominators
CONTAIN the forbidden H values and max H values from smaller n!

n=6: m=10, s=5,  denoms = {1, 5}. The GOLDEN PRIME.
n=7: m=15,       denoms = {1, 3, 5, 15}. 15 = max H at n=5!
n=8: m=21,       denoms = {1, 3, 7, 21}. BOTH FORBIDDEN VALUES!
n=11: m=45,s=45, denoms include 45 = max H at n=6!

The Walsh-Fourier decomposition of H gives the EXACT evolution under the chain:
  E[H(X_t) | X_0] = sum_S hat{H}(S) * (1-2|S|/m)^t * chi_S(X_0)

Session: kind-pasteur-2026-03-16-S116n32
"""
import sys
sys.stdout.reconfigure(line_buffering=True)
from itertools import permutations
from collections import Counter, defaultdict
from math import gcd, comb
from fractions import Fraction

print()
print("  SPECTRAL THEORY OF THE TOURNAMENT FLIP MARKOV CHAIN")
print()
print("=" * 70)
print()

N = 6

# ============================================================
print("  I. THE EIGENVALUE DENOMINATORS")
print("  " + "-" * 50)
print()

print("  Flip chain on tiling space {0,1}^m, m = C(n-1, 2).")
print("  Eigenvalues: lambda_k = (m - 2k)/m for k = 0, ..., m.")
print()

for n in range(3, 14):
    m = (n-1)*(n-2)//2
    # For even m: eigenvalues = (m/2 - k)/(m/2), denoms = divisors of m/2
    # For odd m: eigenvalues = (m - 2k)/m, denoms = divisors of m
    if m % 2 == 0:
        s = m // 2
        label = f"m/2={s}"
    else:
        s = m
        label = f"m={s}"

    # Factor s
    temp = s
    factors = []
    for p in range(2, temp+1):
        while temp % p == 0:
            factors.append(p)
            temp //= p
    if not factors:
        factors = [1]
    factor_str = '*'.join(str(f) for f in factors) if s > 1 else '1'

    # Divisors of s
    divs = sorted(d for d in range(1, s+1) if s % d == 0)

    # Which are prime?
    prime_divs = [d for d in divs if d > 1 and all(d % p != 0 for p in range(2, d))]

    print(f"  n={n:2d}: m={m:3d}, {label}, factor={factor_str:12s}, "
          f"denom divisors={divs}")

print()

# ============================================================
print("  II. THE DENOMINATOR-TOURNAMENT CONNECTION")
print("  " + "-" * 50)
print()

# Key values from the project
forbidden = {7, 21}
max_H = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661}

print("  Forbidden H values: {7, 21}")
print("  Max H values: ", dict(max_H))
print()

for n in range(3, 14):
    m = (n-1)*(n-2)//2
    s = m // 2 if m % 2 == 0 else m
    divs = set(d for d in range(1, s+1) if s % d == 0)

    # Check which forbidden/max values appear as denominators
    forbidden_in = divs & forbidden
    max_in = {k: v for k, v in max_H.items() if v in divs and k < n}

    interesting = ""
    if forbidden_in:
        interesting += f" FORBIDDEN {forbidden_in}"
    if max_in:
        interesting += f" MAX_H {max_in}"

    if interesting:
        print(f"  n={n:2d}: denominators contain{interesting}")

print()

# ============================================================
print("  III. WALSH-FOURIER DECOMPOSITION OF H AT n=6")
print("  " + "-" * 50)
print()

# Compute H for all 2^10 = 1024 tilings
m = (N-1)*(N-2)//2  # 10

tiling_arcs = []
for skip in range(2, N):
    for start in range(N - skip):
        tiling_arcs.append((start, start + skip))

def tiling_to_adj(bits_int, n):
    adj = [[0]*n for _ in range(n)]
    for i in range(n-1):
        adj[i][i+1] = 1
    for idx, (a, b) in enumerate(tiling_arcs):
        if (bits_int >> idx) & 1:
            adj[b][a] = 1  # backward
        else:
            adj[a][b] = 1  # forward
    return adj

def H_dp(adj, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n): dp[(1 << v) * n + v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp[S * n + v]
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                if adj[v][u]: dp[(S | (1 << u)) * n + u] += val
    return sum(dp[((1 << n) - 1) * n + v] for v in range(n))

print(f"  Computing H for all 2^{m} = {1<<m} tilings...")
H_table = [0] * (1 << m)
for bits in range(1 << m):
    H_table[bits] = H_dp(tiling_to_adj(bits, N), N)

mean_H = sum(H_table) / (1 << m)
print(f"  Mean H = {mean_H}")
print()

# Walsh-Fourier transform: hat{H}(S) = (1/2^m) sum_x H(x) (-1)^{<S,x>}
# where <S,x> = sum of x_i for i in S (parity of overlap)

print("  Computing Walsh-Fourier coefficients...")
# Use the definition directly (2^m * 2^m = 2^20 = ~1M operations)
walsh = [0.0] * (1 << m)
for S in range(1 << m):
    total = 0
    for x in range(1 << m):
        # parity of S & x
        overlap = bin(S & x).count('1')
        sign = 1 if overlap % 2 == 0 else -1
        total += sign * H_table[x]
    walsh[S] = total / (1 << m)

print(f"  Done. {sum(1 for w in walsh if abs(w) > 1e-10)} nonzero Walsh coefficients.")
print()

# Verify: reconstruction
errors = 0
for x in range(1 << m):
    val = 0.0
    for S in range(1 << m):
        overlap = bin(S & x).count('1')
        sign = 1 if overlap % 2 == 0 else -1
        val += walsh[S] * sign
    if abs(val - H_table[x]) > 0.01:
        errors += 1
if errors == 0:
    print(f"  VERIFIED: Walsh reconstruction matches all {1<<m} values!")
else:
    print(f"  {errors} reconstruction errors!")
print()

# ============================================================
print("  IV. WALSH SPECTRUM BY DEGREE")
print("  " + "-" * 50)
print()

# Group Walsh coefficients by |S| (the "degree" in Walsh basis)
degree_coeffs = defaultdict(list)
for S in range(1 << m):
    deg = bin(S).count('1')
    if abs(walsh[S]) > 1e-10:
        degree_coeffs[deg].append((S, walsh[S]))

for deg in sorted(degree_coeffs.keys()):
    coeffs = degree_coeffs[deg]
    vals = [c for _, c in coeffs]
    # The eigenvalue for this degree
    lam = Fraction(m - 2*deg, m)
    lam_reduced = Fraction(lam.numerator, lam.denominator)
    denom = lam_reduced.denominator

    # Walsh variance contribution: sum of hat{H}(S)^2 for this degree
    walsh_var = sum(c**2 for c in vals)

    print(f"  Degree {deg:2d}: {len(coeffs):3d} terms, "
          f"eigenvalue = {str(lam_reduced):>6s} (denom {denom}), "
          f"Walsh variance = {walsh_var:8.2f}, "
          f"range = [{min(vals):+.2f}, {max(vals):+.2f}]")
print()

# Total Walsh variance (should equal Var(H))
total_walsh_var = sum(c**2 for S, c in enumerate(walsh) if S > 0 and abs(c) > 1e-10)
# Fix: iterate properly
total_walsh_var = sum(walsh[S]**2 for S in range(1, 1 << m) if abs(walsh[S]) > 1e-10)
var_H = sum((H - mean_H)**2 for H in H_table) / (1 << m)
print(f"  Total Walsh variance: {total_walsh_var:.2f}")
print(f"  Var(H): {var_H:.2f}")
print(f"  Parseval check: {abs(total_walsh_var - var_H) < 0.1}")
print()

# ============================================================
print("  V. THE FIVE EIGENVALUE CLASSES (n=6, denom=5)")
print("  " + "-" * 50)
print()

# For n=6, m=10, s=5: eigenvalues are (5-k)/5 for k=0,...,10
# Denominator classes: 1 (for k=0,5,10) and 5 (for k=1,2,3,4,6,7,8,9)
# But in reduced form: lambda = (5-k)/5
# k=0: 1 (denom 1)
# k=1: 4/5 (denom 5)
# k=2: 3/5 (denom 5)
# k=3: 2/5 (denom 5)
# k=4: 1/5 (denom 5)
# k=5: 0 (denom 1)
# k=6: -1/5 (denom 5)
# k=7: -2/5 (denom 5)
# k=8: -3/5 (denom 5)
# k=9: -4/5 (denom 5)
# k=10: -1 (denom 1)

print("  Eigenvalue spectrum:")
for k in range(m+1):
    lam = Fraction(m - 2*k, m)
    n_coeffs = sum(1 for S in range(1<<m)
                   if bin(S).count('1') == k and abs(walsh[S]) > 1e-10)
    total_power = sum(walsh[S]**2 for S in range(1<<m)
                      if bin(S).count('1') == k and abs(walsh[S]) > 1e-10)
    mult = comb(m, k)
    print(f"  k={k:2d}: lambda={str(lam):>6s}, mult=C({m},{k})={mult:4d}, "
          f"active={n_coeffs:3d}, power={total_power:8.2f}")
print()

# ============================================================
print("  VI. MIXING TIME AND H RELAXATION")
print("  " + "-" * 50)
print()

# Spectral gap = 1 - |lambda_1| = 1 - (1-2/m) = 2/m
gap = Fraction(2, m)
mixing_time = m // 2  # rough estimate: 1/gap

print(f"  Spectral gap: {gap} = {float(gap):.4f}")
print(f"  Approximate mixing time: {mixing_time} steps")
print()

# Simulate: start from transitive (all zeros), track expected H
print("  H relaxation from transitive (all-zero tiling):")
x0 = 0  # transitive = all forward

# E[H(X_t) | X_0 = transitive] = sum_S hat{H}(S) * (1-2|S|/m)^t * chi_S(0)
# chi_S(0) = (-1)^0 = 1 for all S (since all bits of 0 are 0)
# So E[H(X_t)] = sum_S hat{H}(S) * ((m-2|S|)/m)^t

print(f"  {'t':>4s}  {'E[H(Xt)]':>10s}  {'|E[H]-mean|':>12s}")
for t in range(0, 25):
    expected_H = 0
    for S in range(1 << m):
        if abs(walsh[S]) < 1e-12:
            continue
        deg = bin(S).count('1')
        lam = (m - 2*deg) / m
        expected_H += walsh[S] * (lam ** t)
    print(f"  {t:4d}  {expected_H:10.2f}  {abs(expected_H - mean_H):12.4f}")
print()

# Also from the anti-transitive (all ones)
print("  H relaxation from anti-transitive (all-one tiling):")
# chi_S(1...1) = (-1)^|S|
print(f"  {'t':>4s}  {'E[H(Xt)]':>10s}  {'|E[H]-mean|':>12s}")
for t in range(0, 25):
    expected_H = 0
    for S in range(1 << m):
        if abs(walsh[S]) < 1e-12:
            continue
        deg = bin(S).count('1')
        lam = (m - 2*deg) / m
        sign = (-1) ** deg
        expected_H += walsh[S] * (lam ** t) * sign
    print(f"  {t:4d}  {expected_H:10.2f}  {abs(expected_H - mean_H):12.4f}")
print()

# ============================================================
print("  VII. THE DOMINANT WALSH MODES")
print("  " + "-" * 50)
print()

# Sort Walsh coefficients by magnitude
all_walsh = [(S, walsh[S]) for S in range(1 << m) if abs(walsh[S]) > 1e-10]
all_walsh.sort(key=lambda x: abs(x[1]), reverse=True)

print(f"  Top 20 Walsh coefficients (by magnitude):")
print(f"  {'Rank':>4s}  {'|S|':>4s}  {'S (bits)':>12s}  {'hat_H(S)':>10s}  {'eigenvalue':>10s}  {'arcs':>30s}")
for rank, (S, coeff) in enumerate(all_walsh[:20], 1):
    deg = bin(S).count('1')
    bits = [i for i in range(m) if (S >> i) & 1]
    arcs = [tiling_arcs[i] for i in bits]
    lam = Fraction(m - 2*deg, m)
    S_str = format(S, f'0{m}b')[::-1]  # LSB first
    print(f"  {rank:4d}  {deg:4d}  {S_str:>12s}  {coeff:+10.4f}  {str(lam):>10s}  {arcs}")
print()

# ============================================================
print("  VIII. THE DENOMINATOR STRUCTURE")
print("  " + "-" * 50)
print()

# For n=6: all eigenvalues have denominator 1 or 5.
# Group Walsh power by denominator

denom_power = defaultdict(float)
for S in range(1 << m):
    if abs(walsh[S]) < 1e-10:
        continue
    deg = bin(S).count('1')
    lam = Fraction(m - 2*deg, m)
    denom = lam.denominator
    denom_power[denom] += walsh[S]**2

print(f"  Walsh power by eigenvalue denominator (n={N}):")
for d in sorted(denom_power.keys()):
    pct = 100 * denom_power[d] / var_H if var_H > 0 else 0
    print(f"  denom={d}: power={denom_power[d]:.2f} ({pct:.1f}% of variance)")
print()

# The key: what fraction of H's information is in the "denom=5" modes vs "denom=1" modes?
print("  INTERPRETATION:")
print(f"  Denom=1 modes: lambda in {{1, 0, -1}}. These are the 'trivial' modes.")
print(f"    k=0 (constant): hat_H = {walsh[0]:.2f} = mean H")
print(f"    k=5 (lambda=0): these modes are KILLED after 1 step.")
print(f"    k=10 (lambda=-1): these modes flip sign each step.")
print()
print(f"  Denom=5 modes: lambda in {{4/5, 3/5, 2/5, 1/5, -1/5, -2/5, -3/5, -4/5}}.")
print(f"    These are the 'nontrivial' modes that decay geometrically.")
print(f"    They carry {denom_power.get(5,0):.1f} / {var_H:.1f} = "
      f"{100*denom_power.get(5,0)/var_H:.1f}% of H's variance.")
print(f"    The denominator 5 = the GOLDEN PRIME.")
print()

# ============================================================
print("  IX. THE m SEQUENCE AND ITS PRIME STRUCTURE")
print("  " + "-" * 50)
print()

print("  m = C(n-1, 2) = (n-1)(n-2)/2. Odd part s determines denominators.")
print()
print(f"  {'n':>3s}  {'m':>4s}  {'s=odd(m)':>8s}  {'prime_factors(s)':>20s}  {'connection':>30s}")

connections = {
    3: "m=1",
    4: "3 = min cycle length",
    5: "3 = RAMIFIED prime",
    6: "5 = GOLDEN prime",
    7: "15 = max H at n=5",
    8: "21 = FORBIDDEN VALUE",
    9: "7 = FORBIDDEN VALUE",
    10: "9 = H value at n=5",
    11: "45 = max H at n=6",
    12: "55 = 5*11",
    13: "33 = H value at n=6",
}

for n in range(3, 14):
    m = (n-1)*(n-2)//2
    s = m
    while s % 2 == 0:
        s //= 2

    # Factor s
    temp = s
    factors = []
    for p in range(2, temp+1):
        while temp % p == 0:
            factors.append(p)
            temp //= p
    factor_str = '*'.join(str(f) for f in factors) if factors else '1'

    conn = connections.get(n, "")
    print(f"  {n:3d}  {m:4d}  {s:8d}  {factor_str:>20s}  {conn:>30s}")
print()

# ============================================================
print("  X. SYNTHESIS: THE GOLDEN PRIME IN THE MARKOV CHAIN")
print("  " + "-" * 50)
print()
print("  At n=6 (the PIVOT), the flip chain eigenvalue denominator is 5.")
print("  5 = disc(Q(phi)) = the golden prime.")
print("  This means: the MIXING DYNAMICS of tournament H are governed by")
print("  the golden prime. The nontrivial relaxation rates are k/5 for k=1,...,4.")
print()
print("  At n=8, the denominator is 21 = 3*7 = the FORBIDDEN product.")
print("  At n=7, the denominator is 15 = 3*5 = {RAMIFIED}*{GOLDEN}.")
print()
print("  The Markov chain denominator sequence (odd part of C(n-1,2)):")
print("  1, 3, 3, 5, 15, 21, 7, 9, 45, 55, 33, ...")
print("  These are NOT the tournament invariants themselves,")
print("  but the ARENA in which tournament dynamics plays out.")
print("  The prime structure of the arena determines which")
print("  frequencies are available for the spectral decomposition.")
print()
print("  The grand pattern:")
print("  {2,3,5} = arena at n=6 (golden prime controls mixing)")
print("  {2,3,7} = arena at n=8 (forbidden primes control mixing)")
print("  The Markov chain carries the SAME prime structure as the")
print("  tournament theory it describes. The medium IS the message.")
print()
