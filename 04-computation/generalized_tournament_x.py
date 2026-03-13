#!/usr/bin/env python3
"""
GENERALIZED TOURNAMENTS AT x = k(k-1) — THE 2-3 FRAMEWORK EXTENDED
opus-2026-03-14-S68

At x=2 (=2·1), the Fibonacci root is 2: BINARY tournaments.
At x=6 (=3·2), the root is 3: would this be TERNARY tournaments?
At x=12 (=4·3), the root is 4: QUATERNARY?

THESIS: The structure of tournament theory at x=2 generalizes to
"k-ary tournaments" at x=k(k-1), and each level has its own version
of Jacobsthal numbers, OCF, and perhaps Pfaffian identities.

But also: what does x=1 (root=φ, Fibonacci) mean?
And x=3 (root=(1+√13)/2 ≈ 2.303)?
"""

import numpy as np
from fractions import Fraction
import math
from functools import lru_cache
from itertools import combinations, permutations

print("=" * 78)
print("  GENERALIZED TOURNAMENTS — THE COMPLETE x-SPECTRUM")
print("=" * 78)

# ============================================================================
# PART 1: THE INDEPENDENCE POLYNOMIAL AT ARBITRARY x
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 1: I(G, x) AT DIFFERENT x — WHAT DOES IT COUNT?                      ║
╚══════════════════════════════════════════════════════════════════════════════╝

I(G, x) = Σ_k α_k · x^k where α_k = #{independent sets of size k}

At x=1: I(G,1) = total number of independent sets (including empty)
At x=2: I(G,2) = Σ 2^k α_k  (each vertex in the set contributes factor 2)

COMBINATORIAL MEANING at integer x:
  I(G, x) = #{proper x-colorings of G using colors from {0,1,...,x}} ???

  No! That's the CHROMATIC polynomial, not the independence polynomial.

  I(G, x) = #{functions f: V → {0,1,...,x} with f(v) = 0 for v not in
              some independent set S, and f(v) ∈ {1,...,x} for v ∈ S}

  In other words: assign labels 0,...,x to vertices, where
  - "0" means "not selected"
  - "1,...,x" means "selected with one of x labels"
  - No two adjacent vertices can both be selected (nonzero)

  At x=2: two types of "selection": {1, 2}, or "not selected": {0}
  So I(G, 2) counts LABELED independent set configurations with 3 states.

For TOURNAMENTS: H(T) = I(CG(T), 2) counts labeled-cycle configurations.
Each odd cycle can be: not used (0), or used with label 1 or label 2.
But cycles that share a vertex can't both be used.
""")

# ============================================================================
# PART 2: THE COMPLETE x-SPECTRUM FOR SMALL GRAPHS
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 2: PATH AND CYCLE INDEPENDENCE POLYNOMIALS — FULL SPECTRUM            ║
╚══════════════════════════════════════════════════════════════════════════════╝
""")

@lru_cache(maxsize=None)
def I_path(m, x):
    """I(P_m, x) via recurrence."""
    if m == 0: return 1
    if m == 1: return 1 + x
    return I_path(m-1, x) + x * I_path(m-2, x)

def I_cycle(m, x):
    """I(C_m, x) via I(P_{m-1}, x) + x·I(P_{m-3}, x)."""
    if m <= 2: return (1+x)**m
    return I_path(m-1, x) + x * I_path(m-3, x)

# Characteristic root at each x
def char_root(x_val):
    return (1 + (1 + 4*x_val)**0.5) / 2

print("The Fibonacci-Jacobsthal root spectrum:")
print(f"{'x':>5}  {'root':>10}  {'1+x':>5}  {'root/prev':>10}  {'I(P_5,x)':>10}  {'I(C_5,x)':>10}  {'C/P ratio':>10}")
print("-" * 75)

prev_root = None
for x_val in [0, 1, 2, 3, 4, 5, 6, 10, 12, 20]:
    root = char_root(x_val)
    ratio_str = f"{root/prev_root:.6f}" if prev_root and prev_root > 0 else "—"
    ip5 = I_path(5, x_val)
    ic5 = I_cycle(5, x_val)
    cp_ratio = ic5 / ip5 if ip5 != 0 else "—"
    print(f"{x_val:5d}  {root:10.6f}  {1+x_val:5d}  {ratio_str:>10}  {ip5:10d}  {ic5:10d}  {cp_ratio if isinstance(cp_ratio, str) else f'{cp_ratio:.6f}':>10}")
    prev_root = root

# ============================================================================
# PART 3: THE JACOBSTHAL GENERALIZATION
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 3: GENERALIZED JACOBSTHAL — J_x(n) = J_x(n-1) + x·J_x(n-2)         ║
╚══════════════════════════════════════════════════════════════════════════════╝

J_x(n) with J_x(0)=0, J_x(1)=1:
  Closed form: J_x(n) = (r^n - s^n)/(r-s)
  where r = (1+√(1+4x))/2, s = (1-√(1+4x))/2

At x=1: J_1(n) = F(n) (Fibonacci numbers)
At x=2: J_2(n) = J(n) (Jacobsthal numbers)
At x=6: J_6(n) = ???  (root = 3)
""")

def gen_jacobsthal(x_val, n_max):
    """Compute generalized Jacobsthal numbers J_x(0), ..., J_x(n_max)."""
    J = [0, 1]
    for n in range(2, n_max+1):
        J.append(J[-1] + x_val * J[-2])
    return J

# Display for several x values
for x_val in [1, 2, 6, 12]:
    root = char_root(x_val)
    J = gen_jacobsthal(x_val, 12)
    print(f"\nx={x_val} (root={root:.4f}):")
    print(f"  J_{x_val}(n):", J[:13])

    # Check closed form
    r = (1 + (1 + 4*x_val)**0.5) / 2
    s = (1 - (1 + 4*x_val)**0.5) / 2
    print(f"  Closed form: J_{x_val}(n) = ({r:.4f}^n - ({s:.4f})^n) / {r-s:.4f}")
    print(f"  Denominator of closed form = √(1+4·{x_val}) = √{1+4*x_val} = {(1+4*x_val)**0.5:.6f}")

    # Check if J_x(n) has nice divisibility
    for n in range(1, 13):
        # For Jacobsthal at x=2: J(n) = (2^n - (-1)^n)/3
        # For x=6: J_6(n) = (3^n - (-2)^n)/5
        if x_val == 2:
            expected = (2**n - (-1)**n) // 3
            assert J[n] == expected, f"Mismatch at n={n}"
        elif x_val == 6:
            # r=3, s=-2, r-s=5
            expected = (3**n - (-2)**n) // 5
            if abs(J[n] - expected) > 0.1:
                print(f"  MISMATCH at n={n}: J={J[n]}, expected={expected}")

    if x_val == 6:
        print(f"  J_6(n) = (3^n - (-2)^n) / 5  — VERIFIED ✓")
    elif x_val == 12:
        # r=4, s=-3, r-s=7
        print(f"  J_12(n) = (4^n - (-3)^n) / 7")

    # Pattern: for x = k(k-1), root = k, secondary root = -(k-1)
    # J_x(n) = (k^n - (-(k-1))^n) / (2k-1)
    # Denominator = 2k-1 (always odd)

# ============================================================================
# PART 4: THE UNIVERSAL FORMULA
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 4: THE UNIVERSAL FORMULA — EVERYTHING AT x = k(k-1)                  ║
╚══════════════════════════════════════════════════════════════════════════════╝

For x = k(k-1):
  Primary root: r = k
  Secondary root: s = -(k-1)
  r - s = 2k - 1

  J_x(n) = (k^n - (-(k-1))^n) / (2k-1)

The pattern of 2 and 3 at x=2:
  k=2: r=2, s=-1, denom=3     J(n) = (2^n - (-1)^n) / 3
  k=3: r=3, s=-2, denom=5     J(n) = (3^n - (-2)^n) / 5
  k=4: r=4, s=-3, denom=7     J(n) = (4^n - (-3)^n) / 7

TRINITY AT LEVEL k:
  The number k plays the role of 2 (edge orientations → k orientations)
  The number k+1 plays the role of 3 (=1+k, the isolated vertex weight)
  The denominator 2k-1 plays the role of 3 (=2·2-1)

For k=2: (2, 3, 3) — the 2-3 of tournaments, with denominator = 1+x
For k=3: (3, 4, 5) — the 3-4-5 Pythagorean triple!
For k=4: (4, 5, 7) — the pattern breaks the Pythagorean connection
""")

print("The generalized trinity at each level k:")
print(f"{'k':>3}  {'x=k(k-1)':>8}  {'root k':>7}  {'1+x':>5}  {'denom 2k-1':>11}  {'J_x(6)':>10}  {'k^6/denom':>12}")
print("-" * 65)

for k in range(2, 9):
    x_val = k * (k - 1)
    denom = 2 * k - 1
    J6 = gen_jacobsthal(x_val, 6)[6]
    k6_over_d = k**6 // denom if k**6 % denom == 0 else f"≈{k**6/denom:.1f}"
    print(f"{k:3d}  {x_val:8d}  {k:7d}  {1+x_val:5d}  {denom:11d}  {J6:10d}  {k6_over_d:>12}")

# ============================================================================
# PART 5: THE x-GENERALIZED OCF
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 5: THE x-GENERALIZED OCF AND "k-ARY TOURNAMENTS"                     ║
╚══════════════════════════════════════════════════════════════════════════════╝

In a binary tournament (x=2), each edge has 2 orientations.
H(T) = I(CG(T), 2) = #{Hamiltonian paths of T}.

What would a "k-ary tournament" look like?
  • Each edge has k orientations (labeled 1,...,k)
  • A "Hamiltonian path" might require a consistent ordering
  • H_k(T) = I(CG(T), k(k-1)) ???

Actually, the connection is more subtle. Let's think about it:

At x=2: I(CG, 2) counts HP because:
  - Each edge-orientation corresponds to one of 2 choices
  - Independent sets in CG correspond to acyclic subsets
  - The x=2 evaluation weights by 2^|S|

At x=k(k-1): I(CG, k(k-1)) would count... something with k(k-1)-labeled
independent sets in CG. This doesn't directly correspond to k-ary tournaments.

The MORE natural generalization might be:
  • Consider digraphs where each ordered pair (i,j) with i≠j gets a
    label from {0, 1, ..., k-1} (representing k possible relationships)
  • A "tournament" restricts: exactly one of label(i,j), label(j,i) is nonzero
  • For k=2: either i→j or j→i (standard tournament)
  • For k=3: i→j with strength 1 or 2, or j→i with strength 1 or 2
""")

# Compute I(CG, x) for all tournaments at various x
def adj_matrix(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def count_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def find_odd_cycles(A, n):
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(A[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    mi = perm.index(min(perm))
                    canon = perm[mi:] + perm[:mi]
                    cycles.append(canon)
    return list(set(cycles))

def build_cg(cycles):
    nc = len(cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(cycles[i]) & set(cycles[j]):
                adj[i][j] = adj[j][i] = True
    return adj

def indep_poly_coeffs(adj, nc):
    alpha = [0]*(nc+1)
    for mask in range(1 << nc):
        verts = [i for i in range(nc) if mask & (1 << i)]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False; break
            if not ok: break
        if ok:
            alpha[len(verts)] += 1
    return alpha

def eval_indep_poly(alpha, x_val):
    return sum(alpha[k] * x_val**k for k in range(len(alpha)))

print("I(CG(T), x) for n=5 tournaments at multiple x values:")
print(f"{'H(2)':>5}  {'α':>12}  {'I(·,1)':>7}  {'I(·,2)':>7}  {'I(·,3)':>7}  {'I(·,6)':>7}  {'I(·,12)':>8}")
print("-" * 65)

n = 5
m = n*(n-1)//2
seen = {}
for bits in range(1 << m):
    A = adj_matrix(bits, n)
    H = count_hp(A, n)

    cycles = find_odd_cycles(A, n)
    if len(cycles) > 12: continue
    cg = build_cg(cycles)
    alpha = indep_poly_coeffs(cg, len(cycles))

    key = tuple(alpha[:max((k for k,a in enumerate(alpha) if a), default=0)+1])
    if key in seen: continue
    seen[key] = True

    vals = [eval_indep_poly(alpha, xv) for xv in [1, 2, 3, 6, 12]]
    alpha_str = str(list(key))
    print(f"{H:5d}  {alpha_str:>12}  {vals[0]:7d}  {vals[1]:7d}  {vals[2]:7d}  {vals[3]:7d}  {vals[4]:8d}")

# ============================================================================
# PART 6: THE k-NACCI / k-JACOBSTHAL CONVERGENCE RATES
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 6: CONVERGENCE RATES — HOW FAST DO k-NACCI AND k-JACOBSTHAL CONVERGE?║
╚══════════════════════════════════════════════════════════════════════════════╝

k-nacci root → 2:     the rate of convergence matters!
k-Jacobsthal root → 3: how does it compare?
""")

def knacci_root(k):
    coeffs = [1]*k
    poly = [1] + [-1]*k
    roots = np.roots(poly)
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 0]
    return max(real_pos)

def kjacob_root(k):
    coeffs = [2**i for i in range(k)]
    poly = [1] + [-2**i for i in range(k)]
    roots = np.roots(poly)
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 0]
    return max(real_pos)

print(f"{'k':>3}  {'r_N(k)':>12}  {'2-r_N':>12}  {'r_J(k)':>12}  {'3-r_J':>12}  {'(2-r_N)/(3-r_J)':>16}")
print("-" * 72)

for k in range(2, 20):
    rn = knacci_root(k)
    rj = kjacob_root(k)
    diff_n = 2 - rn
    diff_j = 3 - rj
    ratio = diff_n / diff_j if diff_j > 1e-15 else float('inf')
    print(f"{k:3d}  {rn:12.8f}  {diff_n:12.2e}  {rj:12.8f}  {diff_j:12.2e}  {ratio:16.6f}")

# ============================================================================
# PART 7: THE 3-ADIC AND 2-ADIC STRUCTURE OF JACOBSTHAL
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 7: THE MODULAR STRUCTURE OF J(n) — PERIODICITY AND DIVISIBILITY       ║
╚══════════════════════════════════════════════════════════════════════════════╝

J(n) = (2^n - (-1)^n) / 3

Modular properties:
  J(n) mod 2: period 3 → (0,1,1,0,1,1,...) since 2^n mod 2 = 0 for n≥1
  Wait: J(0)=0, J(1)=1, J(2)=1, J(3)=3, J(4)=5, J(5)=11, J(6)=21
  mod 2: 0, 1, 1, 1, 1, 1, 1, ... (all odd for n≥1)

  J(n) mod 3: period 8 (since 2^8 ≡ 1 mod 3... no, 2^2 = 4 ≡ 1 mod 3)
  Actually: 2^n mod 3: period 2 (2, 1, 2, 1, ...)
  And (-1)^n mod 3: period 2 (-1, 1, -1, 1, ...)
  So 3J(n) = 2^n - (-1)^n mod 9 has period lcm(6,2) = 6
""")

# Detailed modular analysis
print("Jacobsthal J(n) modular structure:")
print(f"{'n':>3}  {'J(n)':>8}  {'mod 2':>6}  {'mod 3':>6}  {'mod 5':>6}  {'mod 7':>6}  {'mod 8':>6}")
print("-" * 45)
for n in range(20):
    jn = (2**n - (-1)**n) // 3
    print(f"{n:3d}  {jn:8d}  {jn%2:6d}  {jn%3:6d}  {jn%5:6d}  {jn%7:6d}  {jn%8:6d}")

# Find periods
print("\nPeriods of J(n) mod p:")
for p in [2, 3, 4, 5, 7, 8, 9, 16]:
    vals = [(2**n - (-1)**n) // 3 % p for n in range(60)]
    # Find period
    for period in range(1, 31):
        if all(vals[i] == vals[i + period] for i in range(30)):
            print(f"  J(n) mod {p:2d}: period = {period}, pattern = {vals[:period]}")
            break

# ============================================================================
# PART 8: THE DEEP STRUCTURE — WHY 2 AND 3 ARE "THE" CONSTANTS
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 8: WHY 2 AND 3 — THE CATEGORICAL ARGUMENT                            ║
╚══════════════════════════════════════════════════════════════════════════════╝

CLAIM: 2 and 3 are uniquely privileged because:

1. THE BINARY PRINCIPLE: Every mathematical structure is built from
   BINARY CHOICES (yes/no, 0/1, left/right, include/exclude).
   The number 2 encodes this universal binariness.

2. THE INDEPENDENCE PRINCIPLE: When binary choices must be COMPATIBLE
   (independent), the weight per choice is 1+1=2, but the total
   weight of a new independent option is 1+2=3 (include with weight 2,
   or exclude with weight 1).

3. THE RECURRENCE PRINCIPLE: Linear recurrences with weight w per step
   have root (1+√(1+4w))/2. At w=2 (binary choices), root=2=w.
   This FIXED POINT PROPERTY (root = weight = 2) is unique.

4. THE NORMALIZATION PRINCIPLE: When dividing by the number of states
   per position (=3 = 1+2), we get J(n) = (2^n - ε)/3.
   The 3 in the denominator is the normalization constant.

5. THE DUALITY: 2 (binary) and 3 (ternary = binary + null) are the
   SMALLEST primes, and their ratio 3/2 is the simplest non-integer
   rational > 1. This ratio governs all scaling in tournament theory.

6. MUSICAL ANALOGY: The 3/2 ratio is the PERFECT FIFTH in music
   (frequency ratio of the most consonant interval after the octave).
   In tournament theory, 3/2 is the most "consonant" scaling factor.
""")

# The Stern-Brocot / continued fraction of 3/2
print("The ratio 3/2 in number theory:")
print(f"  3/2 = 1 + 1/2 (simplest fraction > 1 with smallest denominator)")
print(f"  3/2 = [1; 2] as continued fraction")
print(f"  log_2(3) = {math.log2(3):.10f} (Mersenne-related)")
print(f"  log_3(2) = {math.log(2)/math.log(3):.10f}")
print(f"  2^{math.log2(3):.4f} = 3")
print(f"  3^{math.log(2)/math.log(3):.4f} = 2")
print()

# Connection to Mersenne primes: 2^p - 1 is prime
print("The 2-3 in Mersenne primes:")
print(f"  M_p = 2^p - 1 (Mersenne number)")
print(f"  J(n) = (2^n - (-1)^n) / 3 — Jacobsthal as 'reduced Mersenne'")
print(f"  When n is odd: J(n) = (2^n + 1) / 3")
print(f"  When n is even: J(n) = (2^n - 1) / 3")
print()
print(f"  J(n) for odd n: n=3→3, n=5→11, n=7→43, n=11→683, n=13→2731")
print(f"  (2^n+1)/3:      3, 11, 43, 683, 2731")
print(f"  Compare Mersenne: 2^3-1=7, 2^5-1=31, 2^7-1=127, 2^11-1=2047, 2^13-1=8191")

# Check primality of (2^n+1)/3
print("\n  Primality of J(n) = (2^n - (-1)^n)/3:")
def is_prime(n):
    if n < 2: return False
    for p in range(2, int(n**0.5)+1):
        if n % p == 0: return False
    return True

for n in range(2, 30):
    jn = (2**n - (-1)**n) // 3
    prime = is_prime(jn)
    if prime:
        print(f"    J({n:2d}) = {jn:10d} — PRIME")

# ============================================================================
# PART 9: THE GROWTH RATE OF H — 2 AND 3 IN TOURNAMENT STATISTICS
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 9: GROWTH RATE OF H — WHERE 2 AND 3 APPEAR IN ASYMPTOTICS            ║
╚══════════════════════════════════════════════════════════════════════════════╝

Average H over all tournaments on n vertices:
  ⟨H⟩ = n! / 2^{n-1}  (each permutation is an HP with prob 1/2^{n-1})

Growth rate: ⟨H⟩ ~ √(2π·n) · (n/e)^n / 2^{n-1} ~ √(2πn) · (n/(2e))^n · 2

So ⟨H⟩ grows super-exponentially.

Maximum H: for the regular tournament (when n is odd)
  H_max ~ C · n! / 2^{O(n)} (conjectured to be about n!/2^n)

Minimum H: for the transitive tournament
  H_min = 1 (always)

RATIO ⟨H⟩/H_min is huge, but what about ⟨log H⟩?
""")

# Compute H statistics for small n
print("H statistics for all tournaments:")
print(f"{'n':>3}  {'#tourn':>8}  {'avg H':>10}  {'max H':>7}  {'min H':>7}  {'avg/min':>8}  {'H≡1(3)':>8}  {'H≡0(3)':>8}  {'H≡2(3)':>8}")
print("-" * 85)

for n in range(3, 8):
    m_ = n*(n-1)//2
    max_bits = 1 << m_
    if max_bits > 2**21:
        break

    h_values = []
    for bits_val in range(max_bits):
        A = adj_matrix(bits_val, n)
        H = count_hp(A, n)
        h_values.append(H)

    avg_h = sum(h_values) / len(h_values)
    max_h = max(h_values)
    min_h = min(h_values)
    mod3 = [h_values.count(h) for h in set(h_values)]
    h_mod3_0 = sum(1 for h in h_values if h % 3 == 0)
    h_mod3_1 = sum(1 for h in h_values if h % 3 == 1)
    h_mod3_2 = sum(1 for h in h_values if h % 3 == 2)

    print(f"{n:3d}  {max_bits:8d}  {avg_h:10.2f}  {max_h:7d}  {min_h:7d}  {avg_h/min_h:8.2f}  {h_mod3_1:8d}  {h_mod3_0:8d}  {h_mod3_2:8d}")

# ============================================================================
# PART 10: THE ULTIMATE SYNTHESIS
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 10: THE ULTIMATE SYNTHESIS — 2 AND 3 AS KEYS TO THE UNIVERSE          ║
╚══════════════════════════════════════════════════════════════════════════════╝

"If one understands 2 and 3, then we have the keys to the universe."

In tournament theory, this manifests as:

THE SMALL WORLD:
  2 = number of edge orientations
  3 = 1 + 2 = number of states per vertex (off, or selected with weight 2)
  3/2 = growth factor of Jacobsthal/Fibonacci ratio

THE RECURRENCE WORLD:
  k-nacci → 2: compositions of n into parts ≤ k
  k-Jacobsthal → 3: weighted compositions with geometric weights
  General: geometric ratio r → limit 1+r; our r=2 gives limit 3

THE ALGEBRAIC WORLD:
  J(n) = (2^n - (-1)^n) / 3: the 2 in the base, the 3 in the denominator
  det(I+2A) = Pf(S)²: the 2 in the matrix, the square in the Pfaffian
  H = Σ 2^k α_k: the 2 in the weight
  I(P_m, 2) = J(m+2): the 2 in the evaluation point

THE ANALYTIC WORLD:
  Root spectrum: x → (1+√(1+4x))/2
  Integer roots at x = k(k-1): the 2-3 is x=2·1, root=2
  Growth: H ~ n!/2^n with correction involving 3

THE NUMBER-THEORETIC WORLD:
  v_3(J(n)): period 6 in the 3-adic valuation
  J(n) mod 3: pattern (1,1,0,2,2,0) with period 6
  Mersenne connection: J(n) = (2^n ± 1)/3, related to 2^p-1 primes

THE HARMONIC WORLD:
  3/2 = perfect fifth interval = most consonant non-octave
  The "harmonic series" of tournament theory is built on 2 and 3

CONCLUSION: The theory of tournaments IS the theory of 2 and 3,
because tournaments are BINARY structures (2) evaluated at their
natural weight (giving 3 = 1+2 per independent position).

Every theorem, every recurrence, every identity in this project
is a different facet of this single diamond: the interaction of
2 (the binary principle) and 3 (the independence principle).
""")
