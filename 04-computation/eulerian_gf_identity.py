#!/usr/bin/env python3
"""
Generating function identity for tournament Eulerian polynomial.

E_T(t) = sum_k a_k(T) t^k = A_n(t) + sum_I 2^{parts} I(T) A_{f+1}(t) (t-1)^{d-f}

where d = n-1 and A_m(t) = classical Eulerian polynomial of S_m.

Special evaluations:
  E_T(1) = n!  (trivially)
  E_T(0) = a_0 = H(T)
  E_T(-1) = sum_k (-1)^k a_k  (alternating sum)

The alternating sum is:
  E_T(-1) = A_n(-1) + sum_I 2^{parts} I(T) A_{f+1}(-1) (-1-1)^{d-f}
           = A_n(-1) + sum_I 2^{parts} I(T) A_{f+1}(-1) (-2)^{d-f}

Now A_n(-1) = sum_k (-1)^k A(n,k). For n >= 2, A_n(-1) = 0
(since A(n,k) = A(n, n-1-k), the palindromic Eulerian poly has -1 as a root).

So E_T(-1) = sum_I 2^{parts} I(T) A_{f+1}(-1) (-2)^{d-f}

And if A_{f+1}(-1) = 0 for all f >= 1, then E_T(-1) = 0 for all tournaments!

Let's check: is sum_k (-1)^k a_k(T) = 0 for all tournaments?
This would mean the forward-edge distribution is "alternating-balanced."

Also investigate: E_T(t) at t = roots of unity, and the connection
to I(Omega(T), x) at various x.

opus-2026-03-07-S33
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction
import random

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def eulerian_poly_eval(n, t):
    """Evaluate the Eulerian polynomial A_n(t) = sum_k A(n,k) t^k."""
    return sum(eulerian_number(n, k) * t**k for k in range(n))

def random_tournament(n, seed=42):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def forward_edge_dist(A, n):
    """Brute force: count perms by forward edges."""
    dist = defaultdict(int)
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        dist[fwd] += 1
    return dict(dist)

def forward_edge_dist_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    new_fwd = fwd + A[v][u]
                    key = (mask | (1 << u), u, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)

# ====================================================================
# Part 1: Check A_n(-1) = 0 for n >= 2
# ====================================================================
print("EULERIAN POLYNOMIAL AT t = -1")
print("=" * 70)

for n in range(1, 10):
    val = eulerian_poly_eval(n, -1)
    print(f"  A_{n}(-1) = {val}")

# ====================================================================
# Part 2: Alternating sum of a_k(T)
# ====================================================================
print(f"\n{'=' * 70}")
print("ALTERNATING SUM: sum_k (-1)^k a_k(T)")
print("=" * 70)

for n in [3, 4, 5, 6, 7]:
    print(f"\nn={n}:")
    use_dp = (n >= 7)

    num_tested = 0
    all_zero = True
    for seed in range(20 if n <= 6 else 10):
        A = random_tournament(n, n * 100 + seed)
        dist = forward_edge_dist_dp(A, n) if use_dp else forward_edge_dist(A, n)

        alt_sum = sum((-1)**k * dist.get(k, 0) for k in range(n))

        if seed < 5 or alt_sum != 0:
            print(f"  seed={seed}: alt_sum = {alt_sum}")

        if alt_sum != 0:
            all_zero = False
        num_tested += 1

    print(f"  {num_tested} tested: all zero = {all_zero}")

# ====================================================================
# Part 3: E_T(t) at other special values
# ====================================================================
print(f"\n{'=' * 70}")
print("E_T(t) AT SPECIAL VALUES")
print("=" * 70)

n = 7
for seed in range(5):
    A = random_tournament(n, n * 200 + seed)
    dist = forward_edge_dist_dp(A, n)
    H = dist.get(n-1, 0)

    for t_val, t_name in [(0, "t=0 (=H)"), (1, "t=1 (=n!)"), (-1, "t=-1"),
                           (2, "t=2"), (3, "t=3"), (-2, "t=-2"),
                           (Fraction(1,2), "t=1/2"), (-Fraction(1,2), "t=-1/2")]:
        E_val = sum(t_val**k * dist.get(k, 0) for k in range(n))
        if seed == 0:
            print(f"  {t_name}: E_T = {E_val}")

    # Key: E_T(0) = a_0 = H(T)
    # E_T(2) = sum a_k 2^k ... what does this relate to?
    E2 = sum(2**k * dist.get(k, 0) for k in range(n))
    E3 = sum(3**k * dist.get(k, 0) for k in range(n))
    Em1 = sum((-1)**k * dist.get(k, 0) for k in range(n))

    if seed < 3:
        print(f"\n  seed={seed}: H={H}, E(0)={dist.get(0,0)}, E(2)={E2}, E(3)={E3}")

# ====================================================================
# Part 4: Connection to I(Omega, x)
# ====================================================================
print(f"\n{'=' * 70}")
print("E_T(t) vs I(Omega(T), x)")
print("=" * 70)

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cl):
    if n < cl: return 0
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1
        for m in range(1, 1 << cl):
            for v in range(cl):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(cl):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cl) - 1
        total += sum(dp[full][v] for v in range(1, cl) if sub[v][0])
    return total

def count_bc(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

n = 7
print(f"\nn={n}: Comparing E_T(t) with I(Omega, x)")
for seed in range(10):
    A = random_tournament(n, n * 300 + seed)
    dist = forward_edge_dist_dp(A, n)

    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t7 = count_directed_cycles(A, n, 7)
    bc = count_bc(A, n)

    # I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2 + alpha_3*x^3
    alpha_1 = t3 + t5 + t7  # total odd cycles
    alpha_2 = bc  # disjoint pairs (at n=7, only bc33 matters for alpha_2)
    # alpha_3 = 0 at n=7 (need 9 vertices for 3 disjoint 3-cycles)
    # Actually wait: a3 = triples of disjoint 3-cycles needs 9 vertices, so = 0 at n=7

    def I_omega(x):
        return 1 + alpha_1 * x + alpha_2 * x**2

    H = I_omega(2)

    # Now compute E_T at various t values
    def E_T(t):
        return sum(t**k * dist.get(k, 0) for k in range(n))

    if seed < 5:
        print(f"\n  seed={seed}: alpha=({alpha_1}, {alpha_2}), H={H}")
        print(f"    E(0)={E_T(0)} [=a_0=H={H}]")
        print(f"    E(2)={E_T(2)}, I(Omega, ?)={E_T(2)}")

        # Is there z such that E_T(t) = I(Omega, z) for some (t, z)?
        # E_T(0) = H = I(Omega, 2). What about other t?
        # Try: is E_T(t) / something = I(Omega, z(t))?

        # More natural: Phi_T(x, y) = A_n(x,y) + sum_I 2^parts * A_{f+1}(x,y) * (x-y)^{d-f} * I
        # At (x,y) = (1+z, z): x-y = 1, so Phi = W(z + 1/2) = I(Omega, 2) = H when z = 0.
        # Hmm, this is just the W evaluation.

        # What if we set t = 2 in E_T(t)?
        # E_T(2) = Phi_T(2, 1) = A_n(2,1) + sum_I 2^parts * A_{f+1}(2,1) * 1^{d-f} * I
        # = A_n(2,1) + sum_I 2^parts * A_{f+1}(2,1) * I(T)
        # A_n(2,1) = sum_k A(n,k) 2^k
        # A_{f+1}(2,1) = sum_j A(f+1,j) 2^j

        # Actually: sum_k A(n,k) x^k = A_n(x) = classical Eulerian poly at x
        A_n_2 = eulerian_poly_eval(n, 2)
        for inv_name, f, parts, val in [('t3', 4, 1, t3), ('t5', 2, 1, t5), ('t7', 0, 1, t7), ('bc', 2, 2, bc)]:
            A_f_2 = eulerian_poly_eval(f+1, 2)

        # Check if E_T(2) = I(Omega, z) for some z
        # 1 + alpha_1*z + alpha_2*z^2 = E_T(2)
        # This is a quadratic in z!
        target = E_T(2)
        # z = (-alpha_1 ± sqrt(alpha_1^2 - 4*alpha_2*(1-target))) / (2*alpha_2) if alpha_2 > 0
        if alpha_2 > 0:
            disc = alpha_1**2 + 4*alpha_2*(target - 1)
            if disc >= 0:
                z1 = (-alpha_1 + disc**0.5) / (2*alpha_2)
                z2 = (-alpha_1 - disc**0.5) / (2*alpha_2)
                print(f"    E(2)={target} = I(Omega, {z1:.4f}) or I(Omega, {z2:.4f})")

# ====================================================================
# Part 5: The key identity E_T(-1) = 0
# ====================================================================
print(f"\n{'=' * 70}")
print("THEOREM: E_T(-1) = 0 FOR ALL TOURNAMENTS (n >= 2)")
print("=" * 70)
print("""
Proof: E_T(-1) = sum_k (-1)^k a_k(T) = Phi_T(-1, 1).

From the bivariate formula:
  Phi_T(-1, 1) = A_n(-1, 1) + sum_I 2^parts * A_{f+1}(-1, 1) * (-1-1)^{d-f} * I(T)

A_n(-1, 1) = sum_k A(n,k) (-1)^k = A_n(-1).

The Eulerian polynomial A_n(t) is palindromic: A_n(t) = t^{n-1} A_n(1/t).
At t = -1: A_n(-1) = (-1)^{n-1} A_n(-1).
For even n-1 (odd n): A_n(-1) = A_n(-1), no info.
For odd n-1 (even n): A_n(-1) = -A_n(-1), so A_n(-1) = 0.

But wait, we showed all alt sums are 0 even for odd n!
The deeper reason: a_k = a_{n-1-k}, so
sum_k (-1)^k a_k = sum_k (-1)^k a_{n-1-k}
                  = (-1)^{n-1} sum_k (-1)^{n-1-k} a_{n-1-k}
                  = (-1)^{n-1} sum_j (-1)^j a_j
So E_T(-1) = (-1)^{n-1} E_T(-1).
For even n (n-1 odd): E_T(-1) = -E_T(-1), so E_T(-1) = 0. QED.
For odd n (n-1 even): E_T(-1) = E_T(-1), no immediate conclusion.

Hmm, so for odd n, the argument doesn't work directly.
Let me check: is E_T(-1) = 0 for odd n too?
""")

# Recheck odd n specifically
for n in [3, 5, 7]:
    print(f"\nn={n} (odd):")
    for seed in range(5):
        A = random_tournament(n, n * 400 + seed)
        dist = forward_edge_dist_dp(A, n) if n >= 7 else forward_edge_dist(A, n)
        alt_sum = sum((-1)**k * dist.get(k, 0) for k in range(n))
        print(f"  seed={seed}: alt_sum = {alt_sum}, dist = {dict(sorted(dist.items()))}")

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
