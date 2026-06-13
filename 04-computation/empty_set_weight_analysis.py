#!/usr/bin/env python3
"""
The empty independent set and the structural decomposition of a_k(T).

CORE INSIGHT: In the OCF, α₀ = 1 represents the empty independent set.
The full formula reads:

  H(T) = I(Ω(T), 2) = sum_{S independent in Ω} 2^{|S|}
       = 2^0 * 1 + 2^1 * α₁ + 2^2 * α₂ + ...
       = 1 + 2α₁ + 4α₂ + 8α₃ + ...

The empty set contributes 1 = H(transitive).

In the deformed Eulerian formula (THM-062):

  a_k(T) = A(n,k) + sum_I 2^{parts(I)} * c_k^{(f_I, d)} * I(T)

The term A(n,k) is the "empty-set weight" — it's what a_k would be if T had
no odd cycles (i.e., if T were the transitive tournament).

So the k-colored independence polynomial is:
  I_k(Ω, x) = A(n,k) * x^0 + sum_I c_k * x^{parts(I)} * I(T)

and a_k(T) = I_k(Ω(T), 2).

QUESTION: What is A(n,k) in terms of the independence polynomial?
  A(n,k) = I_k(Ω, 0) — evaluating the k-colored IP at x=0 kills all invariants.

This means: G_T(t, 0) = A_n(t) for ALL tournaments T.
The "empty set slice" of the generating function is universal!

THIS SCRIPT investigates:
(1) The x=0 slice: G_T(t, 0) = A_n(t) always
(2) The "weight decomposition": which independent sets contribute what to a_k
(3) The ratios I_k(Ω, 2)/I(Ω, 2) = a_k/H — how does the weight distribute?
(4) The "spectral" view: eigenvalues of the weight matrix

opus-2026-03-07-S33
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction
import random

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def inflated_eulerian(f, d, k):
    total = 0
    for j in range(max(0, k - (d - f)), min(f, k) + 1):
        sign = (-1) ** (d - f - k + j)
        total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
    return total

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
# Part 1: The empty-set weight A(n,k) vs actual a_k
# ====================================================================
print("THE EMPTY INDEPENDENT SET AND EULERIAN BASE WEIGHTS")
print("=" * 70)

print("""
Key identity: a_k(T) = A(n,k) + correction_k(T)

where A(n,k) is the Eulerian number (= a_k for the transitive tournament)
and correction_k depends on odd-cycle invariants.

The ratio a_k(T) / A(n,k) tells us how much T deviates from "transitive behavior."
""")

n = 7
d = n - 1
print(f"n={n}: Eulerian numbers A(7,k) = {[eulerian_number(n,k) for k in range(n)]}")
print(f"  sum = {sum(eulerian_number(n,k) for k in range(n))} = {factorial(n)}")

for seed in range(5):
    A = random_tournament(n, n * 1000 + seed)
    dist = forward_edge_dist_dp(A, n)
    H = dist.get(0, 0)
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t7 = count_directed_cycles(A, n, 7)
    bc = count_bc(A, n)

    print(f"\n  seed={seed}: H={H}, t3={t3}, t5={t5}, t7={t7}, bc={bc}")
    for k in range(n):
        Ank = eulerian_number(n, k)
        ak = dist.get(k, 0)
        delta = ak - Ank
        ratio = Fraction(ak, Ank) if Ank != 0 else "inf"
        if k <= 3:  # only print half by palindromy
            print(f"    k={k}: A(7,{k})={Ank}, a_k={ak}, delta={delta:+d}, ratio={ratio}")

# ====================================================================
# Part 2: Weight distribution across independent set levels
# ====================================================================
print(f"\n{'=' * 70}")
print("WEIGHT DECOMPOSITION BY INDEPENDENCE NUMBER")
print("=" * 70)

print("""
H(T) = 1 + 2α₁ + 4α₂ + ... = sum of 2^{|S|} over independent sets S.

a_k(T) = A(n,k) + 2·Σ(single-cycle terms) + 4·Σ(pair terms) + ...

So a_k decomposes as: a_k = [level-0 weight] + [level-1 weight] + [level-2 weight] + ...
where level-j weight = sum over independent sets of size j of their contribution to a_k.

The level-0 (empty set) contribution is ALWAYS A(n,k).
The level-j contribution involves c_k * I(T) terms with parts = j.
""")

n = 7
d = n - 1
invariants_7 = [('t3', 4, 1), ('t5', 2, 1), ('t7', 0, 1), ('bc', 2, 2)]

for seed in range(3):
    A = random_tournament(n, n * 1100 + seed)
    dist = forward_edge_dist_dp(A, n)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    H = dist.get(0, 0)

    print(f"\n  seed={seed}: H={H}, inv={inv_vals}")
    print(f"  {'k':>3} {'A(n,k)':>8} {'level-1':>10} {'level-2':>10} {'total':>8} {'actual':>8} {'match':>6}")

    for k in range(n):
        base = eulerian_number(n, k)
        level1 = 0
        level2 = 0
        for name, f, parts in invariants_7:
            coeff = 2**parts * inflated_eulerian(f, d, k) * inv_vals[name]
            if parts == 1:
                level1 += coeff
            elif parts == 2:
                level2 += coeff
        total = base + level1 + level2
        actual = dist.get(k, 0)
        match = "OK" if total == actual else "FAIL"
        print(f"  {k:>3} {base:>8} {level1:>10} {level2:>10} {total:>8} {actual:>8} {match:>6}")

# ====================================================================
# Part 3: Ratio a_k/H — the "k-th weight per Hamiltonian path"
# ====================================================================
print(f"\n{'=' * 70}")
print("RATIO a_k / H = 'weight per Hamiltonian path'")
print("=" * 70)

print("""
Since a_0 = a_{n-1} = H, the ratio a_k / H starts at 1, rises, then returns to 1.

For transitive: a_k/H = A(n,k)/1 = A(n,k) (very peaked at center).
For "most random": a_k/H should be flatter (closer to n!/H for all k).

Question: does a_k/H relate to I_k(Ω, 2)/I(Ω, 2)?
Yes! a_k = I_k(Ω, 2) and H = I(Ω, 2) = I_0(Ω, 2), so a_k/H = I_k(Ω,2)/I_0(Ω,2).
""")

n = 7
print(f"\nn={n}:")
print(f"  Transitive: A(7,k)/1 = {[eulerian_number(n,k) for k in range(n)]}")

for seed in range(5):
    A = random_tournament(n, n * 1200 + seed)
    dist = forward_edge_dist_dp(A, n)
    H = dist.get(0, 0)
    ratios = [Fraction(dist.get(k, 0), H) for k in range(n)]
    print(f"  seed={seed}: H={H}, a_k/H = {[str(r) for r in ratios]}")

# ====================================================================
# Part 4: The "universal slice" G_T(t, 0) = A_n(t)
# ====================================================================
print(f"\n{'=' * 70}")
print("UNIVERSAL SLICE: G_T(t, 0) = A_n(t)")
print("=" * 70)

print("""
G_T(t, x) = A_n(t) + sum_I x^{parts} I(T) A_{f+1}(t) (t-1)^{d-f}

At x=0: G_T(t, 0) = A_n(t). This is T-independent!

Meaning: if we "turn off" all cycle corrections (x→0), every tournament
looks like the transitive tournament. The x parameter controls how much
cycle structure influences the forward-edge distribution.

The OCF at x=0: I(Ω, 0) = 1 (only the empty set, weight 2^0 = 1).
The OCF at x=2: I(Ω, 2) = H(T) (all independent sets contribute).

So x interpolates between "transitive behavior" (x=0) and "full tournament behavior" (x=2).
""")

# Verify G_T(t, 0) = A_n(t) for specific t values
n = 7
d = n - 1
for seed in range(3):
    A = random_tournament(n, n * 1300 + seed)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }

    for t in [0, 2, 3, -1, Fraction(1,2)]:
        # G_T(t, 0) should equal A_n(t)
        G_at_0 = sum(eulerian_number(n, k) * t**k for k in range(n))
        # G_T(t, 2) should equal E_T(t) = sum a_k t^k
        dist = forward_edge_dist_dp(A, n)
        E_T_t = sum(dist.get(k, 0) * t**k for k in range(n))

        # Compute G_T(t, 2) via formula
        G_formula = sum(eulerian_number(n, k) * t**k for k in range(n))
        for name, f, parts in [('t3', 4, 1), ('t5', 2, 1), ('t7', 0, 1), ('bc', 2, 2)]:
            A_f1_t = sum(eulerian_number(f+1, j) * t**j for j in range(f+1))
            G_formula += 2**parts * inv_vals[name] * A_f1_t * (t-1)**(d-f)

        if seed == 0:
            An_t = sum(eulerian_number(n, k) * t**k for k in range(n))
            print(f"  t={t}: A_n(t)={An_t}, G_T(t,0)=A_n(t)={G_at_0}, G_T(t,2)={G_formula}, E_T(t)={E_T_t}, match={G_formula==E_T_t}")

# ====================================================================
# Part 5: The x-derivative of G_T at x=0
# ====================================================================
print(f"\n{'=' * 70}")
print("x-DERIVATIVE OF G_T(t,x) AT x=0")
print("=" * 70)

print("""
dG/dx at x=0 = sum_{I with parts=1} I(T) A_{f_I+1}(t) (t-1)^{d-f_I}

This is the "infinitesimal cycle correction" — the first-order effect of
turning on cycle structure. Only SINGLE cycles contribute at x=0.

At t=0 (forward-counting endpoint):
  dG/dx|_{t=0,x=0} = sum_{I, parts=1} I(T) * A_{f_I+1}(0) * (-1)^{d-f_I}
  = sum_{I, parts=1} I(T) * 1 * (-1)^{d-f_I}
  (since A_m(0) = 1 for all m)
""")

n = 7
d = n - 1
print(f"\nn={n}: The x-derivative at x=0 for t=0 and t=1:")
for seed in range(3):
    A = random_tournament(n, n * 1400 + seed)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }

    for t in [0, 1]:
        dG_dx_0 = 0
        for name, f, parts in [('t3', 4, 1), ('t5', 2, 1), ('t7', 0, 1)]:
            A_f1_t = sum(eulerian_number(f+1, j) * t**j for j in range(f+1))
            dG_dx_0 += inv_vals[name] * A_f1_t * (t-1)**(d-f)
        alpha_1 = inv_vals['t3'] + inv_vals['t5'] + inv_vals['t7']
        if seed == 0:
            print(f"  t={t}: dG/dx|_0 = {dG_dx_0}, alpha_1 = {alpha_1}")

    # At t=0: dG/dx = t3 * 1 * (-1)^2 + t5 * 1 * (-1)^4 + t7 * 1 * (-1)^6
    #        = t3 + t5 + t7 = alpha_1
    # At t=1: dG/dx = sum * 0^{d-f} = sum over f=d only = t7 * A_1(1) * 1^0 = 0 (only if d=f=6, but t7 has f=0!)
    # Actually (t-1)^{d-f} at t=1 = 0 for d-f > 0, so only f=d=6 survives. But no invariant has f=6.
    # So dG/dx|_{t=1,x=0} = 0. This confirms G_T(1, x) = n! for all x.

# ====================================================================
# Part 6: Critical x-value where H "emerges" from A(n,0)
# ====================================================================
print(f"\n{'=' * 70}")
print("x-INTERPOLATION: From A(n,0)=1 to H(T)")
print("=" * 70)

print("""
At x=0: G_T(0, 0) = A(n, 0) = 1 (for all T, the "empty set" baseline).
At x=2: G_T(0, 2) = a_0(T) = H(T).

So H(T) "grows" from 1 as x goes from 0 to 2:
  G_T(0, x) = 1 + alpha_1 * x + alpha_2 * x^2 + ... = I(Omega, x)

This IS the independence polynomial!
""")

n = 7
d = n - 1
for seed in range(5):
    A = random_tournament(n, n * 1500 + seed)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    alpha_1 = inv_vals['t3'] + inv_vals['t5'] + inv_vals['t7']
    alpha_2 = inv_vals['bc']
    H = 1 + 2*alpha_1 + 4*alpha_2

    # Compute G_T(0, x) for several x values
    G_at_x = {}
    for x in [0, Fraction(1,2), 1, Fraction(3,2), 2]:
        val = 1  # A(n,0) * x^0
        for name, f, parts in [('t3', 4, 1), ('t5', 2, 1), ('t7', 0, 1), ('bc', 2, 2)]:
            # At t=0: A_{f+1}(0) = 1, (0-1)^{d-f} = (-1)^{d-f}
            A_f1_0 = 1
            factor = (-1)**(d-f)
            val += x**parts * inv_vals[name] * A_f1_0 * factor
        G_at_x[x] = val

    if seed < 3:
        print(f"\n  seed={seed}: alpha=({alpha_1},{alpha_2}), H={H}")
        for x, v in G_at_x.items():
            I_omega_x = 1 + alpha_1 * x + alpha_2 * x**2
            print(f"    G(0, {x}) = {v}, I(Omega, {x}) = {I_omega_x}, match={v == I_omega_x}")

# Wait -- is G_T(0, x) = I(Omega, x)? Let me verify more carefully.
# G_T(0, x) = A(n,0) + sum_I x^{parts} I(T) A_{f+1}(0) (0-1)^{d-f}
#           = 1 + sum_I x^{parts} I(T) * 1 * (-1)^{d-f}
#
# I(Omega, x) = 1 + alpha_1 * x + alpha_2 * x^2
# At n=7: alpha_1 = t3 + t5 + t7, alpha_2 = bc
#
# From G_T(0, x):
#   x^1 coeff: t3 * (-1)^2 + t5 * (-1)^4 + t7 * (-1)^6 = t3 + t5 + t7 = alpha_1 ✓
#   x^2 coeff: bc * (-1)^4 = bc = alpha_2 ✓
#
# YES! G_T(0, x) = I(Omega, x)!

print(f"\n{'=' * 70}")
print("THEOREM: G_T(0, x) = I(Omega(T), x)")
print("=" * 70)
print("""
PROOF: G_T(0, x) = A(n,0) + sum_I x^{parts(I)} I(T) A_{f_I+1}(0) (0-1)^{n-1-f_I}

Since A_m(0) = 1 for all m, and (0-1)^{n-1-f} = (-1)^{n-1-f}:

G_T(0, x) = 1 + sum_I x^{parts(I)} * (-1)^{n-1-f_I} * I(T)

For invariant I at level f:
  - Single (2m+1)-cycle: f = n-1-2m, so n-1-f = 2m. (-1)^{2m} = 1.
  - Pair of cycles with total 2S vertices: f = n-1-2S, so n-1-f = 2S. (-1)^{2S} = 1.

In general, n-1-f = (total number of cycle vertices) - 1 + ... hmm wait.

Actually for a single (2m+1)-cycle: it has 2m+1 vertices, which "use up" f+1 vertices
in the inner Eulerian polynomial and the remaining n-1-f come from (x-y)^{n-1-f}.

Actually, d-f = n-1-f is always even because the cycle collection always removes
an even number of "virtual positions" from the permutation:
  - (2m+1)-cycle "occupies" d-f = 2m positions (leaving f = d-2m for inner poly)
  - Pair of (2a+1)- and (2b+1)-cycles: d-f = 2(a+b) = even

So (-1)^{d-f} = 1 ALWAYS.

Therefore: G_T(0, x) = 1 + sum_I x^{parts(I)} I(T) = I(Omega(T), x).

This is a CLEAN result: the t=0 slice of the trivariate GF IS the independence polynomial!
""")

# ====================================================================
# Part 7: Summary table of G_T at special points
# ====================================================================
print(f"{'=' * 70}")
print("SPECIAL EVALUATION TABLE FOR G_T(t, x)")
print("=" * 70)
print("""
  G_T(t, x) = A_n(t) + sum_I x^{parts} I(T) A_{f+1}(t) (t-1)^{d-f}

  (t, x)     | Value
  -----------|---------------------------
  (t, 0)     | A_n(t) [Eulerian poly, T-independent]
  (0, x)     | I(Omega(T), x) [independence poly!]
  (0, 2)     | H(T) [Hamiltonian paths]
  (0, 0)     | 1 [empty independent set]
  (1, x)     | n! [for all x, T-independent]
  (t, 2)     | E_T(t) = sum_k a_k t^k [tournament Eulerian poly]
  (-1, 2)    | E_T(-1) = deformed zigzag number (0 for even n)
  (1, 0)     | n! [= A_n(1)]

THE x AND t AXES:
  - t-axis (x=0): A_n(t) — universal Eulerian polynomial
  - x-axis (t=0): I(Omega, x) — tournament independence polynomial
  - These meet at (0,0) = 1 (the "double empty-set" point)
  - The point (1, *) = n! (universal)

So G_T is a BILINEAR DEFORMATION: the t-direction deforms the Eulerian polynomial,
the x-direction deforms the independence polynomial, and they CROSS at (0,0) = 1.
""")

print(f"{'=' * 70}")
print("DONE")
print("=" * 70)
