#!/usr/bin/env python3
"""
Deep investigation of u_T(m) root structure.

DISCOVERY (S39): Q_T(w) = u_T(sqrt(w))/sqrt(w) has ALL REAL NON-POSITIVE roots
for every tournament tested (exhaustive n=3,5; sampled n=7).

This means u_T(m) is a "real-rooted odd polynomial" — its roots are
{0, ±i*r_1, ±i*r_2, ...} where r_k > 0.

Key question: WHY does this hold? What structural property guarantees it?

Connections to investigate:
1. I(Omega(T), x) sometimes has complex roots (n=9, THM-025).
   But u_T(m) encodes MORE than just I — it includes the full
   principal specialization ps_1(U_T)(m). Could the extra structure
   from U_T's symmetric function nature force real roots?

2. u_T(m) = sum_k S_k * m^k where S_k are the f-level weighted sums.
   The S_k satisfy constraints from the Eulerian polynomial structure.

3. If Q_T(w) has all real negative roots, then by the theory of
   totally positive matrices / Polya frequency sequences, the
   coefficients of Q_T are log-concave. This would give coefficient
   inequalities on the S_k.

opus-2026-03-07-S39
"""
from itertools import permutations
from collections import defaultdict
import numpy as np
import random


def compute_uT_polynomial(n, A):
    """Compute u_T(m) coefficients."""
    edge_set = set()
    opp_set = set()
    for i in range(n):
        for j in range(n):
            if i != j:
                if A[i][j]:
                    edge_set.add((i, j))
                else:
                    opp_set.add((i, j))

    coeffs = defaultdict(int)
    for sigma in permutations(range(n)):
        visited = [False] * n
        cycles = []
        for s in range(n):
            if visited[s]:
                continue
            cyc = []
            c = s
            while not visited[c]:
                visited[c] = True
                cyc.append(c)
                c = sigma[c]
            cycles.append(tuple(cyc))

        valid = True
        phi = 0
        has_even = False
        for cyc in cycles:
            if len(cyc) == 1:
                continue
            if len(cyc) % 2 == 0:
                has_even = True
                break
            is_T = all((cyc[i], cyc[(i+1)%len(cyc)]) in edge_set for i in range(len(cyc)))
            is_Top = all((cyc[i], cyc[(i+1)%len(cyc)]) in opp_set for i in range(len(cyc)))
            if not is_T and not is_Top:
                valid = False
                break
            if is_T:
                phi += len(cyc) - 1

        if not valid or has_even:
            continue

        sign = (-1) ** phi
        num_cycles = len(cycles)
        coeffs[num_cycles] += sign

    return coeffs


def Q_roots(coeffs, n):
    """Extract Q_T roots from u_T coefficients."""
    q_coeffs = []
    for k in range(1, n + 1, 2):
        q_coeffs.append(coeffs.get(k, 0))
    q_array = np.array(q_coeffs, dtype=float)
    if len(q_array) > 1 and q_array[-1] != 0:
        return np.roots(q_array[::-1]), q_coeffs
    return np.array([]), q_coeffs


def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


# === Exhaustive n=5: Classify root patterns ===
print("=== n=5 exhaustive root analysis ===")
n = 5
root_patterns = defaultdict(int)
min_root = float('inf')
max_root = float('-inf')

for bits in range(1 << 10):
    A = [[0]*n for _ in range(n)]
    edge_list = [(i,j) for i in range(n) for j in range(i+1,n)]
    for k, (i,j) in enumerate(edge_list):
        if bits & (1 << k):
            A[j][i] = 1
        else:
            A[i][j] = 1

    coeffs = compute_uT_polynomial(n, A)
    roots, q = Q_roots(coeffs, n)

    real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-10])
    pattern = tuple(round(r, 4) for r in real_roots)
    root_patterns[pattern] += 1

    for r in real_roots:
        if r < min_root and r != 0:
            min_root = r
        if r > max_root and r != 0:
            max_root = r

print(f"Distinct root patterns: {len(root_patterns)}")
for pattern, count in sorted(root_patterns.items(), key=lambda x: -x[1])[:15]:
    print(f"  {pattern}: {count} tournaments")
print(f"Root range (nonzero): [{min_root:.4f}, {max_root:.4f}]")

# === Connection to independence polynomial ===
print("\n=== Connection: u_T(m) vs I(Omega,x) ===")
# u_T(1) = H(T) = I(Omega, 2)
# u_T(-1) = -H(T) (since odd polynomial)
# What about u_T(i) where i=sqrt(-1)? Then m^2 = -1, so Q_T(-1) = u_T(i)/i
# Q_T(-1) = sum_k S_{2k+1} * (-1)^k
# This is an alternating sum of the S coefficients!

print("For Paley T_7:")
A = [[0]*7 for _ in range(7)]
QR = {1, 2, 4}
for i in range(7):
    for j in range(7):
        if i != j and (j - i) % 7 in QR:
            A[i][j] = 1

coeffs = compute_uT_polynomial(7, A)
print(f"  u_T coefficients: {dict(sorted(coeffs.items()))}")
print(f"  u_T(1) = {sum(coeffs.values())} = H(T)")
print(f"  u_T(-1) = {sum((-1)**k * v for k, v in coeffs.items())}")

# Q_T(w) coefficients
q = []
for k in range(1, 8, 2):
    q.append(coeffs.get(k, 0))
print(f"  Q_T coefficients: {q}")
print(f"  Q_T(1) = {sum(q)} (= u_T(1)/1 = H)")
print(f"  Q_T(-1) = {sum((-1)**i * q[i] for i in range(len(q)))}")

# Relation to independence polynomial
# I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2 + ...
# H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...
# u_T(m) = m^n + 2*alpha_1*m^{n-2} + 4*alpha_2*m^{n-4} + ...
# = m * (m^{n-1} + 2*alpha_1*m^{n-3} + 4*alpha_2*m^{n-5} + ...)
# = m * Q_T(m^2) where Q_T(w) = w^{(n-1)/2} + 2*alpha_1*w^{(n-3)/2} + 4*alpha_2*w^{(n-5)/2} + ...

# So Q_T(w) = sum_k 2^k * alpha_k * w^{(n-1)/2 - k}
# The roots of Q_T are related to the roots of I(Omega, x) by:
# If Q_T(w) = 0, then we need the relationship...

# Actually: Q_T(w) = w^{(n-1)/2} * I(Omega, 2/w) * (something)... let me check

# Q_T(w) = sum_{k=0}^{(n-1)/2} 2^k * alpha_k * w^{(n-1)/2-k}
# = w^{(n-1)/2} * sum_{k=0} alpha_k * (2/w)^k
# = w^{(n-1)/2} * I(Omega, 2/w)

print("\n=== KEY IDENTITY: Q_T(w) = w^m * I(Omega, 2/w) where m=(n-1)/2 ===")
print("If I(Omega, x) = sum alpha_k x^k, then")
print("Q_T(w) = sum_{k=0}^m 2^k * alpha_k * w^{m-k} = w^m * sum alpha_k * (2/w)^k = w^m * I(Omega, 2/w)")
print()
print("This means: Q_T(w) = 0 iff I(Omega, 2/w) = 0 (for w != 0)")
print("i.e., if x_0 is a root of I(Omega, x), then w_0 = 2/x_0 is a root of Q_T(w)")
print()

# Verify for Paley T_7
# I(Omega_T7, x) = 1 + alpha_1*x + alpha_2*x^2 + alpha_3*x^3
# H = I(Omega, 2) = 189
# From earlier work: Paley T_7 has 28 3-cycles, 7*12=84 oriented 5-cycles...
# wait, let me compute alpha_k directly from the Q coefficients

# Q_T = [48, 112, 28, 1] means Q_T(w) = 48 + 112*w + 28*w^2 + 1*w^3
# This should equal w^3 * I(Omega, 2/w) = w^3 * (1 + alpha_1*(2/w) + alpha_2*(4/w^2) + alpha_3*(8/w^3))
# = w^3 + 2*alpha_1*w^2 + 4*alpha_2*w + 8*alpha_3

# So: coefficient of w^3: 1 (check)
#     coefficient of w^2: 2*alpha_1 = 28 => alpha_1 = 14
#     coefficient of w^1: 4*alpha_2 = 112 => alpha_2 = 28
#     coefficient of w^0: 8*alpha_3 = 48 => alpha_3 = 6

alpha_1 = 28 // 2  # 14
alpha_2 = 112 // 4  # 28
alpha_3 = 48 // 8  # 6
print(f"Paley T_7: alpha_1={alpha_1}, alpha_2={alpha_2}, alpha_3={alpha_3}")
print(f"  H = 1 + 2*{alpha_1} + 4*{alpha_2} + 8*{alpha_3} = {1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3}")

# Now the roots of Q_T(w) are the values 2/x_0 where x_0 are roots of I(Omega, x)
# If all roots of I(Omega, x) are real negative (which holds for n<=8 by THM-020),
# then w_0 = 2/x_0 are real negative, so Q_T(w) has all real negative roots!

print("\n=== THEOREM CANDIDATE ===")
print("Q_T(w) = w^m * I(Omega(T), 2/w)  where m = (n-1)/2")
print()
print("COROLLARY: Q_T(w) has all real non-positive roots")
print("  <=> I(Omega(T), x) has all real non-positive roots (for x != 0)")
print()
print("At n <= 8: I(Omega, x) has all real negative roots (Chudnovsky-Seymour,")
print("since Omega is claw-free). So Q_T(w) has all real negative roots.")
print()
print("At n = 9: I(Omega, x) can have complex roots (THM-025).")
print("So Q_T(w) should also have complex roots at n=9!")
print("Let's test this...")

# Build the THM-025 counterexample at n=9
# Score sequence [1,1,3,4,4,4,6,6,7]
# I(Omega, x) = 1 + 94x + 10x^2 + x^3
# Roots of I: -4.995 ± 8.303i, and one real root
# Q_T(w) = w^4 * I(Omega, 2/w) = w^4 + 2*94*w^3 + 4*10*w^2 + 8*w
#         = w*(w^3 + 188*w^2 + 40*w + 8)

print("\n=== THM-025 counterexample (predicted) ===")
I_coeffs = [1, 94, 10, 1]
# Q_T(w) = w^4 + 2*94*w^3 + 4*10*w^2 + 8*1*w = w*(w^3 + 188w^2 + 40w + 8)
q_poly = [0, 8, 40, 188, 1]  # w^0 through w^4
print(f"Q_T(w) = w^4 + 188*w^3 + 40*w^2 + 8*w = w*(w^3 + 188w^2 + 40w + 8)")
inner_roots = np.roots([1, 188, 40, 8])
print(f"Roots of w^3 + 188w^2 + 40w + 8: {inner_roots}")
all_real_inner = all(abs(r.imag) < 1e-8 for r in inner_roots)
print(f"All real: {all_real_inner}")

# Wait - let me reconsider. At n=9, m=(9-1)/2=4
# u_T(m) would have terms m, m^3, m^5, m^7, m^9 (up to m^9)
# Q_T(w) = c_1 + c_3*w + c_5*w^2 + c_7*w^3 + c_9*w^4
# = sum_{k=0}^4 2^k * alpha_k * w^{4-k}
# So Q_T(w) = alpha_0 * w^4 + 2*alpha_1 * w^3 + 4*alpha_2 * w^2 + 8*alpha_3 * w + 16*alpha_4

# For the counterexample: alpha_0=1, alpha_1=94, alpha_2=10, alpha_3=1, alpha_4=0
# Q_T(w) = w^4 + 188*w^3 + 40*w^2 + 8*w + 0
# = w*(w^3 + 188*w^2 + 40*w + 8)

print(f"\nWith Q_T(w) = w^4*I(Omega, 2/w):")
print(f"  Q_T(w) = w^4 + 188w^3 + 40w^2 + 8w")
# Roots: w=0 and roots of w^3 + 188w^2 + 40w + 8

# I(Omega, x) = 1 + 94x + 10x^2 + x^3 has roots at x = -2/w for each nonzero w
# If I has complex roots, then Q_T should also have complex roots
# Let's verify: roots of I(x) = x^3 + 10x^2 + 94x + 1
I_roots = np.roots([1, 10, 94, 1])
print(f"Roots of I(Omega, x): {I_roots}")
print(f"Predicted Q_T roots (2/x_i): {[2/r for r in I_roots]}")

# The inner cubic w^3 + 188w^2 + 40w + 8 should have same reality pattern
# But wait: 2/(-4.995 + 8.303i) is complex, so yes, Q_T has complex roots
print(f"\nAt n=9, when I(Omega) has complex roots, Q_T also has complex roots.")
print(f"At n<=8, I(Omega) always has real roots (claw-free), so Q_T always has real roots.")
print(f"\nCONCLUSION: The real-rootedness of Q_T is EQUIVALENT to real-rootedness of I(Omega, x).")
print(f"This is NOT a new phenomenon — it's the SAME property viewed through a different lens!")
