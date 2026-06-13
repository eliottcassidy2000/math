#!/usr/bin/env python3
"""
Check: does W(r) have only real roots?

If W(r) = w_0 + w_2*r^2 + w_4*r^4 + w_6*r^6 has only real roots,
then the hierarchy has a clean explanation: the roots are "frequencies"
and the coefficients encode their positions.

For odd n, W(r) = r * Q(r^2) where Q is a polynomial in r^2? No...
W(r) has only even powers, so W(r) = P(r^2) where P(t) = w_0 + w_2*t + w_4*t^2 + w_6*t^3.

P(t) having only real NEGATIVE roots => I(x) real-rooted. This would connect
to the independence polynomial real-rootedness!

kind-pasteur-2026-03-06-S25g
"""
import numpy as np
from itertools import permutations, combinations
import random

def count_t3(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]: count += 1
                if A[i][k]*A[k][j]*A[j][i]: count += 1
    return count

n = 7
random.seed(42)

print("W(r) ROOT ANALYSIS at n=7")
print("=" * 60)
print("W(r) = w0 + w2*r^2 + w4*r^4 + w6*r^6")
print("Let P(t) = w0 + w2*t + w4*t^2 + w6*t^3 where t = r^2")
print("W(r) = P(r^2)")
print()

all_real = 0
has_complex = 0

for trial in range(30):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    # Compute W at several r values
    r_vals = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
    W_vals = []
    for r in r_vals:
        W = 0.0
        for p in permutations(range(n)):
            prod = 1.0
            for i in range(n-1):
                prod *= (r + A[p[i]][p[i+1]] - 0.5)
            W += prod
        W_vals.append(W)

    V = np.array([[r**(2*j) for j in range(4)] for r in r_vals])
    coeffs = np.linalg.lstsq(V, np.array(W_vals), rcond=None)[0]
    w0, w2, w4, w6 = coeffs

    # P(t) = w6*t^3 + w4*t^2 + w2*t + w0
    P_coeffs = [w0, w2, w4, w6]  # ascending order
    roots = np.roots([w6, w4, w2, w0])  # descending order for np.roots

    real_roots = [r for r in roots if abs(r.imag) < 1e-6]
    neg_real = [r.real for r in real_roots if r.real < 0]

    t3 = count_t3(A, n)

    if len(real_roots) == 3:
        all_real += 1
        if len(neg_real) == 3:
            tag = "ALL REAL NEGATIVE"
        else:
            tag = f"ALL REAL, {len(neg_real)} negative"
    else:
        has_complex += 1
        tag = "HAS COMPLEX ROOTS"

    if trial < 10:
        root_str = ", ".join(f"{r.real:.4f}{'+' if r.imag>=0 else ''}{r.imag:.4f}i" if abs(r.imag) > 1e-6 else f"{r.real:.4f}" for r in sorted(roots, key=lambda x: x.real))
        print(f"  trial={trial:2d}: t3={t3:2d}, P roots=[{root_str}] {tag}")

print(f"\nSummary: {all_real} all-real, {has_complex} has-complex out of 30")

# Now check: P(t) roots should be NEGATIVE (so that W(r) roots are imaginary)
# This would mean W(r) > 0 for all real r > 0 (since W(1/2) = H > 0)
print(f"\nIf all roots of P(t) are real negative, then W(r) = P(r^2) > 0 for all real r.")
print(f"This would be the 'path' analogue of I(Omega,x) having all real negative roots!")

# Also check at n=5
print(f"\n{'='*60}")
print(f"W(r) ROOT ANALYSIS at n=5")
print(f"{'='*60}")
print(f"P(t) = w0 + w2*t + w4*t^2 = 120*t^2 + (12*t3-30)*t + (-t3+2*t5+1)")

for bits in range(1024):
    A = [[0]*5 for _ in range(5)]
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    t3 = count_t3(A, 5)
    t5_count = 0
    for p in permutations(range(5)):
        if all(A[p[i]][p[(i+1)%5]] for i in range(5)):
            t5_count += 1
    t5 = t5_count // 5

    w0 = -t3 + 2*t5 + 1
    w2 = 12*t3 - 30
    w4 = 120

    roots = np.roots([w4, w2, w0])
    neg_real = all(r.real < 0 and abs(r.imag) < 1e-6 for r in roots)
    if not neg_real:
        print(f"  FAIL at bits={bits}: t3={t3}, t5={t5}, roots={roots}")
        break
else:
    print(f"  ALL 1024 tournaments: P(t) has all real negative roots!")

print("\nDONE")
