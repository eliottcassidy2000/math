"""Is COMPLETE MONOTONICITY of the density the operative hypothesis?

Our measure: dnu/du = (1/s) u^{1/s-1} e^{-u^{1/s}}.  Both factors are completely
monotone (u^{1/s-1} with exponent <= 0; e^{-u^{1/s}} is CM since u^{1/s} is a
Bernstein function for 0 < 1/s <= 1), and CM is closed under products.  So
dnu/du IS completely monotone -- a property the atomic counterexample of H6
manifestly lacks.

By Bernstein, CM density  <=>  dnu(u) = (int_0^inf e^{-a u} dsigma(a)) du,
sigma >= 0.  Then, with Psi_n(y) = int_0^inf (y-w)^n e^{-w} dw the s=1 Appell
polynomial (squarefree, a truncated exponential),

    Phi_n(z) = int_0^inf a^{-n-1} Psi_n(a z) dsigma(a)

-- a positive mixture of DILATES of one fixed squarefree polynomial.

TEST: take sigma atomic (a finite mixture of exponentials, still a CM density)
and search for a multiple root of Phi_n.  If none exist, CM is the right
hypothesis and the conjecture is clean.
"""
import numpy as np
from scipy.optimize import fsolve
from math import comb, factorial

def Psi(n, y):
    return sum(comb(n, j) * (-1) ** j * factorial(j) * y ** (n - j) for j in range(n + 1))

def Phi_mix(n, q, a, z):
    return sum(q[k] * a[k] ** (-n - 1) * Psi(n, a[k] * z) for k in range(len(q)))

print("Sanity: Phi_n for a SINGLE exponential is Psi_n (squarefree, truncated exp).")
for n in (3, 4, 5, 6, 7, 8):
    c = [comb(n, j) * (-1) ** j * factorial(j) for j in range(n + 1)]
    r = np.roots(c)
    gaps = min(abs(r[i] - r[j]) for i in range(len(r)) for j in range(i + 1, len(r)))
    print(f"   n={n}: min root separation of Psi_n = {gaps:.6f}")

print("\nSearch: 2 exponentials, q=(1,q2)>0, a=(1,a2)>0.  4 real eqs, 4 real unknowns.")
rng = np.random.default_rng(11)
for n in (3, 4, 5, 6, 7):
    def F(v):
        q2, a2, x, y = v
        z = complex(x, y); q = [1.0, q2]; a = [1.0, a2]
        A = Phi_mix(n, q, a, z); B = Phi_mix(n - 1, q, a, z)
        return [A.real, A.imag, B.real, B.imag]
    hits = 0
    for _ in range(6000):
        g = [rng.uniform(.01, 20), rng.uniform(.05, 20), rng.uniform(-25, 25), rng.uniform(.01, 20)]
        try:
            s, info, ier, _ = fsolve(F, g, full_output=True, xtol=1e-13)
        except Exception:
            continue
        if ier == 1 and s[0] > 1e-6 and s[1] > 1e-6 and abs(s[3]) > 1e-6 \
           and max(abs(np.array(F(s)))) < 1e-8 * max(1, abs(s[2]) ** n):
            hits += 1
    print(f"   n={n}: positive-weight multiple-root solutions found = {hits}")

print("\nSearch: 3 exponentials, q=(1,q2,q3)>0, a=(1,a2,a3)>0 (underdetermined -> easier).")
for n in (4, 5, 6):
    def F3(v):
        q2, q3, a2, a3, x, y = v
        z = complex(x, y); q = [1.0, q2, q3]; a = [1.0, a2, a3]
        A = Phi_mix(n, q, a, z); B = Phi_mix(n - 1, q, a, z)
        return [A.real, A.imag, B.real, B.imag, 0.0, 0.0]
    hits = 0
    for _ in range(6000):
        g = [rng.uniform(.01, 15), rng.uniform(.01, 15), rng.uniform(.05, 15),
             rng.uniform(.05, 15), rng.uniform(-20, 20), rng.uniform(.01, 15)]
        try:
            s, info, ier, _ = fsolve(F3, g, full_output=True, xtol=1e-13)
        except Exception:
            continue
        if ier == 1 and all(s[i] > 1e-6 for i in range(4)) and abs(s[5]) > 1e-6 \
           and max(abs(np.array(F3(s)))) < 1e-8 * max(1, abs(s[4]) ** n):
            hits += 1
    print(f"   n={n}: positive-weight multiple-root solutions found = {hits}")
