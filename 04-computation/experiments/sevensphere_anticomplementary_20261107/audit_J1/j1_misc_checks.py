#!/usr/bin/env python3
# audit J1 (note section 4, THM-4911 remarks): THM-1950 velocity law and the d^2 = 1/s 'coth' line (real rapidity only for s < 1), s(DT) = s(T), the n = 2 almost-regular case, and equicentric <=> Delta_v constant.
from fractions import Fraction as Fr
import numpy as np
import sympy as sp

def s_exact(A):
    n = len(A); M = sp.Matrix(A) - sp.Matrix(A).T
    return sum((sp.eye(n) + M.T * M).LUsolve(sp.ones(n, 1)))
def rand_T(n, rng):
    U = np.triu(rng.integers(0, 2, (n, n)), 1); return U + np.triu(1 - U, 1).T
def ordered_sum(A1, A2):
    a, b = len(A1), len(A2); A = np.zeros((a + b, a + b), dtype=np.int64); A[:a, :a] = A1; A[a:, a:] = A2; A[:a, a:] = 1; return A
def skew_double(A):
    n = len(A); M = A - A.T; I = np.eye(n, dtype=np.int64); Mp = np.block([[M, M + I], [M - I, -M]]); return (Mp == 1).astype(np.int64)
rng = np.random.default_rng(5)
ok_vel = ok_d2 = 0; s_above_1 = 0; trials = 40
for _ in range(trials):
    A1, A2 = rand_T(int(rng.integers(1, 6)), rng), rand_T(int(rng.integers(1, 6)), rng)
    s1, s2, s12 = s_exact(A1), s_exact(A2), s_exact(ordered_sum(A1, A2))
    ok_vel += sp.simplify(s12 - (s1 + s2) / (1 + s1 * s2)) == 0
    d1, d2, d12 = 1 / s1, 1 / s2, 1 / s12
    ok_d2 += sp.simplify(d12 - (d1 * d2 + 1) / (d1 + d2)) == 0
    s_above_1 += int(bool(s1 > 1) or bool(s2 > 1))
print(f"THM-1950 velocity law s(A=>B) = (sA+sB)/(1+sA sB): {ok_vel}/{trials}; d^2 law d^2(A=>B) = (dA^2 dB^2 + 1)/(dA^2 + dB^2): {ok_d2}/{trials}")
print(f"  pairs with some s > 1 (d^2 < 1, so d^2 = coth(theta) needs complex theta): {s_above_1}/{trials}")
th1, th2 = sp.symbols('t1 t2')
print("  coth addition identity coth(a+b) = (coth a coth b + 1)/(coth a + coth b):",
      sp.simplify(sp.expand_trig(sp.coth(th1 + th2)).rewrite(sp.exp) - ((sp.coth(th1) * sp.coth(th2) + 1) / (sp.coth(th1) + sp.coth(th2))).rewrite(sp.exp)) == 0)
ok_D = 0
for _ in range(10):
    A = rand_T(int(rng.integers(2, 6)), rng); ok_D += sp.simplify(s_exact(skew_double(A)) - s_exact(A)) == 0
print(f"s(DT) = s(T) (THM-4900(3), skew-Sylvester double of THM-447): {ok_D}/10")
A2 = np.array([[0, 1], [0, 0]]); M2 = A2 - A2.T; y = M2 @ (M2 @ np.ones(2, dtype=int))
print("n = 2 single arc: almost regular (scores 1, 0 = n/2, n/2-1), equicentric:", bool((y == y[0]).all()), ", lambda =", int(-y[0]), "(so (3c)'s 'exist exactly when n = 2 mod 4, n >= 6' omits n = 2)")
# equicentric <=> Delta_v constant ((Mx)_v = 2 Delta_v - (n-1)^2), on random tournaments
bad = 0
for _ in range(2000):
    n = int(rng.integers(3, 9)); A = rand_T(n, rng); M = A - A.T; x = M.sum(1); y = M @ x
    Delta = np.array([sum(A[k, v] * A[:, k].sum() for k in range(n)) + sum(A[v, k] * A[k, :].sum() for k in range(n)) for v in range(n)])
    bad += not np.array_equal(y, 2 * Delta - (n - 1) ** 2)
print("identity (M^2 u)_v = 2 Delta_v - (n-1)^2 on 2000 random tournaments: violations", bad, "(so equicentric <=> Delta_v constant, lambda = (n-1)^2 - 2 Delta)")
# Sonneborn-Berger form: SB(v) = sum_{v->k} d+_k; then (M^2 u)_v = 4 SB(v) - 2(n-2) d+_v - (n-1), so equicentric <=> 2 SB(v) - (n-2) d+_v is constant
bad = 0
for _ in range(2000):
    n = int(rng.integers(3, 10)); A = rand_T(n, rng); M = A - A.T; d = A.sum(1); y = M @ (M @ np.ones(n, dtype=np.int64))
    SB = A @ d
    bad += not np.array_equal(y, 4 * SB - 2 * (n - 2) * d - (n - 1))
print("identity (M^2 u)_v = 4 SB(v) - 2(n-2) d+_v - (n-1) (SB = Sonneborn-Berger / second out-score) on 2000 random tournaments: violations", bad)
