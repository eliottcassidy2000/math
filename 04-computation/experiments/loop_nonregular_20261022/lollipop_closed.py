"""Closed form of the lollipop class X_(0,1) (a cycle with one pendant leaf), from the leaf engine's decomposition
weights.  Session opus-2026-10-22-S54.

The lollipop core has a leaf 0, a core vertex 1, the pendant edge (0,1) and the loop (1,1); its value is
  X_(0,1) = (1/|G|) sum_(d1,d2) W[d1,d2] sum_w (1^T P_d1)_w (P_d2)_ww,
with every propagator P_d a rational function f_d of Shat, so (1^T P_d1)_w = (f_d1(-Shat) 1)_w (transpose = Shat -> -Shat).
The kernel T(s,t) = (1/|G|) sum W f_d1(-s) f_d2(t) is a rational function of two commuting symbols; partial fractions
in t over the fourth roots of unity split it into resolvents. Then X_(0,1) = sum_w (a(Shat) delta)_w (B(Shat))_ww.
The script prints T(s,t) in factored form, checks the resulting matrix formula against the engine on random
almost-regular tournaments, and prints the formula.  Usage: python lollipop_closed.py"""
import sys, os
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leaf_engine as LE
from cluster_brute_classes import almost_regular

s, t = sp.symbols('s t')


def fsym(desc, x):
    """the propagator of descriptor desc as a rational function of a commuting symbol x (x stands for Shat)."""
    R = lambda sg: sg * x / (1 - sg * x)
    if desc[0] == 'm1':
        return R(1) + R(-1) if desc[1] == '*' else R(desc[1])
    A = lambda sg: R(1) + R(-1) if sg == '*' else R(sg)
    return sp.together(-A(desc[1]) * A(desc[2]) * (1 - x ** 2) / (1 + x ** 2))


core = LE.get(0, 1)[0]
assert core.V == 2 and core.edges == [(0, 1), (1, 1)]
DESCS = LE.LEF.DESCS
T = 0
it = np.nditer(core.W, flags=['multi_index'])
for w in it:
    if abs(float(w)) == 0:
        continue
    d1 = DESCS[core.used[0][it.multi_index[0]]]; d2 = DESCS[core.used[1][it.multi_index[1]]]
    T += sp.Rational(int(round(float(w)))) * fsym(d1, -s) * fsym(d2, t)
T = sp.factor(sp.together(T / core.G))
print('kernel T(s, t) =', T)
# the vector side: f(-s) 1 has a factor s, so write T = s * a(s) * b(t)
# separate variables if possible
a_s = sp.factor(T.subs(t, sp.Rational(1, 3)) / T.subs({t: sp.Rational(1, 3), s: sp.Rational(1, 5)}))
b_t = sp.factor(T.subs(s, sp.Rational(1, 5)))
sep = sp.simplify(T - a_s * b_t)
print('separable:', sep == 0)
print('a(s) =', a_s, '   b(t) =', b_t)
fa = sp.lambdify(s, a_s, 'numpy'); fb = sp.lambdify(t, b_t, 'numpy')


def mat_func(f_expr, var, Sh):
    """evaluate a rational function of var at the matrix Sh via eigen-decomposition (Sh is normal)."""
    ev, V = np.linalg.eig(Sh)
    fe = sp.lambdify(var, f_expr, 'numpy')
    return (V @ np.diag(fe(ev.astype(complex))) @ np.linalg.inv(V))


rng = np.random.default_rng(3)
for n in (8, 10, 14):
    A = almost_regular(n, rng)
    Sh = (A - A.T).astype(float) / n
    one = np.ones(n)
    # X = sum_w (a(Shat) applied with the s factor: f(-s)1 contains Shat 1 = delta) ... evaluate T bilinearly:
    Am = mat_func(a_s, s, Sh); Bm = mat_func(b_t, t, Sh)
    val = np.sum((Am @ one) * np.diag(Bm)).real
    eng = LE.X(Sh, 0, 1).real
    I = np.eye(n); S2 = Sh @ Sh
    G = -S2 @ np.linalg.inv(I + S2); H = S2 @ np.linalg.inv(I - S2); d = np.diag(H - G)
    mat = -2 * d @ G @ one
    print(f'n={n}: kernel form {val:+.15e}  -2 d^T G 1 = {mat:+.15e}  engine {eng:+.15e}  diffs {val - eng:+.1e} {mat - eng:+.1e}')
