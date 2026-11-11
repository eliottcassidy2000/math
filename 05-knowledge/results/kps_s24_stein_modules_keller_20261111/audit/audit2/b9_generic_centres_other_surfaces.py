"""b9_generic_centres_other_surfaces.py - second audit, Theorem 5.4: sanity runs on two surfaces the author did not
use.  Pseudo-random residues modulo 2^61 - 1 emulate a centre in general position (sanity check, not a proof;
non-membership is relative to the truncation).

  torus : inclusion C^* x C^* -> C^2 (ring of Laurent polynomials); expected dim S_z = b_1(Y) - 2 + b_2(X) = 3 - 2 + 1 = 2
          at a centre with a b != 0, and S_z = 0 if a = 0 or b = 0;
  hyp   : inclusion of C^2 minus {1 + xy = 0} (README 6.3); dimension 3 at generic centres, 2 at (0, 0).
"""
import random
import eng2
from eng2 import RingL, RingH, add, scal, Echelon, power

random.seed(20261111)
P = eng2.PRIME


def torus(av, bv, Np=20, N=10):
    R = RingL()
    U = add(R.mono(1, 0), scal(R.one(), av), -1)
    V = add(R.mono(0, 1), scal(R.one(), bv), -1)
    W = R.mul(U, V)
    Wx, Wy = R.dx(W), R.dy(W)
    wt = lambda m: abs(m[0]) + abs(m[1])
    E = Echelon(lambda m: (wt(m), m))
    for k in range(-Np, Np + 1):
        for l in range(-Np, Np + 1):
            if abs(k) + abs(l) <= Np:
                m = R.mono(k, l)
                E.insert(add(R.dx(m), R.mul(Wx, m), -1))
                E.insert(add(R.dy(m), R.mul(Wy, m), -1))
    dimV = sum(1 for k in range(-N, N + 1) for l in range(-N, N + 1) if abs(k) + abs(l) <= N)
    cod = dimV - sum(1 for m in E.piv if wt(m) <= N)
    mons = [(k, l) for k in range(-4, 5) for l in range(-4, 5) if abs(k) + abs(l) <= 5]
    zero = [m for m in mons if E.contains(R.mono(*m))]
    return cod, len(zero), len(mons)


def hyp(av, bv, Np=30, N=14):
    R = RingH()
    x, y, one = R.x(), R.y(), R.one()
    U = add(x, scal(one, av), -1)
    V = add(y, scal(one, bv), -1)
    W = R.mul(U, V)
    Wx, Wy = R.dx(W), R.dy(W)
    wt = lambda m: abs(m[0]) + 2 * abs(m[1])
    E = Echelon(lambda m: (wt(m), m))
    for c in range(-Np, Np + 1):
        for n in range(-Np, Np + 1):
            if abs(c) + 2 * abs(n) <= Np:
                h = {(c, n): 1}
                E.insert(add(R.dx(h), R.mul(Wx, h), -1))
                E.insert(add(R.dy(h), R.mul(Wy, h), -1))
    dimV = sum(1 for c in range(-N, N + 1) for n in range(-N, N + 1) if abs(c) + 2 * abs(n) <= N)
    cod = dimV - sum(1 for m in E.piv if wt(m) <= N)
    mons = [(c, n) for c in range(-5, 6) for n in range(-3, 4) if abs(c) + 2 * abs(n) <= 8]
    zero = [m for m in mons if E.contains({m: 1})]
    return cod, len(zero), len(mons)


if __name__ == "__main__":
    cents = [(random.randrange(P), random.randrange(P)) for _ in range(3)] + [(2, 3), (0, 3), (2, 0), (0, 0)]
    for (a, b) in cents:
        cod, z, t = torus(a, b)
        print("torus, centre (%d, %d): c(10) = %d; Laurent monomials (weight <= 5) with zero class: %d of %d" % (a, b, cod, z, t))
    for (a, b) in cents:
        cod, z, t = hyp(a, b)
        print("hyp,   centre (%d, %d): c(14) = %d; basis monomials x^c q^n / y^c q^n with zero class: %d of %d" % (a, b, cod, z, t))
