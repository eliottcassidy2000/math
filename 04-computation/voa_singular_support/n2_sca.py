#!/usr/bin/env python3
"""N=2 superconformal algebra (NS sector) as a mode Lie superalgebra for lie_voa_ss.
Types: 0=L (h=2), 1=J (h=1), 2=G+ (h=3/2, odd, charge +1), 3=G- (h=3/2, odd, charge -1).
 [L_m,L_n]=(m-n)L+c/12(m^3-m); [L_m,J_n]=-nJ; [L_m,G_r]=(m/2-r)G; [J_m,J_n]=c/3 m delta
 [J_m,G+-_r]=+-G+-; {G+_r,G-_s}=2L+(r-s)J+c/3(r^2-1/4)delta; {G+,G+}={G-,G-}=0.
Minimal models: c = 3k/(k+2) (chiral algebra of the LG model W = x^(k+2))."""
from fractions import Fraction
from lie_voa_ss import ModeAlgebra
L, J, GP, GM = range(4)

class N2SCA(ModeAlgebra):
    name = "N2"
    def __init__(self, c, P):
        h32 = Fraction(3, 2)
        self.gens = [("L", 2, 0, (0,)), ("J", 1, 0, (0,)), ("G+", h32, 1, (1,)), ("G-", h32, 1, (-1,))]
        self.c = Fraction(c)
        super().__init__(P)

    def _raw(self, a, ta, b, tb):
        c = self.c
        n = a + b
        out, cen = {}, Fraction(0)
        def add(t, co):
            if co: out[(-n, t)] = out.get((-n, t), 0) + Fraction(co)
        if ta == L and tb == L:
            add(L, a - b)
            if n == 0: cen = c * Fraction(a ** 3 - a, 12)
        elif ta == L and tb == J: add(J, -b)
        elif ta == J and tb == L: add(J, a)
        elif ta == L and tb in (GP, GM): add(tb, Fraction(a, 2) - b)
        elif ta in (GP, GM) and tb == L: add(ta, -(Fraction(b, 2) - a))
        elif ta == J and tb == J:
            if n == 0: cen = c / 3 * a
        elif ta == J and tb == GP: add(GP, 1)
        elif ta == J and tb == GM: add(GM, -1)
        elif ta == GP and tb == J: add(GP, -1)
        elif ta == GM and tb == J: add(GM, 1)
        elif ta == GP and tb == GM:
            add(L, 2); add(J, a - b)
            if n == 0: cen = c / 3 * (a * a - Fraction(1, 4))
        elif ta == GM and tb == GP:
            add(L, 2); add(J, b - a)
            if n == 0: cen = c / 3 * (b * b - Fraction(1, 4))
        return out, cen

    def bracket(self, m1, m2):
        a, ta = -m1[0], m1[1]
        b, tb = -m2[0], m2[1]
        out, cen = self._raw(Fraction(a), ta, Fraction(b), tb)
        terms = [((w if not isinstance(w, Fraction) or w.denominator != 1 else int(w), t), self.F(co))
                 for (w, t), co in out.items() if co]
        return terms, self.F(cen)

    def positive_ops(self):
        h = Fraction(1, 2)
        return [(-h, GP), (-h, GM), (-1, L), (-2, L), (-1, J)]
