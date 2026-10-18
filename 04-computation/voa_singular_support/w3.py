#!/usr/bin/env python3
"""W_3 algebra (Zamolodchikov) as a ModeAlgebra with a nonlinear bracket.
Types: 0 = L (h=2), 1 = W (h=3).  Normalization W(z)W(w) ~ (c/3)/(z-w)^6 + ...
 [L_m,L_n] = (m-n)L_{m+n} + c/12 (m^3-m) delta
 [L_m,W_n] = (2m-n) W_{m+n}
 [W_m,W_n] = c/360 m(m^2-1)(m^2-4) delta + (m-n)[(m+n+3)(m+n+2)/15 - (m+2)(n+2)/6] L_{m+n} + beta (m-n) Lambda_{m+n},
 beta = 16/(22+5c),  Lambda_k = sum_{n<=-2} L_n L_{k-n} + sum_{n>=-1} L_{k-n} L_n - (3/10)(k+2)(k+3) L_k."""
from fractions import Fraction
from lie_voa_ss import ModeAlgebra
L, W = 0, 1

class W3(ModeAlgebra):
    name = "W3"
    def __init__(self, c, P):
        self.gens = [("L", 2, 0, ()), ("W", 3, 0, ())]
        self.c = Fraction(c)
        assert 22 + 5 * self.c != 0
        self.beta = Fraction(16) / (22 + 5 * self.c)
        super().__init__(P)
        self.betaP = self.F(self.beta)

    def bracket(self, m1, m2):
        a, ta = -m1[0], m1[1]
        b, tb = -m2[0], m2[1]
        n = a + b
        c = self.c
        terms, cen = [], Fraction(0)
        if ta == L and tb == L:
            if a - b: terms.append(((-n, L), Fraction(a - b)))
            if n == 0: cen = c * Fraction(a ** 3 - a, 12)
        elif ta == L and tb == W:
            if 2 * a - b: terms.append(((-n, W), Fraction(2 * a - b)))
        elif ta == W and tb == L:
            if 2 * b - a: terms.append(((-n, W), -Fraction(2 * b - a)))
        else:  # W, W : linear part
            co = Fraction(a - b) * (Fraction((n + 3) * (n + 2), 15) - Fraction((a + 2) * (b + 2), 6))
            if co: terms.append(((-n, L), co))
            if n == 0: cen = c / 360 * a * (a * a - 1) * (a * a - 4)
        return [(m, self.F(co)) for m, co in terms], self.F(cen)

    def bracket_extra(self, m1, m2, rest):
        if not (m1[1] == W and m2[1] == W):
            return {}
        a, b = -m1[0], -m2[0]
        if a == b:
            return {}
        k = a + b
        coef = self.betaP * self.F(Fraction(a - b)) % self.P
        lam = self.apply_Lambda(k, rest)
        return {mo: co * coef % self.P for mo, co in lam.items()}

    def apply_Lambda(self, k, rest):
        """Lambda_k applied to monomial rest."""
        P = self.P
        wr = sum(w for (w, t) in rest)
        out = {}
        def addvec(vec, f):
            for mo, co in vec.items():
                out[mo] = (out.get(mo, 0) + co * f) % P
        # sum_{n <= -2} L_n L_{k-n} rest : need k - n <= wr  -> n >= k - wr
        for n in range(min(-2, 10**9), k - wr - 1, -1) if False else range(k - wr, -1):
            if n > -2: continue
            inner = self.apply((n * -1, L) if False else (-(k - n), L), rest)  # L_{k-n}: mode w = -(k-n)
            if inner:
                addvec(self.apply_vec((-n, L), inner), 1)   # L_n: w = -n
        # sum_{n >= -1} L_{k-n} L_n rest : need n <= wr
        for n in range(-1, wr + 1):
            inner = self.apply((-n, L), rest)
            if inner:
                addvec(self.apply_vec((-(k - n), L), inner), 1)
        # -(3/10)(k+2)(k+3) L_k rest
        f = self.F(Fraction(-3 * (k + 2) * (k + 3), 10))
        if f:
            addvec(self.apply((-k, L), rest), f)
        return {mo: co for mo, co in out.items() if co}

    def positive_ops(self):
        return [(-1, L), (-2, L), (-1, W), (-2, W)]
