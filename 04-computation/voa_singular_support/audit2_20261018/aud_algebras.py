#!/usr/bin/env python3
"""Independent N=2 (NS) and W_3 mode algebras for the audit2 engine (vamod.py).  Written from the brackets stated in
the audit request, not from the author's code.

N=2: [L_m,L_n]=(m-n)L_{m+n}+c/12(m^3-m)d, [L_m,J_n]=-nJ_{m+n}, [L_m,G^±_r]=(m/2-r)G^±_{m+r}, [J_m,J_n]=(c/3)m d,
     [J_m,G^±_r]=±G^±_{m+r}, {G^+_r,G^-_s}=2L_{r+s}+(r-s)J_{r+s}+(c/3)(r^2-1/4)d, {G^±,G^±}=0.
     anti-involution: L_n->L_{-n}, J_n->J_{-n}, G^±_r->G^∓_{-r}.
W_3: [L_m,W_n]=(2m-n)W_{m+n};
     [W_m,W_n]=c/360 m(m^2-1)(m^2-4)d + (m-n)[(m+n+3)(m+n+2)/15-(m+2)(n+2)/6]L_{m+n} + beta(m-n)Lambda_{m+n},
     beta=16/(22+5c), Lambda_k = sum_{a<=-2} L_a L_{k-a} + sum_{a>=-1} L_{k-a} L_a - (3/10)(k+2)(k+3)L_k.
     anti-involution: L_n->L_{-n}, W_n->W_{-n}.
"""
from fractions import Fraction
from vamod import ModeVA

L_, J_, GP, GM = 0, 1, 2, 3


class N2(ModeVA):
    def __init__(self, c, P=2147483647):
        super().__init__(P)
        self.c = Fraction(c)
        # (name, H2, parity, charge)
        self.gens = [('L', 4, 0, 0), ('J', 2, 0, 0), ('G+', 3, 1, 1), ('G-', 3, 1, -1)]

    def omega(self, X):
        t, n2 = X
        if t == GP:
            return (GM, -n2)
        if t == GM:
            return (GP, -n2)
        return (t, -n2)

    def bracket(self, X, Y):
        (ta, a2), (tb, b2) = X, Y
        m, n = Fraction(a2, 2), Fraction(b2, 2)
        c = self.c
        s2 = a2 + b2
        terms, cen = [], Fraction(0)
        if ta == L_ and tb == L_:
            terms.append(((L_, s2), m - n))
            if s2 == 0:
                cen = c / 12 * (m ** 3 - m)
        elif ta == L_ and tb == J_:
            terms.append(((J_, s2), -n))
        elif ta == J_ and tb == L_:
            terms.append(((J_, s2), m))
        elif ta == L_ and tb in (GP, GM):
            terms.append(((tb, s2), m / 2 - n))
        elif ta in (GP, GM) and tb == L_:
            terms.append(((ta, s2), -(n / 2 - m)))
        elif ta == J_ and tb == J_:
            if s2 == 0:
                cen = c / 3 * m
        elif ta == J_ and tb in (GP, GM):
            terms.append(((tb, s2), 1 if tb == GP else -1))
        elif ta in (GP, GM) and tb == J_:
            terms.append(((ta, s2), -1 if ta == GP else 1))
        elif ta == GP and tb == GM:
            terms.append(((L_, s2), 2))
            terms.append(((J_, s2), m - n))
            if s2 == 0:
                cen = c / 3 * (m * m - Fraction(1, 4))
        elif ta == GM and tb == GP:
            terms.append(((L_, s2), 2))
            terms.append(((J_, s2), n - m))
            if s2 == 0:
                cen = c / 3 * (n * n - Fraction(1, 4))
        else:  # G+G+ or G-G-
            pass
        return [(Z, self.F(co)) for Z, co in terms if co != 0], self.F(cen)

    def bracket_apply(self, X, Y, rest):
        terms, cen = self.bracket(X, Y)
        return self.lie_bracket_apply(terms, cen, rest)


LW, WW = 0, 1


class W3(ModeVA):
    def __init__(self, c, P=2147483647):
        super().__init__(P)
        self.c = Fraction(c)
        assert 22 + 5 * self.c != 0
        self.beta = Fraction(16) / (22 + 5 * self.c)
        self.gens = [('L', 4, 0, 0), ('W', 6, 0, 0)]

    def omega(self, X):
        return (X[0], -X[1])

    def Lambda_apply(self, k, rest):
        P = self.P
        wr = self.wt2(rest) // 2
        out = {}

        def acc(vec, f):
            for mm, cc in vec.items():
                out[mm] = (out.get(mm, 0) + cc * f) % P
        for a in range(k - wr, -1):  # a <= -2 and k-a <= wr
            if a > -2:
                continue
            inner = self.apply((LW, 2 * (k - a)), rest)
            if inner:
                acc(self.apply_vec((LW, 2 * a), inner), 1)
        for a in range(-1, wr + 1):
            inner = self.apply((LW, 2 * a), rest)
            if inner:
                acc(self.apply_vec((LW, 2 * (k - a)), inner), 1)
        f = self.F(Fraction(-3, 10) * (k + 2) * (k + 3))
        if f:
            acc(self.apply((LW, 2 * k), rest), f)
        return {mm: cc for mm, cc in out.items() if cc}

    def bracket_apply(self, X, Y, rest):
        (ta, a2), (tb, b2) = X, Y
        m, n = a2 // 2, b2 // 2
        c = self.c
        s2 = a2 + b2
        terms, cen = [], Fraction(0)
        extra = None
        if ta == LW and tb == LW:
            terms.append(((LW, s2), Fraction(m - n)))
            if s2 == 0:
                cen = c / 12 * (m ** 3 - m)
        elif ta == LW and tb == WW:
            terms.append(((WW, s2), Fraction(2 * m - n)))
        elif ta == WW and tb == LW:
            terms.append(((WW, s2), -Fraction(2 * n - m)))
        else:
            k = m + n
            co = (m - n) * (Fraction((k + 3) * (k + 2), 15) - Fraction((m + 2) * (n + 2), 6))
            terms.append(((LW, s2), co))
            if s2 == 0:
                cen = c / 360 * m * (m * m - 1) * (m * m - 4)
            if m != n:
                extra = (k, self.F(self.beta * (m - n)))
        out = self.lie_bracket_apply([(Z, self.F(co)) for Z, co in terms if co != 0], self.F(cen), rest)
        if extra is not None:
            k, f = extra
            P = self.P
            for mm, cc in self.Lambda_apply(k, rest).items():
                out[mm] = (out.get(mm, 0) + cc * f) % P
            out = {mm: cc for mm, cc in out.items() if cc}
        return out
