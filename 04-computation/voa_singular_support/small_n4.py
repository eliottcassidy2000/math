#!/usr/bin/env python3
"""Small N=4 superconformal algebra (NS sector) as a mode Lie superalgebra for lie_voa_ss.
Types: 0=L (h=2), 1=e, 2=h, 3=f (sl2 currents, h=1, level k=c/6),
       4=G1 (charge -1), 5=G2 (+1), 6=Gb1 (+1), 7=Gb2 (-1)  (h=3/2, odd).
Charges are eigenvalues of h_0 (= 2 J^0).
Brackets (from the standard su(2) form with J^+=e, J^-=f, J^0=h/2):
 [L_m,L_n]=(m-n)L+c/12(m^3-m); [L_m,X_n]=-n X_{m+n} (X current); [L_m,G_r]=(m/2-r)G
 [h_m,h_n]=2k m; [h_m,e_n]=2e; [h_m,f_n]=-2f; [e_m,f_n]=h+k m
 [h_m,G1_r]=-G1; [h_m,G2_r]=G2; [h_m,Gb1_r]=Gb1; [h_m,Gb2_r]=-Gb2
 [e_m,G1_r]=-G2; [f_m,G2_r]=-G1; [e_m,Gb2_r]=Gb1; [f_m,Gb1_r]=Gb2
 {G1_r,Gb1_s}=2L-(r-s)h+c/3(r^2-1/4); {G2_r,Gb2_s}=2L+(r-s)h+c/3(r^2-1/4)
 {G1_r,Gb2_s}=-2(r-s)f; {G2_r,Gb1_s}=-2(r-s)e; {G,G}={Gb,Gb}=0."""
from fractions import Fraction
from lie_voa_ss import ModeAlgebra

L, E, H, Fm, G1, G2, B1, B2 = range(8)
ODD = {G1, G2, B1, B2}

class SmallN4(ModeAlgebra):
    name = "smallN4"
    def __init__(self, c, P):
        h32 = Fraction(3, 2)
        self.gens = [("L", 2, 0, (0,)), ("e", 1, 0, (2,)), ("h", 1, 0, (0,)), ("f", 1, 0, (-2,)),
                     ("G1", h32, 1, (-1,)), ("G2", h32, 1, (1,)), ("Gb1", h32, 1, (1,)), ("Gb2", h32, 1, (-1,))]
        self.c = Fraction(c)
        self.k = self.c / 6
        super().__init__(P)

    def _raw(self, a, ta, b, tb):
        """[X^{ta}_a, Y^{tb}_b} for mode indices a,b (X_a = mode with w=-a). Returns dict, central (Fractions)."""
        c, k = self.c, self.k
        n = a + b
        out = {}
        cen = Fraction(0)
        def add(t, co):
            if co:
                out[(-n, t)] = out.get((-n, t), 0) + Fraction(co)
        cur = {E, H, Fm}
        G = {G1, G2, B1, B2}
        if ta == L and tb == L:
            add(L, a - b)
            if n == 0: cen = c * Fraction(a ** 3 - a, 12)
        elif ta == L and tb in cur:
            add(tb, -b)
        elif ta in cur and tb == L:
            add(ta, a)          # [X_a, L_b] = -[L_b, X_a] = a X_{a+b}
        elif ta == L and tb in G:
            add(tb, Fraction(a, 2) - b)
        elif ta in G and tb == L:
            add(ta, -(Fraction(b, 2) - a))
        elif ta in cur and tb in cur:
            if ta == H and tb == H:
                if n == 0: cen = 2 * k * a
            elif ta == H and tb == E: add(E, 2)
            elif ta == E and tb == H: add(E, -2)
            elif ta == H and tb == Fm: add(Fm, -2)
            elif ta == Fm and tb == H: add(Fm, 2)
            elif ta == E and tb == Fm:
                add(H, 1)
                if n == 0: cen = k * a
            elif ta == Fm and tb == E:
                add(H, -1)
                if n == 0: cen = k * a
        elif ta in cur and tb in G:
            tab = {(H, G1): (G1, -1), (H, G2): (G2, 1), (H, B1): (B1, 1), (H, B2): (B2, -1),
                   (E, G1): (G2, -1), (Fm, G2): (G1, -1), (E, B2): (B1, 1), (Fm, B1): (B2, 1)}
            if (ta, tb) in tab:
                t, co = tab[(ta, tb)]
                add(t, co)
        elif ta in G and tb in cur:
            o, cc = self._raw(b, tb, a, ta)
            for (w, t), co in o.items():
                out[(w, t)] = -co
            cen = -cc
        elif ta in G and tb in G:
            r, s = a, b
            if {ta, tb} <= {G1, G2} or {ta, tb} <= {B1, B2}:
                pass
            else:
                if ta in (B1, B2):   # {Gb_r, G_s} = {G_s, Gb_r}
                    return self._raw(b, tb, a, ta)
                if ta == G1 and tb == B1:
                    add(L, 2); add(H, -(r - s))
                    if n == 0: cen = c / 3 * (r * r - Fraction(1, 4))
                elif ta == G2 and tb == B2:
                    add(L, 2); add(H, (r - s))
                    if n == 0: cen = c / 3 * (r * r - Fraction(1, 4))
                elif ta == G1 and tb == B2:
                    add(Fm, -2 * (r - s))
                elif ta == G2 and tb == B1:
                    add(E, -2 * (r - s))
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
        return [(0, E), (-1, Fm), (-h, G1), (-h, G2), (-h, B1), (-h, B2), (-1, L), (-2, L), (-1, H)]
