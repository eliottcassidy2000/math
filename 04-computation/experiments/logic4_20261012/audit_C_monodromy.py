"""Audit C of THM-4625 (independent auditor, 2026-10-12): statement (5), monodromy S_(d+1).

  (a) Gao's G: resultant_r2 of the difference quotients = 9 (r1^2-4r1-2)(r1^2-4r1+2)^2/64 (re-derived). Which roots are GENUINE coincidence
      pairs r1 != r2 with p(r1) = p(r2), q(r1) = q(r2)?  (r1^2-4r1+2 = p'(r1)/3: the double factor is the DIAGONAL r2 = r1, not a coincidence.)
  (b) the hypothesis 'w -> (p(w), q(w)) generically injective' is automatic: q'/p' = w/2 recovers w from the tangent direction; checked by
      resultant for random integer p of degree 2..6, and the genuine coincidence pairs are finite in every case.
  (c) specialization check of the Galois group: Phi(W) = q(W) + W(bc - p(W))/2 - a c^2 over Q at random rational (a, b, c) has Galois group
      S_(d+1) for THM-1300 (S_3), Gao's G (S_4), and random p of degree 4, 5 (S_5, S_6); since the Galois group of a separable specialization
      is a subgroup of the generic one, this certifies S_(d+1) over K(F) (d = 2..5) independently of the monodromy argument.
  (d) critical values of F1 on the (F2,F3)-fibre are q(r_i)/c^2 at the roots r_i of p = cb; for THM-1300 no two coincide for ANY b (pairs
      with p(r1) = p(r2), q(r1) = q(r2) force r1 = r2).
"""
import sympy as sp, random
from sympy.polys.numberfields.galoisgroups import galois_group
w, s, r1, r2, W, a, b, c = sp.symbols('w s r1 r2 W a b c')
random.seed(4625)

def mkq(p):
    return sp.expand(sp.integrate((s / 2) * sp.diff(p, w).subs(w, s), (s, 0, w)))

def coincidences(p, q, label):
    e1 = sp.cancel((p.subs(w, r1) - p.subs(w, r2)) / (r1 - r2)); e2 = sp.cancel((q.subs(w, r1) - q.subs(w, r2)) / (r1 - r2))
    res = sp.factor(sp.resultant(sp.Poly(e1, r2), sp.Poly(e2, r2)))
    print(label, ": resultant_r2 =", res)
    genuine = []
    for fac, mult in sp.factor_list(res)[1]:
        if sp.Poly(fac, r1).degree() == 0:
            continue
        for rt in sp.Poly(fac, r1).nroots(n=40):
            z1 = sp.Poly(sp.expand(e1.subs(r1, rt)), r2).nroots(n=40); z2 = sp.Poly(sp.expand(e2.subs(r1, rt)), r2).nroots(n=40)
            common = [u for u in z1 if any(abs(u - v) < 1e-20 for v in z2)]
            others = [u for u in common if abs(u - rt) > 1e-20]
            genuine += [(rt, rr) for rr in others]
            print("    factor %s (mult %d), r1 ~ %s: common roots r2 ~ %s -> genuine (r2 != r1): %s" % (fac, mult, sp.N(rt, 8), [sp.N(u, 8) for u in common], [sp.N(u, 8) for u in others]))
    print("    genuine coincidence pairs (r1, r2) ~", [(sp.N(u, 8), sp.N(v, 8)) for u, v in genuine], " count:", len(genuine))
    for u, v in genuine:
        print("      |p(r1) - p(r2)| =", sp.N(abs(p.subs(w, u) - p.subs(w, v)), 3), " |q(r1) - q(r2)| =", sp.N(abs(q.subs(w, u) - q.subs(w, v)), 3), " cb = p(r1) ~", sp.N(p.subs(w, u), 8))
    return genuine

print("(a) Gao's G")
pG = w**3 - 6 * w**2 + 6 * w; qG = mkq(pG)
print("   p'(w)/3 =", sp.factor(sp.diff(pG, w) / 3))
coincidences(pG, qG, "   Gao G")
r=2+sp.sqrt(6); r_=2-sp.sqrt(6)
print("   exact: p(2+sqrt6) - p(2-sqrt6) =", sp.simplify(pG.subs(w,r)-pG.subs(w,r_)), "; q(2+sqrt6) - q(2-sqrt6) =", sp.simplify(qG.subs(w,r)-qG.subs(w,r_)), "; p(2+sqrt6) =", sp.simplify(pG.subs(w,r)), "; critical value q/c^2 with q(2+sqrt6) =", sp.simplify(qG.subs(w,r)))
print("   so for cb = -4 (i.e. b = -4/c) two of the three critical points r = 2 +- sqrt6 of F1 on the (F2,F3)-fibre share the critical value; the third root of p = -4 is w = 2.")
print("\n(d) THM-1300")
pT = -3 * w**2 + 4 * w; qT = mkq(pT)
coincidences(pT, qT, "   THM-1300")

print("\n(b) random integer p of degree 2..6: resultant nonzero (hypothesis of (5) automatic), genuine pairs finite")
for d in range(2, 7):
    p = sum(random.randint(-5, 5) * w**i for i in range(1, d)) + random.choice([1, -1, 2, -2, 3]) * w**d
    q = mkq(p)
    e1 = sp.cancel((p.subs(w, r1) - p.subs(w, r2)) / (r1 - r2)); e2 = sp.cancel((q.subs(w, r1) - q.subs(w, r2)) / (r1 - r2))
    res = sp.resultant(sp.Poly(e1, r2), sp.Poly(e2, r2))
    # the diagonal contributes the factor p'(r1)^2 (double) when d >= 2? check divisibility by p'(r1)
    pp = sp.Poly(sp.diff(p, w).subs(w, r1), r1)
    quo, rem = sp.div(sp.Poly(res, r1), pp**2)
    print("   d=%d p=%s: resultant nonzero %s, degree %d, divisible by p'(r1)^2: %s" % (d, p, sp.expand(res) != 0, sp.Poly(res, r1).degree(), rem.is_zero))

print("\n(c) Galois groups of separable specializations of Phi(W) = q(W) + W(bc - p(W))/2 - a c^2 (subgroup of the generic group)")
def galois_of(p, label):
    q = mkq(p); d = sp.Poly(p, w).degree()
    for trial in range(3):
        av, bv, cv = [sp.Rational(random.randint(1, 9), random.randint(1, 4)) * random.choice([1, -1]) for _ in range(3)]
        Phi = sp.expand(q.subs(w, W) + W * (bv * cv - p.subs(w, W)) / 2 - av * cv**2)
        Pp = sp.Poly(Phi, W)
        Pp = sp.Poly(sp.expand(Phi * sp.lcm([sp.denom(co) for co in Pp.all_coeffs()])), W)
        lc = Pp.LC()
        # make monic integer: substitute W = V/lc and multiply by lc^(deg-1)
        V = sp.Symbol('V')
        mon = sp.Poly(sp.expand(Pp.as_expr().subs(W, V / lc) * lc**(Pp.degree() - 1)), V)
        if not mon.is_irreducible or sp.discriminant(mon) == 0:
            print("   %s trial %d: specialization not irreducible/separable, skipped" % (label, trial)); continue
        G, alt = galois_group(mon)
        print("   %s (d=%d, deg Phi=%d) at (a,b,c)=(%s,%s,%s): Galois group %s, order %d, is S_%d: %s" % (label, d, mon.degree(), av, bv, cv, G, G.order(), d + 1, G.order() == sp.factorial(d + 1)))
        break
galois_of(pT, "THM-1300")
galois_of(pG, "Gao G")
galois_of(w**4 - 2 * w**3 + 3 * w**2 + w, "random d=4")
galois_of(2 * w**5 + w**4 - 3 * w**3 + w**2 - w, "random d=5")
