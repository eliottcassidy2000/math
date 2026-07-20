"""
opus-2026-07-19-S414 (HYP-8070): the JC counterexample's TORUS STRUCTURE, the
derivation of its rationals, and the family probe.

(W1) EQUIVARIANCE: F intertwines (x,y,z)->(lx, y/l, z/l^2) with target weights
     (-2,-1,1): F1(l.p) = l^-2 F1(p), F2 -> l^-1 F2, F3 -> l F3  [symbolic].
(W2) NORMAL FORM: with s = xy, w = x^2 z:
     F = ( z*P + y^2*Q,  y*R + x*z*S,  x*T ),
     P=(1+s)^3, Q=(1+s)(4+3s), R=1+12s+9s^2, S=3(1+s)^2, T=2-3s-w  [symbolic].
(W3) THE COLLISION IS A WHOLE C*-FAMILY: for every l != 0,
     F(l, -3/(2l), 13/(2l^2)) = F(0,0,-1/(4 l^2)) = (-1/(4 l^2), 0, 0)
     -- the published points are the l = 1 slice  [symbolic in l].
(W4) THE QUOTIENT CONTRACTION: q(s,w) := (T*(s R + w S), T^2*(w P + s^2 Q)) is the
     induced map on invariants (verified: F2*F3 and F1*F3^2 pull back to q's
     components); q contracts the whole curve {T=0}; the exotic orbit =
     {T=0} n {sR+wS=0}: substituting w = 2-3s gives 4s+6 = 0 -- s = -3/2,
     w = 13/2 ARE DERIVED, not data.  Also det J(q) computed (expect NON-const:
     the plane map is a contraction, not a plane Keller counterexample).
(W5) FAMILY PROBE: general equivariant ansatz P = a0+a1 s+a2 s^2+a3 s^3,
     Q = b0+b1 s+b2 s^2, R = c0+c1 s+c2 s^2, S = d0+d1 s+d2 s^2,
     T = e0+e1 s - w; impose det J(F) == const: extract the coefficient
     equations, count the solution-variety's dimension near the known point,
     and exhibit at least one NEW explicit counterexample if a deformation
     direction survives (with its own collision re-derived).
"""
import sympy as sp

x, y, z, l, s, w = sp.symbols('x y z l s w')
u = 1 + x*y

F1 = u**3*z + y**2*u*(4+3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
F = sp.Matrix([F1, F2, F3])

print("(W1) torus equivariance:")
sub = {x: l*x, y: y/l, z: z/l**2}
checks = [sp.simplify(F1.subs(sub, simultaneous=True) - F1/l**2),
          sp.simplify(F2.subs(sub, simultaneous=True) - F2/l),
          sp.simplify(F3.subs(sub, simultaneous=True) - F3*l)]
print("  weights (1,-1,-2) -> (-2,-1,1):", [c == 0 for c in checks])

print("(W2) (s,w) normal form:")
P = (1+s)**3; Q = (1+s)*(4+3*s); R = 1+12*s+9*s**2; S = 3*(1+s)**2; T = 2-3*s-w
sub_sw = {s: x*y, w: x**2*z}
nf = [sp.expand(z*P.subs(sub_sw) + y**2*Q.subs(sub_sw) - F1),
      sp.expand(y*R.subs(sub_sw) + x*z*S.subs(sub_sw) - F2),
      sp.expand(x*T.subs(sub_sw) - F3)]
print("  F == (zP + y^2 Q, yR + xzS, xT):", [e == 0 for e in nf])

print("(W3) the collision C*-family (symbolic in l):")
p_orbit = {x: l, y: sp.Rational(-3,2)/l, z: sp.Rational(13,2)/l**2}
p_line  = {x: 0, y: 0, z: sp.Rational(-1,4)/l**2}
im1 = sp.Matrix([sp.simplify(f.subs(p_orbit)) for f in F])
im2 = sp.Matrix([sp.simplify(f.subs(p_line)) for f in F])
print(f"  F(orbit(l)) = {list(im1)};  F(line(l)) = {list(im2)}; equal: "
      f"{sp.simplify(im1-im2) == sp.zeros(3,1)}")

print("(W4) the quotient map and the derived rationals:")
qA = sp.expand(T*(s*R + w*S)); qB = sp.expand(T**2*(w*P + s**2*Q))
pull1 = sp.expand((F2*F3).subs({}) - (qA.subs(sub_sw)))
pull2 = sp.expand((F1*F3**2) - (qB.subs(sub_sw)))
print("  F2*F3 == qA(s,w), F1*F3^2 == qB(s,w):", pull1 == 0, pull2 == 0)
lin = sp.expand((s*R + w*S).subs(w, 2-3*s))
print(f"  (sR+wS)|_(T=0) = {lin}  -> root s = {sp.solve(lin, s)}, "
      f"w = {sp.simplify((2-3*s).subs(s, sp.Rational(-3,2)))}  [(-3/2, 13/2) DERIVED]")
detq = sp.expand(sp.Matrix([[sp.diff(qA,s), sp.diff(qA,w)],
                            [sp.diff(qB,s), sp.diff(qB,w)]]).det())
print(f"  det J(q) constant? {sp.degree(sp.Poly(detq, s, w)) == 0} "
      f"(a contraction, not a plane Keller map -- JC2 untouched)")

print("(W5) family probe (equivariant ansatz, det == const):")
a = sp.symbols('a0:4'); b = sp.symbols('b0:3'); c = sp.symbols('c0:3')
d = sp.symbols('d0:3'); e = sp.symbols('e0:2')
Pg = a[0]+a[1]*s+a[2]*s**2+a[3]*s**3
Qg = b[0]+b[1]*s+b[2]*s**2
Rg = c[0]+c[1]*s+c[2]*s**2
Sg = d[0]+d[1]*s+d[2]*s**2
Tg = e[0]+e[1]*s-w
G = sp.Matrix([z*Pg.subs(sub_sw)+y**2*Qg.subs(sub_sw),
               y*Rg.subs(sub_sw)+x*z*Sg.subs(sub_sw),
               x*Tg.subs(sub_sw)])
detG = sp.expand(sp.Matrix(G).jacobian([x, y, z]).det())
# det must be a constant: all monomial coefficients in (x,y,z) vanish except 1
poly = sp.Poly(detG, x, y, z)
eqs, const = [], None
for mono, coeff in poly.terms():
    if mono == (0,0,0): const = coeff
    else: eqs.append(sp.expand(coeff))
eqs = [e_ for e_ in sp.groebner(eqs, *a, *b, *c, *d, *e).exprs] if False else eqs
known = {a[0]:1,a[1]:3,a[2]:3,a[3]:1, b[0]:4,b[1]:7,b[2]:3,
         c[0]:1,c[1]:12,c[2]:9, d[0]:3,d[1]:6,d[2]:3, e[0]:2,e[1]:-3}
print(f"  #coefficient equations: {len(set(eqs))}; known point satisfies all: "
      f"{all(sp.expand(q_.subs(known)) == 0 for q_ in set(eqs))}; det there = "
      f"{const.subs(known)}")
sols = sp.solve(list(set(eqs)), [*a, *b, *c, *d, *e], dict=True)
print(f"  solve returned {len(sols)} solution branch(es); free parameters per branch:")
for i, so in enumerate(sols[:6]):
    free = [v for v in [*a,*b,*c,*d,*e] if v not in so]
    print(f"   branch {i}: free = {free}")
    if all(sp.simplify(so.get(k, k) - v) is not None for k, v in known.items()):
        pass
# exhibit a NEW example from a branch containing the known point, if any
for so in sols:
    trial = {k: so.get(k, k) for k in [*a,*b,*c,*d,*e]}
    # try to specialize free params to the known values +1 where free
    subst = {}
    ok = True
    for k in [*a,*b,*c,*d,*e]:
        val = so.get(k)
        if val is None: subst[k] = known[k] + 1  # deform the free direction
    inst = {k: sp.simplify(so.get(k, k).subs(subst)) if so.get(k) is not None
            else subst[k] for k in [*a,*b,*c,*d,*e]}
    detI = sp.simplify(detG.subs(inst))
    if detI.free_symbols == set() and detI != 0:
        RgI = Rg.subs(inst); SgI = Sg.subs(inst); TgI = Tg.subs(inst)
        line = sp.expand((s*RgI + w*SgI).subs(w, inst[e[0]]+inst[e[1]]*s))
        roots = sp.solve(line, s)
        print(f"  NEW EXAMPLE candidate: coeffs {inst}; det = {detI}; "
              f"exotic-orbit s-roots on T=0: {roots}")
        break
