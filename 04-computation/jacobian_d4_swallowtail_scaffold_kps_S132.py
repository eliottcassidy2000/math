# -*- coding: utf-8 -*-
# Exact classical scaffolding for the d=4 (swallowtail) layer of THM-2446 (P1).
# Every claim is an exact sympy identity check. Conventions stated inline.
import sympy as sp

p, q, r, m, s, a, u, w, tt, lam = sp.symbols('p q r m s a u w tt lam')
z, y, T, th, c = sp.symbols('z y T theta c')
z1, z2, z3 = sp.symbols('z1 z2 z3')
A, B, C, D = sp.symbols('A B C D')
a1, a2, a3 = sp.symbols('a1 a2 a3')

FAIL = []
def check(name, expr):
    v = sp.expand(expr)
    tag = "PASS" if v == 0 else "FAIL"
    if v != 0:
        FAIL.append(name)
    print(f"[{tag}] {name}" + ("" if v == 0 else f"   RESIDUE={v}"))

def checkbool(name, b):
    tag = "PASS" if b else "FAIL"
    if not b:
        FAIL.append(name)
    print(f"[{tag}] {name}")

print("="*78)
print("PART 1: depressed quartic disc, resolvent cubic, exact disc equality")
print("="*78)
f = z**4 + p*z**2 + q*z + r
Df = sp.expand(sp.discriminant(f, z))
print("disc_z(z^4+p z^2+q z+r) =", Df)
Dtarget = 16*p**4*r - 4*p**3*q**2 - 128*p**2*r**2 + 144*p*q**2*r - 27*q**4 + 256*r**3
check("disc quartic == 16p^4 r -4p^3 q^2 -128 p^2 r^2 +144 p q^2 r -27 q^4 +256 r^3", Df - Dtarget)

g1 = y**3 - 2*p*y**2 + (p**2 - 4*r)*y + q**2      # given convention
g2 = y**3 + 2*p*y**2 + (p**2 - 4*r)*y - q**2      # alternative convention
Dg1 = sp.expand(sp.discriminant(g1, y))
Dg2 = sp.expand(sp.discriminant(g2, y))
check("disc(f) == disc(g1)  [g1 = y^3-2py^2+(p^2-4r)y+q^2]", Df - Dg1)
check("disc(f) == disc(g2)  [g2 = y^3+2py^2+(p^2-4r)y-q^2]", Df - Dg2)

# Root provenance (z1+z2+z3+z4=0): monic depressed quartic prod(z-zi),
# so p=e2, q=-e3, r=e4.
z4 = -(z1 + z2 + z3)
e2 = sp.expand(z1*z2 + z1*z3 + z1*z4 + z2*z3 + z2*z4 + z3*z4)
e3 = sp.expand(z1*z2*z3 + z1*z2*z4 + z1*z3*z4 + z2*z3*z4)
e4 = sp.expand(z1*z2*z3*z4)
sub_sym = {p: e2, q: -e3, r: e4}
res2 = sp.expand(sp.prod([(y - (z1 + zz)**2) for zz in (z2, z3, z4)]))
check("g2 has roots (z1+z2)^2,(z1+z3)^2,(z1+z4)^2", res2 - sp.expand(g2.subs(sub_sym)))
res1 = sp.expand(sp.prod([(y + (z1 + zz)**2) for zz in (z2, z3, z4)]))
check("g1 has roots (zi+zj)(zk+zl) = -(zi+zj)^2", res1 - sp.expand(g1.subs(sub_sym)))
check("(z1+z2)^2-(z1+z3)^2 == (z1-z4)(z2-z3)  [e1=0]",
      (z1+z2)**2 - (z1+z3)**2 - (z1-z4)*(z2-z3))
roots = [z1, z2, z3, z4]
proddiff = sp.expand(sp.prod([(roots[i]-roots[j])**2 for i in range(4) for j in range(i+1, 4)]))
check("disc(f) == prod_{i<j}(zi-zj)^2", sp.expand(Df.subs(sub_sym)) - proddiff)

# General monic quartic and its standard resolvent cubic
F = z**4 + A*z**3 + B*z**2 + C*z + D
R = y**3 - B*y**2 + (A*C - 4*D)*y - (A**2*D - 4*B*D + C**2)
check("general monic quartic: disc(F) == disc(resolvent R)",
      sp.expand(sp.discriminant(F, z)) - sp.expand(sp.discriminant(R, y)))
# resolvent-with-roots z1z2+z3z4 in the depressed case:
R2 = y**3 - p*y**2 - 4*r*y + (4*p*r - q**2)
res3 = sp.expand((y - (z1*z2 + z3*z4))*(y - (z1*z3 + z2*z4))*(y - (z1*z4 + z2*z3)))
check("R2 has roots z1z2+z3z4 etc.", res3 - sp.expand(R2.subs(sub_sym)))
check("disc(f) == disc(R2)", Df - sp.expand(sp.discriminant(R2, y)))

# Invariants I, J and the master relation 27*disc = 4I^3 - J^2
I4 = p**2 + 12*r
J4 = 2*p**3 + 27*q**2 - 72*p*r
check("27*disc(f) == 4*I^3 - J^2, I=p^2+12r, J=2p^3+27q^2-72pr", 27*Df - (4*I4**3 - J4**2))
# depressed form of the resolvent g2: y = w - 2p/3
dep = sp.expand(g2.subs(y, w - 2*p/3))
check("g2(w-2p/3) == w^3 - (I/3) w - J/27", dep - (w**3 - I4*w/sp.Integer(3) - J4/sp.Integer(27)))
dep1 = sp.expand(g1.subs(y, w + 2*p/3))
check("g1(w+2p/3) == w^3 - (I/3) w + J/27", dep1 - (w**3 - I4*w/sp.Integer(3) + J4/sp.Integer(27)))
# quasi-homogeneity: weights (2,3,4)
tq = sp.Symbol('tq')
check("Delta quasi-homog deg 12 for wt(p,q,r)=(2,3,4)",
      sp.expand(Df.subs({p: tq**2*p, q: tq**3*q, r: tq**4*r})) - tq**12*Df)
check("I quasi-homog deg 4", sp.expand(I4.subs({p: tq**2*p, r: tq**4*r})) - tq**4*I4)
check("J quasi-homog deg 6", sp.expand(J4.subs({p: tq**2*p, q: tq**3*q, r: tq**4*r})) - tq**6*J4)

print()
print("="*78)
print("PART 2: swallowtail {disc=0} in (p,q,r): strata, parametrizations, equations")
print("="*78)
gradD = [sp.expand(sp.diff(Df, v)) for v in (p, q, r)]

# --- A2 cuspidal edge: triple root, f=(z-s)^3(z+3s)
fA2 = sp.expand((z - s)**3*(z + 3*s))
print("(z-s)^3(z+3s) =", fA2)
ps, qs, rs = -6*s**2, 8*s**3, -3*s**4
check("edge coefficients (p,q,r)=(-6s^2, 8s^3, -3s^4)",
      fA2 - (z**4 + ps*z**2 + qs*z + rs))
sub_edge = {p: ps, q: qs, r: rs}
check("Delta == 0 on edge", Df.subs(sub_edge))
for i, v in enumerate("pqr"):
    check(f"d Delta/d{v} == 0 on edge", gradD[i].subs(sub_edge))
E2 = 8*p**3 + 27*q**2
check("I == 0 on edge", I4.subs(sub_edge))
check("J == 0 on edge", J4.subs(sub_edge))
check("8p^3+27q^2 == 0 on edge", E2.subs(sub_edge))
check("8p^3+27q^2 == J + 6 p I  (so ideal (I,J) = (I, 8p^3+27q^2))", E2 - (J4 + 6*p*I4))
# elimination: exact ideal of the edge
GB_par = sp.groebner([p - ps, q - qs, r - rs], s, p, q, r, order='lex')
elim = [g for g in GB_par.exprs if not g.has(s)]
print("edge elimination ideal generators (lex, s eliminated):")
for g in elim:
    print("   ", sp.factor(g))
GB_IJ = sp.groebner([I4, E2], p, q, r, order='grevlex')
GB_el = sp.groebner(elim, p, q, r, order='grevlex')
checkbool("every elimination gen lies in ideal (I, 8p^3+27q^2)",
          all(GB_IJ.reduce(g)[1] == 0 for g in elim))
checkbool("I and 8p^3+27q^2 lie in the elimination ideal",
          all(GB_el.reduce(g)[1] == 0 for g in (I4, E2)))

# --- A1A1 self-intersection: two double roots, f=(z^2-a)^2
fA11 = sp.expand((z**2 - a)**2)
sub_node = {p: -2*a, q: sp.Integer(0), r: a**2}
check("(z^2-a)^2 == z^4 -2a z^2 + a^2  [(p,q,r)=(-2a,0,a^2)]",
      fA11 - (z**4 - 2*a*z**2 + a**2))
check("Delta == 0 on node curve", Df.subs(sub_node))
for i, v in enumerate("pqr"):
    check(f"d Delta/d{v} == 0 on node curve", gradD[i].subs(sub_node))
check("node equations: q==0 holds trivially; p^2-4r on param", (p**2 - 4*r).subs(sub_node))

# --- A3 point: quadruple root
quad = sp.expand((z - c)**4)
print("(z-c)^4 =", quad, "   depressed (z^3 coeff -4c = 0) forces c=0 -> f=z^4, (p,q,r)=(0,0,0)")
# the '(z^2+a)^2 shape' is the node sheet, NOT A3, unless a=0:
gtest = (z**2 + a)**2
gg = sp.gcd(sp.expand(gtest), sp.expand(sp.diff(gtest, z)))
print("gcd((z^2+a)^2, d/dz) =", gg, " -> degree 2 with distinct roots +/-sqrt(-a) for a!=0: A1A1, not A3")
check("(z^2+a)^2 lies on node curve: q=0, p^2-4r with p=2a,q=0,r=a^2",
      (p**2 - 4*r).subs({p: 2*a, q: 0, r: a**2}))

# --- Sing(Delta) == edge  UNION  node curve (set-theoretically), exact proof:
Jgens = [Df] + gradD
def in_radical(h, gens):
    G = sp.groebner(gens + [1 - tt*h], p, q, r, tt, order='grevlex')
    return G.contains(sp.Integer(1))
prods = [(I4, 'I'), (E2, '8p^3+27q^2')]
nodeg = [(q, 'q'), (p**2 - 4*r, 'p^2-4r')]
for (gA, nA) in prods:
    for (gB, nB) in nodeg:
        checkbool(f"({nA})*({nB}) in radical of Jacobian ideal (Delta,grad Delta)",
                  in_radical(sp.expand(gA*gB), Jgens))

# --- local transversal structure (orders of vanishing)
# cusp point s=1: P0=(-6,8,-3)
Dser = sp.expand(Df.subs({p: -6 + u*a1, q: 8 + u*a2, r: -3 + u*a3}))
Pu = sp.Poly(Dser, u)
c0 = Pu.coeff_monomial(1); c1 = Pu.coeff_monomial(u); c2 = Pu.coeff_monomial(u**2)
check("edge point P0=(-6,8,-3): order-0 term", c0)
check("edge point P0: order-1 term", c1)
c2f = sp.factor(c2)
print("edge point P0: order-2 term (quadratic form in direction) =", c2f)
# tangent to edge at s=1 is (-12,24,-12) ~ (-1,2,-1); the quadratic form must kill it:
check("edge tangent (-1,2,-1) in kernel of 2-jet", c2.subs({a1: -1, a2: 2, a3: -1}))
# pick kernel direction of the 2-jet independent of the tangent, check cubic term (cusp):
# c2f will be const*(linear)^2 -> kernel plane; report jets on an explicit kernel vector.
lin = sp.factor_list(c2)[1]
print("factor_list of 2-jet:", sp.factor_list(c2))
# node point a=1: P1=(-2,0,1)
DserN = sp.expand(Df.subs({p: -2 + u*a1, q: u*a2, r: 1 + u*a3}))
PuN = sp.Poly(DserN, u)
check("node point P1=(-2,0,1): order-0 term", PuN.coeff_monomial(1))
check("node point P1: order-1 term", PuN.coeff_monomial(u))
c2N = PuN.coeff_monomial(u**2)
print("node point P1: order-2 term =", sp.factor(c2N))
# origin: lowest total-degree part of Delta
polyD = sp.Poly(Df, p, q, r)
mindeg = min(sum(mon) for mon in polyD.monoms())
low = sum(coef*p**e1_*q**e2_*r**e3_ for (e1_, e2_, e3_), coef in
          zip(polyD.monoms(), polyD.coeffs()) if e1_ + e2_ + e3_ == mindeg)
print("origin: mult(Delta) =", mindeg, ", tangent cone =", low)

# restriction of Delta to the symmetric plane q=0 (needed for Part 3):
check("Delta|_{q=0} == 16 r (p^2-4r)^2", Df.subs(q, 0) - sp.expand(16*r*(p**2 - 4*r)**2))

print()
print("="*78)
print("PART 3: quadrisection modulus vs trisection")
print("="*78)
Q4 = 8*T**4 - 8*T**2 + 1
check("8T^4-8T^2+1 == ChebyshevT_4(T)", Q4 - sp.chebyshevt(4, T))
trig = sp.simplify(sp.expand_trig(sp.cos(4*th)) - (8*sp.cos(th)**4 - 8*sp.cos(th)**2 + 1))
checkbool("cos(4 theta) == 8cos^4-8cos^2+1", trig == 0)
D4m = sp.factor(sp.discriminant(Q4 - m, T))
print("disc_T(8T^4-8T^2+1-m) =", D4m)
check("disc_T(8T^4-8T^2+1-m) == 2^17 (1-m)(1+m)^2",
      sp.expand(sp.discriminant(Q4 - m, T)) - sp.expand(131072*(1 - m)*(1 + m)**2))
D3m = sp.factor(sp.discriminant(4*T**3 - 3*T - m, T))
print("disc_T(4T^3-3T-m) =", D3m)
check("disc_T(4T^3-3T-m) == -432(m^2-1)",
      sp.expand(sp.discriminant(4*T**3 - 3*T - m, T)) + 432*(m**2 - 1))
D3n = sp.factor(sp.discriminant(T**3 - 3*T - 2*m, T))
print("disc_T(T^3-3T-2m)  [T=2cos theta norm] =", D3n)
check("disc_T(T^3-3T-2m) == -108(m^2-1)",
      sp.expand(sp.discriminant(T**3 - 3*T - 2*m, T)) + 108*(m**2 - 1))
D4n = sp.factor(sp.discriminant(T**4 - 4*T**2 + 2 - 2*m, T))
print("disc_T(T^4-4T^2+2-2m) [T=2cos theta norm] =", D4n)
check("disc_T(T^4-4T^2+2-2m) == 2^11 (1-m)(1+m)^2",
      sp.expand(sp.discriminant(T**4 - 4*T**2 + 2 - 2*m, T)) - sp.expand(2048*(1 - m)*(1 + m)**2))
check("T=2cos norm: (2T'^... ) chebyshev shift: (T^4-4T^2+2).subs(T,2x)==2*Cheb4(x)",
      sp.expand((T**4 - 4*T**2 + 2).subs(T, 2*T)) - 2*sp.chebyshevt(4, T))
print("fiber at m=+1:", sp.factor(Q4 - 1))
print("fiber at m=-1:", sp.factor(Q4 + 1))
print("d=3 fiber at m=+1:", sp.factor(4*T**3 - 3*T - 1))
print("d=3 fiber at m=-1:", sp.factor(4*T**3 - 3*T + 1))
# inner double cover u=T^2:
Dinner = sp.expand(sp.discriminant(8*u**2 - 8*u + 1 - m, u))
print("disc_u(8u^2-8u+1-m) =", sp.factor(Dinner))
check("inner cover disc == 32(1+m)", Dinner - 32*(1 + m))
# path of the monicized family in (p,q,r): T^4 - T^2 + (1-m)/8
sub_path = {p: sp.Integer(-1), q: sp.Integer(0), r: (1 - m)/8}
check("Delta on Chebyshev path == (1-m)(1+m)^2/2",
      sp.expand(Df.subs(sub_path)) - sp.expand((1 - m)*(1 + m)**2/2))
check("scaling law: disc(8*monic) = 8^6 * disc(monic)",
      sp.expand(sp.discriminant(Q4 - m, T)) - 8**6*sp.expand(Df.subs(sub_path)))
print("cusp-edge test on path: 8p^3+27q^2 =", E2.subs(sub_path), " (never 0 -> path avoids A2 edge)")
print("A3 test on path: p =", sub_path[p], " != 0 -> path never hits swallowtail vertex")
print("node-curve crossing: q=0 always; p^2-4r =", sp.simplify((p**2 - 4*r).subs(sub_path)),
      " -> vanishes only at m=-1")
print("A1-sheet crossing: r =", sub_path[r], " -> vanishes only at m=+1")
# smoothness of swallowtail at the m=+1 crossing point (-1,0,0):
gr = [g.subs({p: -1, q: 0, r: 0}) for g in gradD]
print("grad Delta at (-1,0,0) =", gr, " (nonzero -> smooth A1 sheet, transversal crossing)")

print()
print("="*78)
print("PART 4: exact mod-squares statement")
print("="*78)
# already proved: disc(F) == disc(R) exactly (unit factor 1) for the standard resolvent,
# and disc(f)==disc(g1)==disc(g2)==disc(R2). Normalization robustness:
G4 = sp.expand(g2.subs(y, 4*w)/64)   # Ferrari-style rescale y=4w
print("g2(4w)/64 =", G4)
check("disc(g2(4w)/64) * 4^6 == disc(g2)  [square factor (4^3)^2]",
      sp.expand(sp.discriminant(G4, w))*4096 - Dg2)
# generic affine root change y -> lam*y (lam a unit): disc scales by lam^6 = (lam^3)^2
glam = sp.expand(g2.subs(y, lam*y)/lam**3)
check("disc(g2(lam y)/lam^3) * lam^6 == disc(g2)",
      sp.expand(sp.expand(sp.discriminant(glam, y))*lam**6 - Dg2))

print()
print("FAILED CHECKS:", FAIL if FAIL else "NONE")
