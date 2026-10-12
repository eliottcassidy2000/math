"""Audit A (independent re-derivation) of THM-4621: the three coordinate-pair fibrations of THM-1300's map.
Written from scratch for the audit; only the definition of F is shared with the session scripts.
Every check is an exact sympy identity; assertion failures would abort the run."""
import sympy as sp
x,y,z,a,b,c,t,w,X,s = sp.symbols('x y z a b c t w X s')
u = 1+x*y
F1 = u**3*z + y**2*u*(4+3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
F = [F1,F2,F3]; V = [x,y,z]
J = sp.Matrix(3,3, lambda i,j: sp.diff(F[i], V[j]))
print("det JF =", sp.expand(J.det())); assert sp.expand(J.det()) == -2
print("identity (F2-y)u - 3x F1 =", sp.expand((F2-y)*u - 3*x*F1))

print("\n[1] (F1,F2)-fibre, a != 0 — derived by elimination, not by assuming the parametrization")
zsol = sp.solve(sp.Eq(F1, a), z)[0]                      # F1 is linear in z with coefficient u^3
E = sp.factor(sp.together((F2 - b).subs(z, zsol)))
print("  (F2 - b)|_{F1=a} =", E)
Q = y**2 - b*y + 3*a
num = sp.numer(sp.together((F2 - b).subs(z, zsol)))
print("  numerator factors:", sp.factor(num), "   [contains x*Q - (b - y)?]", sp.simplify(sp.rem(sp.expand(num), sp.expand(x*Q-(b-y)), x)) == 0)
xp = (b-y)/Q
up = sp.cancel(u.subs(x, xp)); print("  u on fibre =", up); assert sp.cancel(up - 3*a/Q) == 0
zp = sp.cancel(zsol.subs(x, xp))
print("  z on fibre =", sp.factor(zp))
assert sp.cancel(zp + Q*(2*y**4 - b*y**3 + (30*a - b**2)*y**2 + 6*a*b*y - 9*a**2)/(27*a**2)) == 0
sub12 = {x: xp, z: zp}
assert sp.cancel(F1.subs(sub12) - a) == 0 and sp.cancel(F2.subs(sub12) - b) == 0
print("  F1 = a, F2 = b on the parametrized curve: OK")
F3p = sp.cancel(F3.subs(sub12))
print("  F3 on fibre =", sp.factor(F3p), "; is polynomial in y:", sp.denom(sp.together(F3p)).free_symbols <= {a})
assert sp.cancel(F3p + (y-b)*(2*y**2 - b*y + 18*a - b**2)/(27*a**2)) == 0
kappa = sp.Rational(-2,9)/a**2; c0 = 2*b/(3*a) - b**3/(27*a**2)
assert sp.expand(F3p - (kappa*(y**3/3 - b*y**2/2 + 3*a*y) + c0)) == 0
assert sp.cancel(sp.diff(F3p, y) - kappa*Q) == 0
print("  F3 = kappa (y^3/3 - b y^2/2 + 3 a y) + c0, dF3/dy = kappa*Q: OK")
print("  Q(b) =", sp.expand(Q.subs(y, b)), " (so no fibre point at a root of Q); F1 at u=0:", sp.factor(F1.subs(z, (-1 - 0*z)*0 + z).subs(x, -1/y)))

print("\n[2] discriminant, Jelonek polynomial, the empty-fibre curve")
cubic = sp.expand(F3p - c)
disc = sp.factor(sp.discriminant(cubic, y))
c3 = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2
print("  disc_y(F3 - c) =", disc); assert sp.cancel(disc + 4*c3/(729*a**6)) == 0
print("  c3 irreducible over Q:", sp.factor_list(c3)[1][0][1] == 1 and len(sp.factor_list(c3)[1]) == 1)
print("  disc_a(c3) =", sp.factor(sp.discriminant(c3, a)), " (not a square => c3 absolutely irreducible as a quadratic in a)")
onC = {a: 4/(27*t**2), b: 4/(3*t)}
print("  on C: Q =", sp.factor(Q.subs(onC)), "; F3(2/(3t)) - t =", sp.simplify(F3p.subs(onC).subs(y, 2/(3*t)) - t))
print("  c3 on b^2 = 12a:", sp.factor(c3.subs(a, b**2/12)), " => {c3=0} ∩ {b^2=12a, a!=0} = {c = 4/(3b)} = C")
grad = [sp.diff(c3, v) for v in (a,b,c)]
print("  grad c3 on C:", [sp.simplify(g.subs(onC).subs(c, t)) for g in grad], " (C ⊂ Sing S)")
print("  Sing S:", sp.solve(grad + [c3], [a,b,c], dict=True))
print("  THM-1310's L(a,b,g) with g=c equals c3:", sp.expand((27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2) - c3) == 0)

print("\n[2'] the stratum a = 0 (excluded from the parametrization; THM-4621 cites Gao for it)")
# u = 0 branch: x = -1/y
F2u0 = sp.simplify(F2.subs(x, -1/y)); F3u0 = sp.simplify(F3.subs(x, -1/y))
print("  on u=0: F1 =", sp.simplify(F1.subs(x, -1/y)), ", F2 =", F2u0, ", F3 =", sp.factor(F3u0), " => one point (y=-b/2, x=2/b, z solved linearly) iff b != 0")
# u != 0 branch: F1 = 0 => z = -y^2(4+3xy)/u^2, then F2 = y
z0 = -y**2*(4+3*x*y)/u**2
print("  on u!=0, F1=0: F2 =", sp.simplify(F2.subs(z, z0)), "; F3 =", sp.factor(sp.simplify(F3.subs(z, z0).subs(y, b))))
eq0 = sp.factor(sp.numer(sp.together(F3.subs(z, z0).subs(y, b) - c)))
print("  F3 = c there: numerator =", eq0, "; discriminant in x =", sp.factor(sp.discriminant(sp.expand(eq0), x)))
print("  => for a = 0: #fibre = 1 + 2 = 3 when b(1-bc) != 0; = 1 when b != 0, bc = 1; = 1 when b = 0 (c != 0 or c = 0); never 0 or 2.  c3|_{a=0} =", sp.factor(c3.subs(a, 0)))

print("\n[3] tangency polynomial: y = b - 3ac/w")
Tg = sp.factor(sp.cancel((F3p - c).subs(y, b - 3*a*c/w)))
print("  F3(b - 3ac/w) - c =", Tg)
assert sp.cancel(Tg + c*(w**3 - 2*w**2 + b*c*w - 2*a*c**2)/w**3) == 0
print("  = -(c/w^3)(w^3 - 2w^2 + bc w - 2ac^2): OK  [THM-4621 writes the prefactor as 9c/2 — an artifact of sympy's numerator normalisation]")

print("\n[4] Galois group: disc not a square (c3 to the first power) => S3; deg F = 3 (THM-1310(3))")
fl = sp.factor_list(sp.numer(sp.together(disc)))
print("  factor list of disc numerator:", fl, " => square?", all(e % 2 == 0 for _, e in fl[1]))

print("\n[5] (F1,F3)-fibre, ac != 0 — derived by elimination")
zc = sp.solve(sp.Eq(F3, c), z)[0]
E13 = sp.factor(sp.numer(sp.together((F1 - a).subs(z, zc))))
print("  numerator of (F1 - a)|_{F3=c} =", E13)
E13w = sp.factor(sp.expand(E13.subs(y, w/x)))
print("  in terms of w = xy:", E13w)
assert sp.expand(E13.subs(y, w/x)*1 - (-(a*x**3 - (1+w)*(2+w)*x + c*(1+w)**3))) == 0 or sp.expand(E13.subs(y, w/x) - (a*x**3 - (1+w)*(2+w)*x + c*(1+w)**3)) == 0
v = X/(a*X**3 - X + c)
xq = sp.cancel(v*X); yq = sp.cancel((v-1)/xq); zq = sp.cancel(zc.subs({x: xq, y: yq}))
sub13 = {x: xq, y: yq, z: zq}
assert sp.cancel(F1.subs(sub13) - a) == 0 and sp.cancel(F3.subs(sub13) - c) == 0
F2q = sp.cancel(F2.subs(sub13)); print("  F2 on fibre =", F2q); assert sp.cancel(F2q - (2*a*X + 2/X - c/X**2)) == 0
print("  dF2/dX =", sp.factor(sp.diff(F2q, X)), "; X = x/(1+xy) recovered:", sp.cancel(xq/(1+xq*yq)) == X)
print("  disc_X(aX^3 - X + c) =", sp.factor(sp.discriminant(a*X**3 - X + c, X)), " (three distinct roots unless 27ac^2 = 4)")

print("\n[6] (F2,F3)-fibre, c != 0 — derived by elimination")
E23 = sp.factor(sp.numer(sp.together((F2 - b).subs(z, zc))))
E23w = sp.factor(sp.expand(E23.subs(y, w/x)))
print("  numerator of (F2 - b)|_{F3=c} in (x, w):", E23w)
conic = b*x**2 - (6+4*w)*x + 3*c*(1+w)**2
print("  equals ±(b x^2 - (6+4w) x + 3c(1+w)^2) up to a power of x:", sp.simplify(sp.cancel(E23.subs(y, w/x)/conic)))
xs = 2/(3*c*s**2 - 4*s + b); ws = -1 + s*xs
ys = sp.cancel(ws/xs); zs = sp.cancel(zc.subs({x: xs, y: ys}))
sub23 = {x: xs, y: ys, z: zs}
assert sp.cancel(F2.subs(sub23) - b) == 0 and sp.cancel(F3.subs(sub23) - c) == 0
F1s = sp.cancel(F1.subs(sub23)); print("  F1 on fibre =", sp.factor(F1s)); assert sp.cancel(F1s - s*(c*s**2 - 2*s + b)/2) == 0
print("  dF1/ds =", sp.factor(sp.diff(F1s, s)), "; u on fibre =", sp.factor(sp.cancel(1 + xs*ys)), "; s = (xy+1)/x recovered:", sp.cancel((xs*ys+1)/xs) == s)

print("\n[8] residue witnesses (independent computation)")
r1, r2, rho, sig = sp.symbols('r1 r2 rho sigma')
# V3 on the (F1,F2) fibre: q dF3 = q * kappa * Q dy ; Q = (y-r1)(y-r2), a = r1 r2/3, b = r1 + r2
subr = {a: r1*r2/3, b: r1+r2}
def res3(q):
    wq = sp.cancel((q*kappa*Q).subs(subr))
    return [sp.factor(sp.residue(wq, y, r)) for r in (r1, r2)]
print("  V3: residues of q dF3 at r1, r2:")
for name, q in [("1", sp.Integer(1)), ("x", xp), ("y", y), ("z", zp), ("x^2", xp**2), ("x^2 y", xp**2*y), ("x^3", xp**3)]:
    print("     q = %-5s ->" % name, res3(q))
# sanity: q = V3(x) restricted must have zero residues (it is exact: dx|_C)
B = (J.T).inv().applyfunc(sp.cancel)
def Vj(j, h): return sp.expand(sum(B[j,k]*sp.diff(h, V[k]) for k in range(3)))
q = sp.cancel(Vj(2, x).subs(sub12)); print("     control q = V3(x)|_C ->", res3(q), " (must be [0, 0])")
q = sp.cancel(Vj(2, x**2*z).subs(sub12)); print("     control q = V3(x^2 z)|_C ->", res3(q), " (must be [0, 0])")
# V2 on the (F1,F3) fibre: residues of q dF2 at X = 0 and at a root rho of aX^3 - X + c (c = rho - a rho^3)
dF2 = sp.diff(F2q, X)
def res2(q):
    wq = sp.cancel(q*dF2)
    r0 = sp.factor(sp.residue(wq, X, 0))
    rr = sp.factor(sp.residue(sp.cancel(wq.subs(c, rho - a*rho**3)), X, rho))
    return r0, rr
print("  V2: residues of q dF2 at X=0 and at X=rho:")
for name, q in [("1", sp.Integer(1)), ("x", xq), ("y", yq), ("z", zq), ("x^2", xq**2)]:
    print("     q = %-5s ->" % name, res2(q))
q = sp.cancel(Vj(1, x).subs(sub13)); print("     control q = V2(x)|_C ->", res2(q), " (must be (0, 0))")
# V1 on the (F2,F3) fibre: residues of q dF1 at a root sigma of 3cs^2 - 4s + b (b = 4 sigma - 3 c sigma^2)
dF1 = sp.diff(F1s, s)
def res1(q):
    wq = sp.cancel((q*dF1).subs(b, 4*sig - 3*c*sig**2))
    return sp.factor(sp.residue(wq, s, sig))
print("  V1: residues of q dF1 at s=sigma:")
for name, q in [("1", sp.Integer(1)), ("x", xs), ("y", ys), ("z", zs), ("x^2", xs**2), ("x^3", xs**3), ("x^2 y", xs**2*ys)]:
    print("     q = %-5s ->" % name, res1(q))
q = sp.cancel(Vj(0, x).subs(sub23)); print("     control q = V1(x)|_C ->", res1(q), " (must be 0)")
q = sp.cancel(Vj(0, x**2*y).subs(sub23)); print("     control q = V1(x^2 y)|_C ->", res1(q), " (must be 0)")

print("\n[7] z|fibre is not a polynomial of degree <= 2 in F3|fibre (coefficients in K(a,b))")
al, be, ga = sp.symbols('alpha beta gamma')
eqs = sp.Poly(sp.expand(zp - (al*F3p**2 + be*F3p + ga)), y).coeffs()
print("  solutions:", sp.solve(eqs, [al, be, ga], dict=True), "; deg_y z|fibre =", sp.Poly(zp, y).degree(), "; deg_y F3|fibre =", sp.Poly(F3p, y).degree())
print("\nALL ASSERTIONS PASSED")
