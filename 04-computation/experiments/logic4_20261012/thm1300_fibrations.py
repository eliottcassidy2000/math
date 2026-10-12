"""THM-1300's Keller map: the three coordinate-pair fibrations, exactly.
Session mac-mini-2026-10-12-logic4.  All checks are exact (sympy rational functions).
"""
import sympy as sp
x,y,z,a,b,c = sp.symbols('x y z a b c')
X,s,r1,r2 = sp.symbols('X s r1 r2')
u = 1+x*y
F1 = u**3*z + y**2*u*(4+3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
F = [F1,F2,F3]
J = sp.Matrix([[sp.diff(f,v) for v in (x,y,z)] for f in F])
print("det JF =", sp.expand(J.det()))
print("deg F1,F2,F3 =", [sp.Poly(f,x,y,z).total_degree() for f in F])
print("identity (F2-y)(1+xy) - 3x F1 == 0:", sp.expand((F2-y)*u-3*x*F1)==0)

print("\n=== fibration G12 = (F1,F2), fibre over (a,b), a != 0 ===")
Q = y**2 - b*y + 3*a
xp = (b-y)/Q
up = sp.cancel(1+xp*y); print("u on fibre =", up)
zp = sp.cancel((a - y**2*up*(4+3*xp*y))/up**3); print("z on fibre =", sp.factor(zp))
sub12 = {x:xp, z:zp}
print("F1,F2 on fibre:", sp.cancel(F1.subs(sub12)), sp.cancel(F2.subs(sub12)))
F3p = sp.cancel(F3.subs(sub12)); print("F3 on fibre =", sp.factor(F3p), " (polynomial in y? denominator:", sp.denom(sp.together(F3p)),")")
dF3 = sp.factor(sp.diff(F3p,y)); print("dF3/dy =", dF3)
kappa = sp.cancel(dF3/Q); print("kappa = dF3/dy / Q =", kappa)
c0 = sp.expand(F3p - kappa*(y**3/3 - b*y**2/2 + 3*a*y)); print("F3 - kappa*(y^3/3 - b y^2/2 + 3 a y) =", c0)
# no point of the fibre has Q(y)=0: (b-y)(1+xy)=3ax at Q(y)=0 forces b=y, then Q(b)=3a != 0
print("Q(b) =", sp.expand(Q.subs(y,b)))
# V3 on the fibre: V3 = (1/F3') d/dy ; cokernel witnesses via residues with Q=(y-r1)(y-r2)
subr = {a: r1*r2/3, b: r1+r2}
def res_class(q):
    """class of q in coker(V3|fibre): residues of kappa*Q*q dy at r1, r2 (q given as rational function of y on the fibre)."""
    w = sp.cancel((kappa*Q*q).subs(subr))
    return [sp.simplify(sp.residue(w, y, r)) for r in (r1,r2)]
for name,q in [("1",sp.Integer(1)),("x",xp),("y",y),("z",zp),("x^2",xp**2),("x*y",xp*y),("x*z",xp*zp),("y^2",y**2),("x^2*y",xp**2*y),("x^3",xp**3)]:
    print("  residues of kappa*Q*(%s) dy at r1,r2:"%name, res_class(q))

print("\n=== fibration G13 = (F1,F3), fibre over (a,c), a,c != 0 ===")
v = X/(a*X**3 - X + c)      # v = 1+xy
w = v-1
xq = sp.cancel(v*X); yq = sp.cancel(w/xq)
zq = sp.cancel((2*xq - 3*xq**2*yq - c)/xq**3)
sub13 = {x:xq, y:yq, z:zq}
print("F1,F3 on fibre:", sp.cancel(F1.subs(sub13)), sp.cancel(F3.subs(sub13)))
F2q = sp.cancel(F2.subs(sub13)); print("F2 on fibre =", F2q, " = 2aX + 2/X - c/X^2 ?", sp.cancel(F2q-(2*a*X+2/X-c/X**2))==0)
print("dF2/dX =", sp.factor(sp.diff(F2q,X)))
print("X recovered as x/(1+xy):", sp.cancel(xq/(1+xq*yq)))

print("\n=== fibration G23 = (F2,F3), fibre over (b,c), c != 0 ===")
# on F3=c: z=(2x-3x^2y-c)/x^3 ; with w=xy: F2 = [w x + 3(1+w)^2 (x(2-3w)-c) + 3w^2(4+3w)x]/x^2 ; bracket = 6+4w
wv = sp.symbols('w')
br = sp.expand(wv + 3*(1+wv)**2*(2-3*wv) + 3*wv**2*(4+3*wv)); print("bracket =", br)
# conic: b x^2 - (6+4w) x + 3c(1+w)^2 = 0, through (x,w)=(0,-1); lines w = -1 + s x
xs = 2/(3*c*s**2 - 4*s + b); ws = -1 + s*xs
ys = sp.cancel(ws/xs); zs = sp.cancel((2*xs - 3*xs**2*ys - c)/xs**3)
sub23 = {x:xs, y:ys, z:zs}
print("F2,F3 on fibre:", sp.cancel(F2.subs(sub23)), sp.cancel(F3.subs(sub23)))
F1s = sp.cancel(F1.subs(sub23)); print("F1 on fibre =", sp.factor(F1s))
print("dF1/ds =", sp.factor(sp.diff(F1s,s)))
print("u = 1+xy on fibre =", sp.factor(sp.cancel(1+xs*ys)))

print("\n=== generic degree and Galois group via the (F1,F2)-fibre cubic ===")
cub = sp.Poly(sp.expand(F3p - c), y)
print("cubic in y:", cub.as_expr())
disc = sp.factor(sp.discriminant(cub.as_expr(), y)); print("discriminant =", disc)
num = sp.numer(sp.together(sp.discriminant(cub.as_expr(),y)))
fl = sp.factor_list(num); print("discriminant numerator factor list:", fl)
print("disc is a perfect square in Q(a,b,c)?", all(e%2==0 for _,e in fl[1]), "=> Galois group of K(x,y,z)/K(F) is", "A3" if all(e%2==0 for _,e in fl[1]) else "S3")
c3 = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2
print("Jelonek polynomial equals Gao's c_3(v):", sp.expand(num + 4*c3)==0 or sp.expand(num - 4*c3)==0, " (num =", sp.expand(num),")")
# image complement: the curve C = {(4/(27t^2), 4/(3t), t)}: double root r = 2/(3t) of Q and F3(r) = t
t = sp.symbols('t')
subC = {a: 4/(27*t**2), b: 4/(3*t)}
print("on C: Q =", sp.factor(Q.subs(subC)), "; F3(r) - t =", sp.simplify(F3p.subs(subC).subs(y, 2/(3*t)) - t))
# fibre over (a,b,c) is {roots of cubic} minus {roots of Q}; over S_F \ C exactly one point, over C empty
# witnesses for V_2 (fibration G13): q in im V2 iff residues of q * F2'(X) dX vanish at X=0 and at the three roots of aX^3-X+c
print("\n=== cokernel witnesses for V_2 on the (F1,F3)-fibre (residues of q*dF2 at X=0 and roots of aX^3-X+c) ===")
rho = sp.symbols('rho')  # a root: a rho^3 - rho + c = 0  => c = rho - a rho^3
dF2 = sp.diff(F2q, X)
for name,q in [("1",sp.Integer(1)),("x",xq),("y",yq),("z",zq),("x^2",xq**2),("x*y",xq*yq),("y^2",yq**2),("x*z",xq*zq),("y*z",yq*zq),("z^2",zq**2),("x^2*y",xq**2*yq),("x*y^2",xq*yq**2)]:
    w = sp.cancel(q*dF2)
    r0 = sp.simplify(sp.residue(w, X, 0))
    # residue at a simple root rho of aX^3-X+c: substitute c = rho - a rho^3, then residue at X=rho
    wr = sp.cancel(w.subs(c, rho - a*rho**3))
    rr = sp.simplify(sp.residue(wr, X, rho))
    print("  q=%-6s res_0 = %s ; res_rho = %s"%(name, r0, sp.factor(rr)))
print("\n=== cokernel witnesses for V_1 on the (F2,F3)-fibre (residues of q*dF1 at the roots of 3cs^2-4s+b) ===")
sig = sp.symbols('sigma')  # root: 3 c sig^2 - 4 sig + b = 0 => b = 4 sig - 3 c sig^2
dF1 = sp.diff(F1s, s)
for name,q in [("1",sp.Integer(1)),("x",xs),("y",ys),("z",zs),("x^2",xs**2),("x*y",xs*ys),("y^2",ys**2),("x*z",xs*zs),("y*z",ys*zs),("z^2",zs**2),("x^3",xs**3),("x^2*y",xs**2*ys)]:
    w = sp.cancel((q*dF1).subs(b, 4*sig - 3*c*sig**2))
    rr = sp.simplify(sp.residue(w, s, sig))
    print("  q=%-6s res_sigma = %s"%(name, sp.factor(rr)))

print("\n=== dual fields and non-nilpotency ===")
B = (J.T).inv()
B = B.applyfunc(sp.cancel)
V = [lambda h, j=j: sp.expand(sum(B[j,k]*sp.diff(h,v) for k,v in enumerate((x,y,z)))) for j in range(3)]
print("V_j F_i =", [[sp.simplify(V[j](F[i])) for i in range(3)] for j in range(3)])
print("coefficients of V_j (polynomial?):", [[sp.denom(sp.together(B[j,k])) for k in range(3)] for j in range(3)])
for j in range(3):
    for name,h0 in [("x",x),("y",y),("z",z)]:
        h=h0; degs=[]
        for m in range(1,7):
            h=V[j](h); degs.append(sp.Poly(h,x,y,z).total_degree() if h!=0 else -1)
        print("  deg V_%d^m(%s), m=1..6:"%(j+1,name), degs)
