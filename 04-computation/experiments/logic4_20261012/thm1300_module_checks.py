"""Hostile controls for THM-4621: (1) Ker V_3 on polynomials of degree <= D equals K[F1,F2] in that range (exact nullspace, flint);
(2) the (F1,F2)-fibre cubic is Gao's tangency polynomial under w = 3ac/(b-y);
(3) x, y, z restricted to the generic (F1,F2)-fibre are not polynomials in F3|fibre (so x,y,z lie outside K[F] = the locally nilpotent part)."""
import sympy as sp, itertools, time
from flint import fmpq_mat, fmpq
x,y,z,a,b,c,w = sp.symbols('x y z a b c w')
u = 1+x*y
F1 = u**3*z + y**2*u*(4+3*x*y); F2 = y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y); F3 = 2*x - 3*x**2*y - x**3*z
F=[F1,F2,F3]
J = sp.Matrix([[sp.diff(f,v) for v in (x,y,z)] for f in F]); B=(J.T).inv().applyfunc(sp.cancel)
def V(j,h): return sp.expand(sum(B[j,k]*sp.diff(h,v) for k,v in enumerate((x,y,z))))
# (2) tangency polynomial
Q = y**2 - b*y + 3*a; F3p = -(y-b)*(2*y**2-b*y+18*a-b**2)/(27*a**2)
eq = sp.cancel((F3p - c).subs(y, b - 3*a*c/w))
T = w**3 - 2*w**2 + b*c*w - 2*a*c**2
print("(2) F3(b - 3ac/w) - c =", sp.factor(eq))
print("    equals -(c/w^3) * (w^3 - 2w^2 + bc w - 2ac^2) [Gao's tangency polynomial at (X,Y)=(bc, ac^2)]:", sp.cancel(eq + (c/w**3)*T)==0)
# (3) x,y,z on the fibre vs polynomials in F3|fibre
xp=(b-y)/Q; up=3*a/Q; zp=sp.cancel((a - y**2*up*(4+3*xp*y))/up**3)
# z = alpha F3^2 + beta F3 + gamma ?  solve for alpha,beta,gamma in K(a,b)
al,be,ga = sp.symbols('alpha beta gamma')
diff = sp.Poly(sp.expand(zp - (al*F3p**2 + be*F3p + ga)), y)
sol = sp.solve(diff.coeffs(), [al,be,ga], dict=True)
print("(3) z|fibre as a quadratic polynomial in F3|fibre: solutions =", sol)
print("    deg_y of z|fibre =", sp.Poly(zp,y).degree(), "; x|fibre has poles at the roots of Q (not a polynomial in y):", sp.denom(sp.together(xp)))
# (1) kernel of V_3 on polynomials of degree <= D
D=12
mons=[(i,j,k) for i in range(D+1) for j in range(D+1) for k in range(D+1) if i+j+k<=D]
imgs=[sp.Poly(V(2, x**i*y**j*z**k), x,y,z) for (i,j,k) in mons]
allm=sorted(set(m for P in imgs for m in P.monoms()))
midx={m:i for i,m in enumerate(allm)}
M=fmpq_mat(len(allm), len(mons))
for col,P in enumerate(imgs):
    for m,cf in zip(P.monoms(),P.coeffs()):
        M[midx[m],col]=fmpq(int(cf.p),int(cf.q))
t=time.time(); rank=M.rank(); print("(1) D=%d: #monomials=%d, image rank=%d, dim Ker V_3 =%d (%.1fs)"%(D,len(mons),rank,len(mons)-rank,time.time()-t))
exp=[(i,j) for i in range(3) for j in range(3) if 7*i+6*j<=D]
print("    monomials F1^i F2^j of degree <= D:", exp, "count =", len(exp))
# check these are in the kernel and independent
polys=[sp.Poly(sp.expand(F1**i*F2**j),x,y,z) for (i,j) in exp]
print("    V_3(F1^i F2^j) all zero:", all(V(2,F1**i*F2**j)==0 for (i,j) in exp))
