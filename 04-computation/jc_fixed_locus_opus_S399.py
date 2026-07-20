from fractions import Fraction as Fr
import sympy as sp
def Fmap(x,y,z):
    u=1+x*y
    return (u**3*z + y**2*u*(4+3*x*y),
            y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y),
            2*x - 3*x**2*y - x**3*z)
x,y,z=sp.symbols('x y z')
F=sp.Matrix(Fmap(x,y,z))
print("(1) sigma-equivariance: F(-x,-y,z) =? tau(F(x,y,z)) with tau(a,b,c)=(a,-b,-c)")
Fs=sp.Matrix(Fmap(-x,-y,z))
tauF=sp.Matrix([F[0],-F[1],-F[2]])
print(f"    F(sigma p) - tau(F p) = {sp.simplify(Fs-tauF).T}")
print()
print("(2) THE FIXED LOCUS.  Fix(sigma) = {x=y=0};  Fix(tau) = {b=c=0}")
Ffix=sp.simplify(sp.Matrix(Fmap(0,0,z)))
print(f"    F(0,0,z) = {Ffix.T}   -> lands in Fix(tau)? {sp.simplify(Ffix[1])==0 and sp.simplify(Ffix[2])==0}")
print(f"    so F|Fix : z |-> {Ffix[0]}  -- a degree-1 map, hence BIJECTIVE")
print()
print("(3) CONSEQUENCE: every fiber over a tau-fixed point has EXACTLY ONE sigma-fixed")
print("    preimage, so |fiber| is ODD -- the Redei-shaped 1+2 (kind-pasteur S128c97).")
print("    The triple collision, checked exactly:")
for p in [(Fr(0),Fr(0),Fr(-1,4)), (Fr(1),Fr(-3,2),Fr(13,2)), (Fr(-1),Fr(3,2),Fr(13,2))]:
    v=Fmap(*p)
    sig=(-p[0],-p[1],p[2])
    print(f"      F{p} = {v}   sigma-fixed: {sig==p}")
print()
print("(4) THE POINT: the non-injectivity is ENTIRELY OFF the fixed locus.")
print("    F restricted to Fix(sigma) is injective; the collision uses the FREE orbit.")
print("    Check: are the two colliding non-fixed points a single sigma-orbit?")
P=(Fr(1),Fr(-3,2),Fr(13,2)); Q=(Fr(-1),Fr(3,2),Fr(13,2))
print(f"      sigma{P} = {(-P[0],-P[1],P[2])}  == Q ? {(-P[0],-P[1],P[2])==Q}")
print()
print("(5) det of the Jacobian, and of the RESTRICTED map")
J=F.jacobian([x,y,z]); print(f"    det JF = {sp.expand(J.det())}")
zz=sp.symbols('z')
print(f"    d/dz of F|Fix = {sp.diff(Ffix[0], z)}  (nonzero constant => JC_1 on the fixed locus)")
