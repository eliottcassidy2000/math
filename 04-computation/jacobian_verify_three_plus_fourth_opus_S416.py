"""opus-2026-07-20-S416: verify the owner's three JC counterexamples exactly +
mint a 4th in the original class. Dets via chain rule on verified factors
(F's det = -2 re-proved symbolically; shears/translations det 1; diag/perm exact)
with E1 and the 4th ALSO verified by direct symbolic determinant; collisions via
exact Fraction pipelines."""
import sympy as sp
from fractions import Fraction
x,y,z = sp.symbols('x y z')
u = 1+x*y
Fs = [u**3*z + y**2*u*(4+3*x*y), y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y), 2*x-3*x**2*y-x**3*z]
detF = sp.expand(sp.Matrix(Fs).jacobian([x,y,z]).det())
print("F: det =", detF)
def Fnum(p):
    X,Y,Z = p; U = 1+X*Y
    return (U**3*Z + Y**2*U*(4+3*X*Y), Y + 3*X*U**2*Z + 3*X*Y**2*(4+3*X*Y), 2*X-3*X**2*Y-X**3*Z)
P = [ (Fraction(0),Fraction(0),Fraction(-1,4)), (Fraction(1),Fraction(-3,2),Fraction(13,2)), (Fraction(-1),Fraction(3,2),Fraction(13,2)) ]
# E1 = (F2, F3, 2*F1) o (x, y+z, z)
E1 = [f.subs({y: y+z}, simultaneous=True) for f in Fs]
E1 = [sp.expand(e) for e in (E1[1], E1[2], 2*E1[0])]
d1 = sp.expand(sp.Matrix(E1).jacobian([x,y,z]).det())
p1 = [(Fraction(0),Fraction(1,4),Fraction(-1,4)),(Fraction(1),Fraction(-8),Fraction(13,2)),(Fraction(-1),Fraction(-5),Fraction(13,2))]
im1 = [tuple(sp.nsimplify(e.subs({x:pp[0],y:pp[1],z:pp[2]})) for e in E1) for pp in p1]
print("E1: det =", d1, "; images:", im1, "; all == (0,0,-1/2):", all(t==(0,0,sp.Rational(-1,2)) for t in im1))
# E2 = F o (X,Y,Z+X) o F  -- dets by chain rule; collisions by pipeline
def E2num(p):
    X,Y,Z = Fnum(p); return Fnum((X,Y,Z+X))
im2 = [E2num(pp) for pp in P]
print("E2: det = (-2)*1*(-2) = 4 [chain rule on verified factors]; images:", im2,
      "; all == (-1/4,3/16,-129/256):", all(t==(Fraction(-1,4),Fraction(3,16),Fraction(-129,256)) for t in im2))
def E3num(p):
    X = Fnum(p); U = Fnum(X); return Fnum((U[0]+1, U[1], U[2]))
im3 = [E3num(pp) for pp in P]
print("E3: det = (-2)^3 * 1 = -8 [chain rule]; images:", im3,
      "; all == (-1/2,-3/2,5/2):", all(t==(Fraction(-1,2),Fraction(-3,2),Fraction(5,2)) for t in im3))
# 4th (original class): G4 = diag(3,1,1) o F o diag(-1,1,-1); lambda = 1/2 slice
G4 = [f.subs({x:-x, z:-z}, simultaneous=True) for f in Fs]
G4 = [sp.expand(g) for g in (3*G4[0], G4[1], G4[2])]
d4 = sp.expand(sp.Matrix(G4).jacobian([x,y,z]).det())
p4 = [(Fraction(0),Fraction(0),Fraction(1)),(Fraction(-1,2),Fraction(-3),Fraction(-26)),(Fraction(1,2),Fraction(3),Fraction(-26))]
im4 = [tuple(sp.nsimplify(g.subs({x:pp[0],y:pp[1],z:pp[2]})) for g in G4) for pp in p4]
print("G4: det =", d4, "; images:", im4, "; all equal:", len(set(im4))==1, im4[0])
print("G4 factored:", [sp.factor(g) for g in G4])
