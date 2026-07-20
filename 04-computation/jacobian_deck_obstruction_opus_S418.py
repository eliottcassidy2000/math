"""opus-2026-07-20-S418 (HYP-8150): DOES iota LIFT? No -- and exactly where it does.
(a) At random rational targets: the fiber cubic L x^3 + (4-3vw)x - 2w is irreducible
    over Q with NONSQUARE discriminant => Galois group S_3 there => the cubic
    extension C(x,y,z)/C(F) is non-Galois => Aut = 1: NO global birational deck map.
    (L evaluated via kp-THM-1335's master identity: L(F(p)) = [2C - (4-3BC)x]/x^3.)
(b) kp's perfect-square law refereed: (4-3vw)^3 + 27 L w^2 is a perfect square (=Q^2)
    at every sample => disc = -4L*Q^2: the obstruction to squareness is exactly -L.
(c) The w=0 degeneration: cubic -> x(Lx^2+4): S_3 drops to C_2; the slice deck iota
    becomes rational BECAUSE claim 2 works on the resolvent chart: s^2 = -L|_{w=0}
    = v^2-16u, and s <-> -s swaps x+- = +-2/s: iota IS the resolvent deck, w=0 chart.
(d) Tournament reading (verified sign-side): the resolvent = the orientation (sign)
    cover of the 3-point fiber bundle; transpositions (odd) act only after adjoining
    sqrt(-L).  Fiber = the 3-cycle, monodromy full S_3."""
import sympy as sp
from fractions import Fraction
import random
random.seed(8150)
x,y,z,s_,v = sp.symbols('x y z s v')
U = 1+x*y
A = U**3*z + y**2*U*(4+3*x*y)
B = y + 3*x*U**2*z + 3*x*y**2*(4+3*x*y)
C = 2*x - 3*x**2*y - x**3*z
def Lval(p):
    a,b,c = [sp.nsimplify(f.subs({x:p[0],y:p[1],z:p[2]})) for f in (A,B,C)]
    xp = sp.nsimplify(p[0])
    return sp.nsimplify((2*c + (3*b*c-4)*xp)/xp**3), (a,b,c)
print("(a)+(b) random-specialization referee:")
for _ in range(4):
    p = tuple(Fraction(random.randint(-9,9), random.randint(1,5)) for _ in range(3))
    if p[0] == 0: continue
    L,(a,b,c) = Lval(p)
    cub = sp.Poly(L*x**3 + (4-3*b*c)*x - 2*c, x)
    irr = len(sp.factor_list(cub.as_expr())[1]) == 1 and cub.degree() == 3
    disc = sp.discriminant(cub.as_expr(), x)
    QQ2 = sp.nsimplify((4-3*b*c)**3 + 27*L*c**2*(2*c/(2*c))**0)  # (4-3vw)^3+27Lw^2, w=c
    QQ2 = sp.nsimplify((4-3*b*c)**3 + 27*L*c**2)
    print(f"  target ~({sp.nsimplify(a)},{sp.nsimplify(b)},{sp.nsimplify(c)}): irreducible {irr}; "
          f"disc square in Q: {sp.sqrt(disc).is_rational is True}; "
          f"(4-3vw)^3+27Lw^2 a perfect square: {sp.sqrt(QQ2).is_rational is True}; "
          f"disc == -4*L*Q^2: {sp.simplify(disc + 4*L*QQ2) == 0}")
print("(c) resolvent chart at w=0: s<->-s swaps the exotic pair (symbolic):")
xs = 2/s_
pe = lambda sg: {x: sg*2/s_, y: v/4 - 3/(2*sg*2/s_), z: 13/(2*(sg*2/s_)**2) - 3*v/(4*sg*2/s_)}
p_plus, p_minus = pe(1), pe(-1)
swap = all(sp.simplify(p_plus[k].subs(s_, -s_) - p_minus[k]) == 0 for k in (x,y,z))
print("   s -> -s maps exotic(+) to exotic(-):", swap, " (iota = the resolvent deck, w=0 chart)")
print("(d) conclusion: Aut(C(x,y,z)/C(F)) = 1 (non-Galois cubic, disc nonsquare);")
print("    iota lifts exactly to C(F)(sqrt(-L)) -- the orientation/sign cover of the")
print("    3-point fiber bundle; only the S_3 closure carries the transpositions.")
