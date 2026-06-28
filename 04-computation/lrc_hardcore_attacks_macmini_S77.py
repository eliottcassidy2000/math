"""
S77: pushes and pulls on the LRC(14) hard core (the CORE construction, measure-zero tight configs).
PUSH: the covering-tight dilations d*{1..13} have the EXPLICIT cyclotomic witness t=1/(14d) (off the 14-grid,
      apex-7 floor blocks t=a/14 -> witness promoted to the finer Phi_{14d} grid). 14=|D_7|=C_7(Legendre/Phi_7) x Z_2(Eisenstein/Phi_2).
PULL: the dip is NOT a Q(sqrt-7) norm a^2+ab+2b^2 (1081=23*47, 47 inert/non-QR mod 7).
"""
from fractions import Fraction as F
from math import gcd

def M_at(S,t): return min(min((s*t)%1,1-(s*t)%1) for s in S)

print("PUSH 2 -- covering-tight dilation witness t=1/(14d) (Phi_{14d}); apex-7 floor promotes off the 14-grid:")
for d in (2,3,4,5,7,14):
    S=[d*s for s in range(1,14)]; apex=[s for s in S if s%14==0]; t=F(1,14*d); m=M_at(S,t)
    print(f"  d={d:>2}: mult-of-14 in set={apex if apex else 'none'}; t=1/{14*d}; M={m} (=1/14? {m==F(1,14)})")
print("  => covering-tight DILATION case CONSTRUCTED: explicit witness t=1/(14d).")

print("\nPULL 1 -- dip Q(sqrt-7)-norm test (a^2+ab+2b^2, disc -7):")
def norm_disc7(n):
    for b in range(0,int((4*abs(n)/7)**0.5)+2):
        for a in range(-2*abs(n)-1,2*abs(n)+2):
            if a*a+a*b+2*b*b==n: return (a,b)
    return None
for lab,n in [('1081=23*47 (dip_8 num)',1081),('23 (2 mod7, QR, split)',23),('47 (5 mod7, non-QR, inert)',47)]:
    print(f"  {lab}: norm? {norm_disc7(n)}")
print("  => 1081 NOT a norm (47 inert at odd power); the dip is not a clean imaginary-quadratic norm. PULL.")

print("\nPUSH 1 -- 14=|D_7|=|C_7 x| Z_2: C_7=Legendre/de Moivre/Phi_7 (rotation); Z_2=Eisenstein/Phi_2 (reflection).")
print("  cap modes = D_7 irreps (trivial+sign from Z_2; three 2-dim from C_7); Borsuk-Ulam (kps) = topological Brouwer (S75f).")
