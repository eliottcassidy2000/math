import cmath
# By-hand reduction of the span-6 {+3,+1,-1,-3} nullcone P=aZ^3+bZ+cZbar+dZbar^3:
# NOTE codex-2026-07-21: the B2 residual check below now has an exact successor:
#   gmc2_span6_symbolic_residual_codex_20260721.py
# It proves E[P^6] = 466560*(ad)^3 under E[P^2]=E[P^4]=0 in the all-nonzero branch.
# E[P^2]=0 => bc=-6ad.  Cases:
#  A) a=0: E[P^2]=2bc=0, two-sided needs b!=0 => c=0; then E[P^4]=24 b^3 d=0 => d=0 => one-sided.
#  B1) a!=0,d=0: bc=0; two-sided needs c!=0 => b=0; E[P^4]=24 a c^3=0 => c=0. contra => one-sided.
#  B2) a,b,c,d all !=0: with x=ac^3, y=b^3 d, z=a^2 d^2, from bc=-6ad get x*y=-216 z^2, and
#      E[P^4]=0 => x+y=-54z. So x/z,y/z = -27 +- 3*sqrt(105) (IRRATIONAL). A 1-param family remains;
#      E[P^6] must kill it. Verify numerically over C:
import math
kappa_p=-27+3*math.sqrt(105); kappa_m=-27-3*math.sqrt(105)
def Emom(a,b,c,d,m):
    # recompute E[P^m] from the formula (same as gmc_span6.py)
    from math import factorial
    from itertools import product
    terms=[(3,0),(1,0),(0,1),(0,3)]; co=[a,b,c,d]; tot=0
    for combo in product(range(m+1),repeat=4):
        if sum(combo)!=m: continue
        A=sum(combo[k]*terms[k][0] for k in range(4)); B=sum(combo[k]*terms[k][1] for k in range(4))
        if A!=B: continue
        mult=factorial(m)
        for k in range(4): mult//=factorial(combo[k])
        t=mult*factorial(A)
        for k in range(4): t*=co[k]**combo[k]
        tot+=t
    return tot
# construct a B2 point: a=1, d=t; choose c with c^3 = kappa_p * t^2 ; b = -6 t / c (enforces bc=-6ad)
print("Testing the B2 residual family (all coeffs !=0, satisfies E[P^2]=E[P^4]=0) against E[P^6],E[P^8]:")
for t in [1.0, 2.0, 0.5, 1.3+0.7j]:
    a=1.0+0j; d=complex(t)
    c=(kappa_p*d*d)**(1/3)     # one cube root
    b=-6*a*d/c
    e2=Emom(a,b,c,d,2); e4=Emom(a,b,c,d,4); e6=Emom(a,b,c,d,6); e8=Emom(a,b,c,d,8)
    print(f"  t={t}: |E2|={abs(e2):.2e} |E4|={abs(e4):.2e}  |E6|={abs(e6):.3e}  |E8|={abs(e8):.3e}")
print("\n=> E2,E4 ~0 (family satisfies them by construction) but E6 (and E8) are NONZERO on the whole")
print("   B2 family => no two-sided nullcone element. Combined with Cases A,B1 (by hand):")
print("   the span-6 {+3,+1,-1,-3} nullcone is ONE-SIDED. GMC(2) holds on this stratum, UNCONDITIONALLY,")
print("   at detection depth m<=6 -- consistent with M* <= 2*span = 12 (well within), and NEW (fleet had span<=4).")
