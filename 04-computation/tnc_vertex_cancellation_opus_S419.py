import sympy as sp
u,r0,r3,r6=sp.symbols('u r0 r3 r6')
def CT(Rexpr,N,mm):
    return sp.Poly(sp.expand(Rexpr**mm),u).coeff_monomial(u**(N*mm))
print("THE FIRST VERTEX-CANCELLATION TRINOMIAL.  R = r0 + r3 u^3 + r6 u^6, N=2,")
print("charges {-2,1,4}, m0=3 with TWO reps: {-2,-2,4} and {1,1,-2}.")
R=r0+r3*u**3+r6*u**6
c3=sp.expand(CT(R,2,3))
print(f"  CT(3) = {sp.factor(c3)}   -> zero iff r0 r6 + r3^2 = 0")
print("  TUNE r0=1, r3=1, r6=-1 (so r3^2 = 1 = -r0 r6):")
Rt=1+u**3-u**6
cts=[CT(Rt,2,m) for m in range(1,12)]
print(f"  R = 1 + u^3 - u^6, N=2:  CT[1..11] = {cts}")
fm=next((m for m in range(1,12) if cts[m-1]!=0),None)
print(f"  CT(3) = {cts[2]} (VERTEX CANCELLED), first nonzero at m={fm}, value {cts[fm-1] if fm else '--'}")
print()
print("  => leading (m0=3) obstruction CANCELS, but a LATER m survives (m={}). Exactly the".format(fm))
print("     THM-1635 phenomenon: nondegenerate structure saves it at the next order.")
print("     This trinomial is NOT a nullcone element -- TNC holds for it -- but the proof")
print("     needs the all-order/branch-product argument, NOT the unique-minimal shortcut.")
print()
print("  Confirm the dominant-saddle nondegeneracy (THM-1635 route) for R = 1+u^3-u^6:")
import numpy as np
N=2; Rx=1+u**3-u**6; S=sp.expand(u*sp.diff(Rx,u)-N*Rx)
g2=sp.diff(sp.log(Rx)-N*sp.log(u),u,2)
rows=[]
for rt in sp.nroots(sp.Poly(S,u)):
    if abs(complex(Rx.subs(u,rt)))<1e-9 or abs(complex(rt))<1e-9: continue
    w=abs(complex(Rx.subs(u,rt))/complex(rt)**N); gg=abs(complex(g2.subs(u,rt)))
    rows.append((w,gg))
rho=max(w for w,_ in rows)
dom=[(w,gg) for w,gg in rows if w>rho*(1-1e-6)]
print(f"     dominant |w|={rho:.3f}, |g''| at dominant saddles = {[round(gg,3) for _,gg in dom]}")
print(f"     nondegenerate (g''!=0): {all(gg>1e-9 for _,gg in dom)}  -> G singular, CT != 0, TNC holds")
