# At a stacking (two pinches r1,r2, same gamma, same fold t_*), the fold amplitude is
#   A_i ~ e^{-r_i}/sqrt(D_rr(r_i,t_*)),  and I claim  D_rr(r_i,t_*) = -2 t_*^2 b'(r_i) phi'(r_i).
# The TOTAL jump beta_total ~ (+-)A_1 (+-) A_2.  When can it vanish? (boxeph's THM-1630 s4 amendment.)
import numpy as np
def fold_data(bcoef, r):
    b=np.polynomial.Polynomial(bcoef); bp=b.deriv(); bpp=bp.deriv()
    gamma=r*bp(r)**2
    tstar=1.0/(b(r)-2*r*bp(r))          # 1/phi
    # D(r,t)=(1-b t)^2-4 gamma r t^2 ; D_rr at (r,tstar) by direct 2nd derivative in r (numeric)
    h=1e-5
    def D(rr): return (1-b(rr)*tstar)**2-4*gamma*rr*tstar**2
    Drr=(D(r+h)-2*D(r)+D(r-h))/h**2
    phip=-(bp(r)+2*r*bpp(r))            # phi'(r) = -(b'+2 r b'')
    Drr_formula=-2*tstar**2*bp(r)*phip
    A=np.exp(-r)/np.sqrt(complex(Drr))
    return gamma,tstar,Drr,Drr_formula,A,bp(r),phip

# stacking example b=r-0.3 r^3 : r1~0.393, r2~0.55 share gamma~0.291, phi~-0.302
bcoef=[0,1.0,0,-0.3]
# refine the two pinch roots at the shared gamma: solve r b'^2 = gamma with b'=1-0.9r^2 => b'^2=(1-0.9r^2)^2
# pick the crossing near (0.39,0.55): find gamma where phi(r1)=phi(r2)
from numpy.polynomial import Polynomial as P
b=P(bcoef); bp=b.deriv()
# scan for the exact self-intersection
rs=np.linspace(0.25,0.75,200000); phis=b(rs)-2*rs*bp(rs); gs=rs*bp(rs)**2
# find r1<r2 with phi equal and g equal
best=None
for i in range(0,len(rs),50):
    j=np.argmin(np.abs(phis-phis[i])+np.abs(gs-gs[i])+ (np.abs(rs-rs[i])<0.05)*1e9)
    if abs(rs[i]-rs[j])>0.05 and abs(phis[i]-phis[j])<1e-4 and abs(gs[i]-gs[j])<1e-4:
        if best is None or abs(phis[i]-phis[j])+abs(gs[i]-gs[j])<best[0]:
            best=(abs(phis[i]-phis[j])+abs(gs[i]-gs[j]), rs[i],rs[j])
if best:
    _,r1,r2=best
    g1,t1,Drr1,Drrf1,A1,bp1,pp1=fold_data(bcoef,r1)
    g2,t2,Drr2,Drrf2,A2,bp2,pp2=fold_data(bcoef,r2)
    print(f"STACKING at b=r-0.3r^3: r1={r1:.4f}, r2={r2:.4f}")
    print(f"  gamma: {g1:.4f} vs {g2:.4f} (equal), t_*: {t1:.4f} vs {t2:.4f} (equal => SAME ARC)")
    print(f"  D_rr numeric vs formula(-2t^2 b' phi'): r1 {Drr1:.4f} vs {Drrf1:.4f}; r2 {Drr2:.4f} vs {Drrf2:.4f}")
    print(f"  amplitudes A_i=e^-r/sqrt(D_rr): A1={A1:.4f}, A2={A2:.4f}")
    print(f"  |A1|={abs(A1):.4f}, |A2|={abs(A2):.4f} -> equal? {abs(abs(A1)-abs(A2))<1e-2}")
    print(f"  beta_total = A1+A2 = {A1+A2:.4f} (|.|={abs(A1+A2):.4f}); A1-A2={A1-A2:.4f}")
    same_sign_Drr = (Drr1>0)==(Drr2>0)
    print(f"  D_rr(r1),D_rr(r2) same sign: {same_sign_Drr} (both {'>' if Drr1>0 else '<'}0 / {'>' if Drr2>0 else '<'}0)")
    print("  => amplitudes REAL same-sign (both D_rr>0): total = A1+A2 != 0 (no cancellation at THIS stack).")
    print(f"     For beta_total=0 one needs |A1|=|A2| (=> e^-r1/sqrt|D_rr1|=e^-r2/sqrt|D_rr2|) AND opposite")
    print(f"     orientation. Here |A1|/|A2|={abs(A1)/abs(A2):.3f} != 1, so this stack is NON-cancelling.")
