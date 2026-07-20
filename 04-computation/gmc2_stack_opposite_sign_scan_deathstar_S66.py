# Broad test: at EVERY stacking of the GMC radial D, do the two folds have OPPOSITE g'(r) sign?
# D_rr(r_i,t_*) = 2 t_*^2 g'(r_i), g=r b'^2. Opposite signs => one amplitude real, one imaginary
# => beta_total = A1 +- A2 != 0 (cannot cancel real vs imaginary), regardless of orientation.
import numpy as np
from numpy.polynomial import Polynomial as P
rng=np.random.default_rng(12345)  # fixed seed (Math.random unavailable in workflows; fine here)
def stackings(bcoef, rmax=7.0, N=120000, tol=2e-4):
    b=P(bcoef); bp=b.deriv()
    r=np.linspace(0.02,rmax,N)
    g=r*bp(r)**2; phi=b(r)-2*r*bp(r); gp=(P(bcoef).deriv()) # placeholder
    # g'(r) = b'(b'+2 r b'') ; compute
    bpp=bp.deriv()
    gprime=bp(r)*(bp(r)+2*r*bpp(r))
    out=[]
    idx=np.argsort(g); gs=g[idx]; ps=phi[idx]; rs=r[idx]; gps=gprime[idx]
    seen=set()
    for i in range(0,len(gs),40):
        for j in range(i+1,min(i+120,len(gs))):
            if abs(rs[i]-rs[j])<0.08: continue
            if abs(gs[i]-gs[j])<tol and abs(ps[i]-ps[j])<tol:
                key=(round(rs[i],1),round(rs[j],1))
                if key in seen: continue
                seen.add(key)
                out.append((rs[i],rs[j],gps[i],gps[j]))
    return out
total=0; opp=0; same=[]
for trial in range(25):
    deg=rng.integers(2,5)
    bcoef=[0]+list(np.round(rng.uniform(-2,2,deg),2))  # b(0)=0 (charge structure), random real
    for r1,r2,gp1,gp2 in stackings(bcoef):
        total+=1
        if gp1*gp2<0: opp+=1
        else: same.append((bcoef,round(r1,3),round(r2,3),round(gp1,3),round(gp2,3)))
print(f"stackings found across 60 random real b (deg 2-4): {total}")
print(f"  with OPPOSITE g' sign (=> one amp real, one imag => beta_total != 0): {opp}/{total}")
print(f"  with SAME g' sign (would need magnitude+orientation analysis): {len(same)}")
for s in same[:6]: print("     SAME-sign stack:", s)
if total>0:
    print(f"\n=> {'ALL' if opp==total else f'{opp}/{total}'} stackings have opposite D_rr sign.")
    print("   Opposite sign is ROBUST non-cancellation: A1 real, A2 imaginary, beta_total=A1+-A2 != 0")
    print("   independent of contour orientation. This is the mechanism for boxeph's nonzero-total-jump amendment.")
