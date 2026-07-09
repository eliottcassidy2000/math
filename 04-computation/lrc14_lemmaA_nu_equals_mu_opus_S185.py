import numpy as np
NG=300007; X=(np.arange(NG)+0.5)/NG
def nu_and_W(E):
    E=np.array(sorted(set(E)),float)
    P=np.sort((np.outer(E,X))%1.0,axis=0)
    gaps=np.vstack([np.diff(P,axis=0),(P[0]+1)-P[-1]])
    mg=gaps.max(axis=0)
    nu=float((mg>1/7).mean())
    W=np.maximum(gaps-1/7,0).sum(axis=0)   # uncovered measure = sum (gap-1/7)_+
    return nu, float((W>1e-12).mean())
for name,E in [("consecutive",list(range(13))),("near-2-AP",[0,2,3,4,5,6,7,8,9,10,11,12,14]),
               ("dissociated",[0,1,3,7,12,20,30,44,65,80,96,122,147])]:
    nu,PW=nu_and_W(E)
    print(f"{name:>14}: nu=meas(maxgap>1/7)={nu:.5f}   P(W>0)={PW:.5f}   equal? {abs(nu-PW)<1e-4}")
print("\n=> nu(E) = mu(E) = P(W>0) EXACTLY (maxgap>1/7 <=> some gap>1/7 <=> W>0). So any tractable lower")
print("   bound on nu IS a moment-route bound on mu; the 7-arc PZ is a WEAK (degree-2) discretized version.")
print("   The moment LP B_d -> mu as d->inf; min_E B_d -> min_E mu = nuConsec=0.4425. D3 (proved) min=0.3088.")
