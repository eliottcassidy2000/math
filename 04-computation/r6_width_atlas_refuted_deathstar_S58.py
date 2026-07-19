from fractions import Fraction as F
import math
def la(x,q): r=x%q; return min(r,q-r)
def test(delta_den, kmax):
    delta=F(1,delta_den)
    KILL=list(range(92,kmax+1))
    cands=[]
    for q in range(15,221):
        for a in range(1,q):
            cm=min(F(la(p*a,q),q)-F(1,14) for p in range(1,13))
            if cm>=6*delta:
                # CORRECT "hit": danger(k) intersects [a/q +- delta/2]  <=> ||k a/q|| < 1/14 + k*delta/2
                cov=[1 if F(la(k*a,q),q) < F(1,14)+k*delta/2 else 0 for k in KILL]
                cands.append(cov)
    nK=len(KILL); w=[1.0]*nK; eta=0.5; vals=[]
    for t in range(3000):
        W=sum(w); pw=[x/W for x in w]
        bc=1e9
        for cov in cands:
            c=sum(pw[i]*cov[i] for i in range(nK))
            if c<bc: bc=c; bcov=cov
        vals.append(bc)
        for i in range(nK):
            if bcov[i]: w[i]*=math.exp(eta)
        if t%40==0: m=max(w); w=[x/m for x in w]
    return sum(vals)/len(vals), len(cands)
for dd in (2300, 2325):
    V,n=test(dd, 332)
    print("delta=1/%d, killers[92,332], CORRECT hit-condition: V*~%.4f (need<0.2), %d candidates -> %s"%(
        dd,V,n,"WORKS" if V<0.2 else "FAILS (finest killers hit too wide a band)"))
# show hit-band fraction for k=332
d=1.0/2300; print("k=332: hit band = ||.||<1/14+332*d/2 = %.4f (fraction ~%.3f of residues)"%(1/14+332*d/2, 2*(1/14+332*d/2)))
