from fractions import Fraction as F
from math import gcd
def nrm(x):  # ||x|| for float
    x=x-int(x); x=x+1 if x<0 else x
    return min(x,1-x)
def good_intervals(W,lvl=1/13.0):
    # G_W = {t: min_w ||wt||>lvl}, as float intervals (complement of danger arcs)
    ivs=[]
    for w in W:
        for k in range(0,w+1):
            lo=(k-lvl)/w; hi=(k+lvl)/w
            ivs.append((max(lo,0.0),min(hi,1.0)))
    ivs=[(a,b) for a,b in ivs if b>a]; ivs.sort()
    m=[]
    for a,b in ivs:
        if m and a<=m[-1][1]: m[-1]=(m[-1][0],max(m[-1][1],b))
        else: m.append((a,b))
    comps=[]; prev=0.0
    for a,b in m:
        if a>prev: comps.append((prev,a))
        prev=max(prev,b)
    if prev<1.0: comps.append((prev,1.0))
    return comps
def avg_norm_over_G(W,vmax,N=200000):
    comps=good_intervals(W); tot=sum(b-a for a,b in comps)
    if tot==0: return 0,0
    # sample proportional to length
    import random; random.seed(0)
    s=0.0; cnt=0
    for a,b in comps:
        n=max(2,int(N*(b-a)/tot))
        for i in range(n):
            t=a+(b-a)*(i+0.5)/n
            s+=nrm(vmax*t); cnt+=1
    return s/cnt, tot
def M_exact(fam):
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best
def optdenom(fam):
    Q=2*max(fam)+2; best=F(0); dens=set()
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            fr=F(m,q)
            if fr>best: best=fr; dens={q}
            elif fr==best: dens.add(q)
    return dens
print("Second-moment / equidistribution: avg of ||182k*t|| over G_W. Need <1/13=%.4f for M(V)<1/13."%(1/13))
print("If avg ~ 1/4=0.25 (equidistributed), the far element CANNOT cover => M(V)>=1/13.\n")
cores=[("{2..12,15} nonaligned",[2,3,4,5,6,7,8,9,10,11,12,15]),
       ("{1..11,24} nonaligned",[1,2,3,4,5,6,7,8,9,10,11,24]),
       ("{1..12} AP (aligned)",list(range(1,13)))]
for name,W in cores:
    dens=optdenom(W); aligned=any(d%13==0 for d in dens)
    for vmax in [182,364]:
        avg,meas=avg_norm_over_G(W,vmax)
        print(f"  {name}: opt-denoms={sorted(dens)} aligned={aligned} meas(G_W)={meas:.4f}")
        print(f"     avg||{vmax}*t|| over G_W = {avg:.4f}  (>1/13? {avg>1/13} => far can't cover => M(V)>=1/13)")
        break
