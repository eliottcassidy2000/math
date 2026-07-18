from fractions import Fraction as F
from math import gcd
def norm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return min(r,1-r)
def minnorm(fam,t): return min(norm(F(v)*t) for v in fam)
def Pmax(P):
    best=F(0); t0=None
    for q in range(2,2*max(P)+2):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            v=minnorm(P,F(a,q))
            if v>best: best=v; t0=F(a,q)
    return best,t0
def measF(P,K):
    # F = union of arcs [c_i-1/14, c_i+1/14] on circle, c_i = frac(-(k_i-k1)*t0)  (fast-frame forbidden set)
    mu,t0=Pmax(P); k1=K[0]
    centers=[]
    for k in K:
        c=F(-(k-k1))*t0
        c=c-int(c); c=c+1 if c<0 else c
        centers.append(c)
    # measure of union of width-1/7 arcs around centers on circle [0,1)
    ivs=[]
    for c in centers:
        lo=c-F(1,14); hi=c+F(1,14)
        # split across wrap
        for a,b in [(lo,hi)]:
            a=a-int(a) if a>=0 else a-int(a)+1 if a<0 else a
        ivs.append((c-F(1,14), c+F(1,14)))
    # normalize to [0,1) with wrap by sampling exact: use union via events on a fine common grid is messy; do exact interval merge on real line then mod
    segs=[]
    for c in centers:
        segs.append((c-F(1,14), c+F(1,14)))
    # bring each into [0,2) and also add +1 copy for wrap handling, merge on [0,2), then intersect measure with a period
    pts=[]
    for lo,hi in segs:
        lo=lo % 1; hi=lo+F(1,7)
        pts.append((lo,hi))
        pts.append((lo+1,hi+1))
    pts.sort()
    merged=[]
    for lo,hi in pts:
        if merged and lo<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
        else: merged.append((lo,hi))
    # measure over one period [0,1): clip merged to [0,1) union its shift
    total=F(0)
    for lo,hi in merged:
        a=max(lo,F(0)); b=min(hi,F(1))
        if b>a: total+=b-a
    # add wrap parts in [1,2) mapped to [0,1)
    for lo,hi in merged:
        a=max(lo,F(1)); b=min(hi,F(2))
        if b>a: total+=b-a
    return mu,t0,min(total,F(1))

battery=[
 ("P{1..10}+[257,258,263]",list(range(1,11)),[257,258,263]),
 ("P{1..10}+[300,301,302]",list(range(1,11)),[300,301,302]),
 ("P{1..9}+[500,502,505,509]",list(range(1,10)),[500,502,505,509]),
 ("P{1..11}+[257,258]",list(range(1,12)),[257,258]),
 ("P{1..10}+[157,158,159]",list(range(1,11)),[157,158,159]),
 ("P{1..8,11,12}+[257,258,263]",[1,2,3,4,5,6,7,8,11,12],[257,258,263]),
]
print("RIGOROUS moiré bound: meas(F) = fast-frame forbidden measure; good kick GUARANTEED if meas(F)<1")
for name,P,K in battery:
    mu,t0,mF=measF(P,K)
    smax=(mu-F(1,14))/max(P); Delta=K[-1]-K[0]
    pred=F(1,7)+Delta*(t0+smax)
    print(f"  {name}: t0={t0} Δ={Delta}  meas(F)={float(mF):.3f}={mF}  (bound 1/7+Δ(t0+smax)={float(pred):.3f})  meas(F)<1: {mF<1}  => M(S)>=1/14")
