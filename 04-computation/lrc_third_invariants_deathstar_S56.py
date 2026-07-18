# death-star-S56: test THREE candidate "third invariants" for the near-tight fragmented residual,
# mined from unrelated threads (disc_v autocorrelation / arc-count / isolation; three-gap; displacement).
# Goal: which invariant vanishes EXACTLY at the AP (deep well) and is bounded away from 0 for non-AP residual cores?
from fractions import Fraction as F
from math import gcd
def M_and_args(fam):
    Q=2*max(fam)+2; best=F(0); args=[]
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam); fr=F(r,q)
            if fr>best: best=fr; args=[(a,q)]
            elif fr==best: args.append((a,q))
    return best,args
def good_components(W,lvl=F(1,13)):
    # danger arcs of W at level lvl; good set = complement; return sorted list of (lo,hi) good comps
    ivs=[]
    for w in W:
        for k in range(0,w+1):
            lo=(F(k)-lvl)/w; hi=(F(k)+lvl)/w
            ivs.append((max(lo,F(0)),min(hi,F(1))))
    ivs=[(a,b) for a,b in ivs if b>a]; ivs.sort()
    mg=[]
    for a,b in ivs:
        if mg and a<=mg[-1][1]: mg[-1]=(mg[-1][0],max(mg[-1][1],b))
        else: mg.append((a,b))
    comps=[]; prev=F(0)
    for a,b in mg:
        if a>prev: comps.append((prev,a))
        prev=max(prev,b)
    if prev<F(1): comps.append((prev,F(1)))
    return comps
def displacement(W):
    # for each good component, distance of its CENTER to nearest a/13; max over comps = max displacement s
    comps=good_components(W)
    if not comps: return 0.0, 0
    smax=0.0
    for (lo,hi) in comps:
        c=(lo+hi)/2
        d13=min(abs(float(c)-a/13.0) for a in range(0,14))
        smax=max(smax,d13)
    return smax, len(comps)
def cover_gap(W,vmax):
    # MAX over good set of ||vmax t||.  If >= 1/13, EXISTS t in G_W with ||vmax t||>=1/13 => M(V)>=1/13
    # (far element cannot cover G_W). If < 1/13, far element covers G_W (potential deep well).
    comps=good_components(W); import math
    best=0.0
    for (lo,hi) in comps:
        a=float(lo); b=float(hi)
        # sample densely + include sawtooth peaks (2j+1)/(2 vmax) lying in [a,b]
        pts=[a+(b-a)*i/2000.0 for i in range(2001)]
        jlo=math.ceil((2*a*vmax-1)/2.0); jhi=math.floor((2*b*vmax-1)/2.0)
        for j in range(jlo,jhi+1): pts.append((2*j+1)/(2.0*vmax))
        for t in pts:
            x=(vmax*t)%1.0; best=max(best,min(x,1.0-x))
    return best
def disc_v(W,vmax,cap=200):
    # crude autocorrelation-discrepancy proxy: sum_{m=1}^{cap} |chat_{m*vmax}|^2 over good set, chat via components
    comps=good_components(W); import math
    s=0.0
    for m in range(1,cap+1):
        N=m*vmax; re=0.0; im=0.0
        for (lo,hi) in comps:
            a=float(lo); b=float(hi)
            # integral of e^{2pi i N t} dt over [a,b] = (e^{2pi i N b}-e^{2pi i N a})/(2pi i N)
            re+=(math.sin(2*math.pi*N*b)-math.sin(2*math.pi*N*a))/(2*math.pi*N)
            im+=-(math.cos(2*math.pi*N*b)-math.cos(2*math.pi*N*a))/(2*math.pi*N)
        s+=re*re+im*im
    return s
AP=[i for i in range(1,13)]
cores={
 "AP {1..12} (deep well)": AP,
 "2*AP {2..24} (dilated, deep well)": [2*i for i in range(1,13)],
 "{1..10,22,24}": [1,2,3,4,5,6,7,8,9,10,22,24],
 "{1..10,24,33}": [1,2,3,4,5,6,7,8,9,10,24,33],
 "{1,2,3,5,7..12,17,19}": [1,2,3,5,7,8,9,10,11,12,17,19],
 "{2..12,15}": [2,3,4,5,6,7,8,9,10,11,12,15],
 "{1..11,24}": [1,2,3,4,5,6,7,8,9,10,11,24],
 "3*{1..11}+{34}": [3,6,9,12,15,18,21,24,27,30,33,34],
}
print(f"{'core':40s} {'M':>9s} {'delta':>8s} {'smax':>7s} {'C':>3s} {'coverGap':>9s} {'covers?':>8s}  (smallest vmax)")
for name,W in cores.items():
    M,args=M_and_args(W); delta=float(M)-1/13.0
    smax,C=displacement(W)
    miss=[q for q in range(2,14) if not any(v%q==0 for v in W)]
    L=1
    for q in miss: L=L*q//gcd(L,q)
    vmax=((max(W)//L)+1)*L if L>0 else 182   # SMALLEST valid far element (binding candidate)
    cg=cover_gap(W,vmax)
    covers = "YES(well)" if cg<1/13.0 else "no"
    print(f"{name:40s} {float(M):9.5f} {delta:8.5f} {smax:7.4f} {C:3d} {cg:9.4f} {covers:>8s}  (vmax={vmax},L={L})")
print()
print("READING: coverGap = max over good set G_W of ||vmax*t||.  coverGap>=1/13=0.0769 => a point of G_W escapes the")
print("far element's danger => M(V)>=1/13 (far element CANNOT cover). Only the AP-type (smax=0) has coverGap<1/13:")
print("its good set sits ON the 1/13-lattice where the far element vanishes (the deep well). Non-AP => smax>0 => coverGap>=1/13.")
