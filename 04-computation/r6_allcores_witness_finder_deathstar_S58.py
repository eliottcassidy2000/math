from fractions import Fraction as F
import itertools
def nd(x): x=x%1; return min(x,1-x)
# G_{2/7}(sigma): largest gap left by 5 arcs of half-width 1/7 at AP step sigma
def G27(s):
    w=F(1,7); arcs=[]
    for m in range(5):
        c=(m*s)%1; lo=(c-w)%1; hi=(c+w)%1
        if lo<hi: arcs.append((lo,hi))
        else: arcs.append((lo,F(1))); arcs.append((F(0),hi))
    arcs.sort(); merged=[]
    for lo,hi in arcs:
        if merged and lo<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
        else: merged.append((lo,hi))
    if len(merged)==1 and merged[0][0]<=0 and merged[0][1]>=1: return F(0)
    mg=F(0)
    for i in range(len(merged)):
        a=merged[i][1]; b=merged[(i+1)%len(merged)][0]+(1 if i+1==len(merged) else 0)
        if b-a>mg: mg=b-a
    return mg
def danger(v):
    w=F(1,14*v); out=[]
    for j in range(v):
        c=F(j,v); lo=(c-w)%1; hi=(c+w)%1
        if lo<hi: out.append((lo,hi))
        else: out.append((lo,F(1))); out.append((F(0),hi))
    return out
def subf(S,arcs):
    for clo,chi in sorted(arcs):
        new=[]
        for lo,hi in S:
            if chi<=lo or clo>=hi: new.append((lo,hi)); continue
            if clo>lo: new.append((lo,clo))
            if chi<hi: new.append((chi,hi))
        S=new
    return S
def csafe(P):
    S=[(F(0),F(1))]
    for v in P: S=subf(S,danger(v))
    return S
# per core: scan a grid of t; find t0 with core margin>0 and G27(2 t0)>0, maximize min(core_margin_scaled, G27)
cores=[list(c) for c in itertools.combinations(range(1,13),7)]
cov=0; uncov=[]
for P in cores:
    S=csafe(P)
    found=False; best=F(0)
    for (a,b) in S:
        # sample interior of this core-safe arc
        NN=24
        for i in range(1,NN):
            t=a+(b-a)*F(i,NN)
            g=G27((2*t)%1)
            if g>0:
                # core margin at t: min_v(||v t||-1/14)
                cm=min(nd(v*t)-F(1,14) for v in P)
                if cm>0 and min(g,cm)>best:
                    best=min(g,cm); found=True
    if found: cov+=1
    else: uncov.append(P)
print("robust finder (sub-arc) covers %d/792 cores"%cov)
if uncov: print("uncovered %d:"%len(uncov), uncov[:10])
else: print("ALL 792 COVERED")
