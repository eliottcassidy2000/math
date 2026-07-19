from fractions import Fraction as F
import itertools
def nd(x): x=x%1; return min(x,1-x)
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
# for each core: find widest usable witness sub-window (t-interval where core-safe-margin>0 and G27(2t)>0)
# width -> B0 ~ ceil(1/width)+8. Report max B0.
cores=[list(c) for c in itertools.combinations(range(1,13),7)]
maxB0=0; worst=None
for P in cores:
    S=csafe(P); bestw=F(0)
    for (a,b) in S:
        # scan for a maximal subinterval where G27(2t)>0 (sampled) -- estimate width via samples
        NN=max(10,int((b-a)*600)); run=0; runstart=None; localbest=F(0)
        prevgood=False
        for i in range(NN+1):
            t=a+(b-a)*F(i,NN) if i<NN else b
            good = (G27((2*t)%1)>0) and (min(nd(v*t)-F(1,14) for v in P)>0)
            if good and not prevgood: runstart=t
            if (not good or i==NN) and prevgood:
                w=t-runstart
                if w>localbest: localbest=w
            prevgood=good
        if localbest>bestw: bestw=localbest
    if bestw>0:
        B0=int(1/float(bestw))+10
        if B0>maxB0: maxB0=B0; worst=(P,float(bestw))
print("max B0 (tail kick-in) over 792 cores ~ %d  at core %s"%(maxB0,worst))
print("=> finite head must reach b=%d ; (my prover uses 600)"%maxB0)
