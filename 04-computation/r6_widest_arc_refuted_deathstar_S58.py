# Does the WIDEST core-safe arc always retain a gap > 1/2331 under 5 fine killers?
# Adversarially MINIMIZE the widest-arc max-gap over 5 killers, per core. If min > 1/2331 for all cores,
# then L > 1/2331 via the widest arc alone -- and the widest arc has ~20 gaps (three-gap APPLIES there).
import itertools, random
random.seed(11)
def csafe(P):
    S=[(0.0,1.0)]
    for v in P:
        w=1.0/(14*v); arcs=[]
        for j in range(v):
            c=j/v; lo=(c-w)%1; hi=(c+w)%1
            if lo<hi: arcs.append((lo,hi))
            else: arcs.append((lo,1.0)); arcs.append((0.0,hi))
        for clo,chi in sorted(arcs):
            nn=[]
            for a,b in S:
                if chi<=a or clo>=b: nn.append((a,b)); continue
                if clo>a: nn.append((a,clo))
                if chi<b: nn.append((chi,b))
            S=nn
    return S
def arc_gap(a,bb,K):
    cuts=[]
    for k in K:
        w=1.0/(14*k); jlo=int(a*k); jhi=int(bb*k)+1
        for j in range(jlo,jhi+1):
            c=j/k; lo=c-w
            if lo>=bb: continue
            hi=c+w
            if hi<=a: continue
            cuts.append((a if lo<a else lo, bb if hi>bb else hi))
    cuts.sort(); cur=a; best=0.0
    for lo,hi in cuts:
        if lo>cur and lo-cur>best: best=lo-cur
        if hi>cur: cur=hi
    if bb-cur>best: best=bb-cur
    return best
THR=1.0/2331
allc=[list(c) for c in itertools.combinations(range(1,13),7)]
gmin=1.0; gcore=None
for P in allc:
    S=csafe(P); A0=max(S,key=lambda x:x[1]-x[0]); a,b=A0; lo=13*max(P)+1
    # adversarial min of widest-arc gap over 5 killers (fine: [lo, 332])
    lomin=1.0
    for _ in range(20):
        K=sorted(random.sample(range(max(lo,283),333),5))
        g=arc_gap(a,b,K); imp=True; it=0
        while imp and it<15:
            imp=False; it+=1
            for idx in range(5):
                for d in (-4,-2,-1,1,2,4):
                    K2=K[:];K2[idx]+=d
                    if len(set(K2))<5 or min(K2)<lo or max(K2)>332: continue
                    K2=sorted(K2); g2=arc_gap(a,b,K2)
                    if g2<g-1e-15: K=K2;g=g2;imp=True
        if g<lomin: lomin=g
    if lomin<gmin: gmin=lomin; gcore=(P,float(b-a),lomin)
print("adversarial MIN of widest-arc gap over all 792 cores (fine killers):")
print("  min = %.7f (1/2331=%.7f, ratio %.3f)"%(gmin,THR,gmin/THR))
print("  worst core %s (widest arc width %.4f, ~%.0f finest-killer gaps)"%(gcore[0],gcore[1],gcore[1]*332))
print("  widest-arc bound holds?", "YES -> L>1/2331 via ONE arc, 3-gap applies (many gaps)" if gmin>THR else "NO -- widest arc CAN be fragmented below 1/2331")
