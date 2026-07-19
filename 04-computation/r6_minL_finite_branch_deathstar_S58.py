# (a) hunt for MIN L over (core, 5 killers in [92,332]); is it > 1/2331?
import itertools, random
random.seed(1)
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
def maxgap(S0,K):
    best=0.0
    for (a,bb) in S0:
        cuts=[]
        for k in K:
            w=1.0/(14*k); jlo=int(a*k); jhi=int(bb*k)+1
            for j in range(jlo,jhi+1):
                c=j/k; lo=c-w; hi=c+w
                if hi>a and lo<bb: cuts.append((max(lo,a),min(hi,bb)))
        cuts.sort(); cur=a
        for lo,hi in cuts:
            if lo>cur and lo-cur>best: best=lo-cur
            if hi>cur: cur=hi
        if bb-cur>best: best=bb-cur
    return best
allc=[list(c) for c in itertools.combinations(range(1,13),7)]
random.shuffle(allc)
thr=1.0/2331
minL=1.0; mcfg=None
for P in allc[:60]:
    S0=csafe(P); lo=13*max(P)+1
    for _ in range(8):
        K=sorted(random.sample(range(max(lo,300),333),5))  # finest killers near 332, spread
        L=maxgap(S0,K)
        # hill-climb to MINIMIZE L
        imp=True
        while imp:
            imp=False
            for idx in range(5):
                for d in (-4,-2,-1,1,2,4):
                    K2=K[:];K2[idx]+=d
                    if len(set(K2))<5 or min(K2)<lo or max(K2)>332: continue
                    K2=sorted(K2); L2=maxgap(S0,K2)
                    if L2<L-1e-15: K=K2;L=L2;imp=True
        if L<minL: minL=L;mcfg=(P,K)
print("min L found over (core, 5 killers<=332) = %.7f  (1/2331=%.7f)"%(minL,thr))
print("  L/threshold = %.4f  at core %s killers %s"%(minL/thr,mcfg[0],mcfg[1]))
print("  min L > 1/2331 ?", minL>thr, "  (=> exactly-one-killer>=333 tail holds if bound proven)")
