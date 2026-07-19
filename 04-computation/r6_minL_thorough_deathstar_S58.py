import itertools, random
random.seed(3)
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
thr=1.0/2331; minL=1.0; mcfg=None
for P in allc[:150]:
    S0=csafe(P); lo=13*max(P)+1
    for _ in range(3):
        K=sorted(random.sample(range(max(lo,280),333),5))
        L=maxgap(S0,K); imp=True; it=0
        while imp and it<25:
            imp=False; it+=1
            for idx in range(5):
                for d in (-6,-3,-1,1,3,6):
                    K2=K[:];K2[idx]+=d
                    if len(set(K2))<5 or min(K2)<lo or max(K2)>332: continue
                    K2=sorted(K2); L2=maxgap(S0,K2)
                    if L2<L-1e-15: K=K2;L=L2;imp=True
        if L<minL: minL=L;mcfg=(P,K)
print("min L (150 cores sampled, killers<=332): %.7f  ratio to 1/2331 = %.3f"%(minL,minL/thr))
print("  extremal: core %s killers %s ; L*k5=%.4f"%(mcfg[0],mcfg[1],minL*max(mcfg[1])))
