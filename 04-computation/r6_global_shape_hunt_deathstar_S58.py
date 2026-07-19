import itertools, random
random.seed(0)
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
def R(S0,K):
    L=maxgap(S0,K); return 1.0/(7*L*max(K)) if L>0 else 9.9
def climb(S0,lo,seeds,krange):
    gb=0.0;gc=None
    for _ in range(seeds):
        k5=random.randint(lo+4,lo+krange)
        K=sorted(random.sample(range(lo,k5),4))+[k5]; r=R(S0,K)
        imp=True
        while imp:
            imp=False
            for idx in range(5):
                for d in (-3,-2,-1,1,2,3):
                    K2=K[:];K2[idx]+=d
                    if len(set(K2))<5 or min(K2)<lo: continue
                    K2=sorted(K2); r2=R(S0,K2)
                    if r2>r+1e-12: K=K2;r=r2;imp=True
        if r>gb:gb=r;gc=K
    return gb,gc
# thorough on P0
P0=[1,2,4,7,9,11,12]; S0=csafe(P0)
gb,gc=climb(S0,157,200,50)
print("P0 hard hill-climb (200 seeds): best R=%.6f at %s"%(gb,gc))
# sample 40 other cores, lighter
allc=[list(c) for c in itertools.combinations(range(1,13),7)]
random.shuffle(allc)
gmax=gb; gcfg=(P0,gc)
for P in allc[:40]:
    S=csafe(P); lo=13*max(P)+1
    r,c=climb(S,lo,8,40)
    if r>gmax: gmax=r;gcfg=(P,c)
print("plus 40 sampled cores: global best R=%.6f at %s"%(gmax,gcfg))
print("consecutive record 0.801133 -> extremal", "HOLDS" if gmax<=0.801134 else "REFUTED")
