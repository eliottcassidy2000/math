import itertools
P=[1,2,4,7,9,11,12]
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
S0=csafe(P)
def maxgap(K):
    best=0.0
    for (a,bb) in S0:
        cuts=[]
        for k in K:
            w=1.0/(14*k); jlo=int(a*k)-1; jhi=int(bb*k)+1
            for j in range(jlo,jhi+1):
                c=j/k; lo=c-w; hi=c+w
                if hi>a and lo<bb: cuts.append((max(lo,a),min(hi,bb)))
        cuts.sort(); cur=a
        for lo,hi in cuts:
            if lo>cur and lo-cur>best: best=lo-cur
            if hi>cur: cur=hi
        if bb-cur>best: best=bb-cur
    return best
def R(K):
    L=maxgap(K); return 1.0/(7*L*max(K)) if L>0 else 9.9
lo=157
gb=0.0; gbc=None; beat=[]
for k5 in range(165,205):
    win=list(range(max(lo,k5-16),k5))
    if len(win)<4: continue
    br=0.0; bk=None
    for combo in itertools.combinations(win,4):
        K=list(combo)+[k5]; r=R(K)
        if r>br: br=r; bk=K
    cons=[k5-8,k5-6,k5-4,k5-2,k5]
    rc=R(cons) if cons[0]>=lo else 0.0
    if br>rc+1e-9: beat.append((k5,br,bk,rc))
    if br>gb: gb=br; gbc=bk
print("k5 in [165,205], window opt (4-subsets of [k5-16,k5)):")
print("  global max R over all shapes = %.6f at %s"%(gb,gbc))
print("  configs where opt beat consecutive:",len(beat))
for k5,br,bk,rc in beat[:12]: print("    k5=%d opt R=%.5f %s vs consec %.5f"%(k5,br,bk,rc))
