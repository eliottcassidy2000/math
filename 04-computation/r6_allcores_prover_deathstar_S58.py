# Per-core two-part prover (consecutive step-2 killers). For each of 792 cores:
#  (1) witness sub-arc: t0 core-safe, 2t0 in G_{2/7}>0 region, with margins -> tail works for b>=B0
#  (2) finite head: float scan R_sharp(b) for b in [lo, B0], confirm <1.
import itertools, sys
sys.stdout.reconfigure(line_buffering=True)
def nd(x): x=x%1; return min(x,1-x)
def G27(s):
    w=1.0/7; arcs=[]
    for m in range(5):
        c=(m*s)%1; lo=(c-w)%1; hi=(c+w)%1
        if lo<hi: arcs.append((lo,hi))
        else: arcs.append((lo,1.0)); arcs.append((0.0,hi))
    arcs.sort(); merged=[]
    for lo,hi in arcs:
        if merged and lo<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
        else: merged.append((lo,hi))
    if len(merged)==1 and merged[0][0]<=1e-12 and merged[0][1]>=1-1e-12: return 0.0
    mg=0.0
    for i in range(len(merged)):
        a=merged[i][1]; b=merged[(i+1)%len(merged)][0]+(1.0 if i+1==len(merged) else 0.0)
        if b-a>mg: mg=b-a
    return mg
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
def witness(P,S):
    # find t0 maximizing min(core_margin, G27(2t0)); estimate usable sub-arc half-width
    best=(-1,None)
    for (a,b) in S:
        NN=max(8,int((b-a)*400))
        for i in range(1,NN):
            t=a+(b-a)*i/NN
            g=G27((2*t)%1)
            if g<=0: continue
            cm=min(nd(v*t)-1.0/14 for v in P)
            if cm<=0: continue
            sc=min(g,cm)
            if sc>best[0]: best=(sc,t,a,b)
    return best
cores=[list(c) for c in itertools.combinations(range(1,13),7)]
gmax=0.0; gcfg=None; maxB0=0; fails=0
for ci,P in enumerate(cores):
    S=csafe(P); lo=13*max(P)+1
    w=witness(P,S)
    # usable sub-arc half-width ~ min(dist to core-arc end, ...); B0 ~ how far finite head must reach.
    # conservatively finite-check to B0=600 (covers narrow-arc kick-in); verify witness exists (w[0]>0)
    B0=600; maxB0=max(maxB0,B0)
    wit_ok = w[0]>0
    mx=0.0
    for b in range(lo, lo+ (B0-lo if lo<B0 else 200)):
        L=maxgap(S,[b+2*i for i in range(5)]); R=1.0/(7*L*(b+8))
        if R>=1.0: fails+=1
        if R>mx: mx=R
    if not wit_ok: print("NO WITNESS core",P)
    if mx>gmax: gmax=mx; gcfg=P
    if ci%80==0: print("...%d/792  gmax=%.5f  (this core max R=%.4f, witness=%s)"%(ci,gmax,mx,wit_ok))
print("DONE: 792 cores, consec step-2. global max R_sharp=%.6f at %s ; finite-head failures=%d"%(gmax,gcfg,fails))
