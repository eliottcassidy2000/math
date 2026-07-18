# death-star-S56: COMPLETE closure of compact case max<=34, using the EXACT cover-gap criterion.
# For V=W+{vmax}: M(V)<1/13 iff far element covers G_W iff coverGap(W,vmax)=max_{t in G_W}||vmax t|| < 1/13.
# coverGap is O(#components of G_W), INDEPENDENT of vmax (no O(vmax^2) M_exact). Exact via Fractions.
# WLOG (proved): only near-tight cores (M<=NEAR=1/13+34/2366) matter (else stability); near-tight=>covers 2..10 &
#  residues mod 13 hit all 6 antipodal pairs (or has a mult of 13). Far element vmax = multiple of L=lcm(missing cap{2..13}),
#  >max(W), <= max(W)/(13 delta). Dilated APs = the conclusion, excluded.
from fractions import Fraction as F
from math import gcd, floor, ceil
from itertools import combinations
import multiprocessing as mp
NEAR=F(1,13)+F(34,2366); TH=F(1,13)
RESBIT=[1<<(v%13) for v in range(35)]
CMASK=[0]*35
for v in range(35):
    m=0
    for i,q in enumerate(range(2,11)):
        if v%q==0: m|=(1<<i)
    CMASK[v]=m
FULLC=(1<<9)-1
PAIRS=[(1,12),(2,11),(3,10),(4,9),(5,8),(6,7)]
DTAB=bytearray(8192)
for rm in range(8192):
    if rm&1: DTAB[rm]=1; continue
    ok=True
    for (r,s) in PAIRS:
        if not ((rm>>r)&1 or (rm>>s)&1): ok=False; break
    DTAB[rm]=1 if ok else 0
def M_and_args(fam):
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            fr=F(r,q)
            if fr>NEAR: return None
            if fr>best: best=fr
    return best
def good_components(W):
    ivs=[]
    for w in W:
        for k in range(0,w+1):
            lo=(F(k)-TH)/w; hi=(F(k)+TH)/w
            a=lo if lo>0 else F(0); b=hi if hi<1 else F(1)
            if b>a: ivs.append((a,b))
    ivs.sort()
    mg=[]
    for a,b in ivs:
        if mg and a<=mg[-1][1]:
            if b>mg[-1][1]: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    comps=[]; prev=F(0)
    for a,b in mg:
        if a>prev: comps.append((prev,a))
        if b>prev: prev=b
    if prev<1: comps.append((prev,F(1)))
    return comps
def ndist(x):
    f=x-floor(x)
    return f if f<=F(1,2) else 1-f
def max_norm_on(a,b,vmax):        # exact max of ||vmax t|| over [a,b]
    va=vmax*a; vb=vmax*b
    if vb-va>=1: return F(1,2)
    m=ndist(va)
    d=ndist(vb)
    if d>m: m=d
    klo=ceil(va-F(1,2))            # smallest half-integer >= va is klo+1/2
    if klo+F(1,2)<=vb: return F(1,2)
    return m
def far_covers(comps,vmax):       # True iff coverGap<1/13 (far element covers G_W => M(V)<1/13 => COUNTEREXAMPLE)
    for (a,b) in comps:
        if max_norm_on(a,b,vmax)>=TH: return False   # escape point found => not covered
    return True
def is_dilated_AP(W):
    W=sorted(W); d=W[0]
    return all(W[i]==d*(i+1) for i in range(12))
def work(ij):
    i,j=ij; counter=[]; near=0; checks=0; scanned=0
    D=DTAB; RB=RESBIT; CM=CMASK
    rm0=RB[i]|RB[j]; cm0=CM[i]|CM[j]
    for tail in combinations(range(j+1,35),10):
        scanned+=1
        rm=rm0; cm=cm0
        for v in tail: rm|=RB[v]; cm|=CM[v]
        if not D[rm]: continue
        if cm!=FULLC: continue
        W=(i,j)+tail
        M=M_and_args(W)
        if M is None: continue
        near+=1
        if is_dilated_AP(W): continue
        delta=M-TH
        if delta<=0: continue
        miss=[q for q in range(2,14) if not any(v%q==0 for v in W)]
        L=1
        for q in miss: L=L*q//gcd(L,q)
        ub=int(max(W)/(13*delta)); start=((max(W)//L)+1)*L
        if start>ub: continue
        comps=good_components(W)
        for vmax in range(start,ub+1,L):
            checks+=1
            if far_covers(comps,vmax):
                counter.append((tuple(W),vmax))
    return counter,near,checks,scanned
if __name__=='__main__':
    TASKS=[(i,j) for i in range(1,24) for j in range(i+1,25)]
    with mp.Pool(4) as p:
        results=[]
        for n,r in enumerate(p.imap_unordered(work,TASKS,chunksize=1),1):
            results.append(r)
            if n%40==0: print("  ...%d/%d tasks done"%(n,len(TASKS)),flush=True)
    C=[]; near=0; checks=0; scanned=0
    for c,nn,ch,sc in results:
        C+=c; near+=nn; checks+=ch; scanned+=sc
    print("DONE: %d subsets scanned; near-tight non-AP covering cores=%d; far-element cover-gap checks=%d"%(scanned,near,checks),flush=True)
    print("COUNTEREXAMPLES (coverGap<1/13, i.e. M(V)<1/13, non-AP core, max<=34): %d"%len(C),flush=True)
    if not C:
        print("=> NO counterexample. COMPACT CASE max<=34 CLOSED: every near-tight non-AP core has coverGap>=1/13 (M(V)>=1/13) for ALL far elements.",flush=True)
    else:
        for w in C[:30]: print("  *** COUNTEREXAMPLE:",list(w[0]),"vmax=",w[1],flush=True)
