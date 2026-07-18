# death-star-S56: COMPLETE closure of compact case max<=34, FAST (4 workers, no per-subset allocation).
# Same math as lrc_residual_parallel; optimized D13 gate via a precomputed 13-bit residue-mask lookup table.
# WLOG (proved this session): only NEAR-TIGHT cores (M<=NEAR=1/13+34/2366) can extend to a counterexample (else
# stability window empty); near-tight => covers 2..10 and (residues mod 13 hit all 6 antipodal pairs OR has a mult of 13);
# far element vmax is a multiple of L=lcm(missing cap {2..13}), >max(W), <= max(W)/(13 delta). Dilated APs excluded.
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import multiprocessing as mp
NEAR=F(1,13)+F(34,2366)
RESBIT=[1<<(v%13) for v in range(35)]                 # residue bit for each element value
CMASK=[0]*35
for v in range(35):
    m=0
    for i,q in enumerate(range(2,11)):
        if v%q==0: m|=(1<<i)
    CMASK[v]=m
FULLC=(1<<9)-1
PAIRS=[(1,12),(2,11),(3,10),(4,9),(5,8),(6,7)]
DTAB=bytearray(8192)                                   # DTAB[resmask]=1 iff D13 gate passes
for rm in range(8192):
    if rm & 1:                                         # has a residue-0 element (mult of 13) => D13=0, pass
        DTAB[rm]=1; continue
    ok=True
    for (r,s) in PAIRS:
        if not ((rm>>r)&1 or (rm>>s)&1): ok=False; break
    DTAB[rm]=1 if ok else 0
def M_and_args(fam):
    Q=2*max(fam)+2; best=F(0); args=[]
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            fr=F(r,q)
            if fr>NEAR: return None,None
            if fr>best: best=fr; args=[(a,q)]
            elif fr==best: args.append((a,q))
    return best,args
TH=F(1,13)
def M_bounded_ge(fam,qmax):
    # returns True as soon as a witness time a/q (q<=qmax) gives min_V ||.|| >= 1/13  => M(V)>=1/13.
    # This is a rigorous LOWER-BOUND certificate (O(qmax^2), independent of the far element's size).
    for q in range(2,qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if 13*r>=q: return True
    return False
def M_exact(fam):   # full max over q<=2*max (used only to CONFIRM a rare small-witness failure)
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(r,q)>best: best=F(r,q)
    return best
def is_dilated_AP(W):
    W=sorted(W); d=W[0]
    return all(W[i]==d*(i+1) for i in range(12))
def safe_at(vmax,args):
    for (a,q) in args:
        r=(vmax*a)%q
        if 13*min(r,q-r)>=q: return True
    return False
def work(ij):
    i,j=ij
    counter=[]; near=0; checks=0; scanned=0
    D=DTAB; RB=RESBIT; CM=CMASK
    rm0=RB[i]|RB[j]; cm0=CM[i]|CM[j]
    for tail in combinations(range(j+1,35),10):
        scanned+=1
        rm=rm0; cm=cm0
        for v in tail: rm|=RB[v]; cm|=CM[v]
        if not D[rm]: continue                          # D13 gate (table lookup, no alloc)
        if cm!=FULLC: continue                          # covers 2..10
        W=(i,j)+tail
        M,args=M_and_args(W)
        if M is None: continue
        near+=1
        if is_dilated_AP(W): continue
        delta=M-F(1,13)
        if delta<=0: continue
        miss=[q for q in range(2,14) if not any(v%q==0 for v in W)]
        L=1
        for q in miss: L=L*q//gcd(L,q)
        ub=int(max(W)/(13*delta)); start=((max(W)//L)+1)*L; qmax=2*max(W)+2
        for vmax in range(start,ub+1,L):
            checks+=1
            V=list(W)+[vmax]
            if not M_bounded_ge(V,qmax):                 # no small-denominator witness
                if M_exact(V)<F(1,13):                    # confirm with full search (rare escalation)
                    counter.append((tuple(W),vmax))
    return counter,near,checks,scanned
if __name__=='__main__':
    TASKS=[(i,j) for i in range(1,24) for j in range(i+1,25)]   # W = {i<j} + 10 from {j+1..34}; need 34-j>=10 => j<=24
    with mp.Pool(4) as p:
        results=[]
        for n,r in enumerate(p.imap_unordered(work, TASKS, chunksize=1),1):
            results.append(r)
            if n%40==0: print("  ...%d/%d (i,j)-tasks done"%(n,len(TASKS)),flush=True)
    C=[]; near=0; checks=0; scanned=0
    for c,n,ch,sc in results:
        C+=c; near+=n; checks+=ch; scanned+=sc
    print("DONE: %d subsets scanned; near-tight non-AP covering cores=%d; far-element checks=%d"%(scanned,near,checks),flush=True)
    print("COUNTEREXAMPLES (M(V)<1/13, non-AP core, max<=34): %d"%len(C),flush=True)
    if not C:
        print("=> NO counterexample. COMPACT CASE max<=34 CLOSED: every near-tight non-AP core has M(V)>=1/13 for ALL far elements.",flush=True)
    else:
        for w in C[:30]: print("  *** COUNTEREXAMPLE:",list(w[0]),"vmax=",w[1],flush=True)
