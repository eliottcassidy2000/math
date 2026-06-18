"""
LEAN adversarial verification of ALL-MULT7-LARGE closure. kind-pasteur S4-wf.
Separates cheap witness-construction (bulk) from expensive Mval (spot checks).
All exact (Fraction).
"""
from fractions import Fraction as F
from math import gcd
import itertools, random

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def is_prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def Mval(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b

def safe_u(W,Vstar):
    R=F(1,2*Vstar); h=F(1,14); forb=[]
    for w in set(W):
        j=0
        while True:
            c=F(j,w); lo=c-h/F(w)
            if lo>R: break
            a=max(F(0),lo); b=min(R,c+h/F(w))
            if a<=b: forb.append((a,b))
            j+=1
    forb.sort(); merged=[]
    for a,b in forb:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    prev=F(0); gaps=[]
    for a,b in merged:
        if a>prev: gaps.append((prev,a))
        prev=max(prev,b)
    if prev<R: gaps.append((prev,R))
    for lo,hi in gaps:
        if hi<=0: continue
        u=(lo+hi)/2
        if u<=0: u=hi/2
        if u>0 and all(nrm(F(w)*u)>=h for w in W): return u
    for lo,hi in gaps:
        for u in (lo,hi):
            if u>0 and u<=R and all(nrm(F(w)*u)>=h for w in W): return u
    return None

def witness_margin(S):
    """Return (tau,min) for ALL-MULT7-LARGE S via tau=1/7+u/7, or None if no safe u."""
    m7=[v for v in S if v%7==0]; nm7=[v for v in S if v%7!=0]
    Vstar=max(nm7); W=[v//7 for v in m7]
    u=safe_u(W,Vstar)
    if u is None: return None
    tau=F(1,7)+u/7
    return tau, min(nrm(F(v)*tau) for v in S)

# ---------- D1: abstract teeth-survival sweep (no Mval) ----------
# Dangerous regime = smallest w (widest teeth). Pool wmin..wmin+12 covers it; the
# measure argument (D2) handles the rest. r<=4 matches claim scope; also probe r=5,6.
print("[D1] abstract: for V*=14..399, smallest-w pool (12 wide), |W|<=4, all 7w>V*")
fails=[]; tot=0
for Vstar in range(14,400):
    wmin=Vstar//7+1; wpool=list(range(wmin,wmin+13))
    for r in range(1,5):
        for combo in itertools.combinations(wpool,r):
            tot+=1
            if safe_u(combo,Vstar) is None:
                fails.append((Vstar,combo))
print(f"  configs={tot}, NO-safe-u failures={len(fails)}")
if fails[:5]: print("  sample:",fails[:5])

# ---------- D2: measure bound (abstract guarantee) ----------
print("\n[D2] measure bound: max r with guaranteed safe point at every V*")
def blocked_ub(w,Vstar):
    R=F(1,2*Vstar); h=F(1,14)
    floorwR=(F(w)*R).numerator//(F(w)*R).denominator
    half0=min(R,h/F(w))
    return half0 + floorwR*min(F(1,7*w),R)
maxr={}
for Vstar in range(14,400):
    R=F(1,2*Vstar); wmin=Vstar//7+1
    for r in range(1,8):
        ws=range(wmin,wmin+r)
        if sum(blocked_ub(w,Vstar) for w in ws)>=R:
            maxr[Vstar]=r-1; break
    else: maxr[Vstar]=7
print("  min over V* of measure-guaranteed r:",min(maxr.values()))
print("  (so for any subsystem with r <= that, a safe point provably exists)")

# ---------- max possible r (#mult-of-7) in an ALL-MULT7-LARGE S3 13-set ----------
# covering needs mults of 2..13 from the WHOLE set. mult-of-7 are all > V*, and
# the non-mult7 part must already cover 2..6,8..13 (everything but 7) AND include
# V*. Minimum non-mult7 count? We just bound r empirically in the bulk hunt below.

# ---------- BULK witness hunt (cheap) ----------
print("\n[C] bulk ALL-MULT7-LARGE hunt: worst witness margin (no Mval in loop)")
random.seed(20260618)
worst=None; ntested=0; nowit=0; rseen=set()
for _ in range(300000):
    Vstar=random.randint(13,200)
    if Vstar%7==0: continue
    pool=[v for v in range(1,Vstar+1) if v%7!=0]
    if len(pool)<8: continue
    core=set([Vstar])
    while len(core)<random.randint(8,11):
        core.add(random.choice(pool))
    if max(core)!=Vstar: continue
    nslots=13-len(core)
    if nslots<2: continue
    m7pool=[7*w for w in range(1,300) if 7*w>Vstar]
    m7_14=[x for x in m7pool if x%14==0]
    if not m7_14: continue
    m7=set([random.choice(m7_14)])
    while len(m7)<nslots: m7.add(random.choice(m7pool))
    if len(m7)!=nslots: continue
    S=sorted(core)+sorted(m7)
    if len(S)!=13: continue
    if not is_cov(S) or not is_prim(S): continue
    if not all(x>Vstar for x in m7): continue
    Vmin=min(S); Vmax=max(S); kk=sum(1 for v in S if v>13)
    if not (kk>=2 and Vmax>=13*Vmin): continue
    ntested+=1; rseen.add(len(m7))
    wm=witness_margin(S)
    if wm is None: nowit+=1; print("  NO WITNESS:",S); continue
    tau,mn=wm
    if worst is None or mn<worst[1]: worst=(S,tau,mn)
print(f"  tested={ntested}, no-witness={nowit}, r(#mult7) values seen={sorted(rseen)}")
if worst:
    print("  WORST witness margin:",worst[2],"~",float(worst[2]),">=1/14?",worst[2]>=F(1,14))
    print("    set:",worst[0]," tau:",worst[1])
    # cross-check Mval on this worst set:
    print("    Mval(worst) =",Mval(worst[0]),"~",float(Mval(worst[0])),">=1/14?",Mval(worst[0])>=F(1,14))
print("\nDONE")
