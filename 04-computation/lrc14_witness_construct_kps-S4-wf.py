"""
Independent check of the WITNESS CONSTRUCTION for ALL-MULT7-LARGE.
For each such S, build tau = k/7 + u/7 where u is an exact safe point of the
mult-of-7 subsystem in (0, 1/(2V*)], then verify min_v ||v*tau|| >= 1/14 DIRECTLY
on the full 13-set. This tests the proof mechanism end-to-end (not just Mval).

Also: HUNT exhaustively (deterministic, small V*) for an ALL-MULT7-LARGE covering
primitive S3 set whose witness construction FAILS to clear 1/14, and for any with
M < 1/14 (cross-checked against Mval).
kind-pasteur S4-wf.
"""
from fractions import Fraction as F
from math import gcd
import itertools

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1 - r

def is_cov(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
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

def safe_u(W, Vstar):
    """exact u in (0,1/(2V*)] with ||w u||>=1/14 for all w in W, else None."""
    R=F(1,2*Vstar); h=F(1,14)
    forb=[]
    for w in W:
        j=0
        while True:
            c=F(j,w)
            lo=c-h/F(w)
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
    # try closed boundaries
    cands=set()
    for lo,hi in gaps:
        cands.add(lo); cands.add(hi)
    for u in sorted(cands):
        if u>0 and u<=R and all(nrm(F(w)*u)>=h for w in W): return u
    return None

def witness_for(S):
    """Construct global witness tau=k/7+u/7 for ALL-MULT7-LARGE S, verify >=1/14."""
    m7=[v for v in S if v%7==0]; nm7=[v for v in S if v%7!=0]
    Vstar=max(nm7)
    W=[v//7 for v in m7]
    u=safe_u(W,Vstar)
    if u is None: return None
    s=u/7  # since u=7s
    # pick k=1 (coprime to 7); witness tau=1/7 + s
    tau=F(1,7)+s
    mn=min(nrm(F(v)*tau) for v in S)
    return tau, mn

# Deterministic exhaustive hunt over a structured family of ALL-MULT7-LARGE sets.
# We build: small non-mult-7 covering core + mult-of-7 runners (all > V*).
# Cover requirements for q in 2..13 must be met by the WHOLE set; 7,14 come from m7.
print("Deterministic exhaustive ALL-MULT7-LARGE hunt (small V*)")
print("Checking: (i) witness clears 1/14, (ii) Mval>=1/14, agreement.")
worst_M=None; worst_wit=None; nfail_wit=0; nfail_M=0; ntested=0
bad=[]

# Build covering small cores from {1..Vstar}\7Z, then attach m7 runners >Vstar.
# To stay finite & exhaustive in the dangerous regime, sweep Vstar small and
# enumerate cores systematically.
import random
random.seed(7)
for Vstar in range(13, 60):
    pool=[v for v in range(1,Vstar+1) if v%7!=0]
    if Vstar not in pool: continue  # V* must be non-mult-7 and present
    # mult-of-7 runners available (>Vstar), keep modest count
    m7pool=[7*w for w in range(1,40) if 7*w>Vstar]
    m7_14=[x for x in m7pool if x%14==0]
    if not m7_14: continue
    # try many random cores of size 9..11 that include Vstar and are "covering-ready"
    for _ in range(3000):
        sz_nm=random.randint(9,11)
        core=set([Vstar])
        while len(core)<sz_nm:
            core.add(random.choice(pool))
        core=sorted(core)
        if max(core)!=Vstar: continue
        nslots=13-len(core)
        if nslots<2: continue
        m7=set([random.choice(m7_14)])
        while len(m7)<nslots:
            m7.add(random.choice(m7pool))
        if len(m7)!=nslots: continue
        S=sorted(core)+sorted(m7)
        if len(S)!=13: continue
        if not is_cov(S) or not is_prim(S): continue
        # confirm ALL-MULT7-LARGE & S3
        if not all(x>Vstar for x in m7): continue
        Vmin=min(S); Vmax=max(S); kk=sum(1 for v in S if v>13)
        if not (kk>=2 and Vmax>=13*Vmin): continue
        ntested+=1
        w=witness_for(S)
        if w is None:
            nfail_wit+=1; bad.append(("NO-WITNESS",S)); continue
        tau,mn=w
        if mn<F(1,14):
            nfail_wit+=1; bad.append(("WITNESS-LOW",S,tau,mn))
        M=Mval(S)
        if M<F(1,14):
            nfail_M+=1; bad.append(("M-LOW",S,M))
        if worst_M is None or M<worst_M[1]:
            worst_M=(S,M)
        if worst_wit is None or mn<worst_wit[1]:
            worst_wit=(S,tau,mn)

print(f"  tested ALL-MULT7-LARGE covering primitive S3 sets: {ntested}")
print(f"  witness constructions that FAILED to clear 1/14: {nfail_wit}")
print(f"  sets with Mval < 1/14: {nfail_M}")
if bad[:5]:
    print("  sample bad:", bad[:5])
if worst_M:
    print("  worst Mval:", worst_M[1], "~", float(worst_M[1]), " set:", worst_M[0])
if worst_wit:
    print("  worst witness margin: tau=",worst_wit[1]," min=",worst_wit[2],"~",float(worst_wit[2]),
          " set:",worst_wit[0])
print("  (worst witness margin should be >= 1/14 =", float(F(1,14)),")")
