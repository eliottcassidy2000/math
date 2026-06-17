"""
DISPROVE-side CRUX: minimize meas(G_C) over 12-sets C via HILL-CLIMB, hunting whether
the infimum is > 0 (prove wins) or -> 0 (disprove wins).

meas(G_C) = L(C) (lonely measure of the 12-set). If inf over 12-sets is 0, then a
13-set C u {w} can have L -> 0 and the singular-series route dies. The PROVE side's
decoupling floor was (6/7)*meas(G_C); the open gap is whether meas(G_C) itself stays
bounded below over ALL 12-sets (esp. large-lcm ones).

We hill-climb to MINIMIZE L(C) from many random/structured starts, allowing large
entries, tracking the lcm of the minimizers. We also do the analogous 13-set climb
(minimize L(S)) to confirm the global inf L.

Exact Fraction L. Float prescreen for neighbor evaluation, exact-confirm moves.
"""
from fractions import Fraction as Fr
from math import gcd, lcm, floor
from functools import reduce
import random

def Lf(S):
    A=[]
    for v in S:
        wv=1.0/(14*v)
        for k in range(v+1):
            c=k/v; lo=c-wv; hi=c+wv
            if lo<0: A.append((0.0,hi)); A.append((1+lo,1.0))
            elif hi>1: A.append((lo,1.0)); A.append((0.0,hi-1))
            else: A.append((lo,hi))
    A.sort()
    tot=0.0; cl=ch=None
    for a,b in A:
        if ch is None: cl,ch=a,b
        elif a<=ch:
            if b>ch: ch=b
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot

def danger_arcs(v):
    w=Fr(1,14*v); A=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: A+=[(Fr(0),hi),(1+lo,Fr(1))]
        elif hi>1: A+=[(lo,Fr(1)),(Fr(0),hi-1)]
        else: A.append((lo,hi))
    return A
def L_exact(S):
    A=[]
    for v in S: A.extend(danger_arcs(v))
    A.sort(key=lambda t:(t[0],t[1]))
    tot=Fr(0); cl=ch=None
    for a,b in A:
        if b<=a: continue
        if ch is None: cl,ch=a,b
        elif a<=ch:
            if b>ch: ch=b
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot

def climb(size, maxstart, seed, allow_big=True):
    random.seed(seed)
    glob=None; globS=None
    steps=(-1,1,-2,2,-3,3,-5,5,-7,7,-12,12)
    for st in range(maxstart):
        if st==0 and size==12:
            cur=tuple(range(1,13))
        elif st==0 and size==13:
            cur=tuple(range(1,14))
        else:
            hi=random.choice([15,25,40,80,200,600] if allow_big else [15,20,30])
            cur=tuple(sorted(random.sample(range(1,hi+1),size)))
        curv=Lf(cur)
        for _ in range(400):
            best=None;bestv=curv
            for i in range(size):
                for d in steps:
                    nx=cur[i]+d
                    if nx<1: continue
                    T=sorted(set(cur[:i]+cur[i+1:]+(nx,)))
                    if len(T)!=size: continue
                    v=Lf(tuple(T))
                    if v<bestv-1e-13:
                        bestv=v; best=tuple(T)
            if best is None: break
            cur,curv=best,bestv
        # exact confirm
        ev=L_exact(cur)
        if glob is None or ev<glob:
            glob=ev; globS=cur
    return glob,globS

print("=== Minimize meas(G_C) over 12-sets (hill-climb, big entries allowed) ===")
g,S=climb(12,120,1)
print(f"   12-set: min L = {g} = {float(g):.8f}  at {S}  lcm={lcm(*S)}  gcd={reduce(gcd,S)}")
print(f"     (compare meas(G_{{1..12}}) = {L_exact(tuple(range(1,13)))} = {float(L_exact(tuple(range(1,13)))):.8f})")

print()
print("=== Minimize L(S) over 13-sets (hill-climb, confirm global inf L) ===")
g13,S13=climb(13,120,2)
print(f"   13-set: min L = {g13} = {float(g13):.8f}  at {S13}  lcm={lcm(*S13)}  gcd={reduce(gcd,S13)}")
print(f"     1/1260 = {float(Fr(1,1260)):.8f};  below 1/1260? {g13<Fr(1,1260) and g13>0}")
if g13==0:
    print("     (hit a TIGHT config; min POSITIVE separately)")

print()
print("=== Targeted: push the large-lcm 12-minimizer further ===")
seed_cfg=(1,3,4,5,7,10,11,13,16,17,18,25)
print(f"   seed {seed_cfg}: L={float(L_exact(seed_cfg)):.8f} lcm={lcm(*seed_cfg)}")
# local exhaustive descent from this seed with finer steps
cur=seed_cfg; curv=Lf(cur)
steps=tuple(range(-8,0))+tuple(range(1,9))
for it in range(2000):
    best=None;bestv=curv
    for i in range(12):
        for d in steps:
            nx=cur[i]+d
            if nx<1: continue
            T=sorted(set(cur[:i]+cur[i+1:]+(nx,)))
            if len(T)!=12: continue
            v=Lf(tuple(T))
            if v<bestv-1e-13: bestv=v;best=tuple(T)
    if best is None: break
    cur,curv=best,bestv
print(f"   descended to: L={L_exact(cur)}={float(L_exact(cur)):.8f} at {cur} lcm={lcm(*cur)}")
