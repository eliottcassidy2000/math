# Wide-net covering-set counterexample hunt: L-filter then exact M.
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.path.insert(0, "04-computation")
# reuse L, is_covering, primitive, exact_M, safe_arcs from the main module by import
import importlib.util
spec = importlib.util.spec_from_file_location("m", "04-computation/lrc14_easy_dominates_hard_kps-S2-wf.py")
# can't run main; instead redefine minimal
def danger_arcs(v):
    w=F(1,14*v); out=[]
    for k in range(v+1):
        c=F(k,v); lo,hi=c-w,c+w
        if lo<0: out.append((F(0),hi)); out.append((1+lo,F(1)))
        elif hi>1: out.append((lo,F(1))); out.append((F(0),hi-1))
        else: out.append((lo,hi))
    return [(a,b) for a,b in out if a<b]
def L(S):
    iv=sorted(x for v in set(S) for x in danger_arcs(v))
    if not iv: return F(1)
    tot=F(0); cl,ch=iv[0]
    for lo,hi in iv[1:]:
        if lo<=ch:
            if hi>ch: ch=hi
        else: tot+=ch-cl; cl,ch=lo,hi
    tot+=ch-cl
    return F(1)-tot
def is_covering(S):
    Ss=set(S)
    return all(any(v%q==0 for v in Ss) for q in range(2,15))
def primitive(S): return reduce(gcd,S)==1
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
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
def M(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(v*t) for v in S)
        if v>b: b=v
    return b

random.seed(12345)
checked=0; l0=0; tight=[]; ctx=[]
# random covering 13-sets with mixed small + large coordinated speeds
for trial in range(60000):
    # build a covering set: ensure a multiple of each q in 2..14 cheaply by seeding
    base=set()
    # seed: take random subset of 1..13 then force covering by adding multiples
    k=random.randint(6,11)
    base=set(random.sample(range(1,26),k))
    # add up to a few large coordinated speeds
    mod=random.choice([14,28,42,84,12,7])
    for _ in range(random.randint(1,3)):
        base.add(mod*random.randint(1,40))
    # pad/trim to 13
    while len(base)<13:
        base.add(random.randint(1,300))
    S=sorted(base)[:13]
    if len(S)!=13: continue
    if not is_covering(S) or not primitive(S): continue
    checked+=1
    if L(S)==0:
        l0+=1
        m=M(S)
        if m<=F(1,14):
            (ctx if m<F(1,14) else tight).append((m,tuple(S)))
print(f"WIDE NET: checked {checked} covering+prim random 13-sets; L=0 survivors {l0}")
print(f"  tight (M=1/14): {len(tight)}; REAL counterexamples (M<1/14): {len(ctx)}")
if tight: print("  tight examples:", tight[:5])
if ctx: print("  COUNTEREXAMPLES:", ctx[:5])
