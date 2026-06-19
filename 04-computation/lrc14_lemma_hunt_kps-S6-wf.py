"""
Adversarial hunt for the ACTUAL lemma counterexample: mu_{1/7}(E) < thr_k, k=8..12.
Engine A (brute-validated). Plus control: same machinery beats consec for mu_{2/7}.
"""
from fractions import Fraction as F
from itertools import combinations
from math import floor
import random

def mu_theta(E,theta):
    E=sorted(set(E)); n=len(E); bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1); total=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2; order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        ks=[floor(E[order[t]]*mid) for t in range(n)]; subs=[]
        for t in range(n):
            o1=order[t];o2=order[(t+1)%n];k1=ks[t];k2=ks[(t+1)%n];wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1];c=F(k1-k2+wrap)
            if s==0:
                if c>theta: subs.append((a,b))
            elif s>0:
                lo=max(a,(theta-c)/s)
                if lo<b: subs.append((lo,b))
            else:
                hi=min(b,(theta-c)/s)
                if a<hi: subs.append((a,hi))
        subs.sort(); cur=cb=None
        for lo,hi in subs:
            if cur is None: cur,cb=lo,hi
            elif lo<=cb: cb=max(cb,hi)
            else: total+=cb-cur; cur,cb=lo,hi
        if cur is not None: total+=cb-cur
    return total

thr = {8:F(3637,5880),9:F(2025,4004),10:F(36,91),11:F(25,91),12:F(1,7)}
th=F(1,7)
consec_mu={8:F(691,735),9:F(247,294),10:F(38,49),11:F(1381,2205),12:F(13823,24255)}

print("=== (D) Exhaustive bounded-spread hunt: min mu_{1/7} over k-subsets of {0..S} ===")
for k in range(8,13):
    S = {8:13,9:13,10:13,11:14,12:14}[k]
    best=None; bestE=None; below=0; ties=0
    rng=list(range(1,S+1))
    for combo in combinations(rng,k-1):
        E=(0,)+combo
        m=mu_theta(list(E),th)
        if best is None or m<best:
            best=m; bestE=E
        if m<thr[k]: below+=1
        if m==consec_mu[k]: ties+=1
    print(f"k={k} spread<={S}: min mu={best}={float(best):.4f} at {bestE} | below thr={below} | #tie-with-consec={ties} | consec is argmin: {bestE==tuple(range(k))}")

print("\n=== (D2) Random/structured large-spread hunt, k=8..12 ===")
random.seed(2024)
glob_hit=None
for trial in range(40000):
    k=random.randint(8,12)
    mode=random.random()
    if mode<0.3:
        spread=random.randint(k-1,60); pts={0}
        while len(pts)<k: pts.add(random.randint(1,spread))
        E=sorted(pts)
    elif mode<0.5:  # AP with step
        s=random.randint(1,9); E=[s*i for i in range(k)]
    elif mode<0.7:  # perforated consecutive
        base=list(range(k+random.randint(1,8)))
        E=sorted(random.sample(base,k));
        if 0 not in E: E[0]=0; E=sorted(set(E))
        while len(E)<k: E.append(max(E)+1)
        E=sorted(set(E))[:k]
    elif mode<0.85: # geometric-ish
        E=[0]; v=1
        while len(E)<k: E.append(v); v+= random.randint(1,4)
    else: # Sidon-like
        E=[0,1];
        while len(E)<k:
            cand=E[-1]+random.randint(1,7); E.append(cand)
    E=sorted(set(E))
    if len(E)!=k: continue
    g=__import__('math').gcd
    from functools import reduce
    dg=reduce(g,[e for e in E if e>0])
    if dg>1: E=[e//dg for e in E]
    m=mu_theta(E,th)
    kk=len(E)
    if m<thr[kk]:
        glob_hit=(E,m,kk);
        print("!!! BELOW THR:",E,m,kk); break
print("global below-thr hit:",glob_hit)

print("\n=== (E) CONTROL: same machinery on mu_{2/7} should beat consec (search-not-impotent) ===")
th27=F(2,7)
for k in [10]:
    consec27=mu_theta(list(range(k)),th27)
    best27=consec27; bestE=tuple(range(k))
    random.seed(99)
    for _ in range(20000):
        spread=random.randint(k-1,40); pts={0}
        while len(pts)<k: pts.add(random.randint(1,spread))
        E=sorted(pts)
        from functools import reduce; from math import gcd
        dg=reduce(gcd,[e for e in E if e>0]) if len(E)>1 else 1
        if dg>1: E=[e//dg for e in E]
        m=mu_theta(E,th27)
        if m<best27: best27=m; bestE=tuple(E)
    print(f"mu_2/7: consec_{k}={float(consec27):.4f}  best_found={float(best27):.4f} at {bestE}  BEATEN:{best27<consec27}")
