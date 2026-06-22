# FAST loop-closure: a few genuine admissible covering+primitive S in LRC(6),(10),
# confirm M(S)>=1/n AND that S has a tau with a free 1/q gap (the witness mechanism).
from fractions import Fraction as Fr
import itertools
from math import gcd
from functools import reduce
def norm(y): r=y%1; return min(r,1-r)
def M_exact(S):
    S=sorted(set(int(v) for v in S)); cand=set()
    for v in S:
        for j in range(0,2*v): cand.add(Fr(2*j+1,2*v))
        for j in range(0,v): cand.add(Fr(j,v))
    for a in S:
        for b in S:
            for sg in (a+b,a-b):
                if sg:
                    d=abs(sg)
                    for j in range(0,d+1): cand.add(Fr(j,d))
    best=Fr(0); arg=None
    for t in cand:
        if 0<=t<=1:
            mn=min(norm(v*t) for v in S)
            if mn>best: best=mn; arg=t
    return best,arg
def is_cov(S,n): return all(any(v%qq==0 for v in S) for qq in range(2,n+1))
def is_prim(S): return reduce(gcd,S)==1

# pick the FIRST few admissible S for n=6 (speeds<=30) and n=10 (speeds<=40)
for n in (6,10):
    q=n//2; found=[]
    rng = 28 if n==6 else 42
    import random; random.seed(7)
    seen=set(); tries=0
    while len(found)<6 and tries<300000:
        tries+=1
        S=tuple(sorted(random.sample(range(1,rng),n-1)))
        if S in seen: continue
        seen.add(S)
        if is_prim(S) and is_cov(S,n): found.append(S)
    print(f"=== LRC({n}) (1/{n}={float(Fr(1,n)):.5f}), q={q}, gap 1/q=1/{q} ===")
    allok=True
    for S in found:
        M,t=M_exact(S)
        # at the optimal tau t, what is the circular maxgap of the runner positions {v t}?
        ph=sorted(set((v*t)%1 for v in S))
        if len(ph)<=1: g=Fr(1)
        else:
            g=max(b-a for a,b in zip(ph,ph[1:])); g=max(g,(ph[0]+1)-ph[-1])
        ok = M>=Fr(1,n)
        allok &= ok
        print(f"  S={S}: M={M}={float(M):.5f} {'>=1/n OK' if ok else 'FAIL'}; "
              f"opt tau={t}, maxgap@tau={g}={float(g):.4f} (>1/q={float(Fr(1,q)):.4f}? {g>Fr(1,q)})")
    print(f"  => all admissible M>=1/{n}: {allok}\n")
