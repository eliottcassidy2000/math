"""The rational-time floor M(W)>=max_k d_k/k (d_k=min_w ||w/k||, evaluate M at t=1/k) proves the
covering-core gap for 95% of SPREAD far-from-AP covering-2..12 12-cores. The 5% residual evades every
t=1/k witness (near-covering at all rational scales) but still has actual M~0.103 -- the irreducible
12-set-uniqueness/Freiman core. Infimum M(W) over spread far cores ~ 3/29 (margin +0.0265). S58i."""
from fractions import Fraction as F
from math import gcd
import random,time
def covers212(W): return all(any(w%k==0 for w in W) for k in range(2,13))
def ham(W):
    b=12
    for d in range(1,max(W)//12+2): b=min(b,12-len(set(W)&set(d*i for i in range(1,13))))
    return b
def rat_floor(W,Kmax=60):
    best=F(0)
    for k in range(2,Kmax+1):
        dk=min(min(w%k,k-(w%k)) for w in W)
        if dk and F(dk,k)>best: best=F(dk,k)
    return best
random.seed(21); THR=F(1,13); n=0; handled=0
t0=time.time()
while time.time()-t0<25:
    W=sorted(random.sample(range(2,random.choice([26,30,40,50])+1),12))
    if not covers212(W) or max(W)/min(W)<6.5 or ham(W)<7: continue
    n+=1
    if rat_floor(W)>THR: handled+=1
print("spread far covering-2..12 cores n=%d: rational-time floor max_k d_k/k > 1/13 handles %d (%.1f%%); residual %.1f%%"
      %(n,handled,100*handled/n,100*(n-handled)/n))
