# THM-5109 (B), p = 1 rows (entry at the gadget source): h' = h - e + R; R' = 2^(m-1)(h + R) at q = 1,
# 2^(m-q) h' at q >= 2.  Brute force on 30 random bases of order <= 6, m <= 8 (all q): 1,080 cases.
import random
from tl import *
random.seed(7)
def rand_tour(n):
    out=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if random.random()<0.5: out[i]|=1<<j
            else: out[j]|=1<<i
    return out
bad=tot=0
for n in range(1,7):
    for _ in range(5):
        T=rand_tour(n); v=random.randrange(n)
        H,G,R,E=run2([T])[0]; h,e,r=H,E[v],R[v]
        jobs=[];keys=[]
        for m in range(1,9):
            for q in range(1,m+1):
                jobs.append(join(T,v,TT(m),0)); keys.append((m,q))
        for (m,q),(H2,G2,R2,E2) in zip(keys,run2(jobs)):
            z=len(T)+q-1; hp=h-e+r
            Rp = 2**(m-q)*hp if q>=2 else 2**(m-1)*(h+r)
            tot+=1
            if H2!=hp or R2[z]!=Rp: bad+=1
print("p=1 rows (h' = h - e + R; R' = 2^(m-1)(h+R) at q=1, 2^(m-q) h' at q>=2): checked",tot,'mismatches',bad)
