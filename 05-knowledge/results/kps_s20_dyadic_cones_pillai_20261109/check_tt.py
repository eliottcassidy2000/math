# verify closed forms for transitive gadgets TT_m, entry position p, exit position q (1-indexed)
import random
from tl import *
random.seed(20)
def rand_tour(n):
    out=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if random.random()<0.5: out[i]|=1<<j
            else: out[j]|=1<<i
    return out
bases=[]
for n in range(1,7):
    for _ in range(4):
        T=rand_tour(n); v=random.randrange(n); bases.append((T,v))
def predict(h,e,R,m,p,q):
    hp = h + (2**(p-1))*R if p>=2 else h - e + R
    if p>=2:
        if q==p: Rp = 2**(m-p)*h + 2**(m-1)*R
        elif q>p: Rp = 2**(m-q)*(h + 2**(p-1)*R)
        elif q>=2: Rp = 2**(m-q)*h + 2**(m-p+q)*3**(p-q-1)*R
        else: Rp = 2**(m-1)*h - 2**(m-p)*e + 2**(m-p+1)*3**(p-2)*R
    else: Rp=None
    ez = 1 if q==m else 0
    P = 1 if (p==1 and q==m) else 0
    if q==m and m>p: Euz=2**(p-1)
    elif q<p and p==m: Euz=2**(q-1)
    else: Euz=0
    ep = ez*h - P*e + Euz*R
    return hp,ep,Rp
bad=0;tot=0
for (T,v) in bases:
    H,G,R,E=run2([T])[0]; h,e,r=H,E[v],R[v]
    jobs=[];keys=[]
    for m in range(1,8):
        for p in range(1,m+1):
            for q in range(1,m+1):
                N=join(T,v,TT(m),p-1); jobs.append(N); keys.append((m,p,q))
    res=run2(jobs)
    for (m,p,q),(H2,G2,R2,E2) in zip(keys,res):
        z=len(T)+q-1
        hp,ep,Rp=predict(h,e,r,m,p,q)
        tot+=1
        if H2!=hp or E2[z]!=ep or (Rp is not None and R2[z]!=Rp):
            bad+=1
            if bad<8: print('MISMATCH',(h,e,r),(m,p,q),'got',(H2,E2[z],R2[z]),'pred',(hp,ep,Rp))
print('checked',tot,'mismatches',bad)
