# Adversarial: is the dense even AP [0,2,..] really the worst base B for p0_inf at ratio 2/1?
# Scan ALL bases B subset {0..14} with 0 in B of the relevant size, compute p0_inf(B,2,1).
from fractions import Fraction as F
from itertools import combinations
import sys
sys.path.insert(0,'/Users/e/Documents/GitHub/math/04-computation')

def sector(fr):
    s=int(fr*7); return 6 if s>=7 else s

def cell_law(p,q):
    bps={F(0),F(1)}
    for c in (q,p):
        for m in range(0,c+1):
            for j in range(7):
                x=F(7*m+j,7*c)
                if 0<=x<=1: bps.add(x)
    bps=sorted(bps); mu=[[F(0)]*7 for _ in range(7)]
    for a,b in zip(bps,bps[1:]):
        mid=(a+b)/2
        fq=(q*mid)-int(q*mid); fq+= (1 if fq<0 else 0)
        fp=(p*mid)-int(p*mid); fp+= (1 if fp<0 else 0)
        mu[sector(fq)][sector(fp)]+=b-a
    return mu

def g_B(B):
    bps={F(0),F(1)}
    for b in B:
        if b==0: continue
        for m in range(0,b+1):
            for j in range(7):
                x=F(7*m+j,7*b)
                if 0<=x<=1: bps.add(x)
    bps=sorted(bps); g=[[F(0)]*7 for _ in range(7)]
    for a,b2 in zip(bps,bps[1:]):
        mid=(a+b2)/2; hit=set()
        for b in B:
            fr=(b*mid)-int(b*mid); fr+=(1 if fr<0 else 0); hit.add(sector(fr))
        length=b2-a; missing=set(range(7))-hit
        for i in range(7):
            for j in range(7):
                if missing<={i,j}: g[i][j]+=length
    return g

def p0_inf(B,p,q):
    mu=cell_law(p,q); g=g_B(B); s=F(0)
    for i in range(7):
        for j in range(7): s+=mu[i][j]*g[i][j]
    return s

# k=10: cluster size 10, r=2 far => base size 8. B subset {0..14}, 0 in B, |B|=8.
cap10=F(55,91)
best=(F(0),None)
cnt=0
for rest in combinations(range(1,15),7):
    B=(0,)+rest
    v=p0_inf(list(B),2,1)
    cnt+=1
    if v>best[0]: best=(v,B)
print(f"k=10 (base size 8, B subset 0..14), ratio 2/1: scanned {cnt} bases")
print(f"  worst base = {best[1]}  p0_inf={best[0]}={float(best[0]):.5f}  cap_10={float(cap10):.5f}  below cap:{best[0]<cap10}")
print(f"  dense even AP [0,2,4,6,8,10,12,14] p0_inf = {float(p0_inf([0,2,4,6,8,10,12,14],2,1)):.5f}")
