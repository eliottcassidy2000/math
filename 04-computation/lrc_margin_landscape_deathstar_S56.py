from fractions import Fraction as F
from math import gcd
import random
random.seed(11)
def covers(fam): return all(any(v%q==0 for v in fam) for q in range(2,15))
def prim(fam):
    from functools import reduce
    return reduce(gcd,fam)==1
def M_arg(fam):
    Q=2*max(fam)+2; best=F(0); arg=None
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q); arg=(a,q)
    return best,arg
def rand_cov(Vmax):
    for _ in range(500):
        fam=sorted(random.sample(range(2,Vmax+1),13))
        if covers(fam): return fam
    return None
def neigh(fam,Vmax):
    for i in range(13):
        for d in (-3,-2,-1,1,2,3):
            nv=fam[i]+d
            if 2<=nv<=Vmax and nv not in fam:
                nf=sorted(fam[:i]+[nv]+fam[i+1:])
                if covers(nf): yield nf

# strong descent for the GLOBAL covering floor, tracking witness denom's prime factors
def pfac(n):
    f=set()
    d=2
    while d*d<=n:
        while n%d==0: f.add(d); n//=d
        d+=1
    if n>1: f.add(n)
    return sorted(f)

best=(F(1),None,None)
for Vmax in (28,34,40):
    for start in range(50):
        cur=rand_cov(Vmax)
        if cur is None: continue
        curM,curA=M_arg(cur)
        for _ in range(80):
            imp=False
            for nb in neigh(cur,Vmax):
                m,a=M_arg(nb)
                if m<curM: cur,curM,curA=nb,m,a; imp=True; break
            if not imp: break
        if curM<best[0]: best=(curM,cur,curA)
M,fam,arg=best
print(f"GLOBAL covering-floor descent: M={M}={float(M):.5f}")
print(f"  family={fam}  primitive={prim(fam)}  mult-of-14? {any(v%14==0 for v in fam)}")
print(f"  witness t={arg[0]}/{arg[1]}, denom={arg[1]} = primes {pfac(arg[1])}")
print(f"  compare: 1/14={1/14:.5f} 1/13={1/13:.5f} 1/11={1/11:.5f} 1/9={1/9:.5f} 1/8={1/8:.5f}")
