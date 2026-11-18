"""Independent cell, component and orbit-bound replay."""
from math import comb
from fractions import Fraction as F
from collections import Counter
from itertools import product
from functools import lru_cache
@lru_cache(None)
def cells(d,cr,h):
    g=d-(2 if cr else 1)-h
    C=Counter()
    for p,r in product(range(h+1),range(g+1)):
        for x,y in (product((0,1),repeat=2) if cr else ((0,0),)):
            a,b=p+r+x,h-p+r+y
            if a!=b:C[tuple(sorted((a,b)))]+=comb(h,p)*comb(g,r)*(1 if cr else 2)
    return C

def forest(d,C):
    par=list(range(d))
    def find(i):
        while par[i]!=i:i=par[i]
        return i
    out=[]
    for (a,b),w in sorted(C.items(),key=lambda x:(-x[1],x[0])):
        x,y=find(a),find(b)
        if x!=y:par[x]=y;out.append((a,b,w))
    return out
def rank(d,cr,h):
    if cr:return d-1
    if h==d-1:return d//2
    return d-(2 if h%2==0 else 1)
def orbits(d,cr):
    return range(0,d-1) if cr else range(1,d)
@lru_cache(None)
def beta(n):return F(comb(n,n//2),2**n)
def count(d,cr,h):
    return comb(d,2)*2**d*comb(d-2,h) if cr else d*2**(d-2)*comb(d-1,h)
def factor(d,cr,h):
    a=F(count(d,cr,h))
    for _,_,n in forest(d,cells(d,cr,h)):a*=beta(n)
    return a
def main():
    p=42*F(35,48)**9*F(3,8)**2
    c=52*F(3,4)**10/4
    a=5*F(35,48)**5*F(3,8)
    assert max(p,c,a)<F(3,4)
    for d in range(3,25):
        for cr in (0,1):
            for h in orbits(d,cr):
                C=cells(d,cr,h)
                assert len(forest(d,C))==rank(d,cr,h)
                assert all(w%(2 if cr else 4)==0 for w in C.values())
                if d>=13:
                    g=d-(2 if cr else 1)-h
                    H=h if g>=h else h-2
                    assert H in orbits(d-2,cr)
                    assert all(C[a+1,b+1]>=2*w for (a,b),w in cells(d-2,cr,H).items())
                    if not cr and h==d-1:assert C[0,d-1]>=4
                    else:
                        t=max(1,h)
                        assert C[0,t]>0 and C[d-1-t,d-1]>0
    for d,cap in ((11,F(4,25)),(12,F(11,1000))):
        A=sum(factor(d,cr,h) for cr in (0,1) for h in orbits(d,cr))
        assert A<cap
        print('exact base d',d,'below',cap,'approx only for presentation',float(A))
    print('uniform contraction factors',*(float(x) for x in (p,c,a)))
if __name__=='__main__':main()
