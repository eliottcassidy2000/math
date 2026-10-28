#!/usr/bin/env python3
from itertools import product
from math import gcd,isqrt,prod
from fractions import Fraction

def need(b,m):
    if not b:raise RuntimeError(m)
def sigma(n):
    total=0
    for d in range(1,isqrt(n)+1):
        if n%d==0:total+=d+(n//d if d*d!=n else 0)
    return total
def span(rows):
    s={0}
    for r in rows:s|={x^r for x in list(s)}
    return s
def dot(a,x):return (a&x).bit_count()%2

def main():
    atoms=[4,2,5,3,7,13,11];cols=[4,2,3,0,1,5,2]
    def syndrome(x):
        v=0
        for i,c in enumerate(cols):
            if x>>i&1:v^=c
        return v
    def legal(x):return not(x&1 and x&2)
    C=[x for x in range(128) if syndrome(x)==0];V=[x for x in C if legal(x)]
    need(len(C)==16 and len(V)==12 and sum(map(legal,range(128)))==96,'counts')
    for x in range(128):
        selected=[a for i,a in enumerate(atoms) if x>>i&1]
        need(legal(x)==all(gcd(a,b)==1 for i,a in enumerate(selected) for b in selected[i+1:]),'conflict exact')
        if legal(x):need((isqrt(sigma(prod(selected)))**2==sigma(prod(selected)))==(syndrome(x)==0),'legal code interpretation')
    E=39;O=127^E
    fibers={o:[x for x in V if x&O==o] for o in sorted({x&O for x in C})}
    need(fibers[8]==[8] and fibers[16]==[22,49],'two contexts')
    need([prod(atoms[i] for i in range(7) if x>>i&1) for x in fibers[16]]==[70,364],'named sources')
    need(sigma(1560)==5040 and prod(sigma(a) for i,a in enumerate(atoms) if 47>>i&1)==7056,'illegal product')
    need(sorted(map(len,fibers.values()))==[1,1,1,1,2,2,2,2],'fiber sizes')
    need(sum(Fraction(len(f),12) for f in fibers.values() if len(f)==2)==Fraction(2,3),'conditional entropy')
    for x in C:
        x0=x&1;x1=(x>>1)&1
        need(x0^x1==((x>>4)^(x>>6))&1,'context parity')
        need(x0==(x0*(x0^x1))^(x0*x1),'Boolean certificate')
    # Exhaustive small linear systems: H constraints, A observations, target ell.
    checks=0
    for h,a,ell in product(range(8),repeat=3):
        for rhs in range(2):
            states=[x for x in range(8) if dot(h,x)==rhs]
            if not states:continue
            rec=all(dot(ell,x)==dot(ell,y) for x in states for y in states if dot(a,x)==dot(a,y))
            need(rec==(ell in span([h,a])),'row span iff');checks+=1
    print('PASS legal code12/formal16; one-bit context repair; 70/364 hostile; linear systems',checks)
if __name__=='__main__':main()
