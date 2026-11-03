#!/usr/bin/env python3
"""Independent exact controls. Every test is live under python -O."""
from fractions import Fraction
from itertools import product
from math import gcd, prod
import json

checks = {}
def need(ok, name):
    checks[name] = checks.get(name, 0) + 1
    if not ok:
        raise RuntimeError(name)

def fibs(n):
    f=[0,1]
    while len(f)<=n: f.append(f[-1]+f[-2])
    return f
F=fibs(100)
def zeck(n):
    out=[]
    while n:
        i=max(i for i in range(2,len(F)) if F[i]<=n)
        out.append(i);n-=F[i]
    return tuple(sorted(out))
def val2(n):
    if not n: raise ValueError('valuation of zero')
    n=abs(n);return (n&-n).bit_length()-1

def encode(I): return 3*prod(1+(1<<(i+2)) for i in I)
def decode(h):
    if h%3: raise ValueError('wrong scalar core')
    p=h//3;out=[]
    while p!=1:
        j=val2(p-1);f=1+(1<<j)
        if j<4 or p%f: raise ValueError('not a code')
        out.append(j-2);p//=f
    return tuple(out)

def one(rows,q,letter):
    n=len(rows);new=list(rows)+[((1<<n)-1)^(1<<q)]
    new[q]|=1<<n
    return tuple(new),n if letter=='C' else q

def move(v,letter):
    s,u,d=v
    return (u,s+u+d,s) if letter=='C' else (s+d,2*u,0)

def strong(rows):
    n=len(rows);full=(1<<n)-1
    for root in range(n):
        seen=1<<root
        while True:
            nxt=seen
            for u in range(n):
                if seen>>u&1:nxt|=rows[u]
            if nxt==seen:break
            seen=nxt
        if seen!=full:return False
    return True

def endpoint(rows,q):
    n=len(rows);full=(1<<n)-1
    # dp[mask][v] counts paths on mask starting at v, independently of C/D.
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,full+1):
        bits=mask
        while bits:
            low=bits&-bits;first=low.bit_length()-1;bits-=low
            rest=mask^low
            if not rest:continue
            nxt=rows[first]&rest;total=0
            while nxt:
                lo=nxt&-nxt;nxt-=lo
                total+=dp[rest][lo.bit_length()-1]
            dp[mask][first]=total
    inc=[i for i in range(n) if rows[i]>>q&1]
    need(len(inc)==1,'literal unique guard')
    s=dp[full][q];h=sum(dp[full]);d=dp[full^(1<<q)][inc[0]]
    return s,h-s,d

def build(I):
    rows=(2,4,1);q=0;state=(1,2,0)
    for i in sorted(I):
        for letter in 'D'*(i+1)+'CD':
            rows,q=one(rows,q,letter);state=move(state,letter)
    return rows,q,state

def cmul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def cpow(a,n):
    z=(1,0)
    for _ in range(n):z=cmul(z,a)
    return z
primes=(5,13,17);pis=((2,1),(3,2),(4,1))
def alpha(v):
    z=(1,0)
    for i,n in enumerate(v):
        p=pis[i] if n>=0 else (pis[i][0],-pis[i][1])
        z=cmul(z,cpow(p,abs(n)))
    return z

def carry(v,w):return prod(p**((abs(a)+abs(b)-abs(a+b))//2) for p,a,b in zip(primes,v,w))
def cochain(v):return prod(p**max(-a,0) for p,a in zip(primes,v))
def normreduce(v,w):
    # Deliberately normalize independent pairs of X/Y exponents.
    coefficient=1;out=[]
    for p,a,b in zip(primes,v,w):
        x=max(a,0)+max(b,0);y=max(-a,0)+max(-b,0)
        pairs=min(x,y);coefficient*=p**pairs;out.append(x-y)
    return coefficient,tuple(out)

def circle(v):
    a,b=alpha(v);c=a*a+b*b
    return Fraction(a*a-b*b,c),Fraction(2*a*b,c)

def main():
    for n in range(10000):
        I=zeck(n);h=encode(I)
        need(all(b-a>=2 for a,b in zip(I,I[1:])),'Zeckendorf gap')
        need(sum(F[i] for i in I)==n,'Zeckendorf sum')
        need(decode(h)==I,'strong-code decoding')
    prefixes=[]
    for mask in range(1<<10):
        I=tuple(i+2 for i in range(10) if mask>>i&1)
        prefixes.append(I)
    for a,I in enumerate(prefixes):
        for J in prefixes[a+1:]:
            need(val2(encode(I)-encode(J))==min(set(I)^set(J))+2,'exact 2-adic distance')
    residue_counts=[]
    for K in range(4,20):
        n=K-4;admissible=[];all_codes=set()
        for mask in range(1<<n):
            I=tuple(i+2 for i in range(n) if mask>>i&1)
            value=encode(I)%(1<<K);all_codes.add(value)
            if mask&(mask<<1)==0:admissible.append(value)
        need(len(set(admissible))==F[K-2],'Fibonacci residue cardinality')
        need(len(all_codes)==1<<n,'full-support ball residues')
        need(all(x%16==3 for x in all_codes),'full-support ball coset')
        residue_counts.append([K,len(admissible)])
    literal=[]
    for I in [(),(2,),(3,),(4,),(5,),(6,),(7,),(2,4)]:
        rows,q,state=build(I);h=encode(I);t=h//3
        need(len(rows)==3+sum(i+3 for i in I),'literal order')
        need(strong(rows),'literal strength')
        need(state==(t,2*t,0),'transfer canonical state')
        need(endpoint(rows,q)==(t,2*t,0),'independent endpoint DP')
        literal.append({'support':I,'order':len(rows),'H':h,'guard':q})
    vectors=list(product(range(-2,3),repeat=3))
    for v in vectors:
        for w in vectors:
            coefficient,out=normreduce(v,w);g=carry(v,w)
            need(out==tuple(a+b for a,b in zip(v,w)) and coefficient==g,'integral normal product')
            need(g*cochain(out)==cochain(v)*cochain(w),'carry coboundary')
            raw=cmul(alpha(v),alpha(w));normal=alpha(out)
            need(raw==tuple(g*x for x in normal),'Gaussian specialization')
            need(gcd(abs(raw[0]),abs(raw[1]))==g,'exact Gaussian content')
            for i,p in enumerate(primes):
                need((g%p==0)==(v[i]*w[i]<0),'prime-local vanishing')
            x,y=circle(v);u,z=circle(w);U=x*x;V=u*u
            center=1-U-V+2*U*V
            branches=((x*u-y*z)**2,(x*u+y*z)**2)
            for Z in branches:
                need(Z*Z-2*center*Z+(U+V-1)**2==0,'quotient correspondence')
            need((branches[0]-branches[1])**2==16*U*V*(1-U)*(1-V),'quotient discriminant')
    small=list(product(range(-1,2),repeat=3))
    for u in small:
        for v in small:
            for w in small:
                uv=tuple(a+b for a,b in zip(u,v));vw=tuple(a+b for a,b in zip(v,w))
                need(carry(u,v)*carry(uv,w)==carry(v,w)*carry(u,vw),'carry reassociation')
    U=Fraction(9,25);V=Fraction(25,169)
    for Z in [Fraction(33*33,65*65),Fraction(63*63,65*65)]:
        need(Z*Z-2*(1-U-V+2*U*V)*Z+(U+V-1)**2==0,'norm-65 example')
    print(json.dumps({'status':'PASS','checks':checks,'total_checks':sum(checks.values()),'literal_tournaments':literal,'residue_counts':residue_counts},indent=2,sort_keys=True))
if __name__=='__main__':main()
