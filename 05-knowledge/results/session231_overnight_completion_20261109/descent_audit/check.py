#!/usr/bin/env python3
"""Independent verifier: no imports from the producer or previous packets."""
import itertools as it
import json
from functools import lru_cache
from math import gcd,isqrt

def need(x, label):
    if not x: raise RuntimeError(label)

def tournaments(n):
    edges=list(it.combinations(range(n),2))
    for bits in it.product((False,True),repeat=len(edges)):
        a=[[False]*n for _ in range(n)]
        for (i,j),b in zip(edges,bits): a[i][j]=b; a[j][i]=not b
        yield tuple(tuple(r) for r in a)

def profile(a):
    n=len(a)
    def path(p): return all(a[i][j] for i,j in zip(p,p[1:]))
    h=e=u=0
    for p in it.permutations(range(n)):
        if path(p): h+=1; e+=p[-1]==0
        cut=p.index(0)
        u+=path(p[:cut]) and path(p[cut:])
    d=sum(path(p) for p in it.permutations(range(1,n)))
    return h,e,u,d

def count(a):
    @lru_cache(None)
    def suffix(v,left):
        if not left: return 1
        return sum(suffix(w,left-(1<<w)) for w in range(len(a)) if left>>w&1 and a[v][w])
    return sum(suffix(v,((1<<len(a))-1)-(1<<v)) for v in range(len(a)))

def join(a,c):
    n,m=len(a),len(c)
    b=[[False]*(n+m) for _ in range(n+m)]
    for i in range(n):
        for j in range(n): b[i][j]=a[i][j]
    for i in range(m):
        for j in range(m): b[n+i][n+j]=c[j][i]
    for i in range(n):
        for j in range(m): b[i][n+j]=True
    b[0][n]=False;b[n][0]=True
    return tuple(tuple(r) for r in b)

def reach(a,v):
    s={v}
    while True:
        t=s|{w for u in s for w in range(len(a)) if a[u][w]}
        if t==s:return s
        s=t

def normalize(a):
    kept=sorted(reach(a,0));removed=[v for v in range(len(a)) if v not in kept]
    small=tuple(tuple(a[i][j] for j in kept) for i in kept)
    prefix=tuple(tuple(a[i][j] for j in removed) for i in removed)
    return small,count(prefix) if prefix else 1

def value(x,y):return x[0]*y[0]-x[1]*y[1]+x[2]*y[2]
def proper(x,y):return x[1]<x[0] or y[1]<y[0]
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def inverse(N,bank):
    result=set();cap=(2*N-1)//3
    for h in range(1,cap+1,2):
      for k in range(1,cap+1,2):
       for e in range(h+1):
        for f in range(k+1):
         if e==h and f==k:continue
         t=N-h*k+e*f
         if not 0<t<N:continue
         for u in range(1,isqrt(t)+1):
          if t%u:continue
          for r,s in {(u,t//u),(t//u,u)}:
           if r<h or s<k:continue
           x=(h,e,r);y=(k,f,s)
           if x in bank and y in bank:result.add((x,y))
    return result

def main():
    graphs=[a for n in range(1,6) for a in tournaments(n)]
    states={a:profile(a) for a in graphs}
    norm={a:normalize(a) for a in graphs}
    for a,x in states.items():
        h,e,u,d=x;b,p=norm[a]
        need(h==count(a),'independent Hamilton count')
        need(h%2==1 and d%2==1 and 1<=d<=h<=u,'basic bounds')
        need((e==h)==(not any(a[0])),'sink equivalence')
        if e<h:need(2*e<=h-1,'endpoint bound')
        need(x==tuple(p*z for z in states[b]),'simultaneous normalization')
        need(bool(any(a[0]))==bool(any(b[0])),'properness preservation')
    proper_count=literal=strong_tests=prime_tests=0
    for a,x in states.items():
      for c,y in states.items():
        N=value(x,y);h,e,u,_=x;k,f,v,_=y
        if proper(x,y):
            proper_count+=1
            need(3*h<=2*N-1 and 3*k<=2*N-1,'descent')
            need(u*v<=N-1,'response bound')
            if e<h:need(2*N>=(3*h+1)*k,'left bound')
            if f<k:need(2*N>=(3*k+1)*h,'right bound')
            if prime(N):
                aa,p=norm[a];cc,q=norm[c]
                need(p==q==1,'prime prefix factors')
                need(gcd(gcd(h,e),u)==gcd(gcd(k,f),v)==1,'primitive prime triple')
                prime_tests+=1
        else:need(N==h*k,'degenerate product')
        if len(a)<=4 and len(c)<=4:
            j=join(a,c);need(count(j)==N,'join law');literal+=1
            if proper(x,y):
                strong=all(len(reach(j,w))==len(j) for w in range(len(j)))
                need(strong==(len(reach(a,0))==len(a) and len(reach(c,0))==len(c)),'strongness iff')
                strong_tests+=1
    bank={x[:3] for x in states.values()};inverse_checks=[]
    for N in range(3,42,2):
        expected={(x,y) for x in bank for y in bank if proper(x,y) and value(x,y)==N}
        actual=inverse(N,bank);need(expected==actual,'divisor inverse');inverse_checks.append([N,len(actual)])
    # Equality witnesses make the leading descent constant and product slack real.
    cycle=((False,True,False),(False,False,True),(True,False,False))
    one=((False,),);tt2=((False,True),(False,False))
    need(value(profile(cycle),profile(one))==5,'descent sharp witness')
    need(value(profile(tt2),profile(one))==3,'response sharp witness')
    return dict(status='PASS',universe='all labelled root-zero tournaments orders 1..5; literal joins and strongness orders 1..4',graphs=len(graphs),physical_triples=len(bank),pair_universe=len(graphs)**2,proper_pairs=proper_count,literal_join_checks=literal,strongness_iff_checks=strong_tests,prime_pair_checks=prime_tests,inverse_bounded_bank=inverse_checks,full_oracle_implemented=False,scope='finite exact independent controls; general claims require proof')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
