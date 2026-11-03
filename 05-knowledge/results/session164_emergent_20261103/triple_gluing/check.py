#!/usr/bin/env python3
"""Exact independent controls for the S164 sign-interface/carry packet."""
from collections import Counter
from itertools import product, permutations
from math import gcd, isqrt
import hashlib, json
from pathlib import Path

checks=Counter()
def eq(x,y,kind):
    if x != y: raise RuntimeError(f'{kind}: {x!r} != {y!r}')
    checks[kind]+=1

def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
def power(a,n):
    z=(1,0)
    for _ in range(n): z=mul(z,a)
    return z

def primitive(z):
    g=gcd(*z)
    return ((z[0]//g,z[1]//g),g)

def triple(z):
    z,g=primitive(z)
    x,y=z
    return (abs(x*x-y*y),2*abs(x*y),x*x+y*y)

def rows(r):
    n=r+2
    out=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    out[0]^=1<<(n-1); out[-1]|=1
    return out

def directed(path,out): return all((out[a]>>b)&1 for a,b in zip(path,path[1:]))
def address_path(r,mask):
    return tuple([i+1 for i in range(r) if mask>>i&1]+[r+1,0]+[i+1 for i in range(r) if not(mask>>i&1)])
def brute(out):
    return {p for p in permutations(range(len(out))) if directed(p,out)}
def dp(out):
    n=len(out); d=[[0]*n for _ in range(1<<n)]
    for v in range(n): d[1<<v][v]=1
    for m in range(1,1<<n):
        for v in range(n):
            if not d[m][v]: continue
            choices=out[v]&~m
            while choices:
                bit=choices&-choices; choices-=bit
                d[m|bit][bit.bit_length()-1]+=d[m][v]
    return sum(d[-1])

for r in range(1,9):
    out=rows(r)
    eq(dp(out),1+2**r,'endpoint_DP')
    paths={address_path(r,m) for m in range(1<<r)}
    eq(len(paths),2**r,'address_injection')
    for p in paths: eq(directed(p,out),True,'literal_path')
    if r<=6: eq(paths|{tuple(range(r+2))},brute(out),'permutation_bijection')

for r,s in product(range(1,5),repeat=2):
    out=rows(r+s); image=set()
    for a,b in product(range(1<<r),range(1<<s)):
        pa,pb=address_path(r,a),address_path(s,b)
        # Literal amalgamation: identify the two endpoint pairs, retain
        # each interior order, and put first-factor interiors before second.
        left_a=pa[:pa.index(r+1)]; right_a=pa[pa.index(0)+1:]
        left_b=pb[:pb.index(s+1)]; right_b=pb[pb.index(0)+1:]
        joined=tuple(left_a)+tuple(r+x for x in left_b)+(r+s+1,0)+tuple(right_a)+tuple(r+x for x in right_b)
        eq(joined,address_path(r+s,a|(b<<r)),'glue_decode')
        eq(directed(joined,out),True,'glue_literal')
        image.add(joined)
    eq(len(image),(dp(rows(r))-1)*(dp(rows(s))-1),'glue_count')

P=(5,13,17); G=((2,1),(3,2),(4,1))
def alpha(v):
    z=(1,0)
    for e,p in zip(v,G): z=mul(z,power(p if e>=0 else conj(p),abs(e)))
    return z

def add(v,w): return tuple(x+y for x,y in zip(v,w))
def neg(v): return tuple(-x for x in v)
def norm(v):
    n=1
    for p,e in zip(P,v): n*=p**abs(e)
    return n

def carry(v,w):
    n=1
    for p,x,y in zip(P,v,w): n*=p**((abs(x)+abs(y)-abs(x+y))//2)
    return n

def canon(v): return min(v,neg(v))
def two(v,w): return [canon(add(v,w)),canon(add(v,neg(w)))]
V=list(product(range(-2,3),repeat=3))
for v,w in product(V,repeat=2):
    z=mul(alpha(v),alpha(w)); u=add(v,w); g=carry(v,w)
    eq(z,tuple(g*t for t in alpha(u)),'gaussian_carry')
    eq(g,gcd(*z),'carry_is_coordinate_gcd')
    eq(norm(v)*norm(w),g*g*norm(u),'norm_loss')
W=list(product(range(-1,2),repeat=3))
for u,v,w in product(W,repeat=3):
    eq(carry(u,v)*carry(add(u,v),w),carry(v,w)*carry(u,add(v,w)),'carry_cocycle')
    lhs=Counter(c for q in two(u,v) for c in two(q,w))
    rhs=Counter(c for q in two(v,w) for c in two(u,q))
    eq(lhs,rhs,'two_valued_associativity')

# Direct Euclid enumeration independently confirms fixed-norm fibers.
def euclid(C):
    ans=set()
    for n in range(2,isqrt(C)+1):
        m2=C-n*n
        if m2<=0: continue
        m=isqrt(m2)
        if m*m==m2 and 0<m<n and gcd(m,n)==1 and (m+n)%2:
            ans.add((n*n-m*m,2*n*m,C))
    return ans

def allocations(exps):
    positions=[i for i,e in enumerate(exps) if e]
    ans=[]
    for signs in product((-1,1),repeat=len(positions)):
        v=[0]*3
        for i,sign in zip(positions,signs): v[i]=exps[i]*sign
        ans.append(tuple(v))
    return ans

fiber_records=[]
for E,F in [((1,0,0),(2,0,0)),((1,1,0),(1,0,1)),((1,1,0),(1,1,0)),((2,1,0),(1,0,2)),((2,0,1),(2,0,1))]:
    seen=set(); got=Counter(); by_norm={}
    for v,w in product(allocations(E),allocations(F)):
        pair=(v,w); twin=(neg(v),neg(w)); key=min(pair,twin)
        if key in seen: continue
        seen.add(key); u=add(v,w); T=triple(alpha(u)); got[T]+=1
        by_norm.setdefault(T[2],set()).add(T)
    for C,ts in by_norm.items():
        if C==1: eq(ts,{(1,0,1)},'degenerate_fiber')
        else: eq(ts,euclid(C),'direct_euclid_fiber')
    for T,mult in got.items():
        # Unique absolute exponent branch is reconstructed from norm.
        z=0; C=T[2]
        for p,e,f in zip(P,E,F):
            if e==f and e>0 and C%p: z+=1
        predicted=2**(z-1) if C==1 else 2**z
        eq(mult,predicted,'cancelled_bit_multiplicity')
    fiber_records.append({'input_exponents':[E,F],'outputs':[(t,m) for t,m in sorted(got.items())]})

eq({t for t,_ in fiber_records[0]['outputs']},{(117,44,125),(3,4,5)},'p_p3_branch')
eq({t for t,m in fiber_records[1]['outputs'] if t[2]==5525},euclid(5525),'p2qr_branch')
eq({m for t,m in fiber_records[1]['outputs'] if t[2]==221},{2},'qr_double_fiber')

fibs=[0,1]
for _ in range(310): fibs.append(sum(fibs[-2:]))
fib_records=[]
for n in range(2,102,3):
    a=(fibs[n+1],fibs[n]); b=(fibs[n+2],fibs[n+1]); K=fibs[2*n+2]
    z=mul(a,conj(b)); C=fibs[2*n+1]; D=fibs[2*n+3]
    eq(z,(K,(-1)**(n+1)),'cassini_product')
    eq(gcd(C,D),1,'fibonacci_coprime_norms')
    eq((C%2,D%2),(1,1),'fibonacci_odd_norms')
    T=(K*K-1,2*K,K*K+1)
    eq(triple(z),T,'fibonacci_gap_two')
    A,B,H=T
    labels=(2*(H*H-4*A*B),2*H*H,2*(H*H+4*A*B))
    eq(labels[0]>0,True,'positive_AP')
    eq(labels[0]+labels[2],2*labels[1],'AP_identity')
    eq(tuple(labels[i]+labels[j] for i,j in [(0,1),(0,2),(1,2)]),(4*(A-B)**2,4*H*H,4*(A+B)**2),'square_sum_edges')
    if n<=8: fib_records.append({'n':n,'input_norms':[C,D],'K':K,'triple':T,'AP_labels':labels})

eq(triple(mul((2,1),(3,2))),(33,56,65),'missing_bit_hostile_plus')
eq(triple(mul((2,1),(3,-2))),(63,16,65),'missing_bit_hostile_minus')
# Canonical positive-even-leg chart denominator |n-m|.
guards=[]
for e in (1,3,4):
    x,y=power((2,1),e); n,m=sorted((abs(x),abs(y)),reverse=True)
    guards.append((5**e,n+m,n-m))
eq(guards,[(5,3,1),(125,13,9),(625,31,17)],'ternary_guard_hostile')
result={'status':'PASS','checks':dict(sorted(checks.items())),'total_checks':sum(checks.values()),'fibers':fiber_records,'fibonacci_samples':fib_records,'ternary_guard_hostile':guards,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(result,indent=2,sort_keys=True))
