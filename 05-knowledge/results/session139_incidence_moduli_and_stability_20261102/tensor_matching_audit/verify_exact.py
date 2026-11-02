#!/usr/bin/env python3
"""Independent exact tensor/Pfaffian audit. Controls are not a proof."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product, combinations
from random import Random
import json

rng=Random(139028)
counts={}
mins={}

def require(c,msg):
    if not c: raise RuntimeError(msg)

def record(key,gap):
    require(gap>=0,(key,gap))
    counts[key]=counts.get(key,0)+1
    mins[key]=min(mins.get(key,gap),gap)

def det(a):
    a=[list(map(F,row)) for row in a]; n=len(a); result=F(1)
    for k in range(n):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None: return F(0)
        if pivot!=k: a[k],a[pivot]=a[pivot],a[k]; result=-result
        p=a[k][k]; result*=p
        for i in range(k+1,n):
            if not a[i][k]: continue
            c=a[i][k]/p
            for j in range(k+1,n): a[i][j]-=c*a[k][j]
    return result

def uv(a):
    n=len(a)
    return det([[a[i][j]+(i==j) for j in range(n)] for i in range(n)]),(-1)**n*det(a)

def k_to_a(k):
    n=len(k)
    return [[F(0) if i==j else (1+k[i][j])/2 for j in range(n)] for i in range(n)]

def pf_factory(k):
    @lru_cache(None)
    def pf(ids):
        if not ids:return F(1)
        require(len(ids)%2==0,'odd Pf')
        return sum(((-1)**(j+1)*k[ids[0]][ids[j]]*pf(ids[1:j]+ids[j+1:]) for j in range(1,len(ids))),F(0))
    return pf

def tensor(b,x):
    r=len(x); f=[]
    for s in range(1<<r):
        value=F(0); j=s
        while True:
            value+=b[s^j]*prod(x[i] for i in range(r) if j>>i&1)
            if not j:break
            j=(j-1)&s
        f.append(value)
    return sum(v*v for v in f)

def prod(xs):
    result=F(1)
    for x in xs:result*=x
    return result

# Full low-dimensional lattice controls, including boundary x=+/-1.
for r in range(3):
    for b0 in [-3,-1,1,3]:
        for tail in product([-4,-2,0,2,4],repeat=(1<<r)-1):
            b=(b0,)+tail
            for x in product([F(-1),F(-1,2),F(0),F(1,2),F(1)],repeat=r):
                record('tensor_exhaustive_r'+str(r),tensor(b,x)-prod(1+y*y for y in x))
for r in [3,4,5]:
    for trial in range(100):
        b=[2*rng.randrange(-10,11)+1]+[2*rng.randrange(-10,11) for _ in range((1<<r)-1)]
        x=[rng.choice([F(-1),F(-7,8),F(-1,3),F(0),F(2,5),F(1)]) for _ in range(r)]
        record('tensor_random_r'+str(r),tensor(b,x)-prod(1+y*y for y in x))

def matrix(n,r,bits,x):
    k=[[F(0) for _ in range(n)] for _ in range(n)]
    pairs={(2*i,2*i+1):x[i] for i in range(r)}
    h=0
    for i in range(n):
        for j in range(i+1,n):
            value=pairs.get((i,j))
            if value is None:value=F(1 if bits>>h&1 else -1);h+=1
            k[i][j]=value;k[j][i]=-value
    return k

def check(k,r,do_reserve=True,do_groups=False):
    n=len(k); t=n-2*r; x=[k[2*i][2*i+1] for i in range(r)]
    q=[(1-y*y)/4 for y in x]; main=prod(1-y for y in q); err=prod(q) if t==0 else 0
    a=k_to_a(k);u,v=uv(a)
    record('det_margin_n'+str(n),u-abs(v)-main+err)
    observer=[[F(0) for _ in range(n+1)] for _ in range(n+1)]
    for i in range(n):
        for j in range(n):observer[i][j]=k[i][j]
        observer[i][n]=F(1);observer[n][i]=F(-1)
    pf=pf_factory(observer)
    even=sum((pf(s)**2 for z in range(0,n+1,2) for s in combinations(range(n),z)),F(0))
    odd=sum((pf(s+(n,))**2 for z in range(1,n+1,2) for s in combinations(range(n),z)),F(0))
    require(u+v==F(2,2**n)*even,('even normalization',n))
    require(u-v==F(2,2**n)*odd,('odd normalization',n))
    counts['pf_normalization']=counts.get('pf_normalization',0)+1
    signed_err=(-1)**r*err
    record('signed_even',u+v-main-signed_err)
    record('signed_odd',u-v-main+signed_err)
    if do_groups:
        # Partition by singleton endpoint choices and isolate choices; even base + observer.
        for choices in product([None,0,1],repeat=r):
            free=[i for i in range(r) if choices[i] is None]
            selected=[2*i+c for i,c in enumerate(choices) if c is not None]
            for isolate_bits in range(1<<t):
                base=selected+[2*r+j for j in range(t) if isolate_bits>>j&1]
                if len(base)%2:base.append(n)
                f=[]; b=[]
                for s in range(1<<len(free)):
                    ids=tuple(sorted(base+[vtx for j,i in enumerate(free) if s>>j&1 for vtx in (2*i,2*i+1)]))
                    f.append(pf(ids))
                    zero=[row[:] for row in observer]
                    for j,i in enumerate(free):
                        if s>>j&1:zero[2*i][2*i+1]=zero[2*i+1][2*i]=0
                    b.append(pf_factory(zero)(ids))
                require(b[0].denominator==1 and b[0].numerator%2==1,('base odd',n,base))
                require(all(v.denominator==1 and v.numerator%2==0 for v in b[1:]),('tail even',n,base,b))
                transformed=[]
                for s in range(1<<len(free)):
                    value=F(0); j=s
                    while True:
                        value+=b[s^j]*prod(x[free[h]] for h in range(len(free)) if j>>h&1)
                        if not j:break
                        j=(j-1)&s
                    transformed.append(value)
                require(f==transformed,('Pf transform',n,base))
                record('group_floor',sum(v*v for v in f)-prod(1+x[i]*x[i] for i in free))
    if do_reserve:
        expect=F(0)
        for bits in range(1<<r):
            kb=[row[:] for row in k];weight=F(1)
            for i in range(r):
                p=(1+x[i])/2;bit=bits>>i&1;weight*=p if bit else 1-p
                kb[2*i][2*i+1]=1 if bit else -1;kb[2*i+1][2*i]=-kb[2*i][2*i+1]
            ub,vb=uv(k_to_a(kb));expect+=weight*(ub+abs(vb))
        reserve=1-main-err
        record('nonlinear_reserve_n'+str(n),expect-(u+abs(v))-reserve)

# Exhaustive arbitrary cross-orientations in order 4, at mixed rational probabilities.
for n,r in [(2,1),(3,1),(4,2)]:
    binary_edges=n*(n-1)//2-r
    for bits in range(1<<binary_edges):
        for x in product([F(-3,4),F(-1,3),F(0),F(1,2)],repeat=r):
            check(matrix(n,r,bits,x),r,do_groups=True)
# Larger rational controls, r=0 pure binary included; no finite extrapolation.
for n in range(1,10):
    for r in range(n//2+1):
        for trial in range(20):
            x=[rng.choice([F(-7,8),F(-2,3),F(0),F(1,5),F(4,5)]) for _ in range(r)]
            k=matrix(n,r,rng.getrandbits(n*(n-1)//2-r),x)
            check(k,r,do_reserve=True,do_groups=(trial==0))
# Ordinal ordered pair blocks are a sharpness control for both determinant floor and reserve.
for n in range(2,10):
    for r in range(1,n//2+1):
        x=[F(2*i+1,2*r+1) for i in range(r)]
        k=matrix(n,r,(1<<(n*(n-1)//2-r))-1,x)
        check(k,r,do_groups=True)
        u,v=uv(k_to_a(k));q=[(1-y*y)/4 for y in x];floor=prod(1-z for z in q)-(prod(q) if n==2*r else 0)
        require(u-abs(v)==floor,('ordinal equality',n,r,u,v,floor))
        counts['ordinal_equalities']=counts.get('ordinal_equalities',0)+1

out={'status':'all exact controls passed; proof audit is separate','counts':counts,'minimum_gaps':{k:str(v) for k,v in mins.items()}}
print(json.dumps(out,indent=2,sort_keys=True))
