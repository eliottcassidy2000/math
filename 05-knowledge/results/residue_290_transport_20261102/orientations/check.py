#!/usr/bin/env python3
"""Fresh exact controls for signed wheel/trace transport; stdlib only."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import gcd, log, expm1
import json

GATES=0
def require(ok, label):
    global GATES
    GATES += 1
    if not ok: raise RuntimeError(label)
def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0: p.pop()
    return p
def add(a,b):
    c=[0]*max(len(a),len(b))
    for i,x in enumerate(a): c[i]+=x
    for i,x in enumerate(b): c[i]+=x
    return trim(c)
def scale(a,b): return trim([b*x for x in a])
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return trim(c)
def divmod_poly(a,b):
    a=trim(a); b=trim(b)
    if len(a)<len(b): return [0],a
    q=[0]*(len(a)-len(b)+1)
    while a!=[0] and len(a)>=len(b):
        k=len(a)-len(b); c=F(a[-1],b[-1]); q[k]=c
        for i,x in enumerate(b): a[k+i]-=c*x
        a=trim(a)
    return trim(q),a
def mod(a,b): return divmod_poly(a,b)[1]
@lru_cache(None)
def cyclo(n):
    p=[-1]+[0]*(n-1)+[1]
    for d in range(1,n):
        if n%d==0:
            p,r=divmod_poly(p,cyclo(d)); require(r==[0],f'cyclotomic division {n},{d}')
    require(all(F(x).denominator==1 for x in p),'integral cyclotomic')
    return tuple(int(x) for x in p)
def trace_poly(n):
    phi=cyclo(n); degree=(len(phi)-1)//2
    require(list(phi)==list(reversed(phi)),'reciprocal cyclotomic')
    s0,s1=[2],[0,1]; r=[phi[degree]]
    for j in range(1,degree+1):
        sj=s1 if j==1 else add(mul([0,1],s1),scale(s0,-1))
        if j>1: s0,s1=s1,sj
        r=add(r,scale(sj,phi[degree+j]))
    require(len(r)==degree+1 and r[-1]==1,'trace degree')
    return r

def orbits(points,step):
    todo=set(points); out=[]
    while todo:
        x=min(todo); cyc=[]; y=x
        while y not in cyc:
            require(y in todo,'disjoint permutation orbit')
            cyc.append(y); y=step(y)
        require(y==x,'permutation cycle closes')
        todo.difference_update(cyc); out.append(cyc)
    return out

def wheel(n):
    units=[r for r in range(1,n) if gcd(r,n)==1]
    cyc=orbits(units,lambda r:2*r%n); d=len(cyc[0])
    require(all(len(c)==d for c in cyc),'equal unit orbit lengths')
    half=next((k for k in range(1,d+1) if pow(2,k,n)==n-1),None)
    folded=orbits({min(r,n-r) for r in units},lambda r:min(2*r%n,(-2*r)%n))
    ell=d//2 if half is not None else d
    require(all(len(c)==ell for c in folded),'folded periods')
    if half is not None: require(half*2==d,'negative halfway')
    for c in folded:
        # Lift each chosen quotient section and multiply the sign on every edge.
        sign=1
        for i,r in enumerate(c):
            y=2*r%n; q=c[(i+1)%len(c)]
            require(y in (q,n-q),'quotient edge lifts')
            sign*=1 if y==q else -1
        require(sign==(-1 if half is not None else 1),'odd monodromy')
    return {'N':n,'phi':len(units),'unit_order':d,'unit_cycles':len(cyc),
            'folded_period':ell,'folded_cycles':len(folded),'odd_sign':-1 if half else 1}

all_wheels=[wheel(n) for n in range(3,302,2)]
selected=[x for x in all_wheels if x['N'] in (5,7,9,29,145)]
u290=[r for r in range(1,290) if gcd(r,290)==1]
require(len(u290)==112,'U290 size')
for r in u290:
    require(gcd(147*r%290,290)==1,'odd-lift unit transport')
    require((147*r%290)%145==2*(r%145)%145,'CRT doubling')
    require(pow(147,14,290)*r%290==(-r)%290,'halfway inversion')
    for s in u290:
        require((r*s%290)%145==(r%145)*(s%145)%145,'CRT multiplicativity')
# Character enumeration in explicit primitive-root coordinates at 5 and29.
require(len({pow(2,i,5) for i in range(4)})==4,'2 primitive mod5')
require(len({pow(2,i,29) for i in range(28)})==28,'2 primitive mod29')
char_exponents={par:Counter() for par in (0,1)}
for a in range(4):
    for b in range(28):
        par=(a+b)%2; e=(7*a+b)%28
        require(e%2==par,'character monodromy parity')
        char_exponents[par][e]+=1
for par in (0,1):
    require(char_exponents[par]==Counter({e:4 for e in range(par,28,2)}),'fourfold signed spectrum')

traces=[]
for n in (5,7,9,29,145):
    w=next(x for x in selected if x['N']==n); r=trace_poly(n)
    x=[0,1]; derivative=[1]
    for k in range(1,w['folded_period']+1):
        derivative=mod(mul(derivative,scale(x,2)),r)
        x=mod(add(mul(x,x),[-2]),r)
        require((x==[0,1])==(k==w['folded_period']),f'exact trace return {n},{k}')
    expected=w['odd_sign']*2**w['folded_period']
    require(derivative==[expected],f'trace differential sign {n}')
    traces.append({'N':n,'degree':len(r)-1,'coefficients':[int(v) for v in r],
                   'period':w['folded_period'],'multiplier':expected})
require(trace_poly(7)==[-1,-2,1,1],'R7')
require(trace_poly(9)==[1,-3,0,1],'R9')
# Exact chosen-cycle reversal and derivative products.
for n,c,multiplier in ((7,F(-7,4),F(1)),(9,F(-11,4),F(19))):
    r=trace_poly(n); l=[F(-1,2),-1]; x=l; derivative=[1]
    for k in range(3):
        derivative=mod(mul(derivative,scale(x,2)),r)
        x=mod(add(mul(x,x),[c]),r)
    require(x==l,'reversed quadratic cycle')
    require(derivative==[multiplier],'quadratic multiplier')
    d2=add(mul([-2,0,1],[-2,0,1]),[-2])
    diff=add(add(mul(l,l),[c]),add(d2,[F(1,2)]))
    require(mod(diff,r)==[0],'L reversal identity')
cycle=[F(-7,4),F(5,4),F(-1,4)]; derivative=F(1)
for i,x in enumerate(cycle):
    require(x*x-F(29,16)==cycle[(i+1)%3],'rational 29 cycle')
    derivative*=2*x
require(derivative==F(35,8),'rational cycle multiplier')
require(-((F(-1,2))**2+7)/4==F(-29,16),'fixed marked trace parameter')

mask=[int(n%290==1)-int(n%290==289) for n in range(1,291)]
require(sum(mask)==0,'signed mean')
require(F(sum(-int(gcd(n,290)==1) for n in range(1,291)),580)==F(-28,145),'unit divergence coefficient')
require(F(sum(-int(n%290 in (1,289)) for n in range(1,291)),580)==F(-1,290),'target divergence coefficient')
prefix=[0]
for x in mask: prefix.append(prefix[-1]+x)
require(set(prefix)=={0,1},'exact bounded prefix range')
for start in range(290):
    total=0; values={0}
    for k in range(580):
        total+=mask[(start+k)%290]; values.add(total)
    require(max(values)-min(values)==1,'every offset tail Abel bound')
# Illustrative floats only, excluded from the exact gate count.
numeric=[]
for K in (1,10,100,1000):
    total=sum(expm1(log(n)/n)*(1 if n%290==1 else -1)
              for j in range(K) for n in (290*j+1,290*j+289))
    numeric.append({'block_cut':290*K,'sum_illustrative':round(total,12),
                    'proved_tail_upper_evaluated_approximately':round(expm1(log(290*K+1)/(290*K+1)),12)})
print(json.dumps({'status':'PASS','gates':GATES,'universe':{'odd_moduli':[3,301],
                 'trace_conductors':[5,7,9,29,145],'all_unit_pairs_290':len(u290)**2,
                 'all_tail_offsets':290},'selected_wheels':selected,'traces':traces,
                 'mrb_numeric_not_exact':numeric},indent=2,sort_keys=True))
