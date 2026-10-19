#!/usr/bin/env python3
"""Exact small controls for the relation-quotient/Golay bridge, no Collatz census."""
from collections import Counter
from itertools import combinations, product
from math import gcd
import json

CHECKS = 0
def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)

# F_243 = F_3[X]/(X^5+X^4+2X^3+X^2+2), the standard Golay factor.
G = (2, 0, 1, 2, 1, 1)
ZERO = (0,)*5
ONE = (1,0,0,0,0)
ALPHA = (0,1,0,0,0)
def add(a,b): return tuple((x+y)%3 for x,y in zip(a,b))
def scale(c,a): return tuple(c*x%3 for x in a)
def mul(a,b):
    v=[0]*9
    for i,x in enumerate(a):
        for j,y in enumerate(b): v[i+j]=(v[i+j]+x*y)%3
    for i in range(8,4,-1):
        c=v[i]
        for j in range(6): v[i-5+j]=(v[i-5+j]-c*G[j])%3
    return tuple(v[:5])
def power(a,n):
    y=ONE
    while n:
        if n&1: y=mul(y,a)
        a=mul(a,a); n//=2
    return y
F=list(product(range(3),repeat=5))
for z in F:
    require(power(z,243)==z,'Frobenius fifth power')
    if z!=ZERO: require(power(z,242)==ONE,'all nonzero elements invertible')
AP=[power(ALPHA,i) for i in range(11)]
require(power(ALPHA,11)==ONE and ALPHA!=ONE,'alpha order eleven')
def syndrome(v,minus=False):
    y=ZERO
    for i,c in enumerate(v): y=add(y,scale(c,AP[-i%11] if minus else AP[i]))
    return y

def act(v,b,j):
    w=[0]*11
    for i,c in enumerate(v): w[(pow(3,j,11)*i+b)%11]=c
    return tuple(w)
def reflect(v): return tuple(v[-i%11] for i in range(11))

errors=[(0,)*11]
for weight in (1,2):
    for support in combinations(range(11),weight):
        for coeffs in product((1,2),repeat=weight):
            v=[0]*11
            for i,c in zip(support,coeffs): v[i]=c
            errors.append(tuple(v))
syndromes={syndrome(v):v for v in errors}
require(len(errors)==len(syndromes)==243,'perfect radius-two syndrome bijection')
for v in errors:
    s=syndrome(v)
    require(syndrome(reflect(v))==syndrome(v,True),'reflection swaps syndromes')
    for b in range(11):
        for j in range(5):
            require(syndrome(act(v,b,j))==mul(AP[b],power(s,3**j)),
                    'Golay syndrome equivariance')

def orbit_sizes(universe, generators):
    unseen=set(universe); sizes=[]
    while unseen:
        start=min(unseen); seen={start}; todo=[start]
        while todo:
            v=todo.pop()
            for gen in generators:
                w=gen(v)
                if w not in seen: seen.add(w);todo.append(w)
        require(seen<=unseen,'orbit partition disjoint')
        unseen-=seen;sizes.append(len(seen))
    return dict(sorted(Counter(sizes).items()))
field_orbits=orbit_sizes(F,[lambda z:mul(ALPHA,z),lambda z:power(z,3)])
error_orbits=orbit_sizes(errors,[lambda v:act(v,1,0),lambda v:act(v,0,1)])
require(field_orbits==error_orbits=={1:1,11:2,55:4},'seven matching orbits')

code=[]
for h in product(range(3),repeat=6):
    v=[0]*11
    for i,x in enumerate(h):
        for j,y in enumerate(G):v[i+j]=(v[i+j]+x*y)%3
    v=tuple(v);code.append(v)
    require(syndrome(v)==ZERO,'code in syndrome kernel')
weights=dict(sorted(Counter(sum(c!=0 for c in v) for v in code).items()))
require(len(set(code))==729 and min(k for k in weights if k)>0,'code dimension six')
require(min(k for k in weights if k)==5,'independent minimum distance five')
reflected_fail=[v for v in code if syndrome(reflect(v))!=ZERO]
require(len(reflected_fail)==726,'only repetition code common to reciprocal codes')
reflection_hostile=list(min(reflected_fail))

# Rank of complete finite Fourier transform: one scalar plus two syndromes.
def rank_mod3(rows):
    a=[list(r) for r in rows]; rank=0
    for col in range(len(a[0])):
        piv=next((i for i in range(rank,len(a)) if a[i][col]%3),None)
        if piv is None:continue
        a[rank],a[piv]=a[piv],a[rank]
        inv=pow(a[rank][col],-1,3)
        a[rank]=[x*inv%3 for x in a[rank]]
        for i in range(len(a)):
            if i!=rank:
                c=a[i][col]
                a[i]=[(x-c*y)%3 for x,y in zip(a[i],a[rank])]
        rank+=1
    return rank
cols=[]
for i in range(11):
    v=tuple(int(j==i) for j in range(11))
    cols.append((1,)+syndrome(v)+syndrome(v,True))
require(rank_mod3(list(zip(*cols)))==11,'two-sheet Fourier map is invertible')
reciprocal=tuple(2*c%3 for c in reversed(G))
require(reciprocal!=G,'inverse multiplication has distinct characteristic polynomial')

# Decoder-selected reflection exists as a nonlinear permutation on the same field.
def jdec(z):return syndrome(reflect(syndromes[z]))
for z in F:
    require(jdec(jdec(z))==z,'decoder reflection involution')
    require(jdec(mul(ALPHA,z))==mul(AP[10],jdec(z)),'decoder reflection reverses coordinate translation')
    require(jdec(power(z,3))==power(jdec(z),3),'decoder reflection commutes with Frobenius')
z=add(scale(2,ONE),AP[2]); w=scale(2,AP[3])
nonlinear_left=jdec(add(z,w)); nonlinear_right=add(jdec(z),jdec(w))
require(nonlinear_left!=nonlinear_right,'explicit decoder nonadditivity')
nonlinear_witness={'z':z,'w':w,'J_sum':nonlinear_left,'sum_J':nonlinear_right}

# General signed period quotients, exact ring relation, and golden common quotients.
def ord_mod(a,m):
    if m==1:return 1
    require(gcd(a,m)==1,'unit for multiplicative order')
    x=1
    for d in range(1,m+1):
        x=x*a%m
        if x==1:return d
    raise RuntimeError('order missing')
periods=[]
for k in range(1,9):
    for eps in (1,-1):
        m=3**k-eps
        require(pow(3,k,m)==eps%m,'signed relation translation annihilator')
        if eps==-1:require(ord_mod(3,m)==2*k,'signed faithful complement')
        else:require(ord_mod(3,m)==k,'unsigned faithful complement')
        periods.append({'k':k,'sign':eps,'translation_modulus':m,'dilation_order':ord_mod(3,m)})
common=[]
for r in range(1,9):
    a=3**r; f=a*a-a-1
    for k in range(1,41):
        d=gcd(f,a**k-1)
        require(f%d==0 and (a**k-1)%d==0,'universal golden common quotient')
        if d>1: common.append({'r':r,'k':k,'modulus':d})
require(gcd(3**8-3**4-1,3**20-1)==11,'eleven common quotient')
require(gcd(3**8-3**4-1,3**36-1)==19,'nineteen common quotient')
require(gcd(3**8-3**4-1,3**60-1)==341,'composite common quotient')

# Hensel lifts and all orbit strata at 11^a. Bounds are finite controls.
def golden_root(a):
    r=4; m=11
    for _ in range(1,a):
        cand=[r+j*m for j in range(11) if ((r+j*m)**2-(r+j*m)-1)%(m*11)==0]
        require(len(cand)==1,'unique golden Hensel root')
        r=cand[0];m*=11
    return r
lift=[]
for a in range(1,5):
    m=11**a;r=golden_root(a)
    ot=orbit_sizes(range(m),[lambda x,m=m:3*x%m])
    og=orbit_sizes(range(m),[lambda x,m=m,r=r:r*x%m])
    nt=sum(ot.values())-1;ng=sum(og.values())-1
    require(nt==22*a-20 and ng==2*a,'tripling/golden orbit multiplicity')
    require(ord_mod(3,m)==5*11**max(0,a-2),'tripling delayed order lift')
    require(ord_mod(r,m)==5*11**(a-1),'golden ordinary order lift')
    lift.append({'a':a,'modulus':m,'golden_root':r,'tripling_orbits':ot,'golden_orbits':og})
require(golden_root(2)==37 and pow(37,5,121)==67 and pow(37,44,121)==3,'121 lift separation')
golden10=orbit_sizes(list(product(range(11),repeat=2)),[lambda xy:(4*xy[0]%11,8*xy[1]%11)])
require(golden10=={1:1,5:2,10:11},'golden R10 orbit partition')

# Literal integer reconstruction controls: neither congruence nor a bad clock suffices.
def T(n):return (3*n+1)//2 if n%2 else n//2
def iterate(n,k):
    for _ in range(k):n=T(n)
    return n
require(iterate(223,7)==iterate(233,15)==425,'retained asynchronous source example')
require(iterate(334,6)==iterate(335,6)==425,'consecutive merge')
require(1%11==12%11 and 1!=12,'congruence alone hostile')
require(all(iterate(1,j)!=iterate(2,j) for j in range(4)) and T(T(1))==1 and T(T(2))==2,
        'closed out-of-phase hostile')

report={'status':'FINITE-EXACT controls; proofs and non-consequences in companion note',
        'checks':CHECKS,'golay_generator_low_to_high':G,'reciprocal_generator_low_to_high':reciprocal,
        'syndrome_count':len(syndromes),'field_group_orbits':field_orbits,'error_group_orbits':error_orbits,
        'code_weight_distribution':weights,'reflection_hostile_codeword':reflection_hostile,'decoder_nonlinearity_witness':nonlinear_witness,
        'signed_period_quotients':periods,'common_golden_quotients':common,
        'prime_power_lifts':lift,'golden_R10_orbits':golden10}
print(json.dumps(report,indent=2,sort_keys=True))
