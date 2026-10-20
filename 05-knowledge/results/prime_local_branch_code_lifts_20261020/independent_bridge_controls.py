#!/usr/bin/env python3
"""Independent linear-algebra, ring-enumeration, and CRT controls."""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import comb,factorial,gcd,lcm
import json

checks={}
def ck(ok,tag):
    checks[tag]=checks.get(tag,0)+1
    if not ok:raise RuntimeError(tag)

def rank(A,p):
    A=[list(row) for row in A]
    if not A:return 0
    r=0
    for j in range(len(A[0])):
        at=next((i for i in range(r,len(A)) if A[i][j]%p),None)
        if at is None:continue
        A[r],A[at]=A[at],A[r];u=pow(A[r][j]%p,-1,p)
        A[r]=[z*u%p for z in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][j]%p:
                u=A[i][j]%p;A[i]=[(a-u*b)%p for a,b in zip(A[i],A[r])]
        r+=1
        if r==len(A):break
    return r

def parts(n,a=1):
    if n==0:yield ()
    for b in range(a,n+1):
        for tail in parts(n-b,b):yield (b,)+tail

def perm(mu):
    out=[];j=0
    for a in mu:out.extend(list(range(j+1,j+a))+[j]);j+=a
    return out

def exponent(mu):return sum((a-1)//2 for a in mu)+sum(gcd(a,b) for i,a in enumerate(mu) for b in mu[i+1:])-len(mu)+1

def z(mu):
    v=1
    for a,c in Counter(mu).items():v*=a**c*factorial(c)
    return v

def ordperm(g):
    seen=set();ans=1
    for i in range(len(g)):
        if i in seen:continue
        j=i;a=0
        while j not in seen:seen.add(j);a+=1;j=g[j]
        ans=lcm(ans,a)
    return ans

def matrix(g):
    n=len(g);pairs=list((i,j) for i in range(1,n) for j in range(i+1,n));index={a:i for i,a in enumerate(pairs)}
    M=[]
    for i,j in pairs:
        row=[0]*len(pairs)
        for a,b,s in [(g[i],g[j],1),(g[0],g[j],-1),(g[0],g[i],1)]:
            if a==0 or b==0:continue
            if a>b:a,b,s=b,a,-s
            row[index[a,b]]+=s
        M.append(row)
    return M

for p in [2,3,5,7,11]:
    for n in range(1,10):
        for mu in parts(n):
            if any(a%p==0 for a in mu):continue
            M=matrix(perm(mu));d=len(M)
            ck(d-rank([[(x-(i==j))%p for j,x in enumerate(row)] for i,row in enumerate(M)],p)==exponent(mu),'regular_fixed_rank')
    for s in [1,2,3]:
        for n in range(1,19):
            val=sum((Fraction(p**(s*exponent(mu)),z(mu)) for mu in parts(n) if all(a%p for a in mu)),Fraction())
            ck(val==Fraction(1,2) if p==2 and n==2 else val.denominator==1,'prime_power_branch_integrality')

def affine_offset(g):
    n=len(g)
    # Zero canonical coordinates: a_ij=0 for i<j and 1 for i>j.
    def a(i,j):return int(i>j)
    return [a(g[i],g[j])-a(g[0],g[j])+a(g[0],g[i]) for i in range(1,n) for j in range(i+1,n)]

ring=[]
for p,s,n in [(2,2,3),(2,3,3),(2,2,4),(2,3,4),(3,2,3),(3,2,4),(5,2,3)]:
    mod=p**s;d=(n-1)*(n-2)//2;states=list(product(range(mod),repeat=d));total=reg=afftotal=affreg=0
    for g in permutations(range(n)):
        M=matrix(g)
        f=sum(all(sum(a*b for a,b in zip(row,x))%mod==x[i] for i,row in enumerate(M)) for x in states)
        total+=f
        offset=affine_offset(g)
        af=sum(all((sum(a*b for a,b in zip(row,x))+offset[i])%mod==x[i] for i,row in enumerate(M)) for x in states)
        afftotal+=af
        if ordperm(g)%p:
            reg+=f;affreg+=af
            ck(af==f,'regular_affine_linear_fixed_equality')
    predicted=sum((Fraction(mod**exponent(mu),z(mu)) for mu in parts(n) if all(a%p for a in mu)),Fraction())
    ck(Fraction(reg,factorial(n))==predicted,'direct_ring_regular_branch')
    ck(total%factorial(n)==0,'direct_ring_burnside')
    ck(afftotal%factorial(n)==0 and Fraction(affreg,factorial(n))==predicted,'direct_ring_affine_branch')
    ring.append({'p':p,'s':s,'n':n,'states':len(states),'regular_branch':str(predicted),'linear_all_orbits':total//factorial(n),'affine_all_orbits':afftotal//factorial(n)})

# Full code-space equality, using cyclic ideal generator instead of digit arguments.
for p,P in [(2,16),(3,27),(5,25),(7,7)]:
    for F in range(P+1):
        pascal=[[comb(j,h)%p for j in range(P)] for h in range(F)]
        a=P-F
        generator=[comb(a,j)*(-1)**(a-j)%p if j<=a else 0 for j in range(P)]
        ideal=[[generator[(j-k)%P] for j in range(P)] for k in range(F)]
        ck(rank(pascal,p)==F and rank(ideal,p)==F and rank(pascal+ideal,p)==F,'pascal_cyclic_ideal_equality')

# Exact direct orbit traversals for the Jordan ring map and mixed moduli.
def orbit_census(mod,F):
    unseen=set(product(range(mod),repeat=F));out=Counter()
    while unseen:
        x=next(iter(unseen));y=x;t=0
        while True:
            unseen.remove(y);t+=1
            y=tuple((y[i]+(y[i+1] if i+1<F else 0))%mod for i in range(F))
            if y==x:break
        out[t]+=1
    return out

crt=[]
for mods,F in [([4,3],3),([2,9],3),([4,5],2),([3,5],3)]:
    local=[orbit_census(mod,F) for mod in mods]
    predicted=Counter()
    for (a,ca),(b,cb) in product(local[0].items(),local[1].items()):predicted[lcm(a,b)]+=ca*cb*gcd(a,b)
    actual=orbit_census(mods[0]*mods[1],F)
    ck(actual==predicted,'crt_cycle_counts')
    crt.append({'moduli':mods,'F':F,'cycle_counts':dict(sorted(actual.items()))})

# Independent fixed-state counts at a slightly larger prime-power lift.
for p in [2,3,5,7]:
    # Enumerate constraints analytically coordinate by coordinate, and check
    # the corresponding integer matrix congruence on a deterministic sample.
    F=p+1;mod=p*p
    M=[[comb(p,j-i)%mod if 0<j-i<=p else 0 for j in range(F)] for i in range(F)]
    samples=list(product(range(mod),repeat=F)) if p<=3 else [tuple((t*(i+1)*(i+2)+i*i)%mod for i in range(F)) for t in range(500)]
    for x in samples:
        kernel=all(sum(a*b for a,b in zip(row,x))%mod==0 for row in M)
        predicted=all(x[i]%p==0 for i in range(2,p)) and (x[p]+p*x[1])%mod==0
        ck(kernel==predicted,'first_lift_fixed_criterion')

print(json.dumps({'status':'PASS','checks':checks,'total_checks':sum(checks.values()),
 'universes':{'fixed_rank':'all p-regular partition types n<=9, p=2,3,5,7,11',
 'branch_integrality':'p=2,3,5,7,11; 1<=s<=3; 1<=n<=18',
 'pascal_ideal':'all F from 0 through P for (p,P)=(2,16),(3,27),(5,25),(7,7)',
 'first_lift':'all states for p=2,3; 500 deterministic states each for p=5,7'},
 'direct_ring_gain_cases':ring,'crt_direct_orbit_cases':crt,
 'field_ring_hostile':{'size':4,'n':3,'regular_branch_both':2,'F4_all_orbits':4,'Z4_all_orbits':3,
 'reason':'The 1-dimensional quotient action is sign. Sign is trivial over F4 but has three orbits over Z4.'}},indent=2,sort_keys=True))
