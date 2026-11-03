#!/usr/bin/env python3
"""Fresh finite controls for reconstructed ROOT/lattice statements; stdlib only."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import gcd, factorial
import json

GATES=0

def check(ok, detail):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(detail)

def val2(n):
    n=abs(n)
    if n==0:
        raise ValueError('valuation zero')
    return (n & -n).bit_length()-1

def norm(x):
    a,b=x
    return a*a-a*b+b*b

def mul(x,y,mod=None):
    a,b=x;c,d=y
    z=(a*c-b*d,a*d+b*c-b*d)
    return tuple(v%mod for v in z) if mod else z

def uword(n):
    seen=set();word=[];states=[]
    while n!=1:
        if n in seen or len(word)>2000:
            raise RuntimeError('finite U sample not authenticated')
        seen.add(n);states.append(n)
        a=val2(3*n+1);word.append(a);n=(3*n+1)>>a
    return states,word

def h(n):
    states,word=uword(n);S=sum(word)
    return F((-1 if S%2 else 1)*2**S,3**len(word)*n)

u_mass=F(0);u_energy=F(1);max_tau=0
for n in range(3,4002,2):
    states,word=uword(n);tau=len(word);S=sum(word)
    L=tau-1;K=sum((a-1)//2 for a in word)
    weight=F(2*factorial(K)*factorial(L+1),factorial(L+K+2))
    hn=F((-1 if S%2 else 1)*2**S,3**tau*n)
    nextn=(3*n+1)>>word[0]
    c=F((-1 if word[0]%2 else 1)*3*n,3*n+1)
    check(h(nextn)==c*hn,('U edge',n))
    check(abs(hn)**6<=tau+1,('U growth',n))
    product=F(1)
    for x in states:product*=F(3*x,3*x+1)
    check(product==F(3**tau*n,2**S),('U telescoping',n))
    u_mass+=weight;u_energy+=weight*hn*hn;max_tau=max(max_tau,tau)
check(u_mass<F(13,3),'finite mass')
check(u_energy<55,'finite norm')
for j in range(1,64):
    r=F(j,64);D=1+r+r*r;kap=1-r*r/D
    A=r*(1+r*r)/D;e=(1-r)*r;B=(r+r*r)/D
    b2=e*B+(A-e)*kap
    target=r*r*(1+2*r+2*r*r-r**4)/D**2
    check(b2==target,('B2 simplification',j))
    check(A<=r and b2<=r*r,('generation comparison',j))

# P words are recorded as (a,j), not just depth or endpoint.
def pstep(s):
    if s<=0 or s%2==0 or s%8==5:raise ValueError('primitive source')
    a=2 if s%8==1 else 1
    y=(3*s+1)//(2**a);j=0
    while y%8==5:y=(y-1)//4;j+=1
    t=a+2*j
    check(2**t*(3*y+1)==3*(3*s+1)+2**a,('P carry',s))
    return y,a,j,t

pair_checks=0
for depth in range(1,5):
    groups=defaultdict(list)
    for s in range(1,8192,2):
        if s%8==5:continue
        y=s;word=[];clock=0
        for _ in range(depth):
            y,a,j,t=pstep(y);word.append((a,j));clock+=t
        groups[tuple(word)].append((s,y,clock))
    for rows in groups.values():
        for left,right in zip(rows,rows[1:]):
            s,y,T=left;z,w,T2=right
            check(T==T2,('common clock',left,right))
            dm=(s-z)//2;dy=(y-w)//2
            check(dy*2**T==3**depth*dm,('P difference',left,right))
            check(val2(dy)==val2(dm)-T,('P owner',left,right))
            pair_checks+=1
r17=pstep(17);r35=pstep(35)
check(r17[0]==r35[0]==3 and r17[3]==4 and r35[3]==5,'cross-word hostile')

# Hensel roots found by literal polynomial search, independent of lift formula.
q_rows=[]
for k in (1,2,3):
    q=7**k
    roots=[t for t in range(q) if (t*t+t+1)%q==0]
    check(len(roots)==2,('root count',q))
    alpha=next(x for x in roots if x%7==2);beta=next(x for x in roots if x%7==4)
    delta=(beta-alpha)%q;di=pow(delta,-1,q);a=2*di%q
    primitives=Counter();seen=set()
    two_obs=set()
    for m in range(q):
        for n in range(q):
            u=(m+alpha*n)%q;v=(m+beta*n)%q
            seen.add((u,v));two_obs.add((u,(m+(alpha+2)*n)%q))
            nn=(v-u)*di%q;mm=(u-alpha*nn)%q
            check((mm,nn)==(m,n),('paired inverse',q,m,n))
            check((norm((m,n))-u*v)%q==0,('paired norm',q,m,n))
            du=(m+(alpha+2)*n)%q;dv=(m+(beta+2)*n)%q
            check((du,dv)==(((1-a)*u+a*v)%q,(-a*u+(1+a)*v)%q),('D matrix',q,m,n))
            check((du==v)==(k==1 or n%(7**(k-1))==0),('exchange iff',q,m,n))
            check(norm((m+2*n,n))-norm((m,n))==2*n*(2*m+n),('norm defect',q,m,n))
            if m%7 or n%7:primitives[u]+=1
    check(len(seen)==q*q and len(two_obs)==q*q,('observer bijection',q))
    check(sum(primitives.values())==q*q-(q//7)**2,('primitive total',q))
    check(all(primitives[u]==(q if u%7 else 6*q//7) for u in range(q)),('fibre mass',q))
    q_rows.append({'q':q,'roots':[alpha,beta],'primitive_states':sum(primitives.values()),
                   'D_matrix':[[int((1-a)%q),a],[int(-a%q),int((1+a)%q)]]})

# Exhaustive product-ring compatibility in the smallest split quotient.
q=7;al=2;be=4
for m in range(q):
 for n in range(q):
  for c in range(q):
   for d in range(q):
    x=mul((m,n),(c,d),q)
    check((x[0]+al*x[1])%q==((m+al*n)*(c+al*d))%q,'ring plus')
    check((x[0]+be*x[1])%q==((m+be*n)*(c+be*d))%q,'ring minus')

# Literal Euclid-word endpoint and determinant, root-first application.
def U(x):m,n=x;return 2*m-n,m
def A(x):m,n=x;return 2*m+n,m
def D(x):m,n=x;return m+2*n,n
word_hosts=[]
for q in (7,49,343):
    x=A((2,1))
    for _ in range(q-1):x=D(x)
    x=D(U(x))
    check(x==(16*q+2,4*q+1),'opposite determinant word endpoint')
    check(tuple(t%q for t in x)==(2,1),'endpoint residue collision')
    check(gcd(*x)==1 and x[0]>x[1]>0 and sum(x)%2==1,'positive source guard')
    word_hosts.append({'q':q,'endpoint':x,'word_det':-1,'empty_det':1})

# Inert local fields: all local elements, complete norm fibres, conjugation.
inert=[]
for p in (2,5,29):
    check(not any((x*x+x+1)%p==0 for x in range(p)),('inert root',p))
    counts=Counter()
    for m in range(p):
      for n in range(p):
        z=(m,n);power=(1,0)
        for _ in range(p):power=mul(power,z,p)
        conj=((m-n)%p,-n%p)
        check(power==conj,('Frobenius',p,z))
        check(mul(z,conj,p)==(norm(z)%p,0),('norm field',p,z))
        counts[norm(z)%p]+=1
    check(counts[0]==1 and all(counts[r]==p+1 for r in range(1,p)),('norm fibres',p))
    inert.append({'p':p,'field_units':p*p-1,'nonzero_norm_fibre':p+1})
# Direct modulus290 counting independently confirms the CRT fibre product.
counts=Counter()
for m in range(290):
 for n in range(290):
    qv=norm((m,n))%290
    if gcd(qv,290)==1:counts[qv]+=1
check(len(counts)==112 and set(counts.values())=={540} and sum(counts.values())==60480,'290 norm fibres')
check(norm((6,2))==28,'E norm28')
check(not any(norm((a,b)) in (145,290) for a in range(-25,26) for b in range(-25,26)),'bounded inert norm hostile')

def gm(x,y):a,b=x;c,d=y;return a*c-b*d,a*d+b*c
check(gm((2,1),(5,-2))==(12,1),'Gaussian145 factors')
check(gm((1,1),(12,1))==(11,13),'Gaussian290 first')
check(gm((1,1),gm((2,1),(5,2)))==(-1,17),'Gaussian290 second')
check(F(112,580)==F(28,145),'MRB density is separate')

# Gaussian split/ramified comparison, exhaustive vectors at145 and290.
roots145=[r for r in range(145) if (r*r+1)%145==0]
check(roots145==[12,17,128,133],'Gaussian roots145')
inv24=pow(24,-1,145)
for m in range(145):
 for n in range(145):
    u=(m+12*n)%145;v=(m-12*n)%145
    nn=(u-v)*inv24%145;mm=(u-12*nn)%145
    check((mm,nn)==(m,n),'Gaussian paired inverse145')
    check((m*m+n*n-u*v)%145==0,'Gaussian paired norm145')
gc=Counter()
for m in range(290):
 for n in range(290):
    qv=(m*m+n*n)%290
    if gcd(qv,290)==1:gc[qv]+=1
check(len(gc)==112 and set(gc.values())=={224} and sum(gc.values())==25088,'Gaussian290 normfibres')
zz=gm((145,145),(145,145))
check((145,145)!=(0,0) and tuple(t%290 for t in zz)==(0,0),'Gaussian290 nilpotent')
check(all((r*r)%290 or r==0 for r in range(290)),'scalar290 reduced')

print(json.dumps({'status':'PASS','gates':GATES,'U_sources':2000,'U_max_tau':max_tau,
 'U_partial_mass':str(u_mass),'U_partial_norm':str(u_energy),'P_sources_per_depth':3072,
 'P_depths':[1,2,3,4],'P_consecutive_common_word_pairs':pair_checks,
 'P_cross_word_hostile':{'sources':[17,35],'endpoint':3,'clocks':[4,5]},
 'split_quotients':q_rows,'opposite_determinant_hosts':word_hosts,'inert_local_fields':inert,
 'mod290_norm_unit_count':60480,'mod290_norm_fibre_size':540,'Gaussian145_roots':roots145,
 'Gaussian290_norm_unit_count':25088,'Gaussian290_norm_fibre_size':224,
 'scope':'Finite controls of proof candidates; no universal ROOT or debt claim'},indent=2,sort_keys=True))
