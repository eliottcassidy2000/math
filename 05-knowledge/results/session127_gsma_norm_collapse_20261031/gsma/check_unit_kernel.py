#!/usr/bin/env python3
"""Exact finite rank-two kernel census for k-nacci inverse powers.
Universe: k in {3,4,5}, integer n=0..500, basis 1,x,...,x^(k-1).
Independent replay uses monic polynomial reduction and binary powering,
not the coefficient update that produces the census. No all-exponent claim.
"""
import hashlib,json,math
from fractions import Fraction
from pathlib import Path
BOUND=500
CHECKS=0

def require(c,m):
    global CHECKS
    CHECKS+=1
    if not c: raise RuntimeError(m)

def canonical(v):
    return json.dumps(v,separators=(',',':')).encode()

def mul(a,b,k,mod=None):
    """Polynomial convolution, then x^k = 1+x+...+x^(k-1)."""
    z=[0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):z[i+j]+=u*v
    for i in range(len(z)-1,k-1,-1):
        t=z[i]
        for j in range(k):z[i-k+j]+=t
        z[i]=0
    z=(z+[0]*k)[:k]
    return [t%mod for t in z] if mod else z

def power(a,n,k,mod=None):
    b=[1]+[0]*(k-1)
    while n:
        if n&1:b=mul(b,a,k,mod)
        a=mul(a,a,k,mod);n>>=1
    return b

def det_bareiss(a):
    a=[r[:] for r in a];n=len(a);den=1;sign=1
    for k in range(n-1):
        if a[k][k]==0:
            r=next((r for r in range(k+1,n) if a[r][k]),None)
            if r is None:return 0
            a[k],a[r]=a[r],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=a[i][j]*pivot-a[i][k]*a[k][j]
                require(num%den==0,'Bareiss exact division')
                a[i][j]=num//den
        for i in range(k+1,n):a[i][k]=0
        den=pivot
    return sign*a[-1][-1]

def norm(a,k):
    cols=[mul(a,[int(i==j) for i in range(k)],k) for j in range(k)]
    return det_bareiss([list(row) for row in zip(*cols)])

def form(p,q,k):return p**k-sum(p**j*q**(k-j) for j in range(k))

def trim(a):
    while a and a[-1]==0:a.pop()
    return a

def rem(a,b,p):
    a=trim([t%p for t in a]);b=trim([t%p for t in b])
    while len(a)>=len(b):
        c=a[-1]*pow(b[-1],-1,p)%p;d=len(a)-len(b)
        for j,t in enumerate(b):a[d+j]=(a[d+j]-c*t)%p
        trim(a)
    return a

def gcd_poly(a,b,p):
    while b:a,b=b,rem(a,b,p)
    return [(t*pow(a[-1],-1,p))%p for t in a] if a else []

rows=[]
for k,p in [(3,3),(4,2),(5,5)]:
    one=[1]+[0]*(k-1);x=[0,1]+[0]*(k-2)
    inv=[-1]*(k-1)+[1]
    require(mul(x,inv,k)==one,'explicit polynomial inverse')
    # Rabin irreducibility criterion: degrees 3,4,5 need one proper test each.
    prime_divisor=2 if k==4 else k
    frob=power(x,p**k,k,p)
    require(frob==x,'Frobenius terminal identity')
    proper=power(x,p**(k//prime_divisor),k,p)
    proper[1]=(proper[1]-1)%p
    g=gcd_poly([-1]*k+[1],proper,p)
    require(g==[1],'proper Frobenius gcd')
    v=one;hits=[];all_rows=[];checkpoints=[]
    for n in range(BOUND+1):
        # Independent path: exponentiate the explicit inverse polynomial.
        independent=power(inv,n,k)
        require(v==independent,f'independent polynomial replay k{k} n{n}')
        require(mul(v,power(x,n,k),k)==one,f'forward cancellation k{k} n{n}')
        all_rows.append([n,v])
        if n in [0,1,k,16,50,100,250,500]:
            checkpoints.append({'n':n,'vector':v,'sha256':hashlib.sha256(canonical(v)).hexdigest()})
        if all(t==0 for t in v[2:]):
            a,b=v[:2];nf=norm(v,k)
            require(nf==(-1)**((k+1)*n),'unit norm sign')
            require(nf==form(a,-b,k),'linear form equals multiplication norm')
            hits.append({'n':n,'coefficients':v,'norm':nf,'binary_pair':[a,-b]})
        c=v[0];v=[t-c for t in v[1:]]+[c]
    expected={3:[0,3,16],4:[0,4],5:[0,5]}[k]
    require([h['n'] for h in hits]==expected,'frozen finite census')
    require(inv[-1]!=0 and not all(t==0 for t in inv[2:]),'hostile norm-one-to-kernel implication')
    require(abs(norm(inv,k))==1,'hostile is genuinely a unit')
    require(power(inv,k,k)==[2,-1]+[0]*(k-2),'uniform baseline unit identity')
    # All primitive pairs in a stated small box check independent determinant norm,
    # telescoping residual, and local gcd restrictions.
    primitive_count=0
    for a in range(-8,9):
        for b in range(-8,9):
            if math.gcd(a,b)!=1:continue
            primitive_count+=1;f=form(a,b,k)
            require(norm([a,-b]+[0]*(k-2),k)==f,'norm form determinant')
            require((a-b)*f==a**(k+1)-2*a**k*b+b**(k+1),'telescoping homogeneous residual')
            require(math.gcd(f,a*b)==1,'primitive norm coprime to pq')
            require((k-1)%math.gcd(f,a-b)==0,'diagonal local gcd gate')
            if a:
                t=Fraction(b,a)
                gf=t**(k-1)/(1-sum(t**j for j in range(1,k+1)))
                require(gf==Fraction(a*b**(k-1),f),'rational GF evaluation')
                require(gf.denominator==abs(f),'exact uncancelled norm denominator')
                require((gf.denominator==1)==(abs(f)==1),'integrality iff unit norm')
    rows.append({'k':k,'irreducibility_modulus':p,'proper_frobenius_gcd':g,'bound':BOUND,
                 'rank_two_hits':hits,'all_vectors_sha256':hashlib.sha256(canonical(all_rows)).hexdigest(),
                 'checkpoints':checkpoints,'primitive_box_count':primitive_count,
                 'hostile_inverse_vector':inv,'hostile_inverse_norm':norm(inv,k)})
result={'status':'FINITE-EXACT','universe':{'orders':[3,4,5],'n_min':0,'n_max':BOUND,
        'primitive_pair_regression_box':[-8,8]},'method':'coefficient recurrence plus independent exact monic-polynomial powering and forward cancellation',
        'checks':CHECKS,'rows':rows,'scope':'No all-exponent classification is asserted for k=4 or k=5; cubic global classification is separate cited work.'}
print(json.dumps(result,indent=2))
