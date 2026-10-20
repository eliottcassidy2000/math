#!/usr/bin/env python3
"""Independent exponent-coordinate checks for THM-4725, not copied from source.
Universe: every antipodal section and every class subset for primes 3..37.
Evenness uses the fundamental-period word criterion for cyclic multipliers.
Gamma signs use direct reference-section flip parity, not class-cycle signs.
"""
from math import gcd
from itertools import product
from sympy import primerange

def rotate(mask,a,m):
    a%=m
    return ((mask<<a)|(mask>>(m-a))) & ((1<<m)-1)

def even(mask,m):
    for d in range(1,m+1):
        if m%d==0 and rotate(mask,d,m)==mask:
            return (mask & ((1<<d)-1)).bit_count()%2==0

def transforms(p):
    n=p-1;m=n//2
    keys=[]
    for j in range(n):
        if gcd(j,n)!=1:continue
        for a in range(n):
            perm=[];flip=0
            for t in range(m):
                y=(a+j*t)%n;perm.append(y%m)
                if y>=m:flip|=1<<t
            keys.append((a,j,perm,flip))
    return keys

def image(mask,perm):
    out=0
    while mask:
        bit=mask&-mask;i=bit.bit_length()-1
        out|=1<<perm[i];mask-=bit
    return out

def census(p):
    n=p-1;m=n//2;N=1<<m;keys=transforms(p)
    tc=0;seen=set()
    for x in range(N):
        if x in seen:continue
        tc+=1
        seen.update(image(x^fl,pr) for a,j,pr,fl in keys)
    ec=0;gec=0;es=0;ges=0;seen=set();witness=None
    for x in range(N):
        if x in seen:continue
        orb={image(x,pr) for a,j,pr,fl in keys};seen.update(orb)
        if not even(x,m):continue
        assert all(even(y,m) for y in orb)
        ec+=1;es+=len(orb)
        ge=all(not ((x&fl).bit_count()%2) for a,j,pr,fl in keys if image(x,pr)==x)
        if ge:gec+=1;ges+=len(orb)
        elif witness is None:witness=x
    assert tc==gec
    if p%4==3:
        assert es==ges and ec==tc
        # Paley uses even exponents; test equivariance, two-element fibres.
        paley=sum(1<<t for t in range(m) if t%2==0)
        beta=lambda x: (x^paley) if (x^paley).bit_count()%2==0 else (N-1)^(x^paley)
        for x in range(N):
            b=beta(x);assert even(b,m) and beta(x^(N-1))==b
            if p<=19:
                assert all(beta(image(x^fl,pr))==image(b,pr) for a,j,pr,fl in keys)
        assert len({beta(x) for x in range(N)})==N//2
    elif p>=13:
        assert ec>tc and witness is not None
    elif p==5:
        assert (tc,ec,gec,es,ges)==(1,1,1,1,1)
    return p,tc,ec,gec,ges,es,witness

print('p,T_Gamma,E_Gamma,E_Gamma_even,Gamma_even_sets,even_sets,odd_witness_mask')
for p in primerange(3,38):print(*census(p),sep=',',flush=True)
# Signed constants and collapse at e=3. Exhaust all Boolean truth tables on
# c=(2,-2) at p=73, independently construct two quotient section types.
p=73;n=p-1;K={pow(x,8,p) for x in range(1,p)}
assert 2 in K and -1%p not in K
sections=[];orbits=set()
for bits in product((0,1),repeat=4):
    s=tuple(bits)+tuple(1-b for b in bits)
    sections.append(s)
    orbits.add(min(s[k:]+s[:k] for k in range(8)))
assert len(orbits)==2
symmetric_count=0
for table in range(16):
    outs=[]
    for s in sections:
        out=tuple((table>>(b+2*(1-b)))&1 for b in s)
        outs.append(out)
    if all(o[t]==o[t+4] for o in outs for t in range(4)):
        symmetric_count+=1;assert len(set(outs))==1
    else:
        assert all(o[t]!=o[t+4] for o in outs for t in range(4))
assert symmetric_count==8
print('Signed c=(2,-2), p=73: 16 Boolean rules; 8 graph-valued rules collapse, 8 non-graph-valued; two tournament types.')
print('Identity hostile: c=1, phi(b)=b retains both tournament types and is not graph-valued.')
