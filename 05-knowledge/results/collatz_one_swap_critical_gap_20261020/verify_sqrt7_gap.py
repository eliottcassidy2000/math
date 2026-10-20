#!/usr/bin/env python3
"""Exact critical-clock and small-weight one-swap audit; stdlib only."""
import json
import math

def check(ok,msg):
    if not ok:
        raise ArithmeticError(msg)

def carry(word):
    B=0
    for i,bit in enumerate(word):
        B=(3 if bit else 1)*B+bit*2**i
    return B

counts={'shapes':0,'swaps':0,'crossings':0,'legal_type_checks':0}
small={'p_min':2,'p_max_exclusive':100,'shapes':0,'swaps':0,'crossings':0,'integer_cycle_hits':0}
maximum=None
for p in range(2,513):
    q=(3**p).bit_length()
    if math.gcd(p,q)!=1:
        continue
    k,r=divmod(q,p)
    check(k==1,'critical k')
    check(4*r<=3*p,'critical residue')
    c=[((i+1)*p)//q-(i*p)//q for i in range(q)]
    B=carry(c)
    D=2**q-3**p
    after=[0]*(q+1)
    for i in range(q-1,-1,-1):
        after[i]=after[i+1]+c[i]
    s=[(j*k+(j*r+p-1)//p-1)%q for j in range(p)]
    labels={position:j for j,position in enumerate(s)}
    counts['shapes']+=1
    if p<100:
        small['shapes']+=1
    for position,j in labels.items():
        a=(-j*r)%p
        for epsilon in (-1,1):
            legal=(c[(position+epsilon)%q]==0)
            expected=(a>=p-r) if epsilon==-1 else (a<r)
            check(legal==expected,('legality',p,q,j,a,epsilon))
            counts['legal_type_checks']+=1
    for i in range(q):
        j=(i+1)%q
        if c[i]==c[j]:
            continue
        crossing=(i==q-1)
        if crossing:
            edited=c.copy()
            edited[i],edited[j]=edited[j],edited[i]
            Bw=carry(edited)
        else:
            epsilon=1 if c[i] else -1
            Bw=B+epsilon*2**i*3**after[i+2]
        g=math.gcd(D,Bw)
        check(g*g<7**p,('sqrt7 content',p,q,i,g))
        check(g<D,('integer cycle',p,q,i,g))
        counts['swaps']+=1
        counts['crossings']+=int(crossing)
        if p<100:
            small['swaps']+=1
            small['crossings']+=int(crossing)
            small['integer_cycle_hits']+=int(g==D)
        if maximum is None or g>maximum['g']:
            maximum=dict(p=p,q=q,swap_start=i,g=g)
check(3**5<2**8 and 3**3<2**5,'rational constants')
print(json.dumps({'status':'PASS','universe':'2<=p<=512; q=bit_length(3^p); gcd(p,q)=1; all unequal cyclic adjacent swaps',
    'counts':counts,'small_weight_certificate':small,'maximum_observed_content':maximum},indent=2,sort_keys=True))
