#!/usr/bin/env python3
"""Read-only exact controls for modular refinement of bank head bounds."""
from fractions import Fraction as F
from itertools import product
import json

checks=0

def check(test,label):
    global checks
    checks+=1
    if not test:
        raise RuntimeError(label)

def ceil(x):
    return -(-x.numerator//x.denominator)

def coeff(bits):
    a,b,m=1,0,1
    for bit in bits:
        if bit=='1': a,b=3*a,3*b+m
        elif bit!='0': raise ValueError('binary word required')
        m*=2
    return a,b,m

K='110110110110110101011011101010'
words={'K':K,'V':K*4+'10','F':'10','H':'110101000'}
centers=[]
for k in 'KV':
    a,b,m=coeff(words[k]);check(a>m,'expansion');centers.append(F(b,a-m))
C=max(centers)
P=max(F(coeff(words[k])[1],coeff(words[k])[2]-coeff(words[k])[0]) for k in 'FH')
E=C+P
classes={(0,1):13,(0,0):7,(1,1):5,(1,0):23}
exponent_pairs=0
for n in range(3,513):
    for h in range(1,n+1):
        delta=2**n-3**h
        if delta>0:
            d=classes[n%2,h%2]
            check(delta%24==d,'mod24 residue')
            check(delta>=d>=5,'positive residue lower bound')
            exponent_pairs+=1
check(2**2-3==1,'N>=3 hypothesis necessary')
check(2**3-1==7,'h=0 outside stated mod24 table')
# The h=0, even-N case gives residue15, not table entry7.
check((2**4-1)%24==15,'h>=1 hypothesis necessary')
sequences=[''.join(z) for n in range(1,6) for z in product('KV',repeat=n)]
sequences+=['KV'*8,'VK'*8,'K'*32,'V'*16]
counts={str(d):0 for d in classes.values()}
rows=[]
for seq in sequences:
    exp=''.join(words[k] for k in seq)
    ea,eb,em=coeff(exp)
    for pattern in ('F','H','FH','HF','FFH','HHF'):
        reset='';r=F(1);j=0
        while r*F(ea,em)>=1:
            letter=pattern[j%len(pattern)]
            aa,bb,mm=coeff(words[letter])
            reset+=words[letter];r*=F(aa,mm);j+=1
        bits=exp+reset
        a,b,m=coeff(bits)
        n,h=len(bits),bits.count('1')
        delta=m-a;d=classes[n%2,h%2]
        rho=(-b*pow(a,-1,m))%m
        tmin=max(0,(b-delta*rho)//(delta*m)+1)
        check(0<rho<m and delta>0,'lawful positive cylinder and contraction')
        check(F(b,m)<E,'inherited envelope direct control')
        check(delta%24==d and delta>=d,'safe class denominator')
        bound=ceil(E/d)
        check(tmin<=bound<=ceil(E/5)==107,'refined head')
        check(delta*(rho+bound*m)>b,'refined tail endpoint inequality')
        counts[str(d)]+=1
        rows.append(tmin)
# Entire small exponent parity table, not just classes seen in bank controls.
check({d:ceil(E/d) for d in (5,7,13,23)}=={5:107,7:77,13:42,23:24},'exact ceilings')
print(json.dumps(dict(status='PROVED modular refinement; no source-entry assertion',checks=checks,
 exponent_pairs=exponent_pairs,bank_schedules=len(rows),C=str(C),P=str(P),uniform_head=107,
 class_controls=counts,class_bounds={str(d):ceil(E/d) for d in sorted(classes.values())},
 max_observed_t_min=max(rows)),indent=2,sort_keys=True))
