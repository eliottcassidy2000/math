#!/usr/bin/env python3
"""Direct divisor readings and exact minors check the incoming Smith corollary."""
from itertools import combinations,product
from math import gcd,isqrt,prod
from pathlib import Path
import json
import sympy as sp

def need(b,msg):
    if not b:raise RuntimeError(msg)
def aliquot(n):
    total=0
    for d in range(1,isqrt(n)+1):
        if n%d==0:total+=d+(n//d if d*d!=n else 0)
    return total-n
reports=[]
for n,p,q in ((1,2,3),(22,2,11),(722,2,19),(722,2,5)):
    cube=list(product(range(2),repeat=2))
    H=sp.Matrix([[aliquot(n*p**(e[0]+h[0])*q**(e[1]+h[1])) for h in cube] for e in cube])
    minors=[1]
    for k in range(1,5):
        g=0
        for rows in combinations(range(4),k):
            for cols in combinations(range(4),k):g=gcd(g,abs(int(H.extract(rows,cols).det())))
        minors.append(g)
    smith=[minors[i+1]//minors[i] for i in range(4)]
    u=n
    for ell in (p,q):
        while u%ell==0:u//=ell
    if u==1:need(smith==[1,1,n,p*q*n*(p+q-1)],'pure two-prime corollary')
    else:need(n==722 and smith[2]!=722,'prime alphabet guard')
    reports.append(dict(n=n,primes=[p,q],outside_core=u,H=[list(map(int,H.row(i))) for i in range(4)],
                        determinantal_divisors=minors,smith=smith))
Path(__file__).with_name('dilation_722.json').write_text(json.dumps(reports,indent=2)+'\n')
print(json.dumps([{k:v for k,v in x.items() if k not in ('H','determinantal_divisors')} for x in reports],indent=2))
