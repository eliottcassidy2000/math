#!/usr/bin/env python3
"""Independent exact monomial-torsion ring audit; no author helper imports."""
from itertools import product
from math import ceil
from collections import Counter
from pathlib import Path
import json

checks=0
def require(x,why):
    global checks
    checks+=1
    if not x: raise RuntimeError(str(why))

def valuation(n,p):
    v=0
    while n%p==0: n//=p; v+=1
    return v

def clog(base,n):
    r=0; t=1
    while t<n: t*=base; r+=1
    return r

def flog(base,n):
    r=0
    while base**(r+1)<=n: r+=1
    return r

def scalar_clock(p,e,a):
    s=valuation(e,p)
    if p>2: return (a-1+s-1)//s
    if a==1: return 0
    units=max(1,(a-2+s-1)//s)
    return max(units,clog(e,a))

def clock(p,e,aa):
    s=valuation(e,p)
    return max([scalar_clock(p,e,aa[0])]+
               [(aj+flog(p,j)+s-1)//s for j,aj in enumerate(aa) if j])

def mul(x,y,mod):
    out=[0]*len(mod)
    for i in range(len(mod)):
        out[i]=sum(x[j]*y[i-j] for j in range(i+1))%mod[i]
    return tuple(out)

def power(x,e,mod):
    result=(1,)+(0,)*(len(mod)-1)
    while e:
        if e%2: result=mul(result,x,mod)
        x=mul(x,x,mod); e//=2
    return result

rows=[]
cases=[(3,3),(3,9),(3,21),(5,5),(5,45),(2,2),(2,4),(2,6)]
for p,e in cases:
    for b in range(1,5):
        for aa in product(range(1,4),repeat=b):
            if any(aa[i]<aa[i+1] for i in range(b-1)): continue
            if sum(aa)>(6 if p==5 else 7): continue
            mod=[p**a for a in aa]; bound=clock(p,e,aa)
            roots={}
            for x in range(mod[0]):
                if pow(x,e,mod[0])==x:
                    require(x%p not in roots,('root not unique',p,e,aa,x))
                    roots[x%p]=x
            require(len(roots)==p,('root number',p,e,aa,roots))
            worst=0; witness=None
            for x in product(*(range(m) for m in mod)):
                target=(roots[x[0]%p],)+(0,)*(b-1)
                y=x; t=0
                while y!=target and t<=bound:
                    y=power(y,e,mod); t+=1
                require(y==target and t<=bound,('upper bound',p,e,aa,x,bound,t,y,target))
                if t>worst: worst=t; witness=x
            require(worst==bound,('sharpness',p,e,aa,worst,bound,witness))
            rows.append({'p':p,'e':e,'torsion_exponents':aa,'predicted':bound,
                         'observed':worst,'witness':witness,'ring_size':__import__('math').prod(mod)})

out={'status':'PASS','checks':checks,'rings':len(rows),'cases':rows,
     'warning':'p=2 arithmetic-only clock needs separate units and nonunits formula; odd-prime formula is not universally sharp.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('status','checks','rings')}))
