#!/usr/bin/env python3
"""Exact signed allocation. No numerical optimization is used."""
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def parts(n,m=2):
    if n==0:yield ()
    for a in range(m,n+1):
        for tail in parts(n-a,a):yield (a,)+tail
def add(*polys):
    r=defaultdict(F)
    for p in polys:
        for k,v in p.items():r[k]+=v
    return {k:v for k,v in r.items() if v}
def scale(a,p):return {k:a*v for k,v in p.items()}
def lift(k,p):
    r=defaultdict(F)
    for t,v in p.items():r[tuple(sorted(t+(k,)))]+=v*(t.count(k)+1)
    return dict(r)
def margin(n,sign):
    r={():F(1)}
    for s in range(2,n+1):
        for t in parts(s):r[t]=F((-1)**sum(x%2==0 for x in t)+(sign*(-1)**len(t) if s==n else 0))
    return r
def exact(n,b):
    r={}
    for s in range(2,n+1):
        for t in parts(s):r[t]=F((2**len(t) if all(x%2 for x in t) else 0)-(b if s==n else 1)*(-1)**sum(x%2==0 for x in t))
    return r
def fmt(p):return {','.join(map(str,k)):str(v) for k,v in sorted(p.items()) if v}
def run():
    P={(2,)*j:F((-1)**(j+1)) for j in range(1,6)}
    T={(3,):1,(3,3):3,(2,2,3):-1,(2,2,3,3):-1,(2,2,2,3):1,(2,2,2,2,3):-2}
    D5={(5,):1,(2,5):1,(3,5):3,(2,2,5):-1}
    D7={(7,):1,(2,7):1,(2,2,7):-2}
    D8={(8,):1,(2,8):-1}
    free={(2,3):1,(9,):1,(10,):1}
    rows=[]
    for b in [0,2]:
        sign=1 if b==0 else -1
        D4=add(lift(4,margin(7,sign)),{(4,4):1,(2,4,4):-1})
        D6=add(lift(6,margin(5,sign)),{(4,6):1})
        groups={'P':P,'T':T,'D4':D4,'D5':D5,'D6':D6,'D7':D7,'D8':D8,'free':free}
        residual=add(exact(11,b),scale(-1,add(*groups.values())))
        need(all(v>=0 for v in residual.values()),'cycle reserve reused')
        rows.append({'full_cover_multiplier':b,'groups':{k:fmt(p) for k,p in groups.items()},'residual':fmt(residual),'exact_ledger':fmt(exact(11,b))})
    need(F(53,128)-F(7,40)==F(153,640),'remaining triangle charge')
    need(F(4,5)*(1-F(3,4)**5)==F(781,1280),'sharp fair-pair coefficient')
    return {'status':'passed','scope':'Exact partition identities on all weighted order-eleven matrices; analytic local bounds and finite graph certificates are additional obligations.','signed_allocations':rows,'triangle_reserve':'53/128','wedge_payment':'7/40','remaining_triangle_reserve':'153/640','candidate_sharp_coefficient':'781/1280'}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
