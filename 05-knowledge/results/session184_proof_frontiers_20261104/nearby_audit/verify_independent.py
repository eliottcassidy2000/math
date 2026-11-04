"""Independent exponent-lattice proof controls, no producer imports."""
from itertools import product
from math import prod
import json

def req(x,m):
    if not x: raise ValueError(m)
def val(v): return prod(p**a for p,a in zip((2,3,5),v))
def pairs(e):
    E=sum(e);out=[]
    for v in product(*(range(a+1) for a in e)):
        if sum(v) in (0,E):continue
        w=tuple(a-b for a,b in zip(e,v))
        if val(v)<=val(w):out.append((v,w))
    return out
count=0
for e in product(range(5),range(4),range(3)):
    E=sum(e)
    if E<2:continue
    pp=pairs(e)
    for K in (1,2,3):
        candidates=list({v for pair in pp for v in pair if sum(v)<=K})
        if len(candidates)>11:continue
        for mask in range(1<<len(candidates)):
            H={v for j,v in enumerate(candidates) if mask>>j&1}
            weak=bool(pp) and all(v in H or w in H for v,w in pp)
            strictpairs=[(v,w) for v,w in pp if v!=w]
            strict=bool(strictpairs) and all(v in H or w in H for v,w in strictpairs)
            req(not weak or E<=2*K+1,('weak cutoff',e,K,H))
            req(not strict or E<=2*K+2,('strict cutoff',e,K,H))
            if E==2*K+1:req(weak==all(v in H for v in candidates),('weak equality',e,K,H))
            if E==2*K+2:
                nonzero=[j for j,a in enumerate(e) if a]
                pure=len(nonzero)==1
                needed=[]
                if pure:
                    for a in range(1,K+1):
                        v=[0,0,0];v[nonzero[0]]=a;needed.append(tuple(v))
                req(strict==(pure and all(v in H for v in needed)),('strict equality',e,K,H))
            count+=1
# Explicit unbounded-Omega counterexample to replacing Omega by omega.
for j in range(1,13):
    e=(j,1,0);H=lambda v:v[0]>0 and v[1:]==(0,0)
    req(all(H(a) or H(b) for a,b in pairs(e)),('distinct-factor hostile',j))
print(json.dumps({'status':'PASS','independent_exponent_lattice_hole_cases':count,'unbounded_degree_family_controls':12},sort_keys=True,indent=2))
