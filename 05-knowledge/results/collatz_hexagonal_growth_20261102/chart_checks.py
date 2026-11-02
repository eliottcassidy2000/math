"""Finite-horizon residue charts for the selected SF5 inverse operation."""
from collections import Counter
from functools import lru_cache
import json

def choices(r):
    candidates=(2,4,6) if r%3==1 else (3,5,7)
    return tuple(a for a in candidates if (2**a*r-1)%9!=0)

@lru_cache(None)
def poly(r,h):
    if not h:return (1,)
    A=choices(r);out=Counter({i:1 for i in range(max(A)+1)})
    for a in A:
        child=((2**a*r-1)//3)%(3**h) if h>1 else 0
        for i,c in enumerate(poly(child,h-1)):out[i+a+1]+=c
    return tuple(out[i] for i in range(max(out)+1))

@lru_cache(None)
def lower_trace(r,h):
    if not h:return ()
    A=choices(r);a=A[0]
    child=((2**a*r-1)//3)%(3**h) if h>1 else 0
    return (A,)+lower_trace(child,h-1)

def main():
    rows=[];collision=None;literal_checked=0
    for h in range(1,8):
        roots=[r for r in range(3**(h+1)) if r%3]
        traces={lower_trace(r,h) for r in roots}
        if len(traces)!=2*3**h:raise RuntimeError('trace chart faithfulness')
        polys={};counts=set()
        for r in roots:
            P=poly(r,h);counts.add(sum(P))
            if h<=4:
                n=r if r%2 else r+3**(h+1)
                if n==1:n+=2*3**(h+1)
                layer=[n];all_sources=[n]
                for _ in range(h):
                    layer=[(2**a*q-1)//3 for q in layer for a in range(2,8)
                        if (2**a*q-1)%3==0 and (2**a*q-1)%9!=0]
                    all_sources.extend(layer)
                seen={}
                for source in all_sources:
                    path=[source]
                    while path[-1]!=n:
                        if len(path)>8*h:raise RuntimeError('literal chart root clock')
                        x=path[-1];path.append(x//2 if x%2==0 else 3*x+1)
                    for i,x in enumerate(path):
                        d=len(path)-1-i
                        if x in seen and seen[x]!=d:raise RuntimeError('literal chart source clock')
                        seen[x]=d
                D=Counter(seen.values());literal=tuple(D[i] for i in range(max(D)+1))
                if literal!=P:raise RuntimeError('chart versus native path union')
                literal_checked+=1
            if collision is None and P in polys:collision={'h':h,'roots':[polys[P],r],'polynomial':P}
            polys[P]=r
        rows.append({'h':h,'modulus':3**(h+1),'branch_response_classes':len(traces),
            'polynomial_classes':len(polys),'count_classes':len(counts),
            'states_with_clock_up_to_h':3**(h+1)-2})
    return {'status':'PASS','rows':rows,'first_polynomial_collision':collision,
        'cache_entries':poly.cache_info().currsize,'literal_root_classes_checked':literal_checked}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
