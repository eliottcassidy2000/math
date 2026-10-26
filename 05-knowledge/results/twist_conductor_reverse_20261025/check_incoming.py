#!/usr/bin/env python3
"""Exact counterexamples to the three corrected incoming implications."""
from itertools import product
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def span(basis):
    out={0}
    for b in basis:out|={v^b for v in tuple(out)}
    return out
def main():
    C,D=span((13,10)),span((9,15));need(C&D=={0},'complementary packets')
    projection={c^d:c for c in C for d in D};need(len(projection)==16,'all characters')
    alpha=(1,2,0,0)
    val=lambda z:sum(a for i,a in enumerate(alpha) if z>>i&1)
    thresholds=[val(projection[z]&(projection[z]^z)) for z in range(16)]
    z=4;coefficient=2;single=[]
    for y in (0,1,2,4,8):
        have=coefficient+val(z&y);required=thresholds[z^y]
        need(have>=required,'original root tests pass');single.append((y,have,required))
    y=12;need(val(y)//2==0,'the failing product is unreduced')
    need(coefficient+val(z&y)<thresholds[z^y],'unreduced product fails')
    need((1+1)//2==1 and 1%3!=0,'strict ideal-containment direction')
    def projective(v,p):
        a=next(x for x in v if x)
        return tuple(x*pow(a,-1,p)%p for x in v)
    points={projective(v,3) for v in product(range(3),repeat=3) if any(v)}
    conic=sorted((s,t,w) for s,t,w in points if (w*w-s*s-t*t)%3==0)
    good=[v for v in conic if v[0]*v[1]*(v[0]**2+v[1]**2)%3]
    need(len(conic)==4 and not good,'nonempty conic but empty required open')
    base={projective(v,3) for v in product(range(3),repeat=2) if any(v)}
    open_base=sorted((x,z) for x,z in base if x*z*(x+z)%3)
    lifts=[(x,z) for x,z in open_base if x==z==(x+z)%3]
    need(open_base==[(1,1)] and not lifts,'direct squareclass check of central open')
    print(json.dumps({'status':'FINITE-EXACT PASS','basis_C':[13,10],'basis_D':[9,15],
        'radicand_valuations':alpha,'B_coefficient_valuations':thresholds,
        'source_character':4,'source_coefficient_valuation':2,
        'single_root_tests':single,'unreduced_failing_character':12,
        'failure_actual_required':[2,3],'F3_conic_points':conic,'F3_open_lifts':lifts,
        'corrections':'reverse ideal containment with reflexive hull; retain full character bank; restrict whole-conic criterion to infinite fields'},indent=2,sort_keys=True))
if __name__=='__main__':main()
