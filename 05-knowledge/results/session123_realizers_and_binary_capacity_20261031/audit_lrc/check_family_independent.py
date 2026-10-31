#!/usr/bin/env python3
"""Full-family implementation audit; independent complete event-cell geometry."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lrc'))
from family_selector import choose_family

def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def safe(t,vs):
    q=t.denominator;p=t.numerator
    return all(14*min((v*p)%q,q-(v*p)%q)>=q for v in vs)
def cover(vs):
    es={F(0),F(1)}
    for v in set(vs):
        for n in range(v):
            es.add(F(14*n+1,14*v));es.add(F(14*n+13,14*v))
    es=sorted(es)
    out=[(x,x) for x in es if safe(x,vs)]
    for x,y in zip(es,es[1:]):
        if safe((x+y)/2,vs):out.append((x,y))
    return out

def check(core,tails,full=False):
    t,st=choose_family(core,tails);vs=(*core,*tails)
    require(0<=t<=1 and safe(t,vs),'unsafe original-coordinate witness')
    if full:
        J=cover(vs);require(any(a<=t<=b for a,b in J),'not in independent full cover')
    require(st['finite_tail_intersections']<=len(tails),'too many finite steps')
    if st['method']=='conductor':
        ts=st['terminal_stats'];require(ts['bisections']<=ts['depth_bound'],'depth')
    return st

def main():
    methods={};full=0;bounded=0
    for core in combinations(range(1,14),8):
        st=check(core,(14,15,16,17,18),True);full+=1
        methods[st['method']]=methods.get(st['method'],0)+1
        check(core,(1,1,1,13,13));bounded+=1
    # Full bounded pool has genuinely point-only geometry; selector must retain it.
    pool=tuple(range(1,14));J=cover(pool)
    require(J and all(a==b for a,b in J),'expected point-only completed cover')
    check(pool,(),True)
    # All larger core sizes, arbitrary repeated tails.
    larger=0
    for m in range(9,14):
        for core in combinations(range(1,14),m):
            check(core,(2**97+13,)*(13-m));larger+=1
    boundary=check(tuple(range(1,13)),(24,),True)
    require(boundary['method']=='one_tail' and boundary['cutoffs']==[24],'one-tail exact boundary')
    huge=[]
    for bits in (64,257,1025,4096):
        for tails in ((14,27,2**bits+1,2**bits+5,3*2**bits+7),(1,1,14,27,2**bits+1)):
            st=check(tuple(range(1,9)),tails)
            huge.append({'bits':bits,'method':st['method'],'max_expanded_speed':st['max_expanded_speed']})
    invalid=[((1,)*8,(14,)*5),(tuple(range(1,9)),(0,14,15,16,17)),(tuple(range(1,9)),(-1,14,15,16,17)),(tuple(range(1,9)),(F(3,2),14,15,16,17)),(tuple(range(1,9)),(14,)*4),(tuple(range(1,8)),(14,)*6),((1,2,3,4,5,6,7,14),(14,)*5)]
    for c,t in invalid:
        try:choose_family(c,t)
        except ValueError:pass
        else:raise RuntimeError('invalid family accepted')
    print(json.dumps({'status':'PASS','independent_full_cover_rows':full,'repeated_bounded_tail_rows':bounded,'larger_core_rows':larger,'point_only_completed_cover_points':len(J),'invalid_controls':len(invalid),'one_tail_exact_boundary':True,'methods_on_full_cover_rows':methods,'huge_controls':huge},indent=2,sort_keys=True))
if __name__=='__main__':main()
