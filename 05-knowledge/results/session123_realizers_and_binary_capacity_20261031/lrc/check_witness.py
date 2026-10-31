#!/usr/bin/env python3
"""Deterministic exact controls; large tests never enumerate periodic teeth."""
from fractions import Fraction as F
from itertools import product
from math import ceil,floor
from pathlib import Path
import json,random
from witness import choose,conductor,bad_measure,norm,need

def literal_bad_measure(lo,hi,v,p,r):
    total=F(0)
    for n in range(floor(v*lo+p-r)-1,ceil(v*hi+p+r)+2):
        total+=max(F(0),min(hi,(n+r-p)/v)-max(lo,(n-r-p)/v))
    return total

def literal_safe(lo,hi,vs,ps,r):
    out=[(lo,hi)]
    for v,p in zip(vs,ps):
        comb=[]
        for n in range(floor(v*lo+p)-1,ceil(v*hi+p)+2):
            a=max(lo,(n+r-p)/v);b=min(hi,(n+1-r-p)/v)
            if a<=b:comb.append((a,b))
        out=[(max(a,c),min(b,d)) for a,b in out for c,d in comb if max(a,c)<=min(b,d)]
    return out

def main():
    measures=0
    for lo,W,v,p,r in product((F(-2,3),F(0),F(3,7)),(F(1,9),F(1),F(7,3)),range(1,7),(F(0),F(2,9),F(-1,5)),(F(1,14),F(1,6),F(1,4))):
        need(bad_measure(lo,lo+W,v,p,r)==literal_bad_measure(lo,lo+W,v,p,r),'measure oracle');measures+=1
    rng=random.Random(123);counts={};cases=0;maxdepth=0
    for k in range(2,7):
        r=F(1,14);a=2*r;N=ceil(a*(k-1)/(1-k*a));c=k*a*N
        for ix in range(180):
            u=rng.randrange(1,15);vs=[u]+[rng.randrange(u,25) for _ in range(k-1)]
            lo=F(rng.randrange(-10,10),rng.randrange(1,10));W=c/u+F(rng.randrange(0,3),19)
            ps=[F(rng.randrange(-8,9),rng.randrange(1,12)) for _ in vs]
            t,st=conductor((lo,lo+W),vs,r,ps)
            need(any(l<=t<=h for l,h in literal_safe(lo,lo+W,vs,ps,r)),'literal common cover')
            counts[st['branch']]=counts.get(st['branch'],0)+1;cases+=1;maxdepth=max(maxdepth,st['bisections'])
    # General budget equality: isolated-point witness and positive-overlap repair.
    specials=[]
    for ps,expected in [([F(0),F(1,2)],'equality_boundary'),([F(0),F(-1,8)],'overlap_repair')]:
        t,st=choose((0,1),[1,1],F(1,4),ps);need(st['branch']==expected,'equality branch control')
        need(all(norm(t+p)>=F(1,4) for p in ps),'equality actual witness');specials.append({'witness':str(t),'stats':st})
    rejected=0
    for args in [((0,F(1,10)),[1,1],F(1,4),[0,F(1,2)]),((0,1),[1,1],F(1,3),[0,0])]:
        try:choose(*args)
        except ValueError:rejected+=1
        else:raise RuntimeError('hostile accepted')
    huge=[]
    for bits in (64,256,1024,4096):
        u=2**bits+57;vs=[u,2*u+1,3*u+7,5*u+3,7*u+11];r=F(1,14);W=F(10,7*u)
        lo=F(1,11);ps=[F(0),F(1,13),F(-2,17),F(3,19),F(4,23)]
        t,st=conductor((lo,lo+W),vs,r,ps);need(all(norm(v*t+p)>=r for v,p in zip(vs,ps)),'huge witness')
        huge.append({'speed_bits':bits,'witness_numerator_bits':t.numerator.bit_length(),'witness_denominator_bits':t.denominator.bit_length(),'stats':st})
    return {'status':'PASS','literal_measure_controls':measures,'literal_common_cover_controls':cases,'branch_counts':counts,'small_max_bisections':maxdepth,'equality_controls':specials,'rejected_budget_hostiles':rejected,'huge_binary_input_controls':huge}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
