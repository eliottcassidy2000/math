#!/usr/bin/env python3
"""Finite-exact hostile checks; no universal integer-halting conclusion."""
from fractions import Fraction
from itertools import product
import hashlib, json
from pathlib import Path

SEQ = [6,16,34,64,114,196,334,564,946,1584,2646,4416,7366,12284]
def require(condition):
    if not condition: raise AssertionError("verification gate failed")
def g(x):
    q,r=divmod(x,3)
    return None if r==2 else 5*q+(6 if r==0 else 9)
def native_cost(x):
    q,r=divmod(x,3)
    return [5*q*q+19*q+15,5*q*q+25*q+27,6*q+12][r]
def word(x,m):
    out=[]
    for _ in range(m):
        r=x%3
        if r==2:return None
        out.append(r); x=g(x)
    return tuple(out)
def address(w):
    x=Fraction(0)
    for r in reversed(w): x=(3*x-18-4*r)/5
    return x
def residue(x,mod): return (x.numerator*pow(x.denominator,-1,mod))%mod
def sign3(x):
    x=Fraction(x); n,d=x.numerator,x.denominator
    require(n)
    while n%3==0:n//=3
    while d%3==0:d//=3
    return (n*pow(d,-1,3))%3
def sigma(n): return sum(d for d in range(1,n+1) if n%d==0)

out={}
x=0; orbit=[]; costs=[]; branches=[]
while x is not None:
    costs.append(native_cost(x)); branches.append(x%3)
    x=g(x)
    if x is not None:orbit.append(x)
require(orbit==SEQ)
require(sum(costs)==47176870)
out['orbit']=orbit;out['branch_word_with_halt']=''.join(map(str,branches))
out['native_costs']=costs;out['native_total']=sum(costs)
out['sigma_ratios']={str(x):str(Fraction(sigma(x),x)) for x in SEQ}
require([x for x in SEQ if sigma(x)%x==0]==[6])
out['sigma_successor_coincidences']=[x for x in [0]+SEQ if x>0 and g(x)==sigma(x)]
require(out['sigma_successor_coincidences']==[946])
require(g(-9)==-9 and g(-11)==-11)

censuses=[];orientation_pairs=0
for m in range(1,10):
    mod=3**m
    inverse={w:residue(address(w),mod) for w in product((0,1),repeat=m)}
    direct={word(x,m):x for x in range(mod) if word(x,m) is not None}
    require(len(inverse)==len(set(inverse.values()))==len(direct)==2**m)
    require(inverse==direct)
    for w,r in inverse.items():
        require(word(r+mod,m)==w)
        x,y=r,r+mod
        for _ in range(m):x,y=g(x),g(y)
        require(y-x==5**m)
    if m<=8:
        ordered=sorted(inverse,key=lambda w:tuple(v if j%2==0 else 1-v for j,v in enumerate(w)))
        for i,w in enumerate(ordered):
            for v in ordered[i+1:]:
                require(sign3(address(v)-address(w))==1)
                orientation_pairs+=1
    censuses.append({'m':m,'modulus':mod,'safe_residues':len(direct)})
out['cylinders']=censuses;out['alternating_lex_orientation_pairs']=orientation_pairs
out['negative_nonhalting_controls']=[-9,-11]
out['scope']='Finite verification only; general cylinder/measure/tournament results have separate elementary proofs.'
print(json.dumps(out,indent=2))
