#!/usr/bin/env python3
"""Exact gcd-fibre transfer, independently checked against literal safe combs."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, ceil
import json

def require(ok, msg):
    if not ok: raise RuntimeError(msg)
def norm(x):
    x=x-x.numerator//x.denominator
    return min(x,1-x)
def grid_value(s,d,v):
    g=gcd(d,v); h=d//g; u=v//g
    return F(h-1,2*h)+norm(u*s-F(h-1,2))/h

def lift(s,d,v):
    g=gcd(d,v); h=d//g; u=v//g
    if h==1: return s/d
    z=F(h,2)-u*s+F(1,2)
    m=(z.numerator//z.denominator)%h
    j=(m*pow(u,-1,h))%h
    return (s+j)/d

def safe(core,r):
    out=[(F(0),F(1))]
    for v in core:
        teeth=[(F(k,v)+r/v,F(k+1,v)-r/v) for k in range(v)]
        nxt=[]; i=j=0
        while i<len(out) and j<len(teeth):
            a,b=out[i]; c,e=teeth[j]
            if max(a,c)<=min(b,e): nxt.append((max(a,c),min(b,e)))
            if b<e:i+=1
            elif e<b:j+=1
            else:i+=1;j+=1
        out=nxt
    return out

def interval_hits(a,b,u,lo,hi):
    # Does u*[a,b] meet any integer translate of the CLOSED phase interval?
    return ceil(u*a-hi)<=((u*b-lo).numerator//(u*b-lo).denominator)
def transfer_accept(intervals,d,v,r):
    g=gcd(d,v); h=d//g; u=v//g
    if not intervals:return False
    a=h*r-F(h-1,2)
    if a<=0:return True
    if h%2:
        return any(interval_hits(l,z,u,a,1-a) for l,z in intervals)
    b=F(1,2)-a
    return any(interval_hits(l,z,u,-b,b) for l,z in intervals)

def run():
    gridchecks=0; transferchecks=0; lawchecks=0
    for d in range(1,25):
      for v in range(1,31):
       for den in range(1,14):
        for num in range(den+1):
         s=F(num,den)
         actual=max(norm(v*(s+j)/d) for j in range(d))
         require(actual==grid_value(s,d,v),'grid formula')
         t=lift(s,d,v)
         require(0<=t<=1 and norm(v*t)==actual,'constructive lift')
         require(norm(d*t-s)==0,'fibre membership')
         gridchecks+=1
    cores=[c for k in (1,2,3) for c in combinations(range(1,6),k)]
    for r in map(F,['1/14','1/4','1/3','2/5','1/2']):
     for c in cores:
      ints=safe(c,r)
      for d in range(1,7):
       for v in range(1,16):
        actual=bool(safe(tuple(d*x for x in c)+(v,),r))
        require(actual==transfer_accept(ints,d,v,r),'whole-core transfer')
        transferchecks+=1
        if r<=F(1,4) and ints:
         expected=True if v%d else bool(safe(c+(v//d,),r))
         require(actual==expected,'rejection dilation law')
         lawchecks+=1
    # Sharp per-fibre threshold: u*s is (h-1)/2 modulo one.
    for h in range(2,40):
      s=F((h-1)%2,2)
      require(grid_value(s,h,1)==F(h-1,2*h),'sharp fibre threshold')
    # Global high-radius hostile: core 1 at r=2/5 safe; {2,1} is not.
    require(bool(safe((1,),F(2,5))), 'hostile core nonempty')
    require(not safe((2,1),F(2,5)), 'above-threshold global hostile')
    # Failed parity erasure at a point: h2, s1/2 -> clearance1/4, not1/2.
    require(grid_value(F(1,2),2,1)==F(1,4),'even parity hostile')
    # Huge input witness without enumerating any dilation fibre.
    d=10**100+151; v=10**150+7; s=F(1,13)
    t=lift(s,d,v)
    require(all(norm(d*c*t)>=F(1,14) for c in range(1,13)), 'large core')
    require(norm(v*t)>=F(1,14),'large attached speed')
    return {'status':'PASS','grid_formula_and_native_witness_checks':gridchecks,
      'independent_literal_comb_transfer_checks':transferchecks,
      'low_radius_rejection_law_checks':lawchecks,
      'sharp_fibre_orders_checked':'2..39','high_radius_hostile':{'C':[1],'d':2,'v':1,'r':'2/5'},
      'huge_input_witness':{'d':str(d),'v':str(v),'t':str(t)},
      'scope':'Exact theorems for every core, dilation, attached speed; no arbitrary-core LRC14 nonemptiness claim.'}
if __name__=='__main__': print(json.dumps(run(),indent=2,sort_keys=True))
