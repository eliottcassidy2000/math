#!/usr/bin/env python3
"""Independent literal closed comb audit; no producer code imports."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, floor, ceil
from pathlib import Path
import json

def require(test, detail):
    if not test: raise RuntimeError(detail)
def norm(x):
    y=x-floor(x)
    return min(y,1-y)
def intersect(a,b):
    out=[]; i=j=0
    while i<len(a) and j<len(b):
        l=max(a[i][0],b[j][0]); r=min(a[i][1],b[j][1])
        if l<=r:
            if out and l<=out[-1][1]: out[-1]=(out[-1][0],max(r,out[-1][1]))
            else: out.append((l,r))
        if a[i][1]<b[j][1]: i+=1
        else: j+=1
    return out
def comb(v,r): return [(F(j+r,v),F(j+1-r,v)) for j in range(v)]
def safe(C,r):
    a=[(F(0),F(1))]
    for v in sorted(set(C)): a=intersect(a,comb(v,r))
    return a
def lift(s,d,v):
    g=gcd(d,v); h=d//g; u=v//g
    z=F(h,2)-u*s; m=floor(z+F(1,2))
    j=0 if h==1 else (pow(u,-1,h)*m)%h
    return (s+j)/d

def quotient(intervals,d,v,r):
    g=gcd(d,v); h=d//g; u=v//g; a=h*r-F(h-1,2)
    if not intervals: return None
    if a<=0: return lift(intervals[0][0],d,v)
    lo,hi=(a,1-a) if h%2 else (-h*(F(1,2)-r),h*(F(1,2)-r))
    for L,R in intervals:
        n=ceil(u*L-hi)
        if n<=floor(u*R-lo): return lift(max(L,(n+lo)/u),d,v)
    return None

counts={"grid_formula_and_native_witness":0,"literal_comb_transfer":0,"rejection_identity":0,"bank_finite_queries":0}
points=sorted({F(j,n) for n in range(1,14) for j in range(n+1)})
for d in range(1,13):
 for v in range(1,18):
  g=gcd(d,v); h=d//g; u=v//g
  for s in points:
   literal=max(norm(v*(s+j)/d) for j in range(d))
   formula=F(h-1,2*h)+norm(u*s-F(h-1,2))/h
   t=lift(s,d,v)
   require(literal==formula==norm(v*t),(d,v,s,literal,formula,t))
   counts['grid_formula_and_native_witness']+=1
radii=[F(1,14),F(1,4),F(1,3),F(2,5),F(1,2)]
cores=[c for k in range(1,4) for c in combinations(range(1,6),k)]
for C in cores:
 for r in radii:
  base=safe(C,r)
  for d in range(1,7):
   p=next((p for p in range(2,d+1) if d%p==0),None)
   threshold=F(1,2) if d==1 else F(p-1,2*p)
   for v in range(1,13):
    actual=safe([d*c for c in C]+[v],r)
    t=quotient(base,d,v,r)
    require(bool(actual)==(t is not None),('transfer',C,r,d,v))
    if t is not None: require(all(norm(w*t)>=r for w in [d*c for c in C]+[v]),('witness',C,r,d,v,t))
    counts['literal_comb_transfer']+=1
    if base and r<=threshold:
     expected=True if v%d else bool(safe(list(C)+[v//d],r))
     require(bool(actual)==expected,('rejection',C,r,d,v))
     counts['rejection_identity']+=1
require(not safe([1,2],F(2,5)), 'global high-radius hostile')
require(norm(F(1,4))==max(norm((F(1,2)+j)/2) for j in range(2)), 'parity hostile')
require(safe([1],F(1,2))==[(F(1,2),F(1,2))], 'singleton')
for h in range(1,31):
 v=1 if h%2==0 else 2
 require(max(norm(v*(F(1,2)+j)/h) for j in range(h))==F(h-1,2*h),('sharp safe fibre',h))
bank=json.loads((Path(__file__).resolve().parents[1]/'lrc/source_bank.json').read_text())
r=F(bank['target_radius']); rows=[]
for entry in bank['cores']:
 C=entry['core']; base=safe(C,r)
 require(bool(base),'empty bank core')
 width=max(R-L for L,R in base)
 require(width>0,'no bank positive interval')
 cutoff=ceil(2*r/width)
 require(cutoff==entry['all_speeds_from'],('bank cutoff',C,cutoff))
 require(len(base)==entry['components'],('bank components',C,len(base)))
 for v in range(1,cutoff):
  require(bool(intersect(base,comb(v,r))),('bank rejection',C,v))
  counts['bank_finite_queries']+=1
 rows.append({'core':C,'components':len(base),'derived_cutoff':cutoff})
require(len(rows)==28 and counts['bank_finite_queries']==1458,'bank counts')
print(json.dumps({'status':'PASS','method':'independent standard-library Fraction literal closed combs; no producer imports','counts':counts,'bank':rows,'scope':'28 inherited core shapes and their dilations only; arbitrary-core LRC14 remains OPEN'},indent=2,sort_keys=True))
