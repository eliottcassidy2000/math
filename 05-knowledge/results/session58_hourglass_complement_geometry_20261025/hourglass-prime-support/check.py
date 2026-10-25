#!/usr/bin/env python3
"""Exact, bounded cross-check; all guards survive python -O."""
from math import isqrt,prod,gcd
from itertools import product,combinations
from pathlib import Path
import json

def require(b,m):
 if not b: raise RuntimeError(m)
def factor(n):
 out={};p=2
 while p*p<=n:
  while n%p==0: out[p]=out.get(p,0)+1;n//=p
  p+=1
 if n>1:out[n]=out.get(n,0)+1
 return out
def mul(z,w):return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def power(z,k):
 r=(1,0)
 for _ in range(k):r=mul(r,z)
 return r
def pi(p):
 for b in range(1,isqrt(p)+1):
  a=isqrt(p-b*b)
  if a*a+b*b==p:return (a,b)
 raise RuntimeError('split representation')
def val(n,p):
 require(n>0,'valuation zero guard');k=0
 while n%p==0:k+=1;n//=p
 return k
def gaussian(e):
 f=factor(e);split=[(p,a) for p,a in f.items() if p%4==1];m=prod(p**a for p,a in f.items() if p%4!=1);out={}
 for js in product(*(range(-a,a+1) for p,a in split)):
  if not any(js):continue
  if next(j for j in js if j)!=abs(next(j for j in js if j)):continue
  z=(m,0)
  for (p,a),j in zip(split,js):
   w=pi(p);w=w if j>=0 else (w[0],-w[1]);z=mul(z,power(w,2*abs(j)));z=(z[0]*p**(a-abs(j)),z[1]*p**(a-abs(j)))
  d=abs(2*z[0]*z[1]);require(d not in out and 0<d<e*e,'Gaussian injectivity/positivity');out[d]=js
 return out,split,m

def brute(e):
 E=e*e;out=set()
 for a in range(1,e):
  d=E-a*a;b=isqrt(E+d)
  if b*b==E+d:out.add(d)
 return out

def main():
 bound=3000;total=0;schur=[];oneprime=0;balanced_extra=None;sample={}
 for e in range(1,bound+1):
  gs,ps,m=gaussian(e);ds=set(gs);require(ds==brute(e),f'independent enumeration e={e}')
  require(len(ds)==(prod(2*a+1 for p,a in ps)-1)//2,'count')
  base,_,_=gaussian(e//m);require(ds=={m*m*d for d in base},'inert/2 stripping')
  for d,js in gs.items():
   for (p,a),j in zip(ps,js):
    v=val(d,p)
    require(v==2*(a-abs(j)) if j else v>=2*a,'valuation law')
    if not j and v>2*a and balanced_extra is None:balanced_extra={'e':e,'gap':d,'prime':p,'valuation':v,'exponent':a,'indices':js}
  if len(ps)==1:
   oneprime+=1;p,a=ps[0];require(sorted(val(d,p) for d in ds)==list(range(0,2*a,2)),'prime power levels')
  for x,y in combinations(sorted(ds),2):
   if x+y in ds:schur.append([e,x,y,x+y])
  # Exact hereditary gap slice for every prime divisor, including inert and 2.
  for p in factor(e):
   child,_,_=gaussian(e//p);require({d for d in ds if d%(p*p)==0}=={p*p*d for d in child},'hereditary slice')
  total+=len(ds)
  if e in (5,25,65,125,1105):sample[e]={'gap_count':len(ds),'gaps':sorted(ds)}
 require(not schur,'bounded Schur census changed')
 require(balanced_extra=={'e':65,'gap':3000,'prime':5,'valuation':3,'exponent':1,'indices':(0,1)},'balanced-index hostile')
 # Positive/negative controls for unique-min filter (not hourglass existence).
 require(6+10==16 and all([val(6,2),val(10,2),val(16,2)].count(min(val(t,2) for t in (6,10,16)))>=2 for _ in [0]),'valuation positive control')
 require([val(t,5) for t in (24,336,15000)]==[0,0,4],'valuation tie control')
 # Projective indexing for squarefree split centers, with exact line incidence.
 def canonical(v):
  w=tuple(x%3 for x in v)
  require(any(w),'projective zero')
  lead=next(x for x in w if x)
  return tuple((x*lead)%3 for x in w)
 pg={}
 for e in (65,1105):
  gs,ps,m=gaussian(e);bypt={canonical(j):d for d,j in gs.items()};pts=set(bypt)
  require(len(pts)==(3**len(ps)-1)//2,'projective point bijection')
  lines=set()
  for u,v in combinations(pts,2):
   line=frozenset(canonical(tuple(a*x+b*y for x,y in zip(u,v))) for a,b in product(range(3),repeat=2) if a or b)
   require(len(line)==4 and line<=pts,'projective line')
   lines.add(line)
  require(all(sum(u in L and v in L for L in lines)==1 for u,v in combinations(pts,2)),'unique projective line')
  for i,(p,a) in enumerate(ps):
   sub={d for pt,d in bypt.items() if pt[i]==0};child,_,_=gaussian(e//p)
   require(sub=={p*p*d for d in child},'coordinate deletion equals hereditary scaling')
  collinear_schur=sum(any(x+y==z for x,y,z in combinations(sorted(bypt[q] for q in L),3)) for L in lines)
  require(collinear_schur==0,'line does not generate a Schur triple')
  pg[e]={'points':len(pts),'lines':len(lines),'points_per_line':4,'lines_containing_schur_triple':collinear_schur}
 gs,ps,m=gaussian(65);inv={j:d for d,j in gs.items()}
 require((inv[(1,0)],inv[(0,1)],inv[(1,1)])==(4056,3000,3696),'projective addition hostile')
 out={'status':'PASS','center_bound':bound,'positive_gaps_checked':total,'one_split_prime_centers':oneprime,'schur_triples':len(schur),'balanced_index_extra_divisibility':balanced_extra,'samples':sample,'projective_checks':pg}
 Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
