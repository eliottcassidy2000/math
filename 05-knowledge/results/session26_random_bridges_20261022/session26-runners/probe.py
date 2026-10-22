from fractions import Fraction as F
from itertools import combinations
from collections import defaultdict
import json

def segments(v):
 walls=sorted({F(k,2*a) for a in v for k in range(2*a+1)})
 out=[]
 for l,r in zip(walls,walls[1:]):
  mid=(l+r)/2; lines=[]
  for a in v:
   q=(a*mid).__floor__(); frac=a*mid-q
   lines.append((a,-q) if frac<F(1,2) else (-a,q+1))
  for i,(s,b) in enumerate(lines):
   L,R=l,r
   for j,(u,c) in enumerate(lines):
    if i==j:continue
    if s==u:
     if b>c or (b==c and i>j):L=R;break
    elif s>u:R=min(R,F(c-b,s-u))
    else:L=max(L,F(c-b,s-u))
   if L<R:out.append((L,R,s,b))
 out.sort()
 merged=[]
 for l,r,s,b in out:
  if merged and merged[-1][1]==l and merged[-1][2:]==(s,b):merged[-1]=(merged[-1][0],r,s,b)
  else:merged.append((l,r,s,b))
 if sum(r-l for l,r,s,b in merged)!=1:raise ValueError('partition')
 return merged

def curve(v):
 ss=segments(v); atoms=defaultdict(F)
 for l,r,s,b in ss:
  lo,hi=sorted((s*l+b,s*r+b));atoms[hi]+=F(1,abs(s));atoms[lo]-=F(1,abs(s))
 return ss,{h:c for h,c in sorted(atoms.items()) if c}

def run():
 rows=[];count=0;hostile=None
 for n in range(1,5):
  for v in combinations(range(1,9),n):
   ss,a=curve(v);count+=1
   if hostile is None and any(h>0 and c<0 for h,c in a.items()):hostile=(v,a)
   M=max(a)
   if a[M]<=0:raise ValueError('terminal sign')
   if sum(c*h for h,c in a.items())!=1:raise ValueError('first moment')
   for p in range(5):
    # Independent time integration of affine segments.
    integ=sum(((s*r+b)**(p+1)-(s*l+b)**(p+1))/F(s*(p+1)) for l,r,s,b in ss)
    mom=sum(c*h**(p+1) for h,c in a.items())/F(p+1)
    if integ!=mom:raise ValueError('moment')
   if n<=3:
    _,a2=curve(tuple(2*x for x in v))
    if a!=a2:raise ValueError('dilation')
 # Positive controls, canonical tight boundary, refuted historic deep-well uniqueness.
 for v in [(1,2),(1,2,3),tuple(range(1,14)),tuple(range(1,13))+(182,)]:
  ss,a=curve(v)
  rows.append({'speeds':v,'segments':len(ss),'knots':len(a),'maximum':str(max(a)),'terminal_residue':str(a[max(a)]),'positive_knots':sum(c>0 for c in a.values()),'negative_positive_height_knots':sum(h>0 and c<0 for h,c in a.items())})
 print(json.dumps({'universe':'all nonempty subsets of {1,...,8} of cardinality <=4','rows':count,'first_signed_hostile':{'speeds':hostile[0],'atoms':{str(k):str(v) for k,v in hostile[1].items()}},'controls':rows},indent=2))
if __name__=='__main__':run()
