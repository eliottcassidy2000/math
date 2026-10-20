#!/usr/bin/env python3
"""All clocks p<=10, p<q<=80, two actual adjacent swaps; exact arithmetic."""
import json,math
from collections import Counter
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise ArithmeticError(msg)
def edits(mask,q):
 for i in range(q):
  j=(i+1)%q
  if ((mask>>i)^(mask>>j))&1:yield i,mask^(1<<i)^(1<<j)
def necklace(mask,q):
 full=(1<<q)-1
 return min(((mask>>i)|(mask<<(q-i)))&full for i in range(q))
def carry(mask,q):
 ones=[i for i in range(q) if (mask>>i)&1]
 return sum(2**i*3**j for j,i in enumerate(reversed(ones)))
def fmt(mask,q):return ''.join(str((mask>>i)&1) for i in range(q))
rows=[];hits=[];totalpaths=totaloutputs=0;offcritical=0;maxstrict={}
for p in range(1,11):
 for q in range(p+1,81):
  D=2**q-3**p;sign='negative' if D<0 else 'positive'
  criticalq=(3**p).bit_length()-(D<0);critical=q==criticalq
  c=sum((((i+1)*p)//q-(i*p)//q)<<i for i in range(q))
  n0=necklace(c,q);one=list(edits(c,q));n1={necklace(v,q) for i,v in one}
  final={};paths=0
  for i,u in one:
   for j,w in edits(u,q):final.setdefault(w,(i,j));paths+=1
  maxg=0;strict=0;nhits=0
  for w,path in final.items():
   B=carry(w,q);g=math.gcd(abs(D),B);maxg=max(maxg,g)
   n=necklace(w,q);dist=0 if n==n0 else 1 if n in n1 else 2;strict+=dist==2
   if g==abs(D):
    nhits+=1;hits.append(dict(p=p,q=q,word=fmt(w,q),path=path,fixed_point=B//D,necklace_distance=dist,critical=critical))
    need(q<=2*p if D>0 else q<(3**p).bit_length(),'universal clock gate')
   if dist==2:
    key=('primitive' if math.gcd(p,q)==1 else 'nonprimitive')+'_'+sign
    if key not in maxstrict or g>maxstrict[key]['content']:
     maxstrict[key]=dict(p=p,q=q,content=g,path=path,word=fmt(w,q),critical=critical)
  totalpaths+=paths;totaloutputs+=len(final);offcritical+=not critical
  rows.append(dict(p=p,q=q,clock=D,critical=critical,paths=paths,distinct_outputs=len(final),exact_necklace_distance_two_outputs=strict,max_content=maxg,integer_hits=nhits))
need(all(h['necklace_distance']==0 for h in hits),'strict distance2 integer hit')
print(json.dumps(dict(status='FINITE_EXACT',universe='1<=p<=10, p<q<=80, all zero-intercept mechanical shapes and two actual cyclic adjacent unequal swaps',checks=checks,shape_count=len(rows),off_critical_shapes=offcritical,paths=totalpaths,distinct_outputs_by_shape=totaloutputs,integer_hits=hits,maximum_exact_distance_two_contents=maxstrict,rows=rows),indent=2))
