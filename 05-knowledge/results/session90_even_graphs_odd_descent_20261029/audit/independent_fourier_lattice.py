#!/usr/bin/env python3
"""Independent explicit Walsh invariant lattice, through n=3, no producer imports."""
from independent_decorated_counts import setup
from math import prod,gcd
import json
C=0
def ck(x,msg):
 global C
 C+=1
 if not x:raise RuntimeError(msg)
def det(a):
 a=[r[:] for r in a];sign=1;last=1
 for k in range(len(a)-1):
  if not a[k][k]:
   r=next(i for i in range(k+1,len(a)) if a[i][k]);a[k],a[r]=a[r],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,len(a)):
   for j in range(k+1,len(a)):
    x=a[i][j]*pivot-a[i][k]*a[k][j]
    ck(x%last==0,'exact Bareiss division');a[i][j]=x//last
   a[i][k]=0
  last=pivot
 return sign*a[-1][-1]
def run(n):
 es,acts=setup(n);m=len(es);q=1<<m;pack=lambda t,r,d:t|(r<<m)|(d<<(2*m));parts=lambda z:(z&(q-1),(z>>m)&(q-1),z>>(2*m))
 src=[];srcidx={};srcaut=[];cols=[];targ=[];taut=[];seen=set()
 for z in range(1<<(2*m+n)):
  t,r,d=parts(z)
  if z not in srcidx:
   orb={pack(ep[t]^c,ep[r],vp[d]) for p,ep,vp,c,rev in acts};i=len(src);src.append(z);srcaut.append(len(acts)//len(orb))
   for w in orb:srcidx[w]=i
  if z not in seen:
   signed={};bad=False
   for p,ep,vp,c,rev in acts:
    w=pack(ep[t],ep[r],vp[d]);sgn=(-1)**((t&rev).bit_count()%2)
    if w in signed and signed[w]!=sgn:bad=True
    signed[w]=sgn
   seen.update(signed)
   if not bad:cols.append(signed);targ.append(z);taut.append(len(acts)//len(signed))
 s=len(src);ck(s==len(cols),'square matrix')
 M=[]
 for z in src:
  t,r,d=parts(z);row=[]
  for col in cols:
   value=0
   for w,sgn in col.items():
    x,rr,dd=parts(w)
    if (rr,dd)==(r,d):value+=sgn*(-1)**((t&x).bit_count()%2)
   row.append(value)
  M.append(row)
 N=[]
 for w in targ:
  x,r,d=parts(w);row=[0]*s
  for t in range(q):row[srcidx[pack(t,r,d)]]+=(-1)**((t&x).bit_count()%2)
  N.append(row)
 for i in range(s):
  for j in range(s):
   ck(sum(M[i][k]*N[k][j] for k in range(s))==q*(i==j),'MN scalarI')
   ck(sum(N[i][k]*M[k][j] for k in range(s))==q*(i==j),'NM scalarI')
 determinant=det(M);ck(determinant*determinant*prod(taut)==(1<<(m*s))*prod(srcaut),'weighted determinant identity')
 ck(abs(determinant)&(abs(determinant)-1)==0,'dyadic determinant')
 ngcd=q
 for row in N:
  for a in row:ngcd=gcd(ngcd,abs(a))
 return dict(n=n,rank=s,det_sign=1 if determinant>0 else -1,det_v2=abs(determinant).bit_length()-1,cokernel_exponent=q//ngcd,source_group_product=prod(srcaut),target_group_product=prod(taut))
def main():print(json.dumps(dict(status='PASS',rows=[run(n) for n in range(4)],checks=C),indent=2))
if __name__=='__main__':main()
