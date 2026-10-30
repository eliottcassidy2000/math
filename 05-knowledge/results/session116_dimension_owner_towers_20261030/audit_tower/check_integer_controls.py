"""Finite checks only; universal conclusions use classical square theorems."""
import itertools,json,math
from collections import Counter
limit=64;B=math.isqrt(limit);r={};witness={};gates=0
def check(ok):
 global gates
 gates+=1
 if not ok:raise RuntimeError(gates)
for k in (3,4):
 counts=Counter()
 for v in itertools.product(range(-B,B+1),repeat=k):
  m=sum(t*t for t in v)
  if m<=limit:
   counts[m]+=1
   if k==4:witness[m]=v
 r[k]=counts
for m in range(1,limit+1):
 n=m
 while n%4==0:n//=4
 check((r[3][m]>0)==(n%8!=7))
 check(r[4][m]==8*sum(d for d in range(1,m+1) if m%d==0 and d%4))
 check(m in witness)
 p,q,t,s=witness[m];a=(p,q,t,s,0,0,0);b=(-q,p,-s,t,0,0,0)
 check(sum(x*x for x in a)==m==sum(x*x for x in b))
 check(sum(x*y for x,y in zip(a,b))==0)
 N=sum(x*x for x in a+b);Delta=m*m
 check(N==2*m and N*N==4*Delta)
check(r[3][1]*r[4][1]==48);check(r[3][7]*r[4][7]==0)
a=(2,1,1,1,0,0,0);b=(-1,2,-1,1,0,0,0)
check(sum(x*x for x in a)==7==sum(x*x for x in b));check(sum(x*y for x,y in zip(a,b))==0)
print(json.dumps({'status':'PASS','gates':gates,'shell_limit':limit,'W_norm2_count':r[3][1]*r[4][1],'W_norm14_count':r[3][7]*r[4][7]},sort_keys=True))
