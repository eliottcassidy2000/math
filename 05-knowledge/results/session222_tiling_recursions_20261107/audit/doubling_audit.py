#!/usr/bin/env python3
"""Independent active (including python -O) oriented-arc doubling audit."""
from itertools import combinations,product
import json
from pathlib import Path

def need(x):
 if not x:raise RuntimeError('audit failure')
def double(E,n):
 # Directed arcs specified by sets, independently of signed matrix code.
 out=set(E)|{(n+j,n+i) for i,j in E}
 for i in range(n):
  out.add((i,n+i))
  for j in range(n):
   if i==j:continue
   out.add((i,n+j) if (i,j) in E else (n+j,i))
 return out

def bits(E,p):return {(i,j):int((p[j],p[i]) in E) for i in range(len(p)) for j in range(i+2,len(p))}
E=set();p=[0];tower=[]
for k in range(6):
 n=len(p);b=bits(E,p);w=sum(b.values());tower.append({'order':n,'backward_weight':w,'forward_weight':len(b)-w})
 if k==5:break
 D=double(E,n);q=p+[n+i for i in reversed(p)];B=bits(D,q)
 need(all((a,z) in D for a,z in zip(q,q[1:])))
 need(sum(B.values())==2*w+n*(n-1)//2)
 need(len(B)-sum(B.values())==2*(len(b)-w)+n-1+n*(n-1)//2)
 for i in range(n-1):need(B[i,2*n-1-i]==0)
 E,p=D,q
cases=0
for n in range(1,7):
 free=[(i,j) for i in range(n) for j in range(i+2,n)]
 for vals in product([0,1],repeat=len(free)):
  E={(i,i+1) for i in range(n-1)}
  E|={(j,i) if v else (i,j) for (i,j),v in zip(free,vals)}
  D=double(E,n);p=list(range(n))+list(range(2*n-1,n-1,-1));b=bits(D,p)
  c={key:int(key[0]<n<=key[1] and key[1]!=2*n-1-key[0]) for key in b}
  need(sum(c.values())==n*(n-1))
  for (i,j),v in b.items():
   sym=(2*n-1-j,2*n-1-i)
   need(v^b[sym]==c[i,j]);need(c[sym]==c[i,j])
  need(sum(b.values())==2*sum(vals)+n*(n-1)//2)
  cases+=1
out={'status':'PASS','framed_cases_orders1_through6':cases,'tower':tower,'syndrome_c_is_reflection_fixed':True,'nonblue_condition':'n>=2 only'}
Path(__file__).with_name('DOUBLING_RESULT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
