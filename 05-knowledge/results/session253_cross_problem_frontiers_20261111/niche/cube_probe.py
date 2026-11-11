from math import isqrt
from collections import defaultdict
import json

def factor(n):
 d={};p=2
 while p*p<=n:
  while n%p==0:d[p]=d.get(p,0)+1;n//=p
  p+=1
 if n>1:d[n]=1
 return d

def D(n):
 fs=factor(n);a=fs.get(3,0);q,r=divmod(a,3)
 if a==1:raise Exception('impossible 3-adic branch')
 out=3**(q+(0 if r==0 else 1 if r==2 else 2))
 for p,a in fs.items():
  if p%3==2:out*=p**(a-2*(a//3))
 return out
reps=defaultdict(list)
for y in range(2,121):
 for x in range(1,y):reps[x**3+y**3].append((x,y))
for n,rs in reps.items():
 d=D(n)
 for x,y in rs:
  if (x+y)%d or not n<(x+y)**3<4*n:raise Exception((n,d,x,y))
 # Even incomplete coordinate boxes must satisfy the upper bound on all reps.
 cap=0;j=1
 while (j*d)**3<4*n:
  if n<(j*d)**3:cap+=1
  j+=1
 if len(rs)>cap:raise Exception((n,d,rs,cap))
print(json.dumps({'positive_distinct_pairs':7140,'target_values':len(reps),'collision_values_in_box':sum(len(v)>1 for v in reps.values()),'all_pass':True},indent=2))
