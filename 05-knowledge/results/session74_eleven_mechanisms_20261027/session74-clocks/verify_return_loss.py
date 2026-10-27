#!/usr/bin/env python3
"""Independent finite additive-quotient control for the Smith return-loss formula."""
from collections import deque,Counter
from math import gcd
import json
count=0;rows=[];multiplication_checks=0
def mul(x,y):
 a,b=x;c,d=y
 return a*c+b*d,a*d+b*c+b*d
def power(n):
 z=(1,0)
 for _ in range(n):z=(z[1],z[0]+z[1])
 return z
def residue(x,n):
 a,b=power(n);a-=1;d=abs(a*(a+b)-b*b)
 u,v=x
 return ((a+b)*u-b*v)%d,(-b*u+a*v)%d
def scale(c,x):return c*x[0],c*x[1]
for m in range(1,19):
 f0,f1=0,1
 for _ in range(m):f0,f1=f1,f0+f1
 a,b=f1-f0-1,f0
 n=abs(a*(a+b)-b*b);d=gcd(a,b)
 steps=(((a+b)%n,(-b)%n),((-b)%n,a%n))
 seen={(0,0)};todo=deque([(0,0)])
 while todo:
  x=todo.popleft()
  for s in steps:
   y=((x[0]+s[0])%n,(x[1]+s[1])%n)
   if y not in seen:seen.add(y);todo.append(y)
 if len(seen)!=n:raise RuntimeError('quotient size')
 tests=[]
 for r in range(1,23):
  image=Counter(((r*x)%n,(r*y)%n) for x,y in seen)
  predicted=gcd(r,d)*gcd(r,n//d)
  if set(image.values())!={predicted} or len(image)*predicted!=n:raise RuntimeError('Smith loss mismatch')
  count+=1;tests.append([r,predicted,len(image)])
  q=(0,0)
  for j in range(r):
   p=power(m*j);q=(q[0]+p[0],q[1]+p[1])
  if residue(mul(q,q),m*r)!=residue(scale(r,q),m*r):raise RuntimeError('Q squared twist')
  multiplication_checks+=1
  for x in ((1,0),(0,1)):
   for y in ((1,0),(0,1)):
    if residue(mul(mul(q,x),mul(q,y)),m*r)!=residue(scale(r,mul(q,mul(x,y))),m*r):raise RuntimeError('bilinear twist')
    multiplication_checks+=1
 rows.append({'m':m,'smith':[d,n//d],'order':n,'loss_tests_r_kernel_image':tests})
q=(0,0)
for j in range(11):
 p=power(5*j);q=(q[0]+p[0],q[1]+p[1])
for a in range(11):
 for b in range(11):
  if residue(mul(scale(a,q),scale(b,q)),55)!=(0,0):raise RuntimeError('eleven-element square-zero ideal')
  multiplication_checks+=1
print(json.dumps({'status':'PASS','tests':count,'multiplication_checks':multiplication_checks,'square_zero_ideal_order':11,'rows':rows},sort_keys=True,indent=2))
