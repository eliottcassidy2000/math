#!/usr/bin/env python3
"""Exact controls for the inverse factorization and divisor-lattice theorems."""
import json
from collections import defaultdict

def need(ok, msg):
 if not ok:raise RuntimeError(msg)
def F(n):
 c=n&-n;return n-n//c+2*c-1
def inverse(m):
 out=[];c=2
 while 3*(c-1)<=m-1:
  if (m-1)%(c-1)==0:out.append(c*((m-1)//(c-1)-2))
  c*=2
 return sorted(out)
forward=defaultdict(list)
for n in range(2,10001,2):forward[F(n)].append(n)
for m in range(2,5001,2):need(inverse(m)==forward[m],('inverse',m))
checks=2500
for A in range(1,301):
 C=1<<A;m=C*(2*C-1)
 want=sorted((1<<a)*((m-1)//((1<<a)-1)-2) for a in range(1,A+1) if A%a==0)
 need(inverse(m)==want,('fixed basin',A));checks+=1
 for n in want:need(F(n)==m,('fixed image',A,n))
 if A>=2:need(m-3==(C+1)*(2*C-3),('composite',A))
for b in range(1,301):
 m=1<<b
 want=sorted((1<<a)*((m-1)//((1<<a)-1)-2) for a in range(1,b) if b%a==0)
 need(inverse(m)==want,('power basin',b));checks+=1
 for n in want:need(F(n)==m,('power image',b,n))
need(inverse(496)==[496,652,986],'496')
print(json.dumps({'status':'PASS','complete_inverse_targets':2500,'fixed_point_exponents':300,'power_exponents':300,'aggregate_cases':checks},sort_keys=True,indent=2))
