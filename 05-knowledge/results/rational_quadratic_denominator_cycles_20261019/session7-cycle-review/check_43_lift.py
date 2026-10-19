#!/usr/bin/env python3
"""Independent cycle decomposition of the full mod43^2/mod43^3 graphs,
and exact local Taylor verification at the rational three-cycle anchors."""
from collections import Counter
from fractions import Fraction
import json

def cycles(mod):
 inv=pow(4,-1,mod)
 image=[((x*x-29)*inv)%mod for x in range(mod)]
 seen=set();out=[]
 for root in range(mod):
  if root in seen:continue
  walk=[];index={};x=root
  while x not in seen and x not in index:
   index[x]=len(walk);walk.append(x);x=image[x]
  if x in index:out.append(tuple(walk[index[x]:]))
  seen.update(walk)
 return out
p=43;c2=cycles(p*p);c3=cycles(p**3)
six2={x for cyc in c2 if len(cyc)==6 for x in cyc}
if len(six2)!=126:raise ValueError('wrong modp2 count')
lifts={x+k*p*p for x in six2 for k in range(p)}
liftcycles=[cyc for cyc in c3 if set(cyc)&lifts]
if set().union(*(set(c) for c in liftcycles))!=lifts:raise ValueError('lift cycles not exhaustive')
if {len(c) for c in liftcycles}!={258}:raise ValueError('unexpected lift period')

def mul(a,b):
 out=[Fraction(0) for _ in range(len(a)+len(b)-1)]
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out

def g_at(alpha):
 a=[Fraction(alpha),Fraction(1)]
 for _ in range(3):
  a=[x/4 for x in mul(a,a)];a[0]-=Fraction(29,4)
 a[0]-=alpha
 return a

def val(x,p):
 if not x:return None
 n=x.numerator;d=x.denominator;r=0
 while n%p==0:n//=p;r+=1
 while d%p==0:d//=p;r-=1
 return r
polys={str(alpha):g_at(alpha) for alpha in (-7,5,-1)}
linear=Fraction(35,8)
for a in polys.values():
 if a[0] or a[1]!=linear:raise ValueError('wrong Taylor center or derivative')
 if any(val(z,p)<0 for z in a if z):raise ValueError('non-integral coefficient')
quad_g2={alpha:a[2]*linear*(linear+1) for alpha,a in polys.items()}
if val(linear**2-1,p)!=1:raise ValueError('wrong linear valuation')
if any(val(z,p)<1 for z in quad_g2.values()):raise ValueError('quadratic not killed')
print(json.dumps(dict(status='PASS',full_graph_mod43_squared_cycle_lengths=dict(sorted(Counter(map(len,c2)).items())),full_graph_mod43_cubed_cycle_lengths=dict(sorted(Counter(map(len,c3)).items())),period6_points_mod43_squared=len(six2),lift_count=len(lifts),lift_cycle_count=len(liftcycles),lift_exact_periods=sorted(set(map(len,liftcycles))),lambda_value=str(linear),lambda_squared_minus1=str(linear**2-1),quadratic_coefficients_g={a:str(v[2]) for a,v in polys.items()},quadratic_coefficients_g_squared= {a:str(v) for a,v in quad_g2.items()}),indent=2))
