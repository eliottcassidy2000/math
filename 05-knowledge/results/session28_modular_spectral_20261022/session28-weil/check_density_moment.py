#!/usr/bin/env python3
from fractions import Fraction as F
import json
# Input table pi(10^k), k=0..11; input is separately sourced by parent lane.
pi=[0,4,25,168,1229,9592,78498,664579,5761455,50847534,455052511,4118054813]
def require(condition):
 if not condition: raise RuntimeError("Exact certificate failed")

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
q=[0,1]
for a in [[-1,1],[-1,1],[1,1],[-1,2]]:q=mul(q,a)
sq=mul(q,q)
delta=[F(pi[k]-pi[k-1],9*10**(k-1)) for k in range(1,12)]
value=sum(c*delta[i] for i,c in enumerate(sq))
require(value==F(-2447633,22500000000))
# Control: uniform measure moments 1/(n+1) must be positive.
control=sum(F(c,i+1) for i,c in enumerate(sq))
require(control>0)
print(json.dumps({'input_pi':pi,'delta':[str(x) for x in delta],'q_coefficients':q,'q_squared_coefficients':sq,'prime_density_L_q2':str(value),'uniform_measure_L_q2':str(control)},indent=2,sort_keys=True))
