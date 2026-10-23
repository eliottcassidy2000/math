#!/usr/bin/env python3
from fractions import Fraction as F
from math import gcd,isqrt
from pathlib import Path
import json
c=lambda t:13-6*(t+1/t)
D=lambda c:(c-1)*(c-25)
square=lambda x:x>=0 and isqrt(x.numerator)**2==x.numerator and isqrt(x.denominator)**2==x.denominator
samples=[]
for a in range(-30,31):
 for b in range(1,31):
  if a==0 or gcd(abs(a),b)!=1:continue
  t=F(a,b);v=c(t)
  assert c(1/t)==v and square(D(v))
  assert v-1==-6*(t-1)**2/t
  assert v-25==-6*(t+1)**2/t
  samples.append(t)
# Quadratic arithmetic (a+b sqrt5); t=(3+sqrt5)/2, reciprocal=(3-sqrt5)/2.
mul=lambda x,y:(x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
t=(F(3,2),F(1,2));u=(F(3,2),F(-1,2))
assert mul(t,u)==(1,0)
assert 13-6*(t[0]+u[0])==-5 and t[1]+u[1]==0
assert not square(F(180))
assert c(F(4,3))==F(1,2)
assert F(4,3)**3==F(64,27) and F(4,3)**4==F(256,81)
assert 5*F(-22,5)+34*F(1,2)==-5
assert 114*F(1,2)==57
# Integer c with rational t: difference of squares k^2-m^2=144,
# enumerate ALL positive factor pairs (k-m)(k+m)=144.
integer_cs=set()
for a in range(1,145):
 if 144%a:continue
 b=144//a
 if a>b or (a+b)%2:continue
 k=(a+b)//2
 integer_cs.update((13-k,13+k))
assert sorted(integer_cs)==[-24,-7,-2,0,1,25,26,28,33,50]
# Fusion-ring control: dimensions encoded as a+b sqrt2.
qd={'1':(1,0),'psi':(1,0),'sigma':(0,1)}
mult=lambda x,y:(x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
assert mult(qd['sigma'],qd['sigma'])==(2,0)
assert mult(qd['psi'],qd['psi'])==qd['1']
result={'rational_parameter_samples':len(samples),'integer_c_with_rational_parameter':sorted(integer_cs),'c_minus5_parameter':'(3 +/- sqrt(5))/2','rational_mirror_at_Ising_c':['64/27','256/81'],'rational_VOA_same_c_minus5':{'Lee_Yang_factors':5,'Ising_factors':34},'rational_VOA_c57_Ising_factors':114,'Ising_sigma_quantum_dimension':'sqrt(2)','scope':'Exact algebraic controls; does not computationally prove Gelfond-Schneider or Virasoro rationality classification.'}
Path(__file__).with_name('coulomb_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
