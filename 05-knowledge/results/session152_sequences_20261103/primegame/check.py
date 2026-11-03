from math import lcm
from fractions import Fraction
from itertools import combinations
fs=list(map(Fraction,['17/91','78/85','19/51','23/38','29/33','77/29','95/23','77/19','1/17','11/13','13/11','15/14','15/2','1/7','55/1']))
# Directly evaluate inclusion-exclusion coefficients by mask, checking every priority branch.
def op(n):
 out={}
 for i,f in enumerate(fs):
  for mask in range(1<<i):
   d=f.denominator
   for j in range(i):
    if mask>>j&1: d=lcm(d,fs[j].denominator)
   if n%d==0:
    v=n*f.numerator//f.denominator
    out[v]=out.get(v,0)+(-1)**mask.bit_count()
 return {k:v for k,v in out.items() if v}
n=2
for t in range(25):
 expected=next((n*f.numerator//f.denominator for f in fs if n%f.denominator==0),None)
 if op(n)!=({expected:1} if expected else {}): raise AssertionError((t,n))
 n=expected
print('25 literal orbit steps: inclusion-exclusion agrees with first-applicable rule')
