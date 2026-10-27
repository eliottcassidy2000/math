"""Independent source-by-source replay of the scalar all-source upper bound."""
from fractions import Fraction as F
from math import factorial
import json
s=F(0); e=F(0); max_steps=0
for n in range(3,5002,2):
 x=n;t=0;k=0
 while x!=1:
  y=3*x+1;a=(y&-y).bit_length()-1;x=y>>a;t+=1;k+=(a-1)//2
  if t>10000: raise RuntimeError('finite replay cap reached')
 w=F(2*factorial(k)*factorial(t),factorial(t+k+1));s+=w;e+=w/F(6*n+1);max_steps=max(max_steps,t)
u=e+(F(13,3)-s)/F(6*5003+1)
if not e<u<F(1,42):raise RuntimeError('budget bound failed')
scale=10**10
lo=(e.numerator*scale)//e.denominator
hi=-((-u.numerator*scale)//u.denominator)
if not F(lo,scale)<=e<u<=F(hi,scale)<F(1,42):raise RuntimeError('rational rounded certificate failed')
print(json.dumps({'status':'PASS','cutoff':5001,'sources':2500,'maximum_odd_steps':max_steps,'paid_lower_rational':str(F(lo,scale)),'paid_upper_rational':str(F(hi,scale)),'strict_upper_1_over_42':True,'scope':'independent direct first-hit replay; inherited total nonroot W<=13/3 used for all larger sources'},indent=2,sort_keys=True))
