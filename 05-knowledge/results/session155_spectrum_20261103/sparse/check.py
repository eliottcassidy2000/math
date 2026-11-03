def require(x):
 if not x: raise AssertionError("sparse check failed")
#!/usr/bin/env python3
from decimal import Decimal, getcontext
from fractions import Fraction
getcontext().prec=60
lo=Decimal('0.28452464');hi=Decimal('0.28452465')
def f(x):return (-x*Decimal(17).ln()).exp()+(-x*Decimal(8).ln()).exp()-1
require(f(lo)>0 and f(hi)<0)
for _ in range(150):
 mid=(lo+hi)/2
 if f(mid)>0:lo=mid
 else:hi=mid
print('delta approx:',(lo+hi)/2)
# Every action sequence through length 15, with exact integer checks.
layer=[(25,34,0,0,0)]
for length in range(16):
 for Y,X,r,s,parity in layer:
  require(r+s==length)
  require(Fraction(Y,1)>=Fraction(193,8)*17**r*8**s)
  if length:require(2*X>=3*Y)
  if parity==0:
   require(3*Y+20>=Fraction(579,8)*17**r*8**s)
   require(3*Y+22>=Fraction(579,8)*17**r*8**s)
 if length<15:
  layer=[u for Y,X,r,s,p in layer for u in ((11*Y+4*X,18*Y+7*X,r+1,s,p),(5*Y+2*X,8*Y+3*X,r,s+1,p^1))]
print('all 65535 action words through length 15: exact bounds pass')
