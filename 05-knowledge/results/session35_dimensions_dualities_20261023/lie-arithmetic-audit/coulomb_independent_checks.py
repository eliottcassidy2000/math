from fractions import Fraction as Q
from math import isqrt,gcd
from pathlib import Path
import json

def sqrtrational(x):
 if x<0:return None
 a,b=isqrt(x.numerator),isqrt(x.denominator)
 return Q(a,b) if a*a==x.numerator and b*b==x.denominator else None

def charge(t):return Q(1)-6*(t-1)**2/t
rows=[]
# Independently traverse central charges, not author parameter enumeration.
for den in range(1,21):
 for num in range(-200,601):
  c=Q(num,den)
  if gcd(abs(num),den)!=1:continue
  d=(c-1)*(c-25);s=sqrtrational(d)
  if s is not None:
   roots=[(13-c+s)/12,(13-c-s)/12]
   assert roots[0]*roots[1]==1
   assert all(charge(t)==c for t in roots)
   rows.append([str(c),*[str(t) for t in roots]])
integer_cs=[c for c in range(-1000,1001) if sqrtrational(Q((c-1)*(c-25))) is not None]
# Exhaustiveness follows |13-c|=(a+b)/2 with positive same-parity factors ab144.
exact=set()
for a in range(1,13):
 if 144%a==0:
  b=144//a
  if (a-b)%2==0:exact.update([13-(a+b)//2,13+(a+b)//2])
assert sorted(exact)==integer_cs==[-24,-7,-2,0,1,25,26,28,33,50]
assert charge(Q(2,5))+charge(Q(3,5))==-5
assert 5*charge(Q(2,5))+34*charge(Q(3,4))==-5
assert 114*charge(Q(3,4))==57
# Rational t=4/3 exact mirror identities without logarithm approximations.
x=Q(64,27);y=Q(256,81)
assert y/x==Q(4,3) and x**4==y**3
assert 2**4==4**2
out={'rational_charge_inversions_checked':len(rows),'all_integer_c_rational_parameters':integer_cs,'small_c_minus5_tensor_charges':[str(charge(Q(2,5))),str(charge(Q(3,5)))],'mirror_Ising_exact':[str(x),str(y)],'scope':'Exact finite algebra controls; algebraicity and Gelfond-Schneider conclusions audited as proofs, not inferred numerically.'}
Path('lie-arithmetic-audit/coulomb_independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
