from fractions import Fraction

def v2(n):
 return (abs(n)&-abs(n)).bit_length()-1

def v2q(x):
 return v2(x.numerator)-v2(x.denominator)

checks=0
for d in range(1,8):
 alpha=1-2*Fraction(4,3)**d
 if alpha.denominator!=3**d: raise AssertionError('center denominator')
 for t in range(1,9):
  H=1+2*3**(d-1)*t
  n=1+Fraction(4,3)**d*(2**H-2)
  if n.denominator!=1 or int(n)%2!=1: raise AssertionError('integral')
  n=int(n)
  x=n
  for j in range(d):
   z=3*x+1
   if v2(z)!=2: raise AssertionError('word')
   x=z//4
  if x!=2**H-1: raise AssertionError('endpoint')
  if v2q(Fraction(n)-alpha)!=H+2*d: raise AssertionError('source budget')
  checks+=1
 print('preperiodic_center',d,str(alpha),'first_H',1+2*3**(d-1))
print('pullback_family_checks',checks)
