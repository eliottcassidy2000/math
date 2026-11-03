from fractions import Fraction
P=[Fraction(17,143),Fraction(133,85),Fraction(17,19),Fraction(23,17),Fraction(1015,69),Fraction(23,29),Fraction(31,23),Fraction(38073,62),Fraction(31,37),Fraction(41,31),Fraction(86,287),Fraction(41,43),Fraction(13,41)]
def bits(n): return tuple(n%f.denominator==0 for f in P)
def step(n):
  for f in P:
    q=n*f
    if q.denominator==1: return int(q)
  return None
first={}
for n in range(1,1000000):
 b=bits(n); k=step(n)
 if k is None: continue
 if b in first:
  x,y,c=first[b]
  if bits(k)!=c:
   print('input pair',x,n,'same bits',''.join(map(lambda b:str(int(b)),b)),'successors',y,k,'output bits',''.join(map(lambda b:str(int(b)),c)),''.join(map(lambda b:str(int(b)),bits(k))));break
 else:first[b]=(n,k,bits(k))

else: raise AssertionError('expected status collision')
