import itertools,json,math
from fractions import Fraction

def req(x,m):
 if not x:raise RuntimeError(m)
records=[]
for b,k in [(2,7),(3,5),(6,4),(10,3)]:
 # Includes negative multipliers and offsets, unlike the original checker.
 coeff=[((-1 if a%2 else b+1),-((-1 if a%2 else b+1)*a)+b*(a-2)) for a in range(b)]
 n=b**k;seen=set()
 for w in itertools.product(range(b),repeat=k):
  # Invert the actual branches modulo the needed remaining precision.
  x=0;mod=1
  for a in reversed(w):
   mod*=b;u,v=coeff[a];x=((b*x-v)*pow(u,-1,mod))%mod
  req(x not in seen,'inverse collisions');seen.add(x)
  z=x;observed=[]
  for _ in range(k):
   a=z%b;observed.append(a);u,v=coeff[a];req((u*z+v)%b==0,'integrality');z=(u*z+v)//b
  req(tuple(observed)==w,'inverse/forward mismatch')
 req(len(seen)==n,'missing residue')
 records.append({'b':b,'k':k,'words':n})
x=Fraction(19,101);orb=[]
for _ in range(7):
 orb.append(str(x));a=x.numerator%2;x=((3*x+1) if a else x)/2
req(x==Fraction(19,101),'rational cycle closure')
print(json.dumps({'status':'PASS','inverse_branch_tests':records,'rational_cycle':orb},indent=2))
