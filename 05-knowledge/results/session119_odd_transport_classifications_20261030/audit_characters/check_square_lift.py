"""Independent square-action, antipodal ANF, and literal N231 controls."""
from math import gcd
from itertools import product
from collections import Counter
import json

def req(c,m):
 if not c:raise RuntimeError(m)
def fac(n):
 out={};p=2
 while p*p<=n:
  while n%p==0:out[p]=out.get(p,0)+1;n//=p
  p+=1
 if n>1:out[n]=1
 return out
def leg(x,p):
 req(x%p!=0,'Legendre unit')
 return 1 if pow(x%p,(p-1)//2,p)==1 else -1
moduli=eligible=0;first_strict=None
for N in range(3,302,2):
 moduli+=1;F=fac(N);U=[u for u in range(N) if gcd(u,N)==1];sq={u*u%N for u in U}
 remaining=set(range(1,N));orbits=[]
 while remaining:
  x=min(remaining);O=frozenset(s*x%N for s in sq);orbits.append(O);remaining-=O
 fixed=sum(O==frozenset((-x)%N for x in O) for O in orbits)
 possible=all(p%4==3 for p in F)
 req((fixed==0)==possible,'full odd square-invariance existence')
 if not possible:continue
 eligible+=1;C=len(orbits)//2
 pred=(__import__('math').prod(1+2*a for a in F.values())-1)//2
 req(C==pred,'square-class exponent')
 L=sum(len(fac(d)) for d in range(2,N+1) if N%d==0)
 closed=sum(a*__import__('math').prod(b+1 for q,b in F.items() if q!=p) for p,a in F.items())
 req(L==closed,'local character exponent')
 req((C==L)==(len(F)<=2),'equality boundary')
 if C>L and first_strict is None:first_strict=N
req(first_strict==231,'first admissible difference')
# Literal complete census of antipodal Boolean functions on three bits.
counts=Counter();tables=[]
for word in range(256):
 f=[word>>x&1 for x in range(8)]
 if not all(f[x]^f[x^7]==1 for x in range(8)):continue
 coeff=f[:]
 for bit in [1,2,4]:
  for x in range(8):
   if x&bit:coeff[x]^=coeff[x^bit]
 req(coeff[7]==0,'no cubic term')
 b=coeff[3];req(coeff[5]==coeff[6]==b,'triangle quadratic part')
 req((coeff[1]^coeff[2]^coeff[4])==(1^b),'linear parity')
 counts['majority_quadratic' if b else 'odd_affine']+=1;tables.append(word)
req(len(tables)==16 and dict(counts)=={'odd_affine':8,'majority_quadratic':8},'8+8 classification')
for signs in product((-1,1),repeat=3):
 total=sum(signs);poly=(total-signs[0]*signs[1]*signs[2])//2
 req(poly==(1 if total>0 else -1),'majority Fourier polynomial')
# Build a full N231 tournament symbol, nonlinear only on the top-order stratum.
N=231;U=[u for u in range(N) if gcd(u,N)==1];squares={u*u%N for u in U};sigma=[0]*N
for d in range(1,N):
 m=N//gcd(d,N);t=d//(N//m);ps=sorted(fac(m));vals=[leg(t,p) for p in ps]
 sigma[d]=(1 if sum(vals)>0 else -1) if m==N else vals[0]
req(all(sigma[d] in (-1,1) and sigma[-d]==-sigma[d] for d in range(1,N)),'full odd symbol')
checks=0
for u in squares:
 for d in range(N):req(sigma[u*d%N]==sigma[d],'literal square covariance');checks+=1
# Show a failure of even the top-stratum signed character rule (sigma(1)=1).
u=next(u for u in U if tuple(leg(u,p) for p in [3,7,11])==(-1,1,1))
v=next(u for u in U if tuple(leg(u,p) for p in [3,7,11])==(1,-1,1))
req(sigma[1]==1 and sigma[u*v%N]!=sigma[u]*sigma[v],'noncharacter witness')
for x in range(N):
 for y in range(x+1,N):req(sigma[(y-x)%N]==-sigma[(x-y)%N],'literal tournament pair')
print(json.dumps({'status':'PASS','odd_moduli_tested':moduli,'eligible_moduli':eligible,'first_strict_modulus':first_strict,'three_bit_antipodal_tables':len(tables),'three_bit_classes':dict(counts),'N231':{'square_class_exponent':13,'local_character_exponent':12,'square_invariant_symbols':8192,'locally_character_symbols':4096,'unit_squares':len(squares),'covariance_checks':checks,'tournament_pairs':N*(N-1)//2,'noncharacter_witness':{'u':u,'v':v,'uv_mod_N':u*v%N,'sigma_u':sigma[u],'sigma_v':sigma[v],'sigma_uv':sigma[u*v%N]}}},indent=2,sort_keys=True))
