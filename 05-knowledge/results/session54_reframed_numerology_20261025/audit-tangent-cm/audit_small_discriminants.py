import json, math
from fractions import Fraction as F
from pathlib import Path
d=json.loads((Path(__file__).resolve().parent/'inputs'/'producer-results.json').read_text())
def req(c,m):
 if not c: raise RuntimeError(m)
forms={}
for D in range(-31,-2):
 if D%4 not in (0,1):continue
 found=[]
 # Independent rectangular enumeration (a,c,b), loose a,c<=|D| bounds.
 for a in range(1,-D+1):
  for c in range(a,-D+1):
   q=D+4*a*c
   if q<0:continue
   r=math.isqrt(q)
   if r*r!=q or r>a:continue
   for b in sorted({r,-r}):
    if math.gcd(a,math.gcd(b,c))!=1:continue
    if b<0 and (r==a or c==a):continue
    found.append([a,b,c])
 req(sorted(found)==sorted(d['reduced_forms_abs_discriminant_below32'][str(D)]),'forms '+str(D))
 forms[D]=found
for row in d['discriminant_bounds']:
 B=row['discriminant_abs_strict_upper']; q=F(row['exp_argument_lower'])
 # Same rational lower argument; only 40 Taylor terms, different count.
 L=sum((q**k/F(math.factorial(k))) for k in range(40))
 req(L>row['trace']+2079,'Taylor '+str(B))
req(forms[-12]==[[1,0,3]],'nonmaximal -12')
req(forms[-16]==[[1,0,4]],'nonmaximal -16')
req(forms[-27]==[[1,1,7]],'nonmaximal -27')
req(forms[-28]==[[1,0,7]],'nonmaximal -28')
print(json.dumps({'status':'PASS','independent_forms':forms,'taylor_terms':40},indent=2))
