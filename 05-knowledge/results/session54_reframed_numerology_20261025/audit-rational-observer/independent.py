import math,json
from pathlib import Path
from itertools import product
r=Path(__file__).resolve().parent
data=json.loads((r/'inputs'/'producer-results.json').read_text())
def req(c,m):
 if not c:raise RuntimeError(m)
def factor(n):
 f={};p=2
 while p*p<=n:
  while n%p==0:f[p]=f.get(p,0)+1;n//=p
  p+=1
 if n>1:f[n]=f.get(n,0)+1
 return f
def direct(n):
 u=[a for a in range(n) if math.gcd(a,n)==1]
 # Explicit quotient equivalence classes, then square their representatives.
 classes={frozenset((a,(-a)%n)) for a in u}
 return len({frozenset((min(c)**2%n,(-min(c)**2)%n)) for c in classes})
def crt_formula(n):
 if n<=2:return 1
 f=factor(n);u2=1;minus_square=True
 for p,a in f.items():
  if p==2:
   u2*=1 if a<3 else 2**(a-3)
   minus_square=minus_square and a<=1
  else:
   u2*=p**(a-1)*(p-1)//2
   minus_square=minus_square and p%4==1
 return u2//(2 if minus_square else 1)
for n in range(3,501):req(crt_formula(n)==direct(n),'global CRT formula '+str(n))
rows=[]
for row in data['bounds']:
 d=row['d']
 # Any admissible odd prime-power q satisfies q<=6d by the local theorem;
 # every admissible two-power q satisfies q<=8d. Enumerate the union bound.
 qs=[q for q in range(2,8*d+1) if len(factor(q))==1 and direct(q)<=d]
 N=math.lcm(*qs)
 req(N==row['N'],'independent LCM '+str(d))
 req(factor(N)=={int(p):a for p,a in row['prime_exponents'].items()},'factorization')
 ps=list(factor(N));es=list(factor(N).values());passing=[]
 for v in product(*[range(e+1) for e in es]):
  n=math.prod(p**a for p,a in zip(ps,v))
  if n>=3 and crt_formula(n)<=d:passing.append(n)
 passing.sort()
 req(sum(n<=500 for n in passing)==row['passing_count'],'producer bounded count')
 rows.append({'d':d,'independent_N':N,'admissible_prime_powers':qs,'all_square_group_survivors':passing,'total_survivors':len(passing),'max_survivor':max(passing)})
# Independent elimination of t from t²+t−1 and J=a*t+b.
# (J-b)²+a*(J-b)−a² = J²+(a−2b)J+b²−ab−a².
a,b=85995,-52515
req(a-2*b==191025 and b*b-a*b-a*a==-121287375,'affine CM elimination')
disc=191025**2+4*121287375
req(math.isqrt(disc)**2!=disc,'irreducibility H15')
# t=2cos(2pi/5), not cos: roots have trace−1 and norm−1.
req(direct(5)==1,'d1 n5')
req(120%15==0 and direct(15)==2,'local-global obstruction')
print(json.dumps({'status':'PASS','independent_bounds':rows,'H15_affine_elimination':True,'H15_irreducible':True},indent=2))
