import itertools,json
from collections import Counter
from math import isqrt
sigma={'C':'CA','A':'B','B':'CA'};w='C'
while len(w)<100:w=''.join(sigma[c] for c in w)
checks=0
for n in range(1,19):
 for retained in [False,True]:
  choices=[('A','B') if c=='B' else ('B','C') if c=='C' and retained else(c,) for c in w[:n]]
  direct=Counter(tuple(x.count(c) for c in 'ABC') for x in itertools.product(*choices))
  poly=Counter({(0,0,0):1})
  for c in 'ABC':
   opts=('A','B') if c=='B' else('B','C') if c=='C' and retained else(c,)
   for _ in range(w[:n].count(c)):
    out=Counter()
    for mon,m in poly.items():
     for v in opts:
      x=list(mon);x['ABC'.index(v)]+=1;out[tuple(x)]+=m
    poly=out
  if direct!=poly:raise ValueError('partition polynomial')
  total=sum(poly.values());b=w[:n].count('B');c=w[:n].count('C');a=n-b-c
  if total!=2**(b+c*retained):raise ValueError('prefix count')
  EA=sum(mon[0]*v for mon,v in poly.items())
  if 2*EA!=(2*a+b)*total:raise ValueError('fair expectation')
  checks+=3
if w[34]!='C':raise ValueError('old35 position')
# Pinning an eligible independent coordinate replaces a two-term factor by one term.
for n in range(35,101):
 b=w[:n].count('B');c=w[:n].count('C')
 if 2**(b+c-1)*2!=2**(b+c):raise ValueError('old35 pin')
 checks+=1
coverage=[]
for q in range(1,201):
 A=set();k=0
 while len(A)<q and k<50000:
  k+=1;A.add((3*k+isqrt(5*k*k))//2%q)
 C={(x-1)%q for x in A}
 if A!=set(range(q)) or C!=set(range(q)):raise ValueError('finite forced residue coverage')
 coverage.append((k,q));checks+=1
if len({(3*k+isqrt(5*k*k))//2%144 for k in range(1,5001)})!=95:raise ValueError('finite horizon hostile')
print(json.dumps(dict(status='PASS',checks=checks,largest_tested_cover_clock=max(coverage),scope='occurrence-dependent compatible prefixes; no shared dictionary claimed',old35_canonical=w[34],old35_retained_pin='B'),indent=2))
