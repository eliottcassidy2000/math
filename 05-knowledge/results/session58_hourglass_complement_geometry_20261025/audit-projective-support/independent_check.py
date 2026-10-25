from itertools import product,combinations
from math import isqrt,prod,gcd
from collections import Counter
from pathlib import Path
import json,sys

def check(x,msg):
 if not x:raise RuntimeError(msg)
if '--negative' in sys.argv:check(False,'injected independent failure')
# Independent enumeration uses u^2+v^2=e^2, not endpoint-square or Gaussian routines.
def gaps(e):
 out=set()
 for u in range(1,isqrt(e*e//2)+1):
  v=isqrt(e*e-u*u)
  if v*v+u*u==e*e and u<v:out.add(2*u*v)
 return out
rows=[]
for ps in [(),(5,),(13,),(5,13),(5,17),(5,13,17),(5,13,17,29)]:
 e=prod(ps);D=gaps(e);k=len(ps);hist=Counter()
 for d in D:
  S=tuple(i for i,p in enumerate(ps) if d%p)
  inactive=prod(p for i,p in enumerate(ps) if i not in S)
  check(gcd(e*e,d)==inactive**2,'gcd support')
  hist[S]+=1
 for size in range(1,k+1):
  for S in combinations(range(k),size):check(hist[S]==2**(size-1),'stratum cardinality')
 slices=0
 for deleted in product((False,True),repeat=k):
  t=prod(p for p,delete in zip(ps,deleted) if delete)
  direct={d//(t*t) for d in D if d%(t*t)==0}
  check(direct==gaps(e//t),'full square division slice equality');slices+=1
 check(len(D)==(3**k-1)//2,'gap count')
 rows.append({'e':e,'count':len(D),'deletion_slices':slices})
# Independent projective orbit construction: remove full scalar orbits, no leading-coordinate normalization.
proj=[]
for k in range(1,6):
 for q in (2,3,5,7):
  remaining=set(product(range(q),repeat=k));remaining.remove((0,)*k);orbits=[]
  while remaining:
   v=min(remaining)
   orbit={tuple(a*x%q for x in v) for a in range(1,q)}
   check(len(orbit)==q-1,'free scalar action')
   remaining.difference_update(orbit);orbits.append(orbit)
  hist=Counter(tuple(i for i,x in enumerate(next(iter(o))) if x) for o in orbits)
  check(len(orbits)==sum(q**i for i in range(k)),'projective total')
  for S,count in hist.items():check(count==(q-1)**(len(S)-1),'projective support')
  proj.append({'k':k,'q':q,'orbits':len(orbits)})
 # Pointed-monoid affine charts: chart i has x_i=1, others 0 or 1.
 # Gluing identifies charts precisely when their nonzero supports agree.
 chart_points=[]
 for i in range(k):
  for bits in product((0,1),repeat=k-1):
   v=bits[:i]+(1,)+bits[i:];chart_points.append((i,tuple(j for j,b in enumerate(v) if b)))
 supports={S for i,S in chart_points}
 check(len(supports)==2**k-1,'literal monoid Hom points')
 check(sum(len(S)==1 for S in supports)==k,'rank-zero Hom points')
 for S in supports:check(sum(T==S for i,T in chart_points)==len(S),'correct overlap multiplicity')
D=gaps(65);check(D=={2016,3000,3696,4056},'hostile exact gaps')
check(not any(x+y in D for x,y in combinations(D,2)),'hostile no Schur')
result={'status':'PASS','gap_enumeration':'positive Pythagorean legs','integer_rows':rows,'projective_orbit_rows':proj,'pointed_monoid_chart_gluing':'PASS k=1..5','negative_flag':'--negative'}
print(json.dumps(result,indent=2));Path(__file__).with_name('independent_results.json').write_text(json.dumps(result,indent=2)+'\n')
