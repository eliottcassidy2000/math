import csv, hashlib, itertools, json, math
from collections import Counter
from fractions import Fraction
from pathlib import Path
root=Path(__file__).resolve().parent.parent
src=root/'collision-sieve'
cert=json.loads((src/'certificate.json').read_text())
def check(x):
 if not x: raise RuntimeError('Independent audit mismatch')
P=tuple(row['p'] for row in cert['controls']['local'])
# Normalize the actual 36 grid-pair differences, retain integral scale.
grid=list(itertools.product((-1,0,1),repeat=2)); directions=Counter(); twos=0
for a,b in itertools.combinations(grid,2):
 x,y=a[0]-b[0],a[1]-b[1];g=math.gcd(x,y);x//=g;y//=g
 if x<0 or (x==0 and y<0):x,y=-x,-y
 directions[x,y]+=1;twos+=g==2
check(twos==8)
check(sorted(directions.values())==[2,2,2,2,5,5,9,9])
# Use linear forms, never producer polynomial/masks.
def passes(u,v,p):return any((a*u+b*v)%p==0 for a,b in directions)
aff=Fraction(1);prim=Fraction(1)
for p in P:
 allowed=sum(passes(u,v,p) for u in range(p) for v in range(p));L=4 if p==3 else 6 if p==5 else 8
 check(allowed==1+L*(p-1));aff*=Fraction(allowed,p*p);prim*=Fraction(allowed-1,p*p-1)
check(str(aff)==cert['controls']['affine_CRT_density']);check(str(prim)==cert['controls']['primitive_projective_CRT_density'])
# Independent Pythagorean-bank construction and all original integer pairs.
B=100000;banks={}
for n in range(1,math.isqrt(B)):
 for m in range(n+1,math.isqrt(B-n*n)+1):
  if (m+n)%2==0 or math.gcd(m,n)>1:continue
  h=m*m+n*n
  for k in range(1,B//h+1):
   e=k*h; a=k*abs(m*m-n*n-2*m*n);d=e*e-a*a
   banks.setdefault(e,set()).add(d)
allpairs=[];survivors=[];difference=0;distinct=0
for e in sorted(banks):
 A=banks[e]
 if len(A)<4:continue
 for v,u in itertools.combinations(sorted(A),2):
  allpairs.append((e,u,v));difference+=u-v in A;distinct+=u!=2*v
  if all(passes(u,v,p) for p in P):survivors.append((e,u,v))
check(len(allpairs)==347561);check(len(survivors)==474);check(difference==0)
remaining=allpairs
for stage in cert['search']['sequential']:
 p=stage['prime'];before=len(remaining);remaining=[z for z in remaining if passes(z[1],z[2],p)]
 check(stage==dict(prime=p,before=before,rejected=before-len(remaining),remaining=len(remaining)))
 check(sum(not passes(u,v,p) for e,u,v in allpairs)==cert['search']['standalone_rejections'][str(p)])
survivors.sort()
hash_survivors=hashlib.sha256(''.join(f'{e},{u},{v}\n' for e,u,v in survivors).encode()).hexdigest()
check(hash_survivors==cert['search']['survivor_sha256'])
rows=list(csv.DictReader((src/'survivors.csv').open()))
check([(int(r['e']),int(r['u']),int(r['v'])) for r in rows]==survivors)
for r in rows:
 e,u,v=int(r['e']),int(r['u']),int(r['v'])
 check(int(r['gcd_uv'])==math.gcd(u,v))
 for key,value in [('root_c_minus_u',e*e-u),('root_c_plus_u',e*e+u),('root_c_minus_v',e*e-v),('root_c_plus_v',e*e+v)]:check(int(r[key])**2==value)
report=dict(pairs=len(allpairs),survivors=len(survivors),difference_gap_passes=difference,distinct_magnitudes=distinct,directions={str(k):v for k,v in sorted(directions.items())},extra_factors_two=twos,affine_density=str(aff),primitive_density=str(prim),survivor_sha256=hash_survivors,status='PASS')
print(json.dumps(report,indent=2))
