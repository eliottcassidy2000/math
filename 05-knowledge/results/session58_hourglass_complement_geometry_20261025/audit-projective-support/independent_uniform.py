from pathlib import Path
from itertools import product
from math import isqrt,prod,gcd
from collections import Counter
import subprocess,sys,json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
rows=[]
for ell in (3,5,7):
 for ps in ((5,),(5,13)):
  a=(ell-1)//2;e=prod(p**a for p in ps);k=len(ps);gaps=set()
  for u in range(1,isqrt(e*e//2)+1):
   v=isqrt(e*e-u*u)
   if v>u and u*u+v*v==e*e:gaps.add(2*u*v)
  actual=Counter(gcd(e*e,d) for d in gaps);expected=Counter()
  for mags in product(range(a+1),repeat=k):
   r=sum(m>0 for m in mags)
   if r:expected[prod(p**(2*(a-m)) for p,m in zip(ps,mags))]+=2**(r-1)
  need(actual==expected,'independent magnitude/gcd stratification')
  remaining=set(product(range(ell),repeat=k));remaining.remove((0,)*k);sizes=[]
  while remaining:
   v=min(remaining);line={tuple(c*x%ell for x in v) for c in range(1,ell)};remaining-=line
   signed={min(w,tuple(-x%ell for x in w)) for w in line};sizes.append(len(signed))
  need(all(s==a for s in sizes),'signed to projective fibers')
  need(len(gaps)==sum(sizes),'independent gap/projective multiplicity')
  rows.append({'ell':ell,'k':k,'center':e,'gaps':len(gaps),'projective_points':len(sizes),'fiber':a})
root=Path(__file__).resolve().parent
source=(root/'replay_uniform.py').read_text();old="need(all(len(v)==a for v in fibers.values()),'uniform projective fiber')"
need(source.count(old)==1,'unique negative site')
f=root/'uniform_negative_tmp.py';f.write_text(source.replace(old,"need(all(len(v)==a+1 for v in fibers.values()),'uniform projective fiber')"))
negative=[]
for flags in ([],['-O']):
 r=subprocess.run([sys.executable,*flags,str(f)],capture_output=True,text=True)
 need(r.returncode!=0 and 'RuntimeError: uniform projective fiber' in r.stderr,'negative mutation rejected')
 negative.append({'flags':flags,'exit':r.returncode})
f.unlink()
out={'status':'PASS','rows':rows,'negative_mutations':negative}
print(json.dumps(out,indent=2));(root/'independent_uniform_results.json').write_text(json.dumps(out,indent=2)+'\n')
