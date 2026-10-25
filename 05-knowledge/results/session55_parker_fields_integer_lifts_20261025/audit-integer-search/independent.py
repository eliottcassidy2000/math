import json,math,importlib.util
from pathlib import Path
base=Path(__file__).resolve().parent.parent
spec=importlib.util.spec_from_file_location('producer',base/'integer-search/verify_and_search.py'); p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
def req(v):
 if not v:raise RuntimeError('audit failure')
def rotate(a):return [a[6],a[3],a[0],a[7],a[4],a[1],a[8],a[5],a[2]]
def reflect(a):return [a[2],a[1],a[0],a[5],a[4],a[3],a[8],a[7],a[6]]
def orbit(a):
 o=[]
 for _ in range(4):o.extend([a,reflect(a)]);a=rotate(a)
 return o
cnt=0
for x in range(-10,11):
 for y in range(-10,11):
  a=p.grid(30,x,y)
  if len(set(a))<9:continue
  canonical=[z for z in orbit(a) if z[0]-30>z[2]-30>0]
  req(len(canonical)==1);req(canonical[0]==p.grid(30,*sorted([abs(x),abs(y)],reverse=True)));cnt+=1
banks,_,_=p.sieve_offsets(2000); checked=0;found=[]
for e in range(1,2001):
 c=e*e; squares={r*r for r in range(math.isqrt(2*c)+1)}
 corners=[a for a in squares if 2*c-a in squares and a!=c]
 req({abs(a-c) for a in corners}==(banks[e] or set()))
 for a in corners:
  for b in corners:
   # Construct other entries directly from opposite pairs and row/column equations.
   z=[a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a]
   checked+=1
   if len(set(z))==9 and all(v in squares for v in z):found.append(z)
req(not found)
out={'D4_distinct_parameter_tests':cnt,'direct_root_corner_center_bound_inclusive':2000,'direct_root_corner_pairs':checked,'solutions':found}
print(json.dumps(out,indent=2))
