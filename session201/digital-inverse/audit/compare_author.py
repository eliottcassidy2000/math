"""Direct physical-digit controls of the frozen author API, including windows."""
from pathlib import Path
from functools import lru_cache
import importlib.util,json,hashlib
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent if(HERE.parent/'verify.py').is_file()else HERE.parent/'session201-digital-inverse'
sp=importlib.util.spec_from_file_location('digital_author',SOURCE/'verify.py');author=importlib.util.module_from_spec(sp);sp.loader.exec_module(author)
def need(x,m):
 if not x:raise RuntimeError(m)
@lru_cache(None)
def g(u):
 if u==1:return 1
 d,e,h,z=1,0,2,1
 for bit in bin(u-1)[3:]:
  if bit=='1':d,e,h,z=d,0,2*h,z
  else:d,e,h,z=d-e+h,z,2*h-z,z
 return d
def read(row,u):a,b,c,d=row;return a*g(u)+b*u+c*(u%2)+d
rows=set()
for a in(1,3,5,7,9,15,21):
 for b in(a,a+1,a+3,2*a,3*a):
  for c in(-a,-a//2,0):
   for d in(a-b,(a-b)//2,0):rows.add((a,b,c,d))
controls=first_controls=0
for row in sorted(rows):
 expected={}
 for v in range(2,202):
  n=read(row,v)
  if n<=201:expected.setdefault(n,[]).append(v)
 for n in range(1,202):
  got,stats=author.fixed_inverse(row,n)
  need(got==expected.get(n,[]),'complete author/direct comparison');controls+=1
 for n in range(1,202,17):
  for lo,hi in((2,2),(3,19),(5,5),(7,103),(40,11)):
   want=[v for v in expected.get(n,[])if lo<=v<=hi]
   got,stats=author.fixed_inverse(row,n,lo,hi);need(got==want,'author window');controls+=1
   first,stats=author.fixed_inverse(row,n,lo,hi,True)
   need((not first)==(not want)and len(first)<=1 and all(v in want for v in first),'first mode valid witness/absence');first_controls+=1
need(author.fixed_inverse((1,1,0,0),3)[0]==[2],'scalar singleton maps to physical index two')
# This control records, but does not authorize, the source helper's out-of-domain call.
try:
 edge=author.fixed_inverse((1,1,0,0),2,vmin=1)[0]
 edge_status='unguarded vmin=1 returns scalar index one'if edge==[1]else'lower domain clamped'
except RuntimeError:edge_status='vmin>=2 API precondition enforced'
print(json.dumps(dict(status='PASS',numerical_rows=len(rows),complete_and_window_controls=controls,
 first_mode_controls=first_controls,scalar_singleton_original_index=2,out_of_domain_boundary=edge_status,
 scope='The physical theorem requires original v>=2; internal h=1 is permitted.',
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
