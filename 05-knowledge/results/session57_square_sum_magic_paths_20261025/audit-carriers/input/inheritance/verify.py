from pathlib import Path
from math import isqrt
import hashlib,json
ROOT=Path(__file__).resolve().parent

def require(b,msg):
 if not b: raise RuntimeError(msg)
def square(n): return n>=0 and isqrt(n)**2==n
def gaps(e):
 c=e*e
 return {c-a*a for a in range(1,e) if square(2*c-a*a)}
def petals(e):
 # Independent edge-root enumeration: p²+r²=2q², q=2e.
 c=e*e; ans=set()
 for p in range(1,2*e):
  if not square(8*c-p*p):continue
  x=p*p-2*c
  if x>0:ans.add((x,2*c,4*c-x))
 return ans
def diamond(ds):
 return [(x,y) for x in sorted(ds) for y in sorted(ds) if x>y and x!=2*y and x-y in ds and x+y in ds]

pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())
for pin in pins:
 require(hashlib.sha256((ROOT/pin['snapshot']).read_bytes()).hexdigest()==pin['sha256'],'source hash '+pin['id'])
count=0
for e in range(1,501):
 ds=gaps(e); short={d for d in ds if 2*d<e*e}
 images={(2*e*e-4*d,2*e*e,2*e*e+4*d) for d in short}
 require(images==petals(e),'bijection at '+str(e));count+=len(images)
 require(not square(2*e*e),'opposite shell')
 for t in images:
  require(t[0]>0 and t[0]+t[2]==2*t[1] and all(square(t[i]+t[j]) for i in range(3) for j in range(i+1,3)),'triangle')
examples={}
for e in [5,41,65,1105]:
 ds=gaps(e);examples[str(e)]={'gaps':sorted(ds),'short_gaps':sorted(d for d in ds if 2*d<e*e),'diamonds':diamond(ds)}
require(examples['5']['gaps']==[24] and not examples['5']['short_gaps'],'long-gap hostile')
require(examples['41']['short_gaps']==[720],'minimum petal')
require(len(examples['1105']['gaps'])==13 and not examples['1105']['diamonds'],'capacity hostile')
print(json.dumps({'status':'FINITE-EXACT','center_roots_checked':[1,500],'positive_petals_checked':count,'source_hashes_checked':len(pins),'examples':examples},sort_keys=True,indent=2))
