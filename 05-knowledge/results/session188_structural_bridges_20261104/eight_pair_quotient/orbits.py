from itertools import permutations
import json,struct
from pathlib import Path
ROOT=Path(__file__).resolve().parent
remaining=set(struct.unpack('<30720I',(ROOT/'blind8.bin').read_bytes()))
records=[];owner={}
def rows(code):
 a=[0]*8;k=0
 for i in range(8):
  for j in range(i+1,8):
   if code>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
   k+=1
 return a
perms=list(permutations(range(8)));pairs=[(i,j) for i in range(8) for j in range(i+1,8)]
while remaining:
 code=min(remaining);a=rows(code);orbit=set()
 for p in perms:
  orbit.add(sum(((a[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(pairs)))
 if not orbit<=remaining:raise RuntimeError('orbit closure')
 for x in orbit:owner[x]=len(records)
 remaining-=orbit
 prime=not any(all((mask>>i&1) or (a[i]&mask) in (0,mask) for i in range(8)) for mask in range(1,255) if mask.bit_count()>1)
 records.append(dict(code=code,rows=a,orbit_size=len(orbit),automorphisms=40320//len(orbit),scores=sorted(x.bit_count() for x in a),prime=prime))
for row in records:row['converse_orbit']=owner[((1<<28)-1)^row['code']]
print(json.dumps(records,indent=2))
