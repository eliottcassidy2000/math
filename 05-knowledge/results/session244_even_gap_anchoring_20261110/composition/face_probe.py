import importlib.util,json
from pathlib import Path
from itertools import product
s=importlib.util.spec_from_file_location('verify',Path(__file__).with_name('verify.py'));v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
def face(w):
 m=v.matrix(w);r=[(a+c,b+d) for a,b,c,d in m]
 return (r[0],r[1]) if r[0]==r[2] and r[1]==r[3] else None
maps={}
for n in range(1,9):
 for x in product('OZJ',repeat=n):
  w=''.join(x);m=face(w)
  if m is not None:maps.setdefault(m,w)
primitive={m:w for m,w in maps.items() if not any(face(w[:i]) is not None for i in range(1,len(w)))}
miss=[];valid=0;covered=0
for H in range(3,502,2):
 for E in range(1,(H-1)//2+1,2):
  valid+=1;yes=False
  for ((a,b),(c,d)),w in maps.items():
   det=a*d-b*c
   if not det:continue
   x,y=d*H-b*E,-c*H+a*E
   if x%det or y%det:continue
   h,e=x//det,y//det
   if h%2==1 and e%2==1 and e>=1 and h>=2*e+1 and h<H:
    yes=True;break
  if yes:covered+=1
  elif H<=51:miss.append((H,E))
report=dict(scope='Exploratory finite face-linear maps, not physical completeness.',max_word_length=8,distinct_maps=len(maps),primitive_return_maps=[dict(word=w,matrix=m,determinant=m[0][0]*m[1][1]-m[0][1]*m[1][0]) for m,w in sorted(primitive.items(),key=lambda x:(len(x[1]),x[1]))],abstract_cone_points=valid,inverse_covered=covered,misses_through_H51=miss)
print(json.dumps(report,sort_keys=True,indent=2))
