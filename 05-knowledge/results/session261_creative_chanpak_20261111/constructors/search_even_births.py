from module_tools import *
import json
seen5={canon(tour(5,b)) for b in range(1024)};seen=set()
for A in seen5:
 for b in range(32):
  T=tuple(A[i]|((1<<5) if b>>i&1 else 0) for i in range(5))+(sum(1<<i for i in range(5) if not b>>i&1),)
  seen.add(canon(T))
seeds=[];hits=[];qset=set()
for T in sorted(seen):
 for O in range(1,64):
  h,w,q=m.data(T,O)
  if h-w!=1:continue
  seeds.append((T,O,q));qset.add(q)
  if q%2:continue
  for v in range(6):
   if O>>v&1:continue
   Z,P=m.inflate(T,O,v,2);hh,ww,qq=m.data(Z,P)
   if hh-ww==1 and qq==2*q+2:
    z={'T':T,'O':O,'q':q,'v':v,'samples':[]}
    for k in range(1,5):
     Z,P=m.inflate(T,O,v,k);h,w,u=m.data(Z,P);z['samples'].append((h,h-w,u))
    hits.append(z);print(z,flush=True)
print('END',len(seeds),len(hits),sorted(qset))
Path(__file__).with_name('even_births.json').write_text(json.dumps({'seeds':seeds,'hits':hits},indent=2))
