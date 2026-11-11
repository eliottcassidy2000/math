from carrier import *
from collections import defaultdict
I=matrix(profile((0,)));S=matrix(profile((2,0)));K=matrix(profile((0,1)))
for n in range(1,19):
 words=[(I,'')]
 # incremental prefix dynamic retain all states
 D={}
 for bits in product((0,1),repeat=n):
  M=I
  for b in bits:M=mul(M,S if b else K)
  p=from_matrix(M);key=tuple(p[k] for k in ('H','h','Q','V','e'))
  u=(p['U'],p['d'])
  if key in D and D[key][0]!=u:
   print(json.dumps(dict(n=n,key=key,old=D[key],new=(u,bits),profile=p)),flush=True);quit()
  D.setdefault(key,(u,bits))
 print('finished',n,flush=True)
