from carrier import *
from math import gcd
D={};found=[]
for T in tournaments(6):
 p=profile(T);hp.cache_clear();key=tuple(p[k] for k in ('H','h','Q','V','e'))
 if key in D:
  B,pb=D[key]
  if (p['U'],p['d'])!=(pb['U'],pb['d']):
   gg=0
   for k in p:gg=gcd(gg,p[k])
   z=dict(A=B,B=T,pA=pb,pB=p,g=gg)
   if gg==1:print(json.dumps(z),flush=True);break
 else:D[key]=(T,p)
else: print('no primitive pair')
