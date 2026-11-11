from carrier import *
from collections import defaultdict
D=defaultdict(dict)
for n in range(2,7):
 for T in tournaments(n):
  p=profile(T); hp.cache_clear()
  if p['h']!=1:continue
  key=(n,)+tuple(p[k] for k in ('H','h','Q','V','e'))
  D[key].setdefault((p['U'],p['d']),T)
 print('done',n,flush=True)
for key,x in D.items():
 if len(x)>1:print(key,x)
print('max',max(map(len,D.values())))
