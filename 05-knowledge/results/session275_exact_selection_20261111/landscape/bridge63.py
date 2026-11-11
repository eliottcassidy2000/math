from experiment import hcount
from itertools import combinations
import json
seed=[16,5,17,7,74,31,175,63];n=8;es=list(combinations(range(n),2));res={}
for k in (1,2):
 for idx in combinations(range(len(es)),k):
  r=seed[:]
  for z in idx:
   a,b=es[z];r[a]^=1<<b;r[b]^=1<<a
  h=hcount(r)
  if h in (61,65):
   res.setdefault(h,[]).append({'rows':r,'flips':[es[z] for z in idx],'neutral_padding':any(x.bit_count() in (0,7) for x in r)})
print(json.dumps({'seed':seed,'H':hcount(seed),'neighbors':res},indent=2))
