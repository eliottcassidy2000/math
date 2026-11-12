import json,itertools
from scout import corr,canon
n=13;bins={}
for k in range(3,7):
 for C in itertools.combinations(range(n),k):
  m=sum(1<<i for i in C)
  if canon(m,n)!=m:continue
  bins.setdefault((k,corr(m,n)),[]).append(m)
best=None
for spec in map(json.loads,open(__file__.replace('search_phase.py','search13.jsonl'))):
 w=spec['weights']
 for (k,cc),ms in bins.items():
  for a,b in itertools.combinations(ms,2):
   for sh in range(n):
    bb=sum(1<<((j+sh)%n)for j in range(n)if b>>j&1)
    common=[v for v in range(n)if(a>>v&1)==(bb>>v&1)]
    gaps=[-(1-2*(a>>v&1))*sum(w[(j-v)%n]*((a>>j&1)-(bb>>j&1))for j in range(n))for v in common]
    for mask in range(1,1<<len(common)):
     gap=sum(gaps[i]for i in range(len(common))if mask>>i&1)
     if not gap:continue
     if best is None or abs(gap)<abs(best['gap']):
      toggles=[common[i]for i in range(len(common))if mask>>i&1]
      best={'gap':gap,'a':a,'b':bb,'toggles':toggles,'spec':spec}
      print(json.dumps(best),flush=True)
