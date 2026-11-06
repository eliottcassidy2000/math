import itertools,json
from functools import lru_cache
n=14;code=292187853813044375413825338
adj=[set() for _ in range(n)]
for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
 if (code>>k)&1:adj[i].add(j)
 else:adj[j].add(i)
for s in range(n):
 seen={s};q=[s]
 for v in q:
  for w in adj[v]:
   if w not in seen:seen.add(w);q.append(w)
 assert len(seen)==n
@lru_cache(None)
def suffix(v,rest):
 if rest==0:return 1
 return sum(suffix(w,rest^(1<<w)) for w in adj[v] if rest>>w&1)
ends=[suffix(v,((1<<n)-1)^(1<<v)) for v in range(n)]
assert sum(ends)==11390625==15**6
print(json.dumps({'n':n,'code':str(code),'strong':True,'count':sum(ends),'start_counts':ends,'rate':'7/3','status':'PASS'},indent=2))
