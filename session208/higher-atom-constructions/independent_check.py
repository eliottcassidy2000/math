import json,itertools
from pathlib import Path
rows=json.load(open(Path(__file__).resolve().parent/'WITNESSES.json'));out=[]
for row in rows:
 n=row['n'];code=int(row['code']);edges=list(itertools.combinations(range(n),2));adj=[set() for _ in range(n)]
 for k,(i,j) in enumerate(edges):
  if code>>k&1:adj[i].add(j)
  else:adj[j].add(i)
 for s in range(n):
  seen={s};q=[s]
  for v in q:
   for w in adj[v]:
    if w not in seen:seen.add(w);q.append(w)
  assert len(seen)==n
 def paths(v,used):
  if len(used)==n:return 1
  return sum(paths(w,used|{w}) for w in adj[v]-used)
 count=sum(paths(v,{v}) for v in range(n));assert count==row['h'];a,b=row['exponents'];assert count==3**a*5**b
 if n==9:
  direct=sum(all(p[i+1] in adj[p[i]] for i in range(n-1)) for p in itertools.permutations(range(n)));assert direct==count
 out.append({'n':n,'h':count,'strong':True,'method':'full permutations and recursive path enumeration' if n==9 else 'recursive path enumeration without memoization'})
assert 1875**2*2025**3==15**14
assert 54675*46875==15**8
print(json.dumps(out,indent=2))
