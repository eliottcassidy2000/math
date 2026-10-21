from itertools import combinations,permutations,product
import json
checks=0
def need(ok):
 global checks
 checks+=1
 if not ok:raise RuntimeError('corollary failure')
for n in range(4):
 edges=list(combinations(range(n),2));index={p:i for i,p in enumerate(edges)};full=(1<<len(edges))-1
 for a in range(1,4):
  for ts in product(range(full+1),repeat=a):
   decorations=(ts[0],)+tuple(t^ts[0] for t in ts[1:])
   inverse=(decorations[0],)+tuple(g^decorations[0] for g in decorations[1:])
   need(inverse==ts)
   for j in range(1,a):
    flipped=list(decorations);flipped[j]^=full
    recovered=(flipped[0],)+tuple(g^flipped[0] for g in flipped[1:])
    expected=list(ts);expected[j]^=full
    need(recovered==tuple(expected))
   reversed_ts=tuple(t^full for t in ts)
   need((reversed_ts[0],)+tuple(t^reversed_ts[0] for t in reversed_ts[1:])==(decorations[0]^full,)+decorations[1:])
   for p in permutations(range(n)):
    def action(x,tournament):
     return sum(((((x>>i)&1)^(tournament and p[u]>p[v]))<<index[tuple(sorted((p[u],p[v])))]) for i,(u,v) in enumerate(edges))
    pts=tuple(action(t,True) for t in ts)
    need((pts[0],)+tuple(t^pts[0] for t in pts[1:])==(action(decorations[0],True),)+tuple(action(g,False) for g in decorations[1:]))
print(json.dumps({'checks':checks,'orders':list(range(4)),'tuple_sizes':[1,2,3],'result':'PASS'},indent=2))
