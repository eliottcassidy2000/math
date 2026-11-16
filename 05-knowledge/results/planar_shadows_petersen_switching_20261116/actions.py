from itertools import combinations,permutations
from collections import Counter
V=tuple(combinations(range(5),2))
P={tuple(sorted((u,v))) for u,v in combinations(V,2) if not set(u)&set(v)}
MP=[]
for m in combinations(sorted(P),5):
 if len({v for e in m for v in e})==10:
  MP.append(frozenset(m))
assert len(MP)==6 and all(len(a&b)==1 for a,b in combinations(MP,2))
E6=tuple(combinations(range(6),2))
A=[p for p in permutations(range(5)) if sum(p[i]>p[j] for i,j in combinations(range(5),2))%2==0]
def av(p,v):return tuple(sorted(p[i] for i in v))
def ae(p,e):return tuple(sorted(av(p,v) for v in e))
def ai(p,d):return tuple(sorted(MP.index(frozenset(ae(p,e) for e in MP[k])) for k in d))
e0=min(P)
best=[d for d in E6 if all(ai(p,d)==d for p in A if ae(p,e0)==e0)]
assert len(best)==3
surfaces=[]
for d in best:
 beta={ae(p,e0):ai(p,d) for p in A}
 assert set(beta)==P and set(beta.values())==set(E6)
 b={v:frozenset(beta[e] for e in P if v in e) for v in V}
 labels=[frozenset(i for de in b[v] for i in de) for v in V]
 if {len(u)for u in labels}=={3}:
  faces=set(labels);surfaces.append(faces)
  assert len(faces)==10
  assert Counter(tuple(sorted(d)) for f in faces for d in combinations(f,2))=={q:2 for q in E6}
  for i in range(6):
   L={tuple(sorted(f-{i})) for f in faces if i in f}
   assert len(L)==5
   c=Counter(v for z in L for v in z);assert c=={j:2 for j in range(6)if j!=i}
  assert all((len(f&g)==2)==(tuple(sorted((u,v))) in P) for (u,f),(v,g) in combinations(zip(V,labels),2))
  print('triangles',sorted(tuple(sorted(f))for f in faces))
 else:
  assert {len(u)for u in labels}=={6}
  print('syntheme', sorted(tuple(sorted(x))for x in b.values()))
assert len(surfaces)==2 and surfaces[0].isdisjoint(surfaces[1])
assert len(surfaces[0]|surfaces[1])==20
print('Exactly three A5 identifications, complementary triangle surfaces.')
