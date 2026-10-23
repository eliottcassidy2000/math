import json,itertools
checks=0
for n in range(3,65):
 C={tuple(sorted((i,(i+1)%n))) for i in range(n)}
 for u,v in itertools.combinations(range(n),2):
  e=(u,v);G=C|{e}
  if e in C:
   start=next(i for i in range(n) if tuple(sorted((i,(i+1)%n)))==e)
   omit=[(start+1+j)%n for j in range(n)]
   include=[(start+2+j)%n for j in range(n)]
  else:
   omit=list(range(n));include=list(range(u+1,v+1));w=u
   while True:
    include.append(w)
    if w==(v+1)%n:break
    w=(w-1)%n
  for P,required in [(include,True),(omit,False)]:
   if sorted(P)!=list(range(n)):raise ValueError('not spanning')
   E={tuple(sorted(a)) for a in zip(P,P[1:])}
   if not E<=G or (e in E)!=required:raise ValueError('edge certificate')
   checks+=1
G={(0,1),(1,2),(2,3),(0,3),(0,2)}
fixed=[p for p in itertools.permutations(range(4)) if p[0]==0 and p[-1]==1 and all(tuple(sorted(e)) in G for e in zip(p,p[1:]))]
if not fixed or any((0,2) in {tuple(sorted(e)) for e in zip(p,p[1:])} for p in fixed):raise ValueError('fixed endpoint firewall')
print(json.dumps(dict(status='PASS',inclusion_omission_checks=checks,universe='canonical cycles n3..64, every cycle edge and every possible added chord',fixed_endpoint_hostile=dict(graph=sorted(G),endpoints=[0,1],paths=fixed,unused_edge=[0,2])),indent=2))
