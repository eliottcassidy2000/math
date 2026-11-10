import importlib.util,json
from pathlib import Path
s=importlib.util.spec_from_file_location('physical',Path(__file__).with_name('physical.py'));p=importlib.util.module_from_spec(s);s.loader.exec_module(p)
cases=0;states=set()
for n in range(1,5):
 for code in range(1<<(n*(n-1)//2)):
  a=p.labelled(n,code)
  if not p.strong(a):continue
  for root in range(n):
   ha,ea,ua,da=p.state(a,root)
   for b in [p.tt(1),p.tt(2),p.tt(3),p.endpoint(1),p.join(p.tt(1),p.endpoint(1))]:
    hb,qb=p.pair(b);ab=p.join(a,b)
    p.need(p.state(ab,root)[2]==ua*qb,'ordinal response factorization')
    t=p.graft(ab,root,[0]);z=p.state(t,len(t)-1)
    expected=(ha*hb+ua*qb,ha*hb,ha*hb+ua*qb,ha*hb)
    p.need(z==expected and p.strong(t),'outdegree one structural identity')
    p.need(z[0]>=3*z[1] and ((z[0]==3*z[1])==(len(b)==1 and ua==ha)),'sharp threefold bound')
    cases+=1;states.add(z[:2])
classification=0
for n in range(3,6):
 for code in range(1<<(n*(n-1)//2)):
  a=p.labelled(n,code)
  if not p.strong(a):continue
  for r in range(n):
   h,e,u,d=p.state(a,r)
   out=[v for v in range(n) if a[r]>>v&1]
   predicted=False
   if len(out)==1:
    keep=[v for v in range(n) if v!=r]
    b=[sum(1<<j for j,v in enumerate(keep) if a[x]>>v&1) for x in keep]
    predicted=p.tables(b)[1][-1][keep.index(out[0])]==0
   p.need((u==h)==predicted,'saturation classification')
   if u==h:p.need(e==d and h>=3*e and h!=13,'saturation consequences')
   classification+=1
print(json.dumps(dict(status='PASS',classified_strong_roots=classification,physical_structural_cases=cases,distinct_he_pairs=sorted(states),scope='Saturation classification on every strong labelled tournament of orders3 through5, and constructive structural tests.'),sort_keys=True,indent=2))
