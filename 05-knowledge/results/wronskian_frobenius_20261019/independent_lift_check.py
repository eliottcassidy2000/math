from itertools import product
import json

def conv(a,b,q):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%q
 return c
rows=[]
for p,m,D in [(2,3,4),(3,2,3),(5,2,2),(2,4,3)]:
 q=p**m; seen=set();count=0
 for a in product(range(q),repeat=D+1):
  b=[1]
  for _ in range(p):b=conv(b,a,q)
  c=[0]*(p*D+1)
  for j,x in enumerate(a):c[p*j]=x
  if b==c:seen.add(tuple(x%p for x in a))
  count+=1
 expected=set()
 for f in product(range(p),repeat=D+1):
  support=[i for i,x in enumerate(f) if x]
  if len({i%(p**(m-1)) for i in support})<=1:expected.add(f)
 if seen!=expected:raise RuntimeError((p,m,D,seen^expected))
 rows.append(dict(p=p,m=m,D=D,all_lifts=count,admissible_mod_p_rules=len(seen)))
print(json.dumps(dict(status='PASS',universes=rows),indent=2))
