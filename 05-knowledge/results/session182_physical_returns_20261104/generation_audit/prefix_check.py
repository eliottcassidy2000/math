import json
from pathlib import Path
P=Path(__file__).resolve().parent

def birth(r,c):return [x|(((c>>i)&1)<<len(r))for i,x in enumerate(r)]+[((1<<len(r))-1)^c]
def response(r):
 n=len(r);N=1<<n;d=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for v in range(n):d[1<<v][v]=1
 for m in range(1,N):
  for v in range(n):
   if m>>v&1:
    for u in range(n):
     if not(m>>u&1) and r[v]>>u&1:d[m|1<<u][u]+=d[m][v]
  H[m]=sum(d[m])
 return H[-1],sum(H[m]*H[N-1-m]for m in range(N))
menus={4:set(),5:set()}
for a in range(8):
 r=birth([6,4,0],a);menus[4].add(response(r))
 for b in range(16):menus[5].add(response(birth(r,b)))
expected={4:{(1,16),(3,24),(5,32)},5:{(1,32),(3,48),(5,64),(9,84),(11,100),(13,104),(15,108),(15,120)}}
if menus!=expected:raise RuntimeError('menu mismatch')
out={str(k):sorted(v)for k,v in menus.items()};(P/'prefix_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS independent subset-DP profiles for8/128 retained-TT3 birth cuts')
