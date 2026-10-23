"""Independent bounded checks of incoming shared-corner lifts and source-set exclusion."""
from itertools import product
from collections import Counter
from pathlib import Path
import argparse,json,hashlib
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);a=p.parse_args();checks=0

def need(c,m):
 global checks
 checks+=1
 if not c:raise RuntimeError(m)

def parity(x):return x.bit_count()%2
counts=[]
for r,c in product(range(1,4),repeat=2):
 matrices=[tuple((bits>>(i*c))&((1<<c)-1) for i in range(r)) for bits in range(1<<(r*c))]
 for u in range(1<<c):
  for v in range(1<<r):
   outputs=Counter()
   for D in matrices:
    aa=sum(parity(row&u)<<i for i,row in enumerate(D));bb=0
    for i,row in enumerate(D):
     if v>>i&1:bb^=row
    outputs[aa,bb]+=1
   for aa in range(1<<r):
    for bb in range(1<<c):
     allowed=parity(v&aa)==parity(u&bb) and (u!=0 or aa==0) and (v!=0 or bb==0)
     rank=r+c-1 if u and v else r if u else c if v else 0
     need(outputs[aa,bb]==(1<<(r*c-rank) if allowed else 0),'complete shared-corner map and exact fibre size')
 counts.append({'rows':r,'columns':c,'matrices':len(matrices)})
# A bounded correction hostile: all coefficients odd, each capacity one.
C=((1,1),(1,1));u=v=(1,1)
need(all(sum(row)%2==0 for row in C),'coarse level valid')
need(all(sum(row)%4!=0 for row in C),'unchanged bounded representative fails next level')
bounded=[M for flat in product(range(2),repeat=4) for M in [(flat[:2],flat[2:])] if all((M[i][j]-C[i][j])%2==0 for i,j in product(range(2),repeat=2))]
need(bounded==[C],'capacity-one coefficient box has no alternate lift')
U=((3,1),(1,3));need(all(sum(row)%4==0 for row in U) and all(sum(col)%4==0 for col in zip(*U)),'unrestricted odd-difference lift exists')
raw=(a.repo/'05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json').read_bytes();D=json.loads(raw)['order8'];Y=D['radical_lifts_mod4'];X=D['higher_layers'][-1]['lifts'][0]
reports=[]
for sources,expected in [([1238,2242,2246,2315],[0]),([1425,1921,1477,2837],[0,15]),([1981,2458,2659],[0,5])]:
 rows=[[1]*len(sources)]+[[Y[j][s]%2 for s in sources] for j in range(4)]+[[X[s]%2 for s in sources]]
 kernel=[bits for bits in range(1<<len(sources)) if all(sum(row[i] for i in range(len(sources)) if bits>>i&1)%2==0 for row in rows)]
 need(kernel==expected,'literal augmented source parity kernel')
 reports.append({'sources':sources,'augmented_rows':rows,'binary_kernel':kernel})
print(json.dumps({'status':'PASS','checks':checks,'shared_corner_tests':counts,'source_parity_reports':reports,'source_sha256':hashlib.sha256(raw).hexdigest(),'scope':'Bounded independent tests of incoming arithmetic mechanisms, not a replay of all routing, assembly or wiring proofs.'},indent=2))
