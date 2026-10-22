from itertools import permutations, combinations
from collections import Counter
import json
from pathlib import Path

def adj(n,m):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(m>>k)&1;a[j][i]=1-a[i][j]
 return a

def code(a,p):return sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(combinations(range(len(p)),2)))
def canonical(a,vs):return min(code(a,p) for p in permutations(vs))
rows=[]
for m in [40,41]:
 a=adj(5,m);paths=[p for p in permutations(range(5)) if all(a[p[i]][p[i+1]] for i in range(4))]
 auto=[p for p in permutations(range(5)) if all(a[i][j]==a[p[i]][p[j]] for i in range(5) for j in range(5))]
 anti=[p for p in permutations(range(5)) if all(a[i][j]==a[p[j]][p[i]] for i in range(5) for j in range(5))]
 deck=Counter(canonical(a,vs) for vs in combinations(range(5),4))
 reach=[sum(a[i][j]<<j for j in range(5)) for i in range(5)]
 for _ in range(5):
  for i in range(5):
   for j in range(5):
    if (reach[i]>>j)&1:reach[i]|=reach[j]
 rows.append({'mask':m,'out_masks':[sum(a[i][j]<<j for j in range(5)) for i in range(5)],'scores':sorted(map(sum,a)),'paths':paths,'H':len(paths),'automorphisms':auto,'anti_automorphisms':anti,'deck4':dict(deck),'strong':all(r==31 for r in reach)})
if not(rows[0]['deck4']==rows[1]['deck4'] and [r['H'] for r in rows]==[15,13] and all(r['strong'] and len(r['automorphisms'])==len(r['anti_automorphisms']) for r in rows) and [len(r['automorphisms']) for r in rows]==[3,1]):raise RuntimeError('seed failed')
print(json.dumps(rows,indent=2));Path(__file__).with_suffix('.json').write_text(json.dumps(rows,indent=2)+'\n')
