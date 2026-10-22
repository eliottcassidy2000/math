"""Independently generate complete small classes; no producer imports."""
from itertools import combinations,permutations,product
from pathlib import Path
import subprocess,json
base=Path(__file__).resolve().parent
result={}
with (base/'small-bank.txt').open('w') as output:
 for n in (2,3,4,5):
  pairs=list(combinations(range(n),2));classes=set();labelled=0
  for bits in product((0,1),repeat=len(pairs)):
   a=[[0]*n for _ in range(n)]
   for (i,j),b in zip(pairs,bits):a[i][j]=b;a[j][i]=1-b
   scores=[sum(r) for r in a]
   if n%2 and any(x!=(n-1)//2 for x in scores):continue
   if not n%2 and any(x not in (n//2-1,n//2) for x in scores):continue
   labelled+=1
   canonical=min(''.join(str(a[p[i]][p[j]]) for i,j in pairs) for p in permutations(range(n)))
   classes.add(canonical)
  result[n]={'classes':len(classes),'labelled':labelled}
  for bits in sorted(classes):output.write(f'{n} {bits}\n')
with (base/'small-bank.txt').open() as inp:
 text=subprocess.check_output([str(base/'audit')],stdin=inp,text=True)
(base/'small-results.txt').write_text(text)
for line in text.splitlines():
 n,bits,aut,anti,inv,orbits,good=line.split()
 result[int(n)]['self_converse']=result[int(n)].get('self_converse',0)+bool(int(anti))
assert all(v['classes']==v['self_converse']==1 for v in result.values())
assert result[2]['labelled']==result[3]['labelled']==2
assert result[4]['labelled']==result[5]['labelled']==24
(base/'small-receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
