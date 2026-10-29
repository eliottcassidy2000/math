from star_probe import build,observer
from pathlib import Path
import json

def require(c,m):
 if not c:raise ArithmeticError(m)
def witness(w):
 j=max(i+2 for i in range(len(w)-1) if w[i:i+2] in ['be','ee'])
 A,cycle,h=build(w[:j]);n=len(A);b,c=n-2,n-1;r=h[-2]['parent'];Q={r,b,c};R=set(range(n))-Q
 cross=[(u,v) for u in R for v in Q if A[u][v]]
 require(len(cross)==1 and cross[0][1]==r,('bridge',w));q=cross[0][0]
 choices=[]
 for a in sorted(R):
  inc=[s for s in R if A[s][a]]
  if len(inc)==1 and q in [a,inc[0]]:choices.append((a,inc[0],r))
 require(choices,('birth',w));return choices[0],j,c
counts=[];W=['']
for depth in range(8):
 count=0
 for w in W:
  if 'be' not in w and 'ee' not in w:continue
  (a,s,r),j,c=witness(w);A,_,h=build(w);ports,M=observer(A);idx={p:i for i,p in enumerate(ports)}
  require(all(row[idx[a,r]]-row[idx[r,a]]-row[idx[s,r]]+row[idx[r,s]]==0 for row in M),('kernel',w,a,s,r))
  require(all(x['parent']!=c for x in h[j:]),('parent',w))
  count+=1
 counts.append([depth,count]);print(counts[-1],flush=True)
 W=[w+x for w in W for x in ('ab' if w.endswith(('a','c')) else 'cde')]
Path(__file__).with_name('last_birth_receipt.json').write_text(json.dumps({'status':'PASS','depth_counts':counts,'total':sum(x[1] for x in counts)},indent=2))
