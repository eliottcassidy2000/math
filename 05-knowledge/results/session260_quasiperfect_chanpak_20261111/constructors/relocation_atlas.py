#!/usr/bin/env python3
"""Two explicit five-core inflation boxes, selected-pair relocation, and heads."""
from itertools import product
from pathlib import Path
import json,random
from module_ladders import endpoints,data
HERE=Path(__file__).resolve().parent
CORES=[(20,1,2,7,14),(24,1,3,6,14)]
def modules(T,lens):
 starts=[];n=0
 for l in lens:starts.append(n);n+=l
 rows=[0]*n
 for v,l in enumerate(lens):
  for i in range(l):
   z=starts[v]+i
   rows[z]|=sum(1<<(starts[v]+j) for j in range(i+1,l))
   for w,k in enumerate(lens):
    if T[v]>>w&1:rows[z]|=((1<<k)-1)<<starts[w]
 return tuple(rows)
def marked_pairs(B):
 n=len(B);f=(1<<n)-1;ee=endpoints(B)
 rev=tuple(sum(1<<u for u in range(n) if B[u]>>v&1) for v in range(n));ss=endpoints(rev)
 h=sum(ee[f]);K=[[0]*n for _ in range(n)];U=ss[f][:]
 for A in range(1,f):
  es=[(i,x) for i,x in enumerate(ee[A]) if x];bs=[(i,x) for i,x in enumerate(ss[f^A]) if x];ha=sum(x for _,x in es)
  for b,y in bs:
   U[b]+=ha*y
   for a,x in es:K[a][b]+=x*y
 return [(h,a,b,U[a]+U[b]) for a in range(n) for b in range(a+1,n) if K[a][b]+K[b][a]==h-1]
def oddpart(q):
 while q%2==0:q//=2
 return q

def heads(targets):
 rng=random.Random(260);found={};attempts=0
 for n in range(3,9):
  for _ in range(20000):
   rows=[0]*n
   for a in range(n):
    for b in range(a+1,n):
     x,y=(a,b) if rng.randrange(2) else (b,a);rows[x]|=1<<y
   h=sum(endpoints(tuple(rows))[-1]);attempts+=1
   if h in targets and h not in found:found[h]=rows
   if len(found)==len(targets):return found,attempts
 return found,attempts

def main():
 found=[];graphs=0
 for c,T in enumerate(CORES):
  for lens in product(range(1,6),repeat=5):
   if sum(lens)>12:continue
   graphs+=1;B=modules(T,lens)
   for h,a,b,u in marked_pairs(B):found.append(dict(core=c,lengths=lens,H=h,selected=[a,b],U=u,odd_slope=oddpart(u)))
 best={}
 for r in found:
  c=r['odd_slope']
  if c not in best or (sum(r['lengths']),r['U'])<(sum(best[c]['lengths']),best[c]['U']):best[c]=r
 showcase={c:best[c] for c in [7,13,21]};targets=set()
 for c,r in showcase.items():
  q=r['U'];j=(q//c).bit_length()-1
  r['rows']=modules(CORES[r['core']],r['lengths']);r['ray_start_exponent']=j+1
  targets.update(1+c*2**s for s in range(1,j+1))
 small,attempts=heads(targets)
 if set(small)!=targets:raise RuntimeError(('missing heads',targets-set(small)))
 for c,r in showcase.items():
  B=r['rows'];O=sum(1<<v for v in r['selected']);h,w,q=data(B,O)
  if (h-w,q)!=(1,r['U']):raise RuntimeError('showcase direct selection')
 # The atlas is bounded, not a universal grammar or canonical class census.
 out={'status':'PASS','cores':CORES,'universe':{'each_module_length':[1,5],'total_order_at_most':12,'selected_set_size':2,'template_graphs':graphs},'seed_presentations':len(found),'distinct_odd_slopes':len(best),'missing_odd_slopes_below100':[c for c in range(1,100,2) if c not in best],'showcase':showcase,'below_start_witnesses':small,'head_search_candidates':attempts,'seeds':found}
 (HERE/'relocation_atlas_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k!='seeds'},indent=2))
if __name__=='__main__':main()
