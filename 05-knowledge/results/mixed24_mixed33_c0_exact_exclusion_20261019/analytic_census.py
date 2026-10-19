from pathlib import Path
from functools import lru_cache
from collections import Counter
import json,gzip,time
import cone_engine as F
ROOT=Path(__file__).resolve().parent
SCOPES=[(2,4,0),(3,3,0)]
@lru_cache(None)
def natural(n):
 if n==0:return [()]
 return [p+(m,) for p in natural(n-1) for m in range(1<<(n-1)) if all(not(m>>i&1) or not(p[i]&~m) for i in range(n-1))]
@lru_cache(maxsize=1)
def bounds(p):
 n=len(p);up=tuple(sum(1<<j for j in range(n) if p[j]>>i&1) for i in range(n));deg=F.degrees(p);seed=3|sum(1<<i for i,d in enumerate(deg) if d==7)
 alpha=[0 if ((1<<v)|p[v])&seed else 7-deg[v] for v in range(n)];beta=[0 if ((1<<v)|up[v])&seed else 7-deg[v] for v in range(n)]
 return {(u,v):(p[v]&~p[u]).bit_count()+(up[u]&~up[v]).bit_count()+alpha[u]+beta[v] for u in range(n) for v in range(n) if u!=v and not(p[v]>>u&1 or p[u]>>v&1)}
F.defect_bound=lambda p,u,v:bounds(p)[u,v]
def prepare(scope):
 l,u,c=scope;path=ROOT/('cover_%s%s%s.json.gz'%scope)
 if path.exists():return json.loads(gzip.decompress(path.read_bytes()))
 t=time.time();records=[];hist=Counter();qs=natural(l+u+c)
 for q in qs:
  p=[0,sum(1<<i for i in range(2,2+l))]+[v<<2 for v in q]
  for a in range(2+l+c,len(p)):p[a]|=1
  p=tuple(p);F.require(F.closure(p)==p,'closure')
  if max(F.degrees(p))>7:hist['invalid']+=1;continue
  h=bounds(p);w=next(((a,b) for (a,b),v in h.items() if a<b and v<=1 and h[b,a]<=1),None)
  if w:rec={'pred':p,'balanced_pair':w};hist['balanced']+=1
  else:
   cyc=F.forced_cycle(p,h)
   if cyc:rec={'pred':p,'cycle_arcs':cyc};hist['cycle']+=1
   else:rec={'pred':p};hist['residual']+=1
  records.append(rec)
 records.sort(key=lambda r:r['pred']);out={'scope':scope,'natural':len(qs),'valid':len(records),'hist':dict(hist),'seconds':time.time()-t,'cores':records};path.write_bytes(gzip.compress(json.dumps(out,separators=(',',':')).encode(),mtime=0));return out
if __name__=='__main__':
 summaries=[]
 for scope in SCOPES:
  d=prepare(scope);s={k:v for k,v in d.items() if k!='cores'};summaries.append(s);print(json.dumps(s),flush=True)
 (ROOT/'analytic_summary.json').write_text(json.dumps(summaries,indent=2))
