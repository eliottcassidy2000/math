from itertools import combinations
import json,time
from pathlib import Path
BASE=Path(__file__).resolve().parent

def chain(m):
 n=3*m+1;edges=set()
 for k in range(m):
  t,a,b=1+3*k,2+3*k,3+3*k
  edges.update(tuple(sorted(e))for e in [(0,a),(0,b),(a,b),(a,t),(b,t)])
  if k:edges.add((t-3,t))
 return n,edges

def induced(host,vertices,n,edges):
 index={v:i for i,v in enumerate(vertices)};h=len(vertices);ha=[0]*h;pa=[0]*n
 for u,v in host['edges']:
  if u in index and v in index:
   u,v=index[u],index[v];ha[u]|=1<<v;ha[v]|=1<<u
 for u,v in edges:pa[u]|=1<<v;pa[v]|=1<<u
 domains=[sum(1<<j for j in range(h)if ha[j].bit_count()>=pa[i].bit_count() and h-1-ha[j].bit_count()>=n-1-pa[i].bit_count())for i in range(n)]
 mapping=[-1]*n;calls=0
 def dfs(dom,used,count):
  nonlocal calls
  calls+=1
  if count==n:return mapping.copy()
  u=min((i for i in range(n)if mapping[i]<0),key=lambda i:(dom[i]&~used).bit_count())
  candidates=dom[u]&~used
  while candidates:
   bit=candidates&-candidates;candidates-=bit;v=bit.bit_length()-1;mapping[u]=v;dd=dom.copy();good=True
   for i in range(n):
    if mapping[i]<0:
     dd[i]&=(ha[v]if pa[u]>>i&1 else ~ha[v]);dd[i]&=~(used|bit)
     if not dd[i]:good=False;break
   if good:
    result=dfs(dd,used|bit,count+1)
    if result is not None:return result
   mapping[u]=-1
  return None
 result=dfs(domains,0,0)
 if result is not None:result=[vertices[v]for v in result]
 return result,calls

if __name__=='__main__':
 reports=[]
 for order,strong,ms in [(5,False,[2]),(6,False,[3]),(6,True,[3]),(7,True,[3,4,5]),(7,False,[4,5])]:
  host=json.loads((BASE/f'metagraph_n{order}.json').read_text());vertices=host['strong']if strong else list(range(len(host['representatives'])))
  for m in ms:
   start=time.time();n,es=chain(m);found,calls=induced(host,vertices,n,es)
   r={'tournament_order':order,'strong_only':strong,'diamonds':m,'mapping_class_indices':found,'mapping_masks':None if found is None else [host['representatives'][v]for v in found],'search_nodes':calls,'seconds':round(time.time()-start,3)}
   print(json.dumps(r),flush=True);reports.append(r)
   (BASE/'chain_embedding_search.json').write_text(json.dumps(reports,indent=2)+'\n')
