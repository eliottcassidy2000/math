#!/usr/bin/env python3
"""Complete finite certificate for the all-order H=67 classification.
Standard library only. All checks survive python -O. No repository data is read.
"""
from itertools import combinations, permutations, product
from collections import Counter
from pathlib import Path
import hashlib,json

def require(x,msg):
 if not x:raise RuntimeError(msg)

def canonical(rows):
 n=len(rows);colors=[r.bit_count() for r in rows]
 while True:
  sig=[(colors[u],tuple(sorted(colors[v] for v in range(n) if rows[u]>>v&1))) for u in range(n)]
  pal={v:i for i,v in enumerate(sorted(set(sig)))};nc=[pal[s] for s in sig]
  if len(set(nc))==len(set(colors)):break
  colors=nc
 cells=[[u for u in range(n) if nc[u]==c] for c in sorted(set(nc))]
 best=None
 for ps in product(*(permutations(c) for c in cells)):
  p=sum(ps,());inv=[0]*n
  for j,u in enumerate(p):inv[u]=j
  out=tuple(sum(1<<inv[v] for v in range(n) if rows[u]>>v&1) for u in p)
  if best is None or out<best:best=out
 return best

def extend_classes(ts):
 ans=set();m=len(next(iter(ts)))
 for rows in ts:
  for mask in range(1<<m):
   r=[x|(((mask>>u)&1)<<m) for u,x in enumerate(rows)]+[((1<<m)-1)^mask]
   ans.add(canonical(r))
 return ans

def count_paths(rows):
 n=len(rows)
 if not n:return 1
 dp=[[0]*n for _ in range(1<<n)]
 for u in range(n):dp[1<<u][u]=1
 for mask in range(1,1<<n):
  for u,count in enumerate(dp[mask]):
   if not count:continue
   todo=rows[u]&~mask
   while todo:
    bit=todo&-todo;todo-=bit
    dp[mask|bit][bit.bit_length()-1]+=count
 return sum(dp[-1])

def count_cycles(rows):
 n=len(rows)
 if n<3:return 0
 full=(1<<n)-1
 def visit(mask,u):
  if mask==full:return int(bool(rows[u]&1))
  todo=rows[u]&~mask;total=0
  while todo:
   bit=todo&-todo;todo-=bit;total+=visit(mask|bit,bit.bit_length()-1)
  return total
 return visit(1,0)

def literal_counts(rows):
 n=len(rows);h=c=0
 for p in permutations(range(n)):
  if all(rows[a]>>b&1 for a,b in zip(p,p[1:])):
   h+=1
   if n>=3 and p[0]==0 and rows[p[-1]]&1:c+=1
 return h,c

def all_labelled(n):
 pairs=list(combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  rows=[0]*n
  for j,(u,v) in enumerate(pairs):
   if mask>>j&1:rows[u]|=1<<v
   else:rows[v]|=1<<u
  yield rows

def strong(rows):
 n=len(rows)
 if not n:return False
 reach=[r|(1<<u) for u,r in enumerate(rows)]
 for k in range(n):
  for u in range(n):
   if reach[u]>>k&1:reach[u]|=reach[k]
 return all(x==(1<<n)-1 for x in reach)

def transitive_order(rows,vertices):
 p=sorted(vertices,key=lambda u:-sum(rows[u]>>v&1 for v in vertices))
 return p if all(rows[u]>>v&1 for i,u in enumerate(p) for v in p[i+1:]) else None

def apex_presentations(rows):
 n=len(rows);ans=[]
 for a,b in permutations(range(n),2):
  if not rows[a]>>b&1:continue
  order=transitive_order(rows,[v for v in range(n) if v not in (a,b)])
  if order is not None:
   word=[((rows[a]>>v)&1)+2*((rows[b]>>v)&1) for v in order]
   ans.append({'roots':[a,b],'base_order':order,'word':word})
 return sorted(ans,key=lambda x:x['word'])

def automorphisms(rows):
 n=len(rows);total=0
 for p in permutations(range(n)):
  if all(((rows[u]>>v)&1)==((rows[p[u]]>>p[v])&1) for u in range(n) for v in range(n)):
   total+=1
 return total

def main():
 require(67**3==300763 and 5**8==390625 and 67**3<5**8,'Busch cutoff arithmetic')
 literal_checks=0
 for n in range(1,6):
  for r in all_labelled(n):
   require((count_paths(r),count_cycles(r))==literal_counts(r),f'literal counter mismatch {r}')
   literal_checks+=1
 ts={(0,)};inventory=[];hits=[];totals=[1,1,2,4,12,56,456,6880]
 for n in range(1,9):
  if n>1:ts=extend_classes(ts)
  require(len(ts)==totals[n-1],f'unlabelled inventory at {n}')
  hist=Counter();target=[]
  for r in sorted(ts):
   h=count_paths(r);c=count_cycles(r);is_strong=strong(r)
   require((c>0)==(is_strong and n>=3),'cycle/strongness mismatch')
   if c in (1,2):hist[c]+=1
   if h==67:
    target.append({'rows':list(r),'C':c,'strong':is_strong})
    if is_strong:hits.append((n,r,c))
  inventory.append({'order':n,'classes':len(ts),'C1_classes':hist[1],'C2_classes':hist[2],'H67':target})
 require(len(hits)==6 and all(n==7 for n,r,c in hits),'strong H67 classification')
 require(Counter(c for n,r,c in hits)=={3:3,4:3},'cycle distribution')
 require(len(inventory[-1]['H67'])==12 and all(not x['strong'] for x in inventory[-1]['H67']),'weak order8 target classes')
 cores=[];rr=[r for n,r,c in hits]
 for i,(n,r,c) in enumerate(hits,1):
  ps=apex_presentations(r)
  require(ps,'target core lacks two-apex presentation')
  require(all(transitive_order(r,[v for v in range(n) if v!=u]) is None for u in range(n)),'target core is one-apex')
  require(automorphisms(r)==1,'target core is not rigid')
  rev=[sum(1<<v for v in range(n) if r[v]>>u&1) for u in range(n)]
  cores.append({'id':i,'rows':list(r),'H':67,'C':c,'automorphism_order':1,'feedback_vertex_number':2,'converse_id':rr.index(canonical(rev))+1,'two_apex_presentations':ps})
 payload={'bound':67,'maximum_strong_order':8,'literal_labelled_crosschecks':literal_checks,'inventory':inventory,'strong_cores':cores,'all_order_unlabelled_count':'0 for n<7; 6*(n-6) for n>=7'}
 out=Path(__file__).resolve().parent
 data=json.dumps(payload,sort_keys=True,indent=2)+'\n';(out/'certificate.json').write_text(data)
 receipt={'literal_labelled_crosschecks':literal_checks,'classes_by_order':totals,'strong_H67_cores':6,'cycle_distribution':{'3':3,'4':3},'minimum_strong_cycle_count':3,'all_cores_order':7,'all_cores_rigid':True,'all_cores_feedback_vertex_number':2,'certificate_sha256':hashlib.sha256(data.encode()).hexdigest()}
 text=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(out/'receipt.json').write_text(text);print(text,end='')

if __name__=='__main__':main()
