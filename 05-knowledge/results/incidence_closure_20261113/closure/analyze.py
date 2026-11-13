#!/usr/bin/env python3
"""Exact anchor-column optimization and finite closure consequences from data.json."""
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(str(msg))
def exact_cover(masks):
 total=0
 for x in masks:total|=x
 if not total:return {'minimum':0,'columns':[],'failed_budgets':[],'search_nodes':0,'packing':[]}
 # Retain a representative of each coverage mask; dominated columns can be
 # replaced by a covering column without increasing cardinality.
 unique={}
 for i,m in enumerate(masks):
  if m and m not in unique:unique[m]=i
 candidates=[(m,i) for m,i in unique.items() if not any(m!=z and m&z==m for z in unique)]
 covers=[m for m,i in candidates];original=[i for m,i in candidates]
 index={};nb=total.bit_length()
 for bit in range(nb):
  if total>>bit&1:index[bit]=tuple(i for i,m in enumerate(covers) if m>>bit&1)
 order=sorted(index,key=lambda x:len(index[x]));candidate_bits={x:sum(1<<i for i in index[x]) for x in index}
 def packing(U):
  used=0;result=[]
  for bit in order:
   if U>>bit&1 and not used&candidate_bits[bit]:result.append(bit);used|=candidate_bits[bit]
  return result
 need(all(index.values()),'every pair has a separating anchor')
 U=total;greedy=[]
 while U:
  j=max(range(len(covers)),key=lambda j:(U&covers[j]).bit_count());greedy.append(j);U&=~covers[j]
 search_nodes=0
 @lru_cache(None)
 def solve(U,k):
  nonlocal search_nodes
  search_nodes+=1
  if not U:return ()
  if k==0 or len(packing(U))>k:return None
  maxcover=max((U&m).bit_count() for m in covers)
  if maxcover*k<U.bit_count():return None
  bit=next(b for b in order if U>>b&1)
  for j in sorted(index[bit],key=lambda j:-(U&covers[j]).bit_count()):
   ans=solve(U&~covers[j],k-1)
   if ans is not None:return (j,)+ans
  return None
 failed=[];lower=len(packing(total))
 for k in range(lower,len(greedy)+1):
  before=search_nodes;ans=solve(total,k)
  if ans is not None:break
  failed.append({'budget':k,'search_nodes':search_nodes-before})
 need(ans is not None,'greedy bounds exact search')
 cols=[original[i] for i in ans];got=0
 for i in cols:got|=masks[i]
 need(got==total and len(cols)==k,'exact anchor certificate')
 return {'minimum':k,'columns':cols,'failed_budgets':failed,'search_nodes':search_nodes,
         'packing':packing(total),'undominated_columns':len(covers),'initial_lower_bound':lower}

def main():
 data=json.loads((HERE/'data.json').read_text());results=[]
 for row in data['rows']:
  nodes=row['merged_nodes'];codes=[x['code'] for x in nodes];ix={c:i for i,c in enumerate(codes)};N=len(nodes);A=[[0]*N for _ in nodes];B=[[0]*N for _ in nodes]
  for a,b,v in row['merged_A']:A[ix[a]][ix[b]]=v
  for a,b,v in row['merged_blue']:B[ix[a]][ix[b]]=v
  groups=defaultdict(list)
  for i,x in enumerate(nodes):groups[tuple(x['Q'])].append(i)
  pairs=[p for g in groups.values() for p in combinations(g,2)];need(all(i<j for i,j in pairs),'pair order')
  item={'n':row['n'],'Q_collision_pairs':[[codes[i],codes[j]] for i,j in pairs],'anchor_optimization':{}}
  for normalized in (False,True):
   for colored in (False,True):
    matrices=[B,[[A[i][j]-B[i][j] for j in range(N)] for i in range(N)]] if colored else [A]
    masks=[sum(any(M[i][c]*(nodes[j]['degree'] if normalized else 1)!=M[j][c]*(nodes[i]['degree'] if normalized else 1) for M in matrices)<<k for k,(i,j) in enumerate(pairs)) for c in range(N)]
    need((0 if not masks else __import__('functools').reduce(int.__or__,masks,0))==(1<<len(pairs))-1,'full literal columns separate Q pairs')
    result=exact_cover(masks);result['columns']=[codes[i] for i in result['columns']]
    result['packing']=[item['Q_collision_pairs'][i] for i in result['packing']]
    item['anchor_optimization'][('P' if normalized else 'A')+('_colored' if colored else '')]=result
  item['Q_blue_eq_one_classes']=len({(tuple(x['Q']),x['blue_degree']==1) for x in nodes})
  item['Q_blue_and_degree_classes']=len({(tuple(x['Q']),x['blue_degree'],x['degree']) for x in nodes})
  if row['n']==5:
   need(item['Q_blue_eq_one_classes']==N,'one intrinsic Boolean closes n5')
   need(all(v['minimum']==1 for v in item['anchor_optimization'].values()),'one physical anchor closes n5')
  if row['n']==7:
   collision=[x for x in row['closure']['Q']['P']['levels'][1] if len(x)>1]
   need(collision==[[5,11]],'unique n7 delayed pair')
   need(A[ix[5]]==[2*(j==ix[301]) for j in range(N)] and A[ix[11]]==[j==ix[291] for j in range(N)],'delayed deterministic rows')
   need(nodes[ix[291]]['Q']==nodes[ix[301]]['Q'],'equal-Q intermediate targets')
   item['delayed_pair']={'sources':[nodes[ix[c]] for c in [5,11]],'targets':[nodes[ix[c]] for c in [291,301]]}
  results.append(item)
  print('analyzed',row['n'],{k:v['minimum'] for k,v in item['anchor_optimization'].items()},flush=True)
 out={'status':'FINITE-EXACT; candidate pending audit','rows':results,'checks':CHECKS,'data_sha256':hashlib.sha256((HERE/'data.json').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'analysis.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps({'checks':CHECKS,'analysis_sha256':hashlib.sha256((HERE/'analysis.json').read_bytes()).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
