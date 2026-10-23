#!/usr/bin/env python3
"""Exact adaptive observers and frozen electrical path reconfiguration.

Independent of producer algorithms. Uses original pinned atlas bytes and keeps
all actual labelled paths; no quotient or numerical eigenvalue calculations.
"""
from array import array
from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import argparse, csv, gzip, json, sys
HERE=Path(__file__).resolve().parent
CHECKS=0

def need(ok,message):
 global CHECKS
 CHECKS+=1
 if not ok: raise RuntimeError(message)

def distance(p,q): return sum(a!=b for a,b in zip(p,q))

def adjacency(word,n=8):
 A=[[0]*n for _ in range(n)]
 for bit,(i,j) in enumerate(combinations(range(n),2)):
  A[i][j]=(word>>(n*(n-1)//2-1-bit))&1; A[j][i]=1-A[i][j]
 return A

def isomorphism(A,B):
 """Backtracking with degree colors and already assigned literal arc checks."""
 n=len(A); da=list(map(sum,A)); db=list(map(sum,B)); order=sorted(range(n),key=lambda v:(db.count(da[v]),v)); p=[-1]*n
 def search(k,unused):
  if k==n:return tuple(p)
  u=order[k]
  for v in sorted(unused):
   if db[v]!=da[u] or any(A[u][w]!=B[v][p[w]] for w in order[:k]):continue
   p[u]=v; result=search(k+1,unused-{v})
   if result is not None:return result
  p[u]=-1
  return None
 return search(0,set(range(n)))

def all_paths_dfs(A):
 n=len(A); answer=[]
 def visit(p,left):
  if not left:answer.append(tuple(p));return
  for v in sorted(left):
   if not p or A[p[-1]][v]:visit(p+[v],left-{v})
 visit([],set(range(n)));return answer

def adaptive(paths,values):
 """D is minimax. E is mean-optimal total cost, without depth restriction."""
 N=len(paths); n=len(paths[0]); columns=[[sum(1<<i for i,p in enumerate(paths) if p[j]==v) for v in range(n)] for j in range(n)]
 colors=[sum(1<<i for i,v in enumerate(values) if v==value) for value in sorted(set(values))]
 @lru_cache(None)
 def dp(S):
  if any(S&C==S for C in colors):return (0,0,0,None,None)
  candidates=[]
  for j in range(n):
   parts=[S&C for C in columns[j] if S&C]
   if len(parts)<2:continue
   children=[dp(C) for C in parts]
   candidates.append((1+max(x[0] for x in children),S.bit_count()+sum(x[1] for x in children),1+max(x[2] for x in children),j))
  D=min(x[0] for x in candidates); E=min(candidates,key=lambda x:(x[1],x[2],x[3]))
  return D,E[1],E[2],min(x[3] for x in candidates if x[0]==D),E[3]
 @lru_cache(None)
 def cap(S,budget):
  if dp(S)[0]==0:return (0,None)
  if dp(S)[0]>budget:return (10**9,None)
  choices=[]
  for j in range(n):
   parts=[S&C for C in columns[j] if S&C]
   if len(parts)<2 or any(dp(C)[0]>budget-1 for C in parts):continue
   choices.append((S.bit_count()+sum(cap(C,budget-1)[0] for C in parts),j))
  return min(choices)
 def tree(S,budget):
  if dp(S)[0]==0:return {'answer':values[(S&-S).bit_length()-1], 'mass':S.bit_count()}
  j=cap(S,budget)[1]
  children={str(v):tree(S&C,budget-1) for v,C in enumerate(columns[j]) if S&C}
  need(len(children)>=2,'nonconstant adaptive query')
  return {'clock':j,'mass':S.bit_count(),'children':children}
 full=(1<<N)-1; result=dp(full); tree_root=tree(full,result[0]); capped=cap(full,result[0])[0]
 arities=Counter()
 def evaluate(node,indices,depth):
  need(node['mass']==len(indices),'decision tree mass')
  if 'answer' in node:
   need(all(values[i]==node['answer'] for i in indices),'every adaptive terminal has one correct output')
   return depth,len(indices)*depth
  j=node['clock']; arities[len(node['children'])]+=1; depths=[];cost=0
  for owner,child in node['children'].items():
   ids=[i for i in indices if paths[i][j]==int(owner)];d,e=evaluate(child,ids,depth+1);depths.append(d);cost+=e
  return max(depths),cost
 depth,cost=evaluate(tree_root,list(range(N)),0)
 need((depth,cost)==(result[0],capped),'explicit optimal worst-case tree and depth-capped mean')
 return dict(worst_depth=result[0],optimal_mean=str(Fraction(result[1],N)),mean_policy_worst_depth=result[2],
             best_mean_at_optimal_worst_depth=str(Fraction(capped,N)),state_count=dp.cache_info().currsize,
             minimax_first_clock=result[3],mean_first_clock=result[4],tree_arity_histogram=dict(sorted(arities.items())),tree=tree_root)

def reconfiguration(records):
 groups=defaultdict(list)
 for i,r in enumerate(records):groups[(r['path'][0],r['path'][-1],r['phase'])].append(i)
 hist=Counter(); out=[]
 for key,ids in sorted(groups.items()):
  if len(ids)==1:budget=0;tree=[]
  else:
   parent={i:i for i in ids}
   def find(i):
    while i!=parent[i]:parent[i]=parent[parent[i]];i=parent[i]
    return i
   edges=sorted((distance(records[i]['path'],records[j]['path']),i,j) for i,j in combinations(ids,2));tree=[]
   for d,i,j in edges:
    a,b=find(i),find(j)
    if a!=b:parent[a]=b;tree.append([d,i,j])
    if len(tree)==len(ids)-1:break
   budget=max(e[0] for e in tree)
   # Independent threshold flood fill; proves lower threshold fails and upper succeeds.
   for bound,expect in [(budget-1,False),(budget,True)]:
    seen={ids[0]};queue=[ids[0]]
    while queue:
     i=queue.pop()
     for j in ids:
      if j not in seen and distance(records[i]['path'],records[j]['path'])<=bound:seen.add(j);queue.append(j)
    need((len(seen)==len(ids))==expect,'MST bottleneck and direct threshold connectivity agree')
  hist[budget]+=1;out.append(dict(endpoint_phase=list(key),indices=ids,budget=budget,spanning_tree=tree))
 return dict(budget_histogram=dict(sorted(hist.items())),fibres=out)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path);parser.add_argument('--write',action='store_true');args=parser.parse_args()
 pins=json.loads((HERE/'SOURCE_PINS.json').read_text());raws={}
 for row in pins['files']:
  p=HERE/'inputs'/row['local']
  if not p.exists():
   candidates=([args.repo] if args.repo else [])+list(HERE.parents)
   p=next((base/row['path'] for base in candidates if (base/row['path']).exists()),p)
  raw=p.read_bytes();need(sha256(raw).hexdigest()==row['sha256'],'pinned original input '+row['path']);raws[row['path']]=raw
 def raw(suffix):return next(v for k,v in raws.items() if k.endswith(suffix))
 nodes=list(csv.DictReader(raw('metagraph_order8_20261023/merged_nodes.tsv').decode().splitlines(),delimiter='\t'))
 original=json.loads(raw('metagraph_order8_electrical_20261023/receipt.json'))['order8'];x=original['higher_layers'][-1]['lifts'][0]
 need(original['root']==3527 and x[3527]==0 and len(x)==3528,'original grounded mod256 character')
 need([x[i] for i in (1238,2242,2246,2315)]==[97,18,65,212],'original switch entries')
 owners=array('H'); owners.frombytes(gzip.decompress(raw('owner_map.u16.gz')))
 if sys.byteorder!='little':owners.byteswap()
 need(len(owners)==1<<21,'complete native forward owner map')
 merged={};canonical={}
 for row in nodes:
  for letter in ('a','b'):
   ci=int(row['class_'+letter]);merged[ci]=int(row['merged_id']);canonical[ci]=adjacency(int(row['canon_'+letter]))
 chords=[(i,j) for i in range(8) for j in range(i+2,8)];full=(1<<21)-1; reports=[];allrecords={}
 for source,H in zip((1238,2242,2246,2315),(409,407,523,527)):
  ci=int(nodes[source]['class_a']);A=canonical[ci]
  paths=[p for p in permutations(range(8)) if all(A[u][v] for u,v in zip(p,p[1:]))]
  need(paths==all_paths_dfs(A) and len(paths)==H,'two independent complete Hamiltonian enumerations')
  aut=[g for g in permutations(range(8)) if all(A[i][j]==A[g[i]][g[j]] for i,j in combinations(range(8),2))]
  need(len(aut)==1,'these four actual source automorphism groups are trivial')
  records=[]
  for p in paths:
   mask=sum(A[p[i]][p[j]]<<b for b,(i,j) in enumerate(chords));need(owners[mask]==ci,'source forward-bit owner')
   target_ci=owners[full^mask];target=merged[target_ci];phase=(x[source]-x[target])%8
   kept={frozenset(e) for e in zip(p,p[1:])}
   K=[[A[i][j] if i==j or frozenset((i,j)) in kept else 1-A[i][j] for j in range(8)] for i in range(8)]
   need(all(K[u][v] for u,v in zip(p,p[1:])),'literal partner retains actual Hamiltonian path')
   need(sum(K[p[i]][p[j]]<<b for b,(i,j) in enumerate(chords))==full^mask,'literal partner equals complement mask')
   f=isomorphism(K,canonical[target_ci]);need(f is not None,'independent literal partner-to-canonical isomorphism')
   need(all(K[i][j]==canonical[target_ci][f[i]][f[j]] for i,j in combinations(range(8),2)),'target isomorphism certificate')
   reverse=tuple(reversed(p));reverse_mask=sum(A[reverse[j]][reverse[i]]<<b for b,(i,j) in enumerate(chords))
   need(merged[owners[reverse_mask]]==source and merged[owners[full^reverse_mask]]==target,'literal converse reversal retains merged source and target')
   records.append(dict(path=list(p),phase=phase,target=target,target_ordinary=target_ci,mask=mask,target_isomorphism=list(f)))
  allrecords[str(source)]=records
  objectives={key:adaptive(paths,[tuple(r['path']) if key=='path' else r[key] for r in records]) for key in ('phase','target','path')}
  reports.append(dict(source=source,ordinary_source=ci,canonical_word=int(nodes[source]['canon_a']),adjacency=A,H=H,automorphisms=1,
                      objectives=objectives,reconfiguration=reconfiguration(records)))
 # Explicit frozen-character hostile and an unguarded route, selected before final replay.
 source=1238;recs=allrecords[str(source)];lookup={tuple(r['path']):r for r in recs}
 route=[(0,4,2,1,6,7,5,3),(0,4,2,6,7,1,5,3),(0,6,7,4,2,1,5,3),(0,6,7,2,1,5,4,3)]
 witness=[lookup[p] for p in route];fibre=[r for r in recs if (r['path'][0],r['path'][-1],r['phase'])==(0,3,2)]
 need([tuple(r['path']) for r in fibre]==[route[0],route[-1]],'complete pinned frozen-phase fibre has precisely two survivors')
 need(distance(route[0],route[-1])==6 and all(distance(a,b)<=4 for a,b in zip(route,route[1:])),'sharp six-clock move versus unguarded four-clock route')
 need([r['phase'] for r in witness]==[2,0,7,2] and [r['target'] for r in witness]==[2289,2501,2684,2009],'native hostile target and phase records')
 result=dict(status='FINITE-EXACT',source_commit=pins['commit'],sources=reports,hostile=dict(source=1238,pins=[[0,0],[7,3]],phase=2,guarded_minimum_budget=6,unguarded_route_maximum_budget=4,route=witness),checks=CHECKS)
 text=json.dumps(result,sort_keys=True,indent=2)+'\n';records_text=json.dumps(allrecords,sort_keys=True,indent=2)+'\n'
 if args.write:(HERE/'receipt.json').write_text(text);(HERE/'path_records.json').write_text(records_text)
 else:
  need(json.loads((HERE/'receipt.json').read_text())==json.loads(text),'frozen exact receipt')
  need(json.loads((HERE/'path_records.json').read_text())==allrecords,'frozen exact actual path records')
 print(json.dumps(dict(status='PASS',checks=CHECKS,depths={str(r['source']):{k:v['worst_depth'] for k,v in r['objectives'].items()} for r in reports})))

if __name__=='__main__':main()
