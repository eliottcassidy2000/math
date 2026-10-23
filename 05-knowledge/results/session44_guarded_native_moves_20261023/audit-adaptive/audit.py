"""Independent adaptive-observer audit: no producer imports; explicit checks."""
from array import array
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
import argparse,csv,gzip,json,sys
HERE=Path(__file__).resolve().parent
CHECKS=0

def need(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(msg)

def objective(paths,values,published):
 N=len(paths);full=tuple(range(N))
 @lru_cache(None)
 def terminal(S):return len({values[i] for i in S})==1
 @lru_cache(None)
 def choices(S):
  result=[]
  for j in range(8):
   buckets=defaultdict(list)
   for i in S:buckets[paths[i][j]].append(i)
   if len(buckets)>1:result.append((j,tuple(tuple(v) for k,v in sorted(buckets.items()))))
  return tuple(result)
 @lru_cache(None)
 def possible(S,b):
  if terminal(S):return True
  if b==0:return False
  return any(all(possible(C,b-1) for C in children) for j,children in choices(S))
 D=next(b for b in range(9) if possible(full,b))
 need(D==published['worst_depth'],'independent exact minimax depth')
 need(D==0 or not possible(full,D-1),'exhaustive smaller budget infeasibility')
 @lru_cache(None)
 def mean(S):
  if terminal(S):return 0
  return len(S)+min(sum(mean(C) for C in children) for j,children in choices(S))
 @lru_cache(None)
 def cap(S,b):
  if terminal(S):return 0
  if b==0:return None
  costs=[]
  for j,children in choices(S):
   vals=[cap(C,b-1) for C in children]
   if all(v is not None for v in vals):costs.append(len(S)+sum(vals))
  return min(costs) if costs else None
 E=mean(full);C=cap(full,D)
 need(Fraction(E,N)==Fraction(published['optimal_mean']),'independent unrestricted mean optimum')
 need(Fraction(C,N)==Fraction(published['best_mean_at_optimal_worst_depth']),'independent depth-capped mean optimum')
 need(cap(full,D-1) is None if D else cap(full,0)==0,'separate capped lower bound')
 mean_depth=next(b for b in range(D,9) if cap(full,b)==E)
 need(mean_depth==published['mean_policy_worst_depth'],'minimum worst depth of all mean-optimal policies')
 def replay(node,S,history):
  need(node['mass']==len(S),'tree exact state mass')
  if 'answer' in node:
   v=tuple(node['answer']) if isinstance(node['answer'],list) else node['answer']
   need(terminal(S) and all(values[i]==v for i in S),'tree exact terminal observable')
   return 0,0
  j=node['clock'];need(0<=j<8 and j not in history,'tree actual fresh clock')
  buckets=defaultdict(list)
  for i in S:buckets[paths[i][j]].append(i)
  need(len(buckets)>1 and set(node['children'])=={str(v) for v in buckets},'tree exhaustive answer partition')
  vals=[replay(node['children'][str(v)],tuple(ids),history+(j,)) for v,ids in buckets.items()]
  return 1+max(x[0] for x in vals),len(S)+sum(x[1] for x in vals)
 need(replay(published['tree'],full,())==(D,C),'published policy attains both certified objectives')
 sufficient=[]
 for b in range(9):
  for clocks in combinations(range(8),b):
   lookup={};ok=True
   for i,p in enumerate(paths):
    key=tuple(p[j] for j in clocks)
    if key in lookup and lookup[key]!=values[i]:ok=False;break
    lookup[key]=values[i]
   if ok:sufficient.append(clocks)
  if sufficient:break
 need(b==6,'exhaustive nonadaptive clock lower and upper bounds')
 return dict(paths=N,worst_depth=D,smaller_budget_impossible=True,optimal_mean=str(Fraction(E,N)),
   best_mean_at_optimal_worst_depth=str(Fraction(C,N)),least_worst_depth_among_mean_optimal=mean_depth,
   static_optimum=b,all_minimum_static_clock_sets=sufficient,mean_states=mean.cache_info().currsize,
   feasibility_states=possible.cache_info().currsize,capped_states=cap.cache_info().currsize)

def main():
 p=argparse.ArgumentParser();p.add_argument('--packet',type=Path,default=HERE.parent/'session44-actual-guards');p.add_argument('--repo',type=Path);p.add_argument('--write',action='store_true');args=p.parse_args();packet=args.packet
 files=('ADAPTIVE_OBSERVERS.md','check.py','independent.py','receipt.json','path_records.json','SOURCE_PINS.json','MANIFEST.json')
 hashes={f:sha256((packet/f).read_bytes()).hexdigest() for f in files}
 need(hashes==json.loads((HERE/'SOURCE_HASHES.json').read_text()),'producer snapshot frozen')
 R=json.loads((packet/'receipt.json').read_text());records=json.loads((packet/'path_records.json').read_text());pins=json.loads((packet/'SOURCE_PINS.json').read_text());raws={}
 for row in pins['files']:
  p=packet/'inputs'/row['local']
  if not p.exists() and args.repo:p=args.repo/row['path']
  raw=p.read_bytes();need(sha256(raw).hexdigest()==row['sha256'],'original native input hash');raws[row['path']]=raw
 def raw(suffix):return next(v for k,v in raws.items() if k.endswith(suffix))
 nodes=list(csv.DictReader(raw('merged_nodes.tsv').decode().splitlines(),delimiter='\t'))
 E=json.loads(raw('metagraph_order8_electrical_20261023/receipt.json'))['order8'];x=E['higher_layers'][-1]['lifts'][0]
 owner=array('H');owner.frombytes(gzip.decompress(raw('owner_map.u16.gz')))
 if sys.byteorder!='little':owner.byteswap()
 need(len(owner)==1<<21 and len(x)==3528 and E['root']==3527 and x[3527]==0,'native domain and original grounding')
 def matrix(word):
  A=[[0]*8 for _ in range(8)]
  for k,(i,j) in enumerate(combinations(range(8),2)):A[i][j]=(word>>(27-k))&1;A[j][i]=1-A[i][j]
  return A
 merged={};canons={}
 for row in nodes:
  for letter in ('a','b'):
   ci=int(row['class_'+letter]);merged[ci]=int(row['merged_id']);canons[ci]=matrix(int(row['canon_'+letter]))
 chords=[(i,j) for i in range(8) for j in range(i+2,8)];full=(1<<21)-1;answers=[]
 for rr in R['sources']:
  source=rr['source'];ci=int(nodes[source]['class_a']);A=canons[ci];saved=records[str(source)]
  paths=[(i,) for i in range(8)]
  for length in range(1,8):paths=[p+(v,) for p in paths for v in range(8) if v not in p and A[p[-1]][v]]
  paths=sorted(paths)
  need(paths==[tuple(r['path']) for r in saved] and len(paths)==rr['H'],'independent complete breadth-by-length path enumeration')
  aut=[g for g in permutations(range(8)) if all(A[i][j]==A[g[i]][g[j]] for i,j in combinations(range(8),2))]
  need(aut==[tuple(range(8))] and rr['automorphisms']==1,'full actual kernel retained and group trivial')
  for r in saved:
   p=tuple(r['path']);q=sum(A[p[i]][p[j]]<<bit for bit,(i,j) in enumerate(chords));ti=owner[full^q]
   need(q==r['mask'] and owner[q]==ci and merged[ti]==r['target'] and ti==r['target_ordinary'],'actual source target and mask')
   need(r['phase']==(x[source]-x[r['target']])%8,'original character value')
   K=[v.copy() for v in A];kept={frozenset((i,j)) for i,j in zip(p,p[1:])}
   for i,j in combinations(range(8),2):
    if frozenset((i,j)) not in kept:K[i][j]^=1;K[j][i]^=1
   iso=r['target_isomorphism'];B=canons[ti]
   need(sorted(iso)==list(range(8)) and all(K[i][j]==B[iso[i]][iso[j]] for i in range(8) for j in range(8)),'literal partner isomorphism verified every ordered pair')
   reverse=p[::-1];rq=sum(A[reverse[j]][reverse[i]]<<bit for bit,(i,j) in enumerate(chords))
   need(merged[owner[rq]]==source and merged[owner[full^rq]]==r['target'],'actual converse reversal and owners')
  obs={}
  for key in ('phase','target','path'):
   vals=[tuple(r['path']) if key=='path' else r[key] for r in saved]
   obs[key]=objective(paths,vals,rr['objectives'][key])
  answers.append(dict(source=source,H=len(paths),ordinary_source=ci,uniform_measure='each actual labelled Hamiltonian path has mass1/H',automorphism_order=1,objectives=obs))
 result=dict(status='PASS independent adaptive lower and upper certificates',source_sha256=hashes,sources=answers,checks=CHECKS,
   universal_bank_corollary='Depends on independently audited original-network theorem identifying fixed-source output partition exactly with target partition; no extra optimisation needed.')
 text=json.dumps(result,indent=2,sort_keys=True)+'\n'
 if args.write:(HERE/'receipt.json').write_text(text)
 else:need((HERE/'receipt.json').read_text()==text,'frozen independent audit receipt')
 print(json.dumps(dict(status='PASS',checks=CHECKS,source_depths={r['source']:{k:v['worst_depth'] for k,v in r['objectives'].items()} for r in answers},receipt_sha256=sha256(text.encode()).hexdigest()),sort_keys=True))

if __name__=='__main__':main()
