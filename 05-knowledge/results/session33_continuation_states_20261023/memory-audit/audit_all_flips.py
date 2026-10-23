"""Independent full labelled-flip audit; writes new artifacts only."""
import itertools as it,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
checks=0
def ck(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)
source=ROOT/'metagraph-memory'
cert=json.loads((source/'full_flips_receipt.json').read_text());out=[]
for n in (3,4,5):
 free=[(u,v) for u in range(n) for v in range(u+2,n)];pairs=list(it.combinations(range(n),2));N=2**len(free);a=N-1;keys=[];edges=[]
 for x in range(N):
  es={(u,u+1) for u in range(n-1)}|{(u,v) if x>>k&1 else(v,u) for k,(u,v) in enumerate(free)};edges.append(es)
  keys.append(min(tuple(int((p[u],p[v]) in es)^flip for u,v in pairs) for p in it.permutations(range(n)) for flip in (0,1)))
 ck([x for x,k in enumerate(keys) if k==keys[a]]==[a],'singleton merged transitive fibre')
 for x in range(a):
  ck(any((u,v) in edges[x] and (v,w) in edges[x] and (w,u) in edges[x] for u,v,w in it.permutations(range(n),3)),'independent directed triangle search')
 # Radius-d observation signatures, not producer's iterative partition algorithm.
 sizes=[]
 for d in range(len(free)+2):
  shifts=[b for b in range(N) if b.bit_count()<=d]
  sizes.append(len({tuple(keys[x^b] for b in shifts) for x in range(N)}))
  if len(sizes)>1 and sizes[-1]==sizes[-2]:break
 expected=next(r for r in cert['rows'] if r['n']==n)
 ck(sizes==expected['full_flip_refinement_sizes'],'refinement depth')
 tau=[len(set(keys)),len({(keys[x],keys[x^a]) for x in range(N)})]
 if tau[-1]!=tau[-2]:tau.append(tau[-1])
 ck(tau==expected['tau_refinement_sizes'],'tau depth')
 tested=0
 for x,y in it.combinations(range(N),2):
  word=x^a
  ck(keys[x^word]!=keys[y^word],'anchor separating word');tested+=1
 ck(tested==expected['pairwise_separations'],'pair count')
 # Each translated singleton has a distinct single support: exact independent rank certificate.
 supports=[tuple(x for x in range(N) if keys[x^shift]==keys[a]) for shift in range(N)]
 ck(set(supports)=={(x,) for x in range(N)},'translated singleton basis')
 ck(expected['linear_observer_dimension']==expected['positive_hidden_state_minimum']==N,'claimed minima')
 out.append(dict(n=n,tilings=N,full_refinement=sizes,tau_refinement=tau,pair_checks=tested,translated_indicator_rank=N))
result={'status':'PASS','checks':checks,'rows':out,'source_pins':{p:hashlib.sha256((source/p).read_bytes()).hexdigest() for p in ('ADDENDUM_FULL_FLIPS.md','check_full_flips.py','full_flips_receipt.json')}}
Path(__file__).with_name('all_flips_audit_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
