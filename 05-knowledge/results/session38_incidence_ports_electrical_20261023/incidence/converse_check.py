"""Exact converse-paired incidence menu; frozen main packet is not modified."""
from pathlib import Path
from itertools import combinations
from collections import defaultdict,Counter
import argparse,csv,json,hashlib
HERE=Path(__file__).resolve().parent
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise AssertionError(msg)
def incidence(code):
 n=8;A=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):A[i][j]=code>>(27-k)&1;A[j][i]=1-A[i][j]
 out=[sum(A[u][v]<<v for v in range(n)) for u in range(n)];inc=[255^x^(1<<u) for u,x in enumerate(out)];E=[[0]*n for _ in range(n)]
 for u in range(n):
  for v in range(n):
   if A[u][v]:E[u][v]=E[v][u]=(out[v]&inc[u]).bit_count()
 c=[sum(E[u][v] for v in range(n) if A[u][v]) for u in range(n)]
 need(all(c[u]==sum(E[u][v] for v in range(n) if A[v][u]) for u in range(n)),'Eulerian incidence balance')
 return A,c,E
def sig(A,c,E,p,q,converse=False):
 radix=7 if q is None else q
 return b''.join(sorted(bytes([sum(A[v][u] if converse else A[u][v] for v in range(8))]+sorted((c[v] if p is None else c[v]%p)*radix+(E[u][v] if q is None else E[u][v]%q) for v in range(8) if (A[v][u] if converse else A[u][v]))) for u in range(8)))
def forest_info(colors,convs,canons):
 G=defaultdict(list)
 for cid,color in enumerate(colors):G[color].append(cid)
 lab={s:i for i,s in enumerate(G)};m=len(G);par=list(range(2*m));size=[1]*(2*m);edges=[];cycles=[]
 def find(v):
  while par[v]!=v:par[v]=par[par[v]];v=par[v]
  return v
 for cid,s in enumerate(colors):
  a=lab[s];b=m+lab[colors[convs[cid]]];edges.append((a,b,cid));aa,bb=find(a),find(b)
  if aa==bb:cycles.append(cid)
  else:
   if size[aa]<size[bb]:aa,bb=bb,aa
   par[bb]=aa;size[aa]+=size[bb]
 bad=[i for i,s in enumerate(colors) if len(G[s])>1 and len(G[colors[convs[i]]])>1]
 paired=defaultdict(list);unordered=defaultdict(list)
 for i,s in enumerate(colors):
  t=colors[convs[i]];paired[s,t].append(i);unordered[tuple(sorted((s,t)))].append(i)
 pc=[[canons[i] for i in v] for v in paired.values() if len(v)>1]
 component_edges=Counter(find(a) for a,b,cid in edges);component_sizes=Counter(component_edges.values())
 return {'ordinary_signatures':m,'paired_signatures':len(paired),'unordered_paired_signatures':len(unordered),'singleton_or_converse_singleton':not bad,'bad_classes':len(bad),'forest':not cycles,'cycle_rank':len(cycles),'vertices':2*m,'components':len(component_edges),'component_edge_count_histogram':[[k,v] for k,v in sorted(component_sizes.items())],'first_cycle_closing_class':canons[cycles[0]] if cycles else None,'first_paired_collision':pc[0] if pc else None},G,pc

def main():
 p=argparse.ArgumentParser();p.add_argument('--classes',type=Path,default=HERE.parent/'session38-wiring/classes_input_n8.tsv');p.add_argument('--write',action='store_true');args=p.parse_args()
 R=sorted(csv.DictReader(args.classes.open(),delimiter='\t'),key=lambda r:int(r['canon']));need(len(R)==6880,'class count');lookup={int(r['class']):i for i,r in enumerate(R)};convs=[lookup[int(r['sigmaClass'])] for r in R];need(all(convs[convs[i]]==i for i in range(6880)),'converse involution');merged=len({min(i,convs[i]) for i in range(6880)})
 canons=[int(r['canon']) for r in R];keys=[(p,q) for p in range(1,9) for q in range(1,8)]+[(None,2),(2,None),(3,None),(None,None)];by={key:[] for key in keys}
 for cid,row in enumerate(R):
  A,c,E=incidence(int(row['canon']))
  for key,L in by.items():L.append(sig(A,c,E,*key))
 for cid,row in enumerate(R):
  A,c,E=incidence(int(row['canon']))
  for key,L in by.items():need(sig(A,c,E,*key,converse=True)==L[convs[cid]],'independent converse signature mismatch')
 O={'scope':'FINITE-EXACT all6880classes; paired class IDs versus separately retained converse marginals explicitly distinguished','classes_sha256':hashlib.sha256(args.classes.read_bytes()).hexdigest(),'classes':6880,'merged_nodes':merged,'menu':{}}
 ids={}
 for key,L in by.items():
  name=('exact' if key[0] is None else str(key[0]))+'_'+('exact' if key[1] is None else str(key[1]));info,G,pc=forest_info(L,convs,canons);O['menu'][name]=info
  if key in [(3,2),(3,3),(2,3),(2,2),(None,2),(2,None),(3,None)]:
   O['menu'][name]['paired_collision_groups']=pc
  ids[name]={cid:gid for gid,vs in enumerate(G.values()) for cid in vs}
 G33=defaultdict(list);G3e=defaultdict(list)
 for i in range(6880):G33[by[3,3][i]].append(i);G3e[by[3,None][i]].append(i)
 same=sorted(map(tuple,G33.values()))==sorted(map(tuple,G3e.values()));need(same,'E33 and E3exact partitions differ');O['E33_E3exact_same_partition']=same
 O['E33_collision_converses']=[{'classes':[canons[i] for i in ids],'converses':[canons[convs[i]] for i in ids],'converse_fiber_sizes':[len(G33[by[3,3][convs[i]]]) for i in ids]} for ids in G33.values() if len(ids)>1]
 need(all(v==1 for g in O['E33_collision_converses'] for v in g['converse_fiber_sizes']),'E33 converse singleton test')
 O['checks']=checks
 with (HERE/'converse_signature_ids.tsv').open('w') as f:
  W=csv.writer(f,delimiter='\t');W.writerow(['class_id','canon','converse_class']+list(ids))
  for i in range(6880):W.writerow([i,canons[i],convs[i]]+[D[i] for D in ids.values()])
 dest=HERE/'converse_receipt.json'
 if args.write:dest.write_text(json.dumps(O,indent=2)+'\n')
 else:need(json.loads(dest.read_text())==O,'frozen converse receipt differs')
 print(json.dumps(O,indent=2))
if __name__=='__main__':main()
