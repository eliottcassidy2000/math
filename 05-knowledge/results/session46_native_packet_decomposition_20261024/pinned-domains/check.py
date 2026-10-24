"""Authenticate complete common-clock domains for the 28-source supply."""
import argparse,csv,gzip,hashlib,json,sys
from array import array
from itertools import combinations,permutations
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--supply',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
R=a.repo/'05-knowledge/results'; f=R/'metagraph_order8_20261023';checks=0

def need(c,m):
 global checks
 checks+=1
 if not c:raise RuntimeError(m)

sources=[f/'merged_nodes.tsv',f/'owner_map.u16.gz',R/'metagraph_order8_electrical_20261023/receipt.json',a.supply]
nodes=list(csv.DictReader(sources[0].open(),delimiter='\t'));owner=array('H');owner.frombytes(gzip.decompress(sources[1].read_bytes()))
if sys.byteorder!='little':owner.byteswap()
merged={int(n[k]):int(n['merged_id']) for n in nodes for k in ('class_a','class_b')}
D=json.loads(sources[2].read_text())['order8'];Y=D['radical_lifts_mod4'];X=D['higher_layers'][-1]['lifts'][0]
rows=json.loads(a.supply.read_text())['records'];need(len(rows)==28,'full supplied28')
pairs=list(combinations(range(8),2));chords=[(i,j) for i,j in pairs if j>=i+2];indices={e:i for i,e in enumerate(chords)};full=(1<<21)-1

def label(s):return [Y[j][s]%2 for j in range(4)]+[X[s]%128]
def reflect(w):return sum(((w>>i)&1)<<indices[7-v,7-u] for i,(u,v) in enumerate(chords))
def orbit(w):return {w,full^w,reflect(w),full^reflect(w)}
def partner(B,path):
 E={frozenset(e) for e in zip(path,path[1:])}
 return [[0 if u==v else B[u][v] if frozenset((u,v)) in E else B[v][u] for v in range(8)] for u in range(8)]
def transpose(B):return list(map(list,zip(*B)))
def autos(B):
 deg=list(map(sum,B));groups=[[j for j,d in enumerate(deg) if d==deg[i]] for i in range(8)];assignment={};out=[]
 def rec():
  if len(assignment)==8:out.append(tuple(assignment[i] for i in range(8)));return
  i=next(i for i in range(8) if i not in assignment)
  for j in groups[i]:
   if j in assignment.values() or any(B[i][u]!=B[j][v] for u,v in assignment.items()):continue
   assignment[i]=j;rec();del assignment[i]
 rec();return out
reports=[];globalresources={};globalports=set()
for r in rows:
 s=r['source'];code=int(nodes[s]['canon_a']);B=[[0]*8 for _ in range(8)]
 for bit,(u,v) in enumerate(pairs):B[u][v]=(code>>(27-bit))&1;B[v][u]=1-B[u][v]
 need(autos(B)==[tuple(range(8))],'actual source has trivial automorphism; no extra hidden lift capacity')
 P,Q=r['before']['path'],r['after']['path'];fixed={j:P[j] for j in range(8) if P[j]==Q[j]};free=[j for j in range(8) if j not in fixed]
 need(len(free)==r['cost'],'actualchangedclocks')
 start,end,j,c,d=r['pins'];need(len({start,end,c,d})==4,'four distinct declaredsharedroles');need(P[0]==Q[0]==start and P[-1]==Q[-1]==end and P[j:j+2]==Q[j:j+2]==[c,d],'literal sharedpinrecords')
 owner_map={v:({start:'shared:first',end:'shared:last',c:'shared:arc0',d:'shared:arc1'}.get(v,f'private:{s}:{v}')) for v in range(8)}
 need(len(set(owner_map.values()))==8,'ownerframeinjective')
 domain=[]
 for perm in permutations([P[j] for j in free]):
  path=P[:]
  for j,x in zip(free,perm):path[j]=x
  if not all(B[u][v] for u,v in zip(path,path[1:])):continue
  w=sum(B[path[u]][path[v]]<<i for i,(u,v) in enumerate(chords));t=merged[owner[full^w]];O=orbit(w)
  need(merged[owner[w]]==s,'literal source owner');need(len(O)==4,'black normalized resource')
  need(all(path[j]==v for j,v in fixed.items()),'allfixedclocks')
  need(min(O) not in globalresources,'normalizedresourcegloballydistinct across fullpinmenus')
  globalresources[min(O)]=[s,t,path];need(not globalports.intersection(O),'closedportdisjointness');globalports.update(O)
  K=partner(B,path);C=transpose(B);rev=path[::-1]
  need(partner(K,path)==B,'partnerinvolution');need(partner(C,rev)==transpose(K),'actualconversesquare');need(all(C[u][v] for u,v in zip(rev,rev[1:])),'converseHamiltonianpath')
  need(sum(C[rev[u]][rev[v]]<<i for i,(u,v) in enumerate(chords))==reflect(w),'literalreflectedword')
  need(all(owner_map[path[j]]==owner_map[P[j]] for j in fixed),'everyfixedsharedorprivatecoordinate')
  domain.append({'path':path,'word':w,'target':t,'resource':min(O),'closed_ports':sorted(O),'target_label':label(t),'target_in_original_pair':t in(1965,2234),'source_phase_mod8':(X[s]-X[t])%8})
 need(any(d['path']==P and d['word']==r['before']['mask'] for d in domain),'beforepresent');need(any(d['path']==Q and d['word']==r['after']['mask'] for d in domain),'afterpresent')
 selected=[d for d in domain if d['target_in_original_pair']]
 need(all(d['target_label']==label(1965) for d in selected),'selectedpaircode')
 reports.append({'source':s,'matrix':B,'owner_map':owner_map,'shared_pin_record':r['pins'],'fixed_clock_owners':sorted(map(list,fixed.items())),'changed_clocks':free,'complete_pinned_paths':domain,'within_pair_paths':len(selected),'within_pair_target_counts':{str(t):sum(d['target']==t for d in selected) for t in(1965,2234)}})
result={'status':'PASS','checks':checks,'scope':'All common-clock pinned domains of the supplied28 options. Extra targets are excluded only when original U/V target margins remain fixed; no assertion that raw pins alone define binary choices.','counts':{'sources':len(rows),'complete_pinned_paths':sum(len(r['complete_pinned_paths']) for r in reports),'within_pair_paths':sum(r['within_pair_paths'] for r in reports),'normalized_resources':len(globalresources),'closed_ports':len(globalports)},'domains':reports,'pins':{'inherited_supply:'+p.name if p==a.supply else str(p.relative_to(a.repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'counts':result['counts'],'domain_counts':[[r['source'],len(r['complete_pinned_paths']),r['within_pair_paths'],r['within_pair_target_counts']] for r in reports]},indent=2))
