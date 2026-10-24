#!/usr/bin/env python3
"""Independent native common-clock audit. No producer imports; -O-safe checks."""
import argparse,csv,gzip,hashlib,itertools,json,struct
from collections import Counter,defaultdict
from pathlib import Path

P=argparse.ArgumentParser();P.add_argument('--repo',type=Path,required=True);P.add_argument('--supply',type=Path,required=True);P.add_argument('--out',type=Path,required=True);A=P.parse_args()
checks=0
def require(ok,message):
 global checks
 checks+=1
 if not ok:raise RuntimeError(message)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
R=A.repo/'05-knowledge/results';G=R/'metagraph_order8_20261023'
inputs={'nodes':G/'merged_nodes.tsv','owners':G/'owner_map.u16.gz','edges':G/'merged_edges.tsv.gz','electrical':R/'metagraph_order8_electrical_20261023/receipt.json','supply':A.supply}
hashes={k:sha(v)for k,v in inputs.items()}
nodes=[{k:int(v)for k,v in r.items()}for r in csv.DictReader(inputs['nodes'].open(),delimiter='\t')]
owners=struct.unpack('<'+'H'*(1<<21),gzip.decompress(inputs['owners'].read_bytes()))
merge={r[k]:r['merged_id']for r in nodes for k in ('class_a','class_b')}
edges=[tuple(int(r[k])for k in ('source','target','blue','black'))for r in csv.DictReader(gzip.open(inputs['edges'],'rt'),delimiter='\t')]
edge_colors={(i,j):(b,k)for i,j,b,k in edges}
E=json.loads(inputs['electrical'].read_text())['order8'];Y=E['radical_lifts_mod4'];X=E['higher_layers'][-1]['lifts'][0]
records=json.loads(A.supply.read_text())['records'];U,V=1965,2234
require(len(nodes)==3528 and len(edges)==438216,'original full graph dimensions')
require(len(records)==len({r['source']for r in records})==28,'exact distinct supply')
require(X[U]==X[V]==146,'actual equal primitive target potential')
require(all((y[U]-y[V])%2==0 for y in Y),'all target paired parity phases')
pairs=tuple(itertools.combinations(range(8),2));chords=tuple((i,j)for i,j in pairs if j-i>1);chord_id={e:k for k,e in enumerate(chords)};full=(1<<21)-1
all_permutations=list(itertools.permutations(range(8)))
def matrix(code):
 M=[[0]*8 for _ in range(8)]
 for bit,(i,j)in zip(range(27,-1,-1),pairs):M[i][j]=(code>>bit)&1;M[j][i]=1-M[i][j]
 return M
def transpose(M):return [list(z)for z in zip(*M)]
def ham(M,p):return all(M[i][j]for i,j in zip(p,p[1:]))
def word(M,p):return sum(M[p[i]][p[j]]<<k for k,(i,j)in enumerate(chords))
def rev_word(w):return sum(((w>>chord_id[7-j,7-i])&1)<<k for k,(i,j)in enumerate(chords))
def tile(w):return {w,w^full,rev_word(w),rev_word(w)^full}
def switch(M,p):
 keep={tuple(sorted(e))for e in zip(p,p[1:])}
 K=[[0]*8 for _ in range(8)]
 for i,j in pairs:K[i][j]=M[i][j]if(i,j)in keep else M[j][i];K[j][i]=1-K[i][j]
 return K
def code(s):return tuple(y[s]%2 for y in Y[:4])+(X[s]%128,)
def lap(v):
 out=[0]*len(nodes)
 for i,j,b,k in edges:
  z=(b+k)*(v[i]-v[j]);out[i]+=z;out[j]-=z
 return out
for y in Y:require(all(v%4==0 for v in lap(y)),'full original mod-four harmonic')
require(all(v%256==0 for v in lap(X)),'full original primitive mod-256 harmonic')
chars=[[64*y for y in yy]for yy in Y[:4]]+[X]
for chi in chars:require(all(v%256==0 for v in lap(chi)),'actual full graph separator')
reports=[];resources={};all_ports=set();labels={};menus={}
for r in records:
 s=r['source'];M=matrix(nodes[s]['canon_a']);degrees=tuple(map(sum,M));actual=[];automorphisms=[]
 for p in all_permutations:
  if ham(M,p):actual.append(p)
  if all(degrees[i]==degrees[p[i]]for i in range(8))and all(M[i][j]==M[p[i]][p[j]]for i,j in pairs):automorphisms.append(p)
 require(automorphisms==[tuple(range(8))],f'full permutation automorphism census {s}')
 require(nodes[s]['class_a']!=nodes[s]['class_b'],'distinct ordinary/converse source classes')
 require(2*len(actual)==nodes[s]['fibre'],'complete actual source fibre')
 capacities=Counter(merge[owners[word(M,p)^full]]for p in actual)
 require(all(capacities[t]==r['complete_target_capacities'][str(t)]for t in (U,V)),'unrestricted native target capacities')
 before,after=tuple(r['before']['path']),tuple(r['after']['path']);fixed=[j for j in range(8)if before[j]==after[j]]
 require(8-len(fixed)==r['cost'],'true changed-clock cost')
 start,end,j,c,d=r['pins'];require(1<=j<=5 and len({start,end,c,d})==4,'four distinct shared roles')
 require(before[0]==after[0]==start and before[-1]==after[-1]==end and before[j:j+2]==after[j:j+2]==(c,d),'actual indexed endpoint and arc pins')
 shared={start:'shared:first',end:'shared:last',c:'shared:arc0',d:'shared:arc1'};frame={i:shared.get(i,f'private:{s}:{i}')for i in range(8)}
 require(len(set(frame.values()))==8,'actual injective occurrence-owner frame')
 complete=[p for p in actual if all(p[j]==before[j]for j in fixed)];domain=[]
 for p in complete:
  w=word(M,p);target=merge[owners[w^full]];T=tile(w);K=switch(M,p);C=transpose(M);KC=transpose(K);rev=p[::-1]
  require(merge[owners[w]]==s,'literal source owner from original atlas')
  require(len(T)==4,'pure black normalized resource')
  require(all(ham(N,z)for N,z in ((M,p),(K,p),(C,rev),(KC,rev))),'four literal Hamiltonian presentations')
  require(switch(K,p)==M and switch(C,rev)==KC,'literal kernel involution and converse square')
  require([word(N,z)for N,z in ((M,p),(K,p),(C,rev),(KC,rev))]==[w,w^full,rev_word(w),rev_word(w)^full],'four actual words equal normalized orbit')
  require(all(frame[p[j]]==frame[before[j]]for j in fixed),'all fixed physical owner clocks, not just four common roles')
  blue,black=edge_colors[tuple(sorted((s,target)))];require(blue==0 and black>=2,'authentic pure-black colored edge')
  require(min(T)not in resources and not all_ports.intersection(T),'global normalized resource and closed-port injection')
  resources[min(T)]=(s,target);all_ports.update(T)
  domain.append({'path':list(p),'word':w,'target':target,'normalized_resource':min(T),'closed_ports':sorted(T),'label':list(code(target)),'chosen_X_phase_mod8':(X[s]-X[target])%8})
 for key in ('before','after'):
  require(any(z['path']==r[key]['path']and z['word']==r[key]['mask']and z['target']==r[key]['target']for z in domain),'actual inventory endpoint recovered without trusting mask')
 selected=[z for z in domain if z['target']in(U,V)]
 require(Counter(z['target']for z in selected)==Counter({U:1,V:1}),'complete common-clock menu has exactly one U and V alternative')
 require(all(z['label']==list(code(U))for z in selected),'complete paired target code retained')
 require(list(code(s))==r['source_label'],'input label independently reconstructed')
 labels[s]=code(s);menus[s]={z['target']:z for z in selected}
 reports.append({'source':s,'actual_hamiltonian_paths':len(actual),'tested_owner_permutations':len(all_permutations),'automorphisms':len(automorphisms),'unrestricted_target_capacities':{str(t):capacities[t]for t in(U,V)},'fixed_clocks':fixed,'fixed_owner_values':[before[j]for j in fixed],'owner_frame':frame,'matrix':M,'complete_pinned_domain':domain})
extras=[{'source':r['source'],**p}for r in reports for p in r['complete_pinned_domain']if p['target']not in(U,V)]
require(sum(len(r['complete_pinned_domain'])for r in reports)==57 and len(resources)==57 and len(all_ports)==228,'complete original raw pin domain sizes')
require(len(extras)==1 and extras[0]['source']==1923 and extras[0]['target']==1559,'unique extra pinned target, excluded by zero target margin')
# Four-slot zero-sum trades, independent pair-bucket enumeration.
labels2={s:(sum(q[j]<<j for j in range(4)),q[4])for s,q in labels.items()};buckets=defaultdict(list)
for a,b in itertools.combinations(sorted(labels),2):buckets[(labels2[a][0]^labels2[b][0],(labels2[a][1]+labels2[b][1])%128)].append((a,b))
trades=[]
for pairlist in buckets.values():
 for a,b in itertools.combinations(pairlist,2):
  if set(a).isdisjoint(b):trades.append((a,b))
trades=sorted(trades)
require(len(set(labels.values()))==28,'no two-slot balanced equal-charge moves')
require(len(trades)==27,'entire 28-source supply has exactly 27 four-slot balanced charge trades')
require(all(sum(s in(1638,1974)for s in a)==sum(s in(1638,1974)for s in b)for a,b in trades),'every four-slot trade preserves b1638+b1974')
require([t for t in trades if any(s in(1638,1974)for side in t for s in side)]==[((1638,1841),(1974,2176)),((1638,2176),(1974,2945))],'exact four-trade boundary touching critical pair')
out={'status':'PASS','checks':checks,'scope':'Complete common-clock domains reconstructed from all 40320 owner permutations per source; full 28 supplied source menus, not a full-atlas replacement. All-dyadic completeness uses independently audited Session45 dyadic filtration. Shared physical owner records do not assert common-host induced-arc agreement.','input_sha256':hashes,'counts':{'sources':28,'owner_permutations_tested':28*40320,'raw_pinned_paths':57,'UV_pinned_paths':56,'normalized_resources':57,'closed_ports':228},'sources':reports,'extra_pinned_paths':extras,'four_slot_trades':[[list(a),list(b)]for a,b in trades],'four_move_obstruction':'All zero-through-four-slot moves preserve b1638+b1974. Original endpoint values are 0 and 2. Target margins exclude the sole target1559 path, so this also proves a universal six-slot lower bound in the complete 28-source common-clock universe.'}
A.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'counts':out['counts'],'four_slot_trades':len(trades),'output':str(A.out)},sort_keys=True))
