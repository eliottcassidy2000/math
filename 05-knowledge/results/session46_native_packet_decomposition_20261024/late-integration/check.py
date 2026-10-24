"""Bounded independent arithmetic comparison of late concurrent repo work."""
import argparse,collections,hashlib,itertools,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--supply',type=Path,required=True);p.add_argument('--domains',type=Path,required=True);p.add_argument('--incoming-dir',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();H=Path(__file__).resolve().parent;checks=0

def need(c,m):
 global checks
 checks+=1
 if not c:raise RuntimeError(m)
P=json.loads((H/'SOURCE_PINS.json').read_text())
for path,sha in P['files'].items():
 q=(a.incoming_dir/Path(path).name) if a.incoming_dir else a.repo/path
 need(hashlib.sha256(q.read_bytes()).hexdigest()==sha,'incoming source pin')
for path,sha in P.get('inherited_files',{}).items():need(hashlib.sha256((a.repo/path).read_bytes()).hexdigest()==sha,'inherited finite census pin')
raw=(a.repo/'05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json').read_bytes();E=json.loads(raw)['order8'];X=E['higher_layers'][-1]['lifts'][0];Y=E['radical_lifts_mod4'];N=len(X)
def lab(s):return tuple(Y[j][s]%2 for j in range(4))+(X[s]%128,)
targets=[s for s in range(N) if lab(s)==lab(1965)]
need(targets==[570,1965,2234],'complete target charge fiber')
need(lab(1965)==(1,1,0,1,18),'original target label')
dom=json.loads(a.domains.read_text());need(dom['counts']['complete_pinned_paths']==57,'complete original common-clock supply')
actual=[(d['source'],r['target']) for d in dom['domains'] for r in d['complete_pinned_paths']]
need(len(actual)==57 and sum(t in (1965,2234) for s,t in actual)==56,'two-target domain count')
need(not any(t==570 for s,t in actual),'arithmetic third target has no pinned path')
need([(s,t) for s,t in actual if t not in (1965,2234)]==[(1923,1559)],'only extra actual target')
rows=json.loads(a.supply.read_text())['records'];options={r['source']:r for r in rows};ids=sorted(options)
prefix=a.repo/'05-knowledge/results/session45_injective_multisource_moves_20261023/source-circuits'
menus=json.loads((prefix/'full-search/all_menus.json').read_text());prior=json.loads((prefix/'receipt.normal.json').read_text())['universe']
need(len(menus)==700,'complete inherited normalized black menu census')
need(not any(570 in r['targets'] for r in menus),'no target570 transition in this inherited grammar')
need({r['source'] for r in menus if r['targets']==[1965,2234]}==set(ids),'same28source pair supply')
need((prior['canonical_sources'],prior['actual_source_paths'],prior['black_endpoint_arc_menus'])==(3528,1102392,700),'inherited census scope not new enumeration')
A=[176,1394,1457,1477,1532,1923,2087,2176,3369];B=[1425,1549,1841,1921,2185,2553,2835,2905,3000]
need(len(A)==len(B)==9 and set(A).isdisjoint(B),'distinct incoming18')
need(sum(options[s]['cost'] for s in A+B)==45,'incoming clock cost45')
sources=A+B;sign=[1]*9+[-1]*9;size=1<<18;count=[0]*size;mass=[0]*size;parity=[0]*size;labels=[sum((Y[j][s]%2)<<j for j in range(4)) for s in sources]
for bits in range(1,size):
 low=bits&-bits;j=low.bit_length()-1;old=bits^low;count[bits]=count[old]+sign[j];mass[bits]=mass[old]+sign[j]*X[sources[j]];parity[bits]=parity[old]^labels[j]
balanced=[x for x in range(size) if count[x]==0];coarse=[x for x in balanced if parity[x]==0 and mass[x]%4==0];fine=[x for x in balanced if parity[x]==0 and mass[x]%128==0]
need(len(balanced)==48620 and len(coarse)==696 and fine==[0,size-1],'incoming complete guard pattern')
need(2*mass[-1]==-256,'incoming raw integer debt')
# Independently verify all four-move invariants, keeping the full invariant space distinction.
c={1691:-3,3000:-3,3369:-3,1841:-2,2176:-2,2581:-2,2945:-2,1480:-1,2553:-1,3360:-1,176:1,1532:1,2837:1,1477:2,1425:3,1638:3,1974:3,3015:3}
pairs=collections.defaultdict(list)
for u,v in itertools.combinations(ids,2):pairs[tuple((lab(u)[j]+lab(v)[j])%2 for j in range(4))+((X[u]+X[v])%128,)].append((u,v))
trades=[]
for pp in pairs.values():
 for aa,bb in itertools.combinations(pp,2):
  if set(aa).isdisjoint(bb):
   need(sum(c.get(s,0) for s in aa)==sum(c.get(s,0) for s in bb),'new exact invariant on every four-trade');trades.append([aa,bb])
need(len(trades)==27,'complete four-trade inventory')
critical=lambda ss:sum(s in(1638,1974) for s in ss)
need(critical(A)==critical(B)==0,'our small invariant is blind to incoming endpoints')
need((sum(c.get(s,0) for s in B),sum(c.get(s,0) for s in A))==(-3,-1),'new invariant separates')
# Exact group-coordinate reconstruction for zero-sum target-charge tuples.
gamma_controls=[]
for modulus in (2,3,4):
 for t in range(2,6):
  tuples=[x for x in itertools.product(range(modulus),repeat=t) if sum(x)%modulus==0]
  prefixes={x[:-1] for x in tuples};need(len(tuples)==len(prefixes)==modulus**(t-1),'target-charge kernel is product on t-minus-one coordinates')
  need(all(x[-1]==(-sum(x[:-1]))%modulus for x in tuples),'last physical target reconstructed')
  gamma_controls.append([modulus,t,len(tuples)])
need([t*(1+(t-1)*131) for t in (2,3,4)]==[264,789,1576],'p-group cycle bounds')
result={'status':'PASS','checks':checks,'incoming_commit':P['commit'],'complete_equal_charge_target_fiber':targets,'missing_native_target':570,'inherited_transition_census':{'canonical_merged_sources':3528,'actual_paths':1102392,'menus':700,'menus_involving570':0,'scope':'Inherited audited black nonloop equal-charge endpoint-and-indexed-internal-arc grammar only; no full census rerun.'},'incoming18':{'plus':A,'minus':B,'cost':45,'balanced':len(balanced),'H8':len(coarse),'H256':len(fine),'integer_X_debt':-256,'small_invariant_endpoint_values':[0,0],'new_invariant_endpoint_values':[-3,-1]},'gamma_controls':gamma_controls,'pins':{'original_receipt':hashlib.sha256(raw).hexdigest(),'supply':hashlib.sha256(a.supply.read_bytes()).hexdigest(),'complete_domains':hashlib.sha256(a.domains.read_bytes()).hexdigest()},'scope':'Bounded incoming comparison and original-label guard reconstruction. Not a fresh literal native-path/image audit of the incoming45-clock circuit or a replay of all incoming routes, transport or spectral results.'}
a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2))
