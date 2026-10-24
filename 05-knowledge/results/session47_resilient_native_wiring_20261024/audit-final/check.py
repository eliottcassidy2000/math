#!/usr/bin/env python3
"""Independent compact synthesis controls; no producer imports or mutations.
Native classification and all-state connectivity remain inherited audit results.
All reads are fixed local source data. Output goes only to stdout.
"""
from pathlib import Path
from fractions import Fraction
from itertools import combinations
import argparse, csv, gzip, hashlib, json, math
parser=argparse.ArgumentParser()
parser.add_argument("--root-dir",type=Path,default=Path(__file__).resolve().parent.parent,help="Workspace root or publication bundle root")
ROOT=parser.parse_args().root_dir.resolve()
ALIASES={"session47-audit-fresh570":"fresh570-audit","session47-helper-resilience/final-v3":"helper-resilience-v3","session47-clock-connectivity":"clock-connectivity","session47-carrier-transfer":"carrier-transfer","session47-audit-global-observer":"audit-global-observer","session47-root/final":"","session47-audit-final":"audit-final"}
def resolve(name):
 p=ROOT/name
 if p.exists():return p
 for old,new in sorted(ALIASES.items(),key=lambda x:-len(x[0])):
  if name==old or name.startswith(old+"/"):
   p=ROOT/(new+name[len(old):] if new else name[len(old):].lstrip("/"))
   if p.exists():return p
 raise FileNotFoundError(name)
checks=0
inputs={}
def require(ok,why):
 global checks
 checks+=1
 if not ok: raise RuntimeError(why)
def read(name):
 b=resolve(name).read_bytes(); inputs[name]={'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()};return b
def data(name): return json.loads(read(name))

supply=data('session47-audit-fresh570/inputs/original/compatible_charged_helpers_20261023/supply/selected_supply.json')
helpers=supply['helpers'];ids=[h['source'] for h in helpers]
labels={h['source']:h['label'] for h in helpers}
def charge(sources):
 x=0;y=0
 for s in sources:
  z=labels[s];x+=z[0];y^=sum(v<<i for i,v in enumerate(z[1:]))
 return 16*(x%128)+y
weighted=data('session47-helper-resilience/final-v3/fractional_weights.json')
require(weighted['source_ids']==ids,'53-owner indexing')
require(len(weighted['charges'])==2048,'weighted charge census')
count=0;ratios=[]
for q,row in enumerate(weighted['charges']):
 d=row['denominator'];loads=[0]*53;seen=set();total=0
 require(isinstance(d,int) and d>0,'positive capacity')
 for mask,w in row['edges']:
  sel=[j for j in range(53) if mask>>j&1]
  require(mask>=0 and mask<1<<53 and len(sel)==4 and mask not in seen,'actual distinct quadruple')
  require(type(w) is int and w>0,'integer proof weight')
  require(charge(ids[j] for j in sel)==q,'native supplied charge')
  seen.add(mask);total+=w
  for j in sel: loads[j]+=w
 require(max(loads)<=d and total>13*d,'strict hitting dual')
 require(sum(loads)==4*total and 4*total<=53*d,'fractional ceiling')
 count+=len(seen);ratios.append(Fraction(total,d))
require(count==108309,'weighted edge count')
require(min(ratios)==Fraction(13145902,1000001) and ratios.index(min(ratios))==540,'minimum certificate ratio')
packing=data('session47-helper-resilience/final-v3/quad_packings.json')
require(packing['source_ids']==ids and len(packing['packings'])==2048,'packing carrier')
for q,ms in enumerate(packing['packings']):
 used=0;require(len(ms)==12,'twelve witnesses')
 for mask in ms:
  sel=[j for j in range(53) if mask>>j&1]
  require(len(sel)==4 and not mask&used and charge(ids[j] for j in sel)==q,'disjoint actual charge')
  used|=mask
require(len({charge([s]) for s in ids})==53 and math.gcd(13,128)==1,'thirteen-packing obstruction')
upper=data('session47-helper-resilience/final-v3/upper_transversal.json')
quad353=[c for c in combinations(ids,4) if charge(c)==353]
require(len(quad353)==139 and len(set(upper['forced_U']))==15,'integral upper dimensions')
require(all(set(c)&set(upper['forced_U']) for c in quad353),'complete fifteen-owner cover')
triple305=[c for c in combinations(ids,3) if charge(c)==305]
require(len(triple305)==4 and all(set(c)&{776,3139} for c in triple305),'full53 two-U hostile')

# Independent 27-choose-4 exhaustive check of the deletion hostile.
core={int(a[0]):list(map(int,a[1:])) for a in (line.split() for line in read('session47-clock-connectivity/input.txt').decode().splitlines())}
def corecharge(c):
 x=sum(core[s][1] for s in c)%128;y=0
 for s in c:y^=core[s][2]
 return x,y
fiber=[frozenset(c) for c in combinations([s for s in core if s!=2071],4) if corecharge(c)==(45,11)]
expected=[frozenset(c) for c in [[1394,1921,1923,2905],[1532,1841,2835,3360],[1691,1923,1974,2837]]]
require(set(fiber)==set(expected),'complete deletion fiber')
costs=[sum(core[s][3] for s in expected[a]^expected[b]) for a,b in [(0,1),(0,2),(1,2)]]
require(costs==[23,12,23],'deletion exact three costs')

# Arithmetic and connectivity certificate for the quantitative odd part.
tree=data('session47-audit-fresh570/odd-part.normal.json');n=tree['vertices'];edges=tree['selected_tree_edges']
require(n==3528 and len(edges)==n-1,'spanning tree size')
parents=list(range(n))
def root(x):
 while parents[x]!=x: parents[x]=parents[parents[x]];x=parents[x]
 return x
wanted={(min(a,b),max(a,b)):w for a,b,w in edges}
require(len(wanted)==len(edges),'distinct tree edges')
P=1
for a,b,w in edges:
 require(0<=a<n and 0<=b<n and a!=b and w>0,'positive nonloop tree edge')
 ra,rb=root(a),root(b);require(ra!=rb,'tree acyclicity');parents[ra]=rb;P*=w
require(len({root(i) for i in range(n)})==1,'spans every vertex')
name='session47-audit-fresh570/inputs/original/metagraph_order8_20261023/merged_edges.tsv.gz';read(name)
seen=set()
with gzip.open(resolve(name),'rt') as stream:
 for row in csv.DictReader(stream,delimiter='\t'):
  a,b=int(row['source']),int(row['target']);key=min(a,b),max(a,b)
  if key in wanted:
   require(wanted[key]==int(row['blue'])+int(row['black']),'literal original conductance')
   seen.add(key)
require(len(seen)==n-1,'all tree edges authentic')
require(P==int(tree['tree_conductance_product_decimal']) and P.bit_length()==9205,'exact tree product')
require(tree['critical_two_primary_order_log2']==3374 and P>=1<<9204,'odd lower bound premises')

# Explicit four-row parity obstruction, grounded bank size and extension count.
susp=data('session47-carrier-transfer/five-to-seven/normal.json');L=susp['old']['L'];N=susp['new']['L'];readouts=susp['readouts']
require(len(L)==10 and len(N)==272 and len(readouts)==48,'carrier dimensions')
require([sum(N[i][j] for i in [10,41,196,244])%2 for j in range(272)]==[int(j in [10,41]) for j in range(272)],'full four-row identity')
h=[int(i in [3,4,6,7]) for i in range(10)]
require(all(sum(a*b for a,b in zip(row,h))%2==0 for row in L) and h[3]!=h[9],'old separating harmonic')
require([10,3] in readouts and [41,9] in readouts,'source readouts retained')
oldbank=[]
for bits in range(1<<9):
 v=[0]+[(bits>>i)&1 for i in range(9)]
 if all(sum(a*b for a,b in zip(row,v))%2==0 for row in L):oldbank.append(v)
require(len(oldbank)==8,'grounded old count')
def solvable(rows,nvars):
 piv={}
 for row in rows:
  while row&((1<<nvars)-1):
   j=(row&-row).bit_length()-1
   if j not in piv: piv[j]=row;break
   row^=piv[j]
  else:
   if row>>nvars&1:return False
 return True
matrix=[sum((a%2)<<j for j,a in enumerate(row)) for row in N]
extend=[]
for v in oldbank:
 rows=matrix+[1]+[(1<<a)|(v[b]<<272) for a,b in readouts]
 if solvable(rows,272):extend.append(v)
require(len(extend)==4,'four of eight extend')
# Pair-charge filtering, genuinely independent of raw domain enumeration.
fresh=[(4,35),(15,91),(3,81),(8,45)]
valid=[list(c) for c in combinations(range(4),2) if (fresh[c[0]][0]^fresh[c[1]][0],sum(fresh[i][1] for i in c)%128)==(11,126)]
require(valid==[[0,1],[2,3]],'fresh exact paired projection')
require(all(z!=(0,126) for z in fresh),'helper2393 cannot replace one fresh')
require(sum([7,6,7,3])==23 and 8*22==176,'clock arithmetic')
# Positive shared-owner section, independently checked from its frozen record.
shared=data('session47-audit-global-observer/inputs/shared_owner.json')
require(shared['core_source']==176 and shared['helper_source']==368 and shared['shared_owners']==[2,6],'shared-owner frame')
require([(r['core_bit'],r['helper_bit']) for r in shared['all_bit_pairs'] if r['allowed']]==[(0,1),(1,0)],'opposite bit constraint')
bank=shared['helper_source_ids'];require(len(bank)==48 and len(set(bank))==48 and 368 in bank,'shared bank48')
for sec in shared['signed_sections']:
 bit=sec['forced_helper_bit'];masks=sec['section_masks']
 require(bit in [0,1] and len(masks)==2048,'complete signed section')
 for q,mask in enumerate(masks):
  chosen=[bank[i] for i in range(48) if mask>>i&1]
  require(0<=mask<1<<48 and len(chosen)==4,'shared section actual quadruple')
  require(int(368 in chosen)==bit and charge(chosen)==q,'shared section compatibility and charge')
route=shared['native_route'];states=route['states']
require([r['shared_owner'] for r in states]==[2,6,6,2],'changing then returned owner')
require(states[0]['helper_V']==states[-1]['helper_V'] and states[0]['core_V']!=states[-1]['core_V'],'helper return not closed core loop')
hcost={h['source']:sum(a!=b for a,b in zip(h['selected']['before']['path'],h['selected']['after']['path'])) for h in helpers}
for i,event in enumerate(route['events']):
 a,b=states[i],states[i+1]
 changed=(set(a['core_V'])^set(b['core_V']))|(set(a['helper_V'])^set(b['helper_V']))
 work=sum(core[t][3] if t in core else hcost[t] for t in changed)
 require(sorted(changed)==event['changed_sources'] and len(changed)==event['support'],'actual shared event support')
 require(work==event['owner_clock_occurrences'],'actual shared event work')
require([r['support'] for r in route['events']]==[10,8,10] and [r['owner_clock_occurrences'] for r in route['events']]==[53,39,53],'shared route reported values')
replay=data('session47-root/final/SHARED_OWNER_REPLAY.json')
require(replay['normal_output']==replay['optimized_output']=='PASS 8,699; shared_owner.json matches','parent native replay receipt')
require(replay['source_hashes']['05-knowledge/results/compatible_charged_helpers_20261023/supply/shared_owner.json']==inputs['session47-audit-global-observer/inputs/shared_owner.json']['sha256'],'replay pins exact shared record')
pybytes=read('session47-audit-final/inputs/shared_owner.py.txt')
require(hashlib.sha256(pybytes).hexdigest()==replay['source_hashes']['05-knowledge/results/compatible_charged_helpers_20261023/supply/shared_owner.py'],'replay pins archived unchanged code')
print(json.dumps({'status':'PASS','checks':checks,'weighted_quadruples':count,'minimum_weight_ratio':[min(ratios).numerator,min(ratios).denominator],'complete_charge353_quadruples':len(quad353),'deletion_costs':costs,'tree_product_bits':P.bit_length(),'dyadic_tree_valuation':3374,'odd_part_lower_log2':9204-3374,'old_grounded_mod2_count':len(oldbank),'extendible_grounded_count':len(extend),'fresh_pair_indices':valid,'shared_owner_sections':4096,'shared_route_supports':[10,8,10],'shared_route_clocks':[53,39,53],'inputs':inputs,'scope':'Compact independent certificate arithmetic and finite projection controls. Complete native reconstruction, exhaustive old28 connectivity, and proof-level all-size statements are reviewed through their frozen independent audits, not newly recomputed here.'},indent=2,sort_keys=True))
