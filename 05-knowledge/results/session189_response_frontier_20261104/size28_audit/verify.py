#!/usr/bin/env python3
"""Independent size-28 replay. Explicit checks survive python -O."""
from pathlib import Path
import collections,csv,gzip,hashlib,itertools,json,random,subprocess,sys
from prepare_physical import build,compress,TYPES
D=Path(__file__).resolve().parent;P=D/'input'
def need(v,m):
 if not v:raise RuntimeError(m)
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def tab(p):return [{k:int(v) for k,v in r.items()} for r in csv.DictReader(p.open(),delimiter='\t')]
def path_count(rows):
 return sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in itertools.permutations(range(len(rows))))
def profile(rows):
 n=len(rows);paths=[0]*(1<<n)
 for mask in range(1,1<<n):
  nodes=[i for i in range(n) if mask>>i&1]
  paths[mask]=sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in itertools.permutations(nodes))
 prev=[0]*(1<<n);prev[0]=1;answer=[]
 for h in range(1,n+1):
  nxt=[0]*(1<<n)
  for mask in range(1,1<<n):
   part=mask
   while part:
    nxt[mask]+=paths[part]*prev[mask^part];part=(part-1)&mask
  answer.append(nxt[-1]);prev=nxt
 return answer

def local_checks():
 records=load(P/'four_child_records.json');need(len(records)==64,'all labeled four sources')
 hist=collections.Counter();families=0;physical_types={}
 for code,r in enumerate(records):
  rows=[0]*4;bit=0
  for i in range(4):
   for j in range(i+1,4):
    if code>>bit&1:rows[i]|=1<<j
    else:rows[j]|=1<<i
    bit+=1
  need(r['code']==code and r['rows']==rows,'four source orientation')
  q=profile(rows);need(q==r['profile'],'independent subset convolution cover profile')
  mods=[]
  for mask in range(1,16):
   if mask.bit_count()>1 and all(len({rows[v]>>w&1 for w in range(4) if mask>>w&1})==1 for v in range(4) if not(mask>>v&1)):mods.append(mask)
  need(mods==r['modules'],'all original four-child modules')
  need(len(r['targets'])==1<<len(mods),'all original four-child families')
  for f,t in enumerate(r['targets']):
   selected=[m for i,m in enumerate(mods) if f>>i&1];target=rows[:]
   for i in range(4):
    for j in range(i+1,4):
     if sum((m>>i&1)*(m>>j&1) for m in selected)%2:target[i]^=1<<j;target[j]^=1<<i
   tq=profile(target);need(t['family']==selected and t['rows']==target and t['profile']==tq,'original family target')
   if q[0]>1:need(tq==q,'nontransitive family profile invariance');families+=1
  if q[0]>1:hist[tuple(q)]+=1
 need(hist=={(3,18,36,24):16,(5,22,36,24):24} and families==160,'four profile classification')
 # Explicit isomorphisms, not merely equal scalar profiles, cover all 40 sources.
 for typ,(rows,_) in TYPES.items():
  keys=set()
  for p in itertools.permutations(range(4)):
   permuted=[0]*4
   for i in range(4):
    for j in range(4):
     if rows[i]>>j&1:permuted[p[i]]|=1<<p[j]
   keys.add(tuple(permuted))
  physical_types[typ]=keys
 need([len(v) for v in physical_types.values()]==[8,8,24],'four labeled orbit sizes')
 need(len(set().union(*physical_types.values()))==40,'four type disjointness')
 need(all(tuple(r['rows']) in set().union(*physical_types.values()) for r in records if r['nontransitive']),'nontransitive orbit coverage')
 for typ,(rows,factor) in TYPES.items():
  q=profile(rows);need([x*pow(factor,-1,8)%8 for x in q]==[1,6,4,0],'normalized complete profile')
 named=load(P/'named_controls.json');cubic=quartic=0;exterior_controls=0
 for ctl in named:
  source=TYPES[ctl['type']][0]
  for label,rows in [('source',source),('converse',[15^r^(1<<i) for i,r in enumerate(source)])]:
   vals=[]
   for m in range(256):
    a=m&15;c=m>>4;aug=[r|(((a>>i)&1)<<4)|(((c>>i)&1)<<5) for i,r in enumerate(rows)]+[(15^a)|32,15^c]
    base=path_count(aug);vals.append(base)
    if m%17==0:
     for dominates in (False,True):
      extended=[r|(0 if dominates else 1<<6) for r in aug]+[(1<<6)-1 if dominates else 0]
      need(path_count(extended)==base,'exact uniform singleton exterior deletion bijection');exterior_controls+=1
   coefficients=vals[:]
   for bit in range(8):
    for mask in range(256):
     if mask>>bit&1:coefficients[mask]-=coefficients[mask^(1<<bit)]
   need(ctl[label]['rows']==rows and ctl[label]['values']==vals and ctl[label]['coefficients']==coefficients,'literal complete named values and coefficients')
   need(all(not v for m,v in enumerate(coefficients) if m.bit_count()>4),'local named degree bound')
   cubic+=sum(v!=0 and m.bit_count()==3 for m,v in enumerate(coefficients));quartic+=sum(v!=0 and m.bit_count()==4 for m,v in enumerate(coefficients))
 return dict(nontransitive_four_sources=40,family_outputs=families,nonzero_cubics=cubic,nonzero_quartics=quartic,exact_singleton_exterior_controls=exterior_controls)

def main():
 pins=load(D/'INPUT_PINS.json')
 for name,meta in pins.items():need(sha(P/name)==meta['sha256'] and (P/name).stat().st_size==meta['bytes'],'frozen input '+name)
 candidates=load(P/'candidates.json')['candidates'];need(len(candidates)==2348,'complete inherited residual domain')
 by={x['node']:x['rows'] for x in candidates};need(len(by)==2348,'unique inherited nodes')
 witness=tab(P/'witnesses.tsv');need(len(witness)==21132,'all nine placements')
 need([(x['node'],x['position']) for x in witness]==[(c['node'],j) for c in candidates for j in range(9)],'exact complete placement order')
 battery=load(P/'battery.json');queries=[(q['a'],q['c']) for q in battery]
 need(queries==[(0,1),(1,0),(1,2),(1,3),(3,1),(3,3),(3,5),(0,7),(0,18),(40,256),(4,3),(130,128),(65,9),(3,32)],'fixed battery')
 need(len(set(queries+[(1<<i,0) for i in range(9)]))==22,'endpoint union cardinality')
 chosen=sorted(set((w['node'],w['query_index'],w['a'],w['c']) for w in witness))
 expected_input=str(len(chosen))+'\n'+''.join(' '.join(map(str,[n,q,a,c]+by[n]))+'\n' for n,q,a,c in chosen)
 need((D/'all_recount.in').read_text()==expected_input,'independent replay domain')
 if '--recount' in sys.argv:
  for src in ['split_words','physical_subsets']:
   subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra','-Werror',str(D/(src+'.cpp')),'-o',str(D/src)],check=True)
  for program,stem in [('split_words','all_recount'),('physical_subsets','large_physical'),('physical_subsets','small_physical')]:
   with (D/(stem+'.in')).open('rb') as inp,(D/(stem+'.out')).open('wb') as out,(D/(stem+'.log')).open('wb') as err:subprocess.run([str(D/program)],stdin=inp,stdout=out,stderr=err,check=True)
 output={};rows_count=0
 for line in (D/'all_recount.out').read_text().splitlines():
  node,query,position,s,t=map(int,line.split());key=node,query,position;need(key not in output,'no duplicate independent recount')
  need(0<=s<8 and 0<=t<8,'independent residue range');output[key]=(s,t);rows_count+=1
 need(rows_count==6283*9,'complete independent recount output')
 need(set(output)=={(n,q,j) for n,q,_,_ in chosen for j in range(9)},'independent output exact domain')
 initial=tab(P/'residues.tsv');extra=tab(P/'selected_residues.tsv');bank={}
 for row in initial:
  for qi in range(12):bank[row['node'],qi,row['position']]=(row['S'+str(qi)],row['T'+str(qi)])
 for row in extra:
  for qi in range(2):bank[row['node'],12+qi,row['position']]=(row['S'+str(qi)],row['T'+str(qi)])
 for key,value in output.items():need(value==bank[key],'every independently recounted retained bank entry')
 hist=collections.Counter();miss=[]
 for row in initial:
  if not any(row['S'+str(q)]!=row['T'+str(q)] for q in range(12)):miss.append(dict(node=row['node'],position=row['position']))
 need(miss==load(P/'initial_misses.json') and len(miss)==13,'exact original battery misses')
 for w in witness:
  n=w['node'];j=w['position'];q=w['query_index'];a,c=queries[q]
  need((a,c)==(w['a'],w['c']),'query index and same named mask')
  need(output[n,q,j]==(w['normalized_source'],w['normalized_target']),'all21132 independent placement witnesses')
  need((w['normalized_source']-w['normalized_target'])%8==4,'every independent gap four')
  earlier=[qq for qq in range(q) if (n,qq,j) in bank and bank[n,qq,j][0]!=bank[n,qq,j][1]];need(not earlier,'first selected witness')
  for bitname,m in [('physical_a',a),('physical_c',c)]:
   mask=sum(1<<v for i in range(9) if m>>i&1 for v in range(3*i,3*i+3))
   if m>>j&1:mask|=1<<27
   need(w[bitname]==mask,'physical mask expansion')
  hist[q]+=1
 # Check every retained physical graph, using separately written direct construction.
 count=0;type_profiles={typ:profile(rows) for typ,(rows,_) in TYPES.items()}
 with gzip.open(P/'physical_witnesses.jsonl.gz','rt') as stream:
  for w in witness:
   q=by[w['node']];op=[511^r^(1<<i) for i,r in enumerate(q)]
   for typ,(_,factor) in TYPES.items():
    record=json.loads(stream.readline());need(all(record[k]==v for k,v in w.items()),'physical record witness metadata')
    need(record['type']==typ and record['profile']==type_profiles[typ],'physical four-child type')
    for target,qq in [(False,q),(True,op)]:
     side='target' if target else 'source';rows,cells=build(qq,w['position'],typ,w['a'],w['c'])
     need(record['augmented_'+side+'_rows']==rows,'every complete physical graph, common orientation and port convention')
     need(record[side+'_residue']==factor*w['normalized_'+side]%8,'physical residue odd-factor restoration')
    need((record['source_residue']-record['target_residue'])%8==4,'physical nonzero gap');count+=1
  need(stream.read()=='','no extra physical records')
 need(count==63396,'all physical query pairs')
 large=load(D/'large_physical_records.json');largeout=[list(map(int,l.split())) for l in (D/'large_physical.out').read_text().splitlines()]
 need(len(large)==len(largeout)==186,'all focused physical graph counts')
 actual_large_input=[]
 for r,out in zip(large,largeout):
  w=r['witness'];qq=by[w['node']]
  if r['opposite']:qq=[511^z^(1<<i) for i,z in enumerate(qq)]
  rows,cells=build(qq,w['position'],r['type'],w['a'],w['c']);cycles,groups=compress(rows,cells,w['position'],r['type'])
  need(rows==r['rows'] and cycles==r['cycles'] and groups==r['group_rows'],'physical-subset input is actual graph')
  need(out==[r['index'],r['expected']],'independent physical subset agrees with normalized counting')
  actual_large_input.append(' '.join(map(str,[r['index'],cycles,len(groups),0]+groups)))
 need((D/'large_physical.in').read_text()==str(len(large))+'\n'+'\n'.join(actual_large_input)+'\n','large physical replay input')
 chosen_placements={(r['witness']['node'],r['witness']['position']) for r in large}
 need(len(chosen_placements)==31 and {(r['node'],r['position']) for r in miss}<=chosen_placements,'large controls contain every initial miss')
 need({r['witness']['position'] for r in large}==set(range(9)),'large controls every four-child position')
 need({r['witness']['query_index'] for r in large}==set(hist),'large controls every used query')
 small=load(D/'small_physical_records.json');smallout=[list(map(int,l.split())) for l in (D/'small_physical.out').read_text().splitlines()]
 need(len(small)==1644 and len(smallout)==3288,'all physical orbit versus literal controls')
 lines=[]
 for r in small:
  idx=r['index'];k=r['k'];code=r['code'];qq=[0]*k;b=0
  for i in range(k):
   for j in range(i+1,k):
    if code>>b&1:qq[i]|=1<<j
    else:qq[j]|=1<<i
    b+=1
  actual,cells=build(qq,r['position'],r['type'],r['a'],r['c']);cycles,group=compress(actual,cells,r['position'],r['type'])
  if idx%7==0:
   perm=list(range(len(actual)));random.Random(idx).shuffle(perm);rel=[0]*len(actual)
   for i,row in enumerate(actual):
    for j in range(len(actual)):
     if row>>j&1:rel[perm[i]]|=1<<perm[j]
   actual=rel
  need(actual==r['physical_rows'] and cycles==r['cycles'] and group==r['group_rows'],'small physical common relabeling and compression')
  need(smallout[2*idx][0]==2*idx and smallout[2*idx+1][0]==2*idx+1,'small output identity')
  need(smallout[2*idx][1]==smallout[2*idx+1][1]%8,'literal physical count equals compressed recurrence')
  lines.append(' '.join(map(str,[2*idx,r['cycles'],len(r['group_rows']),0]+r['group_rows'])))
  lines.append(' '.join(map(str,[2*idx+1,0,len(r['physical_rows']),1]+r['physical_rows'])))
 need((D/'small_physical.in').read_text()==str(len(lines))+'\n'+'\n'.join(lines)+'\n','literal physical replay input')
 local=local_checks()
 receipt=dict(status='PASS',scope='Independent finite size28 response audit and source-relative lower bound29; no general injectivity or arbitrary ambient lift',residual_quotients=2348,placements=21132,unique_recount_queries=6283,independent_placement_query_records=rows_count,independent_residues=2*rows_count,retained_physical_pairs=count,retained_physical_graphs=2*count,focused_physical_placements=31,focused_physical_graphs=186,small_literal_graphs=1644,first_query_histogram=[hist[i] for i in range(14)],local_controls=local,producer_proof_sha256=sha(P/'PROOF.md'))
 receipt['audit_sha256']={name:sha(D/name) for name in ['split_words.cpp','physical_subsets.cpp','prepare_physical.py','verify.py','INPUT_PINS.json','DEPENDENCIES.json','PRODUCER_PINS.json','all_recount.in','all_recount.out','large_physical.in','large_physical.out','large_physical_records.json','small_physical.in','small_physical.out','small_physical_records.json']}
 (D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
