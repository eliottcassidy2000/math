#!/usr/bin/env python3
"""Independent size29 audit: explicit checks remain active under python -O."""
import collections,csv,gzip,hashlib,itertools,json,random,subprocess,sys
from model import D,P,FIVE_A,FIVE_B,FOURS,CYCLE,inflate,ambient,compress,profile,need
from prepare_physical import readw
from local_checks import main as local_checks,decode

def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def integers(p):return [list(map(int,line.split())) for line in p.read_text().splitlines()]
def tab(p):return [{k:int(v) for k,v in row.items()} for row in csv.DictReader(p.open(),delimiter='\t')]
def rootlocals(w):return {j:FOURS[0] for j in w['positions']} if w['type']=='two_four' else {w['positions'][0]:FIVE_A if w['type']=='five_6' else FIVE_B}
def qop(q):return [511^r^(1<<i) for i,r in enumerate(q)]
def replay():
 for program in ['split_words','ambient_words','cycle_covers','physical_subsets']:
  subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra','-Werror',str(D/(program+'.cpp')),'-o',str(D/program)],check=True)
 for program,stem in [('split_words','root_recount'),('ambient_words','ambient_recount'),('cycle_covers','cycle_recount'),('physical_subsets','large_physical'),('physical_subsets','small_physical'),('physical_subsets','paired_physical')]:
  with (D/(stem+'.in')).open('rb') as inp,(D/(stem+'.out')).open('wb') as out,(D/(stem+'.log')).open('wb') as err:subprocess.run([str(D/program)],stdin=inp,stdout=out,stderr=err,check=True)

def main():
 for name,meta in load(D/'INPUT_PINS.json').items():need(sha(P/name)==meta['sha256'] and (P/name).stat().st_size==meta['bytes'],'frozen input '+name)
 if '--recount' in sys.argv:replay()
 cs=load(P/'candidates.json')['candidates'];by={r['node']:r['rows'] for r in cs};need(len(cs)==len(by)==2348,'inherited residual domain')
 for r in cs:
  q=r['rows'];need(len(q)==9 and all(0<=v<512 and not(v>>i&1) for i,v in enumerate(q)) and all(((q[i]>>j)&1)+((q[j]>>i)&1)==1 for i in range(9) for j in range(i)),'actual quotient tournament')
  need(len(r['starts'])==len(r['ends'])==9 and all((a+b)%2==0 for a,b in zip(r['starts'],r['ends'])),'inherited endpoint-zero classification')
 W=readw();need(len(W)==126792,'complete root witness table')
 need([(w['node'],w['configuration']) for w in W]==[(r['node'],i) for r in cs for i in range(1,55)],'every root configuration')
 queries=[(q['a'],q['c']) for q in load(P/'root_battery.json')]
 need(queries==[(0,1),(1,0),(1,2),(1,3),(3,1),(3,3),(3,5),(0,7),(0,18),(40,256),(4,3),(130,128),(65,9),(3,32),(0,2),(1,17)],'fixed root battery')
 need(len(set(queries+[(1<<i,0) for i in range(9)]))==24,'endpoint union')
 bank={};initial=tab(P/'root_residues.tsv');extra=tab(P/'root_selected_residues.tsv');need(len(initial)==32872 and len(extra)==22,'full producer bank shape')
 for r in initial+extra:
  key=r['node'],r['a'],r['c'];need(key not in bank,'unique retained query');v=[r[f'{s}{i}'] for i in range(46) for s in 'ST'];need(all(0<=x<8 and x%2 for x in v),'normalized bank range');bank[key]=v
 chosen=sorted(set((w['node'],w['query_index'],w['a'],w['c']) for w in W));need(len(chosen)==8769,'independent selected-query domain')
 expected=str(len(chosen))+'\n'+''.join(' '.join(map(str,[n,qi,a,c]+by[n]))+'\n' for n,qi,a,c in chosen)
 need((D/'root_recount.in').read_text()==expected,'root replay input exact domain')
 independent={}
 for row in integers(D/'root_recount.out'):
  need(len(row)==95,'all46 outputs');key=tuple(row[:3]);need(key not in independent,'no duplicated independent query');need(row[3:]==bank[key],'independent root recurrence agrees');independent[key]=row[3:]
 need(set(independent)=={(n,a,c) for n,qi,a,c in chosen},'complete independent root recount')
 pairs=list(itertools.combinations(range(9),2));hist=collections.Counter();miss=[]
 for w in W:
  n,cfg=w['node'],w['configuration'];a,c=w['a'],w['c'];col=cfg if cfg<46 else 0;qi=w['query_index'];need(queries[qi]==(a,c),'same coarse query')
  need(independent[n,a,c][2*col:2*col+2]==[w['normalized_source'],w['normalized_target']],'every selected root separator')
  need((w['normalized_source']-w['normalized_target'])%8==4,'root gap four')
  typ='five_6' if cfg<10 else 'two_four' if cfg<46 else 'five_2';positions=[cfg-1] if cfg<10 else list(pairs[cfg-10]) if cfg<46 else [cfg-46]
  need(w['type']==typ and w['positions']==positions,'configuration interpretation')
  factor=1 if typ=='five_2' else 3;need(w['physical_factor']==factor,'canonical physical odd factor')
  if not any(bank[n,*queries[j]][2*col]!=bank[n,*queries[j]][2*col+1] for j in range(14)):
   miss.append(dict(node=n,configuration=cfg))
  need(not any(bank[n,*queries[j]][2*col]!=bank[n,*queries[j]][2*col+1] for j in range(qi)),'first witness')
  hist[qi]+=1
 need(miss==load(P/'root_initial_misses.json') and len(miss)==11,'exact initial misses')
 need({(r['node'],r['a'],r['c']) for r in initial}=={(r['node'],a,c) for r in cs for a,c in queries[:14]},'complete initial bank universe')
 need({(r['node'],r['a'],r['c']) for r in extra}=={(r['node'],a,c) for r in miss for a,c in queries[14:]},'complete selected extra bank universe')
 physical=0
 with gzip.open(P/'root_physical_witnesses.jsonl.gz','rt') as f:
  for w in W:
   record=json.loads(f.readline());need(all(record[k]==v for k,v in w.items()),'root physical metadata');q=by[w['node']];local=rootlocals(w)
   for name,qq in [('source',q),('target',qop(q))]:
    rows,cells=inflate(qq,local,w['a'],w['c']);need(record['augmented_'+name+'_rows']==rows,'every retained root physical graph')
    need(record[name+'_residue']==w['physical_factor']*w['normalized_'+name]%8,'physical residue restoration')
   for name,mask in [('physical_a',w['a']),('physical_c',w['c'])]:need(w[name]==sum(1<<v for i,cell in enumerate(cells) if mask>>i&1 for v in cell),'root physical mask')
   physical+=1
  need(f.read()=='','no extra physical root records')
 need(physical==126792,'all canonical physical root pairs')
 # Entire independent cycle-cover calculation, not merely selected parities.
 cyctab=tab(P/'cycle_parity.tsv');cmap={}
 for r in cyctab:
  for i,(a,c) in enumerate(queries[:12]):cmap[r['node'],a,c]=(r[f'S{i}'],r[f'T{i}'])
 cycout=integers(D/'cycle_recount.out');need(len(cycout)==28177,'complete determinant cycle bank')
 expected_rows=[(r['node'],a,c,r['rows']) for r in cs for a,c in queries[:12]]+[(56186,2,4,by[56186])]
 need((D/'cycle_recount.in').read_text()==str(len(expected_rows))+'\n'+''.join(' '.join(map(str,[n,a,c]+q))+'\n' for n,a,c,q in expected_rows),'cycle determinant replay inputs')
 cind={}
 for row,e in zip(cycout,expected_rows):
  n,a,c,s,t=row;need((n,a,c)==e[:3] and s in [0,1] and t in [0,1],'cycle output domain')
  if (n,a,c) in cmap:need((s,t)==cmap[n,a,c],'cycle-cover determinant equals anchored cycle bank')
  else:need((n,a,c,s,t)==(56186,2,4,1,1),'extra cycle parity')
  cind[n,a,c]=(s,t)
 need(sum(s!=t for s,t in cmap.values())==7624,'retained cancellation/parity changes')
 # Direct13-symbol ambient path recount does not use q2 or the cycle identity.
 ambient_w=tab(P/'witnesses.tsv');need(len(ambient_w)==2348 and [w['node'] for w in ambient_w]==[r['node'] for r in cs],'all ambient witness classes')
 expected=str(len(ambient_w))+'\n'+''.join(' '.join(map(str,[w['node'],w['a'],w['c']]+by[w['node']]))+'\n' for w in ambient_w)
 need((D/'ambient_recount.in').read_text()==expected,'direct ambient replay input')
 ambient_out=integers(D/'ambient_recount.out');need(len(ambient_out)==2348,'all direct ambient counts')
 with gzip.open(P/'physical_witnesses.jsonl.gz','rt') as f:
  for w,row in zip(ambient_w,ambient_out):
   n,a,c=w['node'],w['a'],w['c'];need(row==[n,a,c,w['ambient_source'],w['ambient_target']],'every independently recounted ambient separator')
   need((row[-2]-row[-1])%8==4,'ambient gap four')
   need(cind[n,a,c]==(w['cycle_source'],w['cycle_target']),'ambient selected cycle parities')
   for side in ['source','target']:need((5*w['core_'+side]+4-4*w['cycle_'+side])%8==w['ambient_'+side],'q2/cycle formula agrees with direct ambient count')
   r=json.loads(f.readline());need(r['node']==n and r['coarse_a']==a and r['coarse_c']==c,'ambient same query metadata')
   for side,qq in [('source',by[n]),('target',qop(by[n]))]:
    rows,cells=ambient(qq,a,c);need(r[side+'_rows']==rows and r[side+'_residue']==w['ambient_'+side],'every retained physical cyclic ambient graph')
   for name,mask in [('physical_a',a),('physical_c',c)]:need(r[name]==sum(7<<(3*i) for i in range(9) if mask>>i&1)+(1<<28),'ambient mask incidences')
  need(f.read()=='','no extra ambient physical records')
 # Third counting method: literal physical-subset recurrence with rotational
 # compression only on actual physical triangle modules.
 wmap={(w['node'],w['configuration']):w for w in W};amap={w['node']:w for w in ambient_w}
 large=load(D/'large_physical_records.json');out=integers(D/'large_physical.out');need(len(large)==len(out)==226,'all focused physical graphs')
 lines=[]
 for r,o in zip(large,out):
  w=r['witness'];q=by[w['node']];qq=qop(q) if r['opposite'] else q
  if r['kind']=='root':
   need(w==wmap[w['node'],w['configuration']],'focused root witness identity');locals_={int(k):v for k,v in r['locals'].items()};rows,cells=inflate(qq,locals_,w['a'],w['c']);factor=1
   for i in range(9):factor=factor*profile(locals_.get(i,CYCLE))[0]%8
   need(r['expected']==factor*w['normalized_target' if r['opposite'] else 'normalized_source']%8,'focused physical expected residue independently restored')
  else:
   need(w==amap[w['node']],'focused ambient witness identity');rows,cells=ambient(qq,w['a'],w['c'])
   need(r['expected']==w['ambient_target' if r['opposite'] else 'ambient_source'],'focused ambient direct-word expected residue')
  cycles,groups=compress(rows,cells);need(rows==r['physical_rows'] and cycles==r['cycles'] and groups==r['group_rows'],'focused graph and true rotation blocks')
  need(o==[r['index'],r['expected']],'physical-subset count matches direct word counts')
  lines.append(' '.join(map(str,[r['index'],cycles,len(groups),0]+groups)))
 need((D/'large_physical.in').read_text()==str(len(lines))+'\n'+'\n'.join(lines)+'\n','large physical input')
 need({r['witness']['configuration'] for r in large if r['kind']=='root'}==set(range(1,55)),'all54 configuration positions in focused physical controls')
 small=load(D/'small_physical_records.json');out=integers(D/'small_physical.out');need(len(small)==11472 and len(out)==22944,'all small physical controls');lines=[]
 for r in small:
  idx=r['index'];rows,cells=inflate(decode(r['code'],r['k']),{int(k):v for k,v in r['locals'].items()},r['a'],r['c']);cycles,groups=compress(rows,cells)
  if idx%17==0:
   perm=list(range(len(rows)));random.Random(idx).shuffle(perm);rel=[0]*len(rows)
   for i,row in enumerate(rows):
    for j in range(len(rows)):
     if row>>j&1:rel[perm[i]]|=1<<perm[j]
   rows=rel
  need(rows==r['physical_rows'] and cycles==r['cycles'] and groups==r['group_rows'],'small actual physical graph and common permutation')
  need(out[2*idx][0]==2*idx and out[2*idx+1][0]==2*idx+1 and out[2*idx][1]==out[2*idx+1][1]%8,'literal physical integer count agrees with orbit recurrence')
  lines += [' '.join(map(str,[2*idx,cycles,len(groups),0]+groups)),' '.join(map(str,[2*idx+1,0,len(rows),1]+rows))]
 need((D/'small_physical.in').read_text()==str(len(lines))+'\n'+'\n'.join(lines)+'\n','small physical input')
 # Profile-changing families challenge the optional paired-battery statement.
 paired=load(D/'paired_records.json');vals=dict(integers(D/'paired_physical.out'));need(len(vals)==768 and len(paired)==256,'paired integer controls')
 diffs={};pairedhits=collections.Counter();lines=[]
 for r in paired:
  for idx,qq,local in zip(r['indices'],[[2,4,1],[4,1,2],[4,1,2]],[decode(r['code'],5),decode(r['code'],5),decode(r['target_code'],5)]):
   rows,_=inflate(qq,{0:local},r['a'],r['c']);lines.append(' '.join(map(str,[idx,0,len(rows),1]+rows)))
  s,t,u=[vals[i] for i in r['indices']];need(u>=t,'profile-changing target monotonicity');diffs[r['code'],r['a'],r['c']]=(s-t,s-u)
 need((D/'paired_physical.in').read_text()==str(len(lines))+'\n'+'\n'.join(lines)+'\n','paired integer replay input')
 for (code,a,c),(d,actual) in diffs.items():
  need(diffs[code,7^c,7^a][0]==-d,'exact converse-complement sign reversal')
  if d<0:need(actual<0,'negative baseline contrast cannot cancel')
  if d:pairedhits[code]+=1
 need(pairedhits=={2:44,5:44,8:44,11:44},'all profile-changing classes exercise nonzero paired contrasts')
 local_checks();local=load(D/'local_receipt.json')
 receipt=dict(status='PASS',scope='Independent size29 finite root and cyclic-exterior audit; original-source-module collision lower bound30; no general injectivity',residual_quotients=2348,new_root_placements=105660,inherited_profile_five_placements=21132,retained_root_pairs=126792,independent_root_queries=8769,independent_root_configuration_records=8769*46,independent_root_residues=2*8769*46,independent_cycle_query_pairs=28177,independent_ambient_pairs=2348,retained_physical_graphs=2*(126792+2348),focused_large_physical_graphs=226,small_literal_graphs=11472,paired_literal_graphs=768,local=local,first_query_histogram=[hist[i] for i in range(16)],producer_proof_sha256=sha(P/'root_PROOF.md'))
 names=['split_words.cpp','ambient_words.cpp','cycle_covers.cpp','physical_subsets.cpp','model.py','prepare_physical.py','paired_controls.py','local_checks.py','verify.py','INPUT_PINS.json','DEPENDENCIES.json','PRODUCER_PINS.json']
 names += [f'{stem}.{suffix}' for stem in ['root_recount','ambient_recount','cycle_recount','large_physical','small_physical','paired_physical'] for suffix in ['in','out']]
 names += ['large_physical_records.json','small_physical_records.json','paired_records.json','local_identity_controls.json','local_receipt.json']
 receipt['audit_sha256']={name:sha(D/name) for name in names}
 (D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
