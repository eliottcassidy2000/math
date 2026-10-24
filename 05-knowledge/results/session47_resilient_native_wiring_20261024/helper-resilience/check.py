#!/usr/bin/env python3
"""Exact deletion/packing certificate replay; all checks remain active under -O."""
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse,csv,gzip,json,sys
HERE=Path(__file__).resolve().parent
COUNT=0
U,V=1965,2234
FULL=(1<<21)-1
CHORDS=[(i,j) for i in range(8) for j in range(i+2,8)]
PAIRS=list(combinations(range(8),2))
REV=[CHORDS.index((7-j,7-i)) for i,j in CHORDS]
def need(ok,why):
 global COUNT
 COUNT+=1
 if not ok:raise RuntimeError(why)
def plus(a,b):return ((((a>>4)+(b>>4))%128)<<4)|((a^b)&15)
def minus(a,b):return ((((a>>4)-(b>>4))%128)<<4)|((a^b)&15)
def digest(p):
 b=p.read_bytes()
 if p.suffix in ('.json','.py','.md','.tsv'):b=b.replace(b'\r\n',b'\n')
 return sha256(b).hexdigest()
def reflect(m):return sum(((m>>i)&1)<<REV[i] for i in range(21))
def orbit(m):return {m,FULL^m,reflect(m),FULL^reflect(m)}
def main():
 p=argparse.ArgumentParser();p.add_argument('--results',type=Path);p.add_argument('--write',action='store_true');args=p.parse_args()
 R=args.results
 if R is None:
  candidates=[HERE.parent,HERE.parents[1],Path.cwd()/'session35-publication/clean-repo/05-knowledge/results']
  R=next((z for z in candidates if (z/'compatible_charged_helpers_20261023').exists()),None)
 need(R is not None,'locate inherited results, or pass --results')
 supply_path=R/'compatible_charged_helpers_20261023/supply/selected_supply.json'
 supply=json.loads(supply_path.read_text());hs=supply['helpers'];sources=[z['source'] for z in hs];hset=set(sources)
 electrical_path=R/'metagraph_order8_electrical_20261023/receipt.json'
 E=json.loads(electrical_path.read_text())['order8'];X=E['higher_layers'][-1]['lifts'][0];Y=E['radical_lifts_mod4'];pot=[X]+Y[:5];mods=[256]+[4]*5
 def label(s):return ((X[s]%128)<<4)|sum((Y[j][s]%2)<<(j-1) for j in range(1,5))
 def charge(ids):
  q=0
  for s in ids:q=plus(q,label(s))
  return q
 nodes=list(csv.DictReader((R/'metagraph_order8_20261023/merged_nodes.tsv').open(),delimiter='\t'))
 merged={int(z[k]):int(z['merged_id']) for z in nodes for k in ('class_a','class_b')}
 owner=array('H');owner.frombytes(gzip.decompress((R/'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
 if sys.byteorder!='little':owner.byteswap()
 need(len(sources)==53 and len(hset)==53,'exact53-helper source inventory')
 need(len({label(s) for s in sources})==53,'selected53 labels distinct')
 old_path=R/supply['old_domain'];old=json.loads(old_path.read_text())['sources'];core=sorted(z['source'] for z in old)
 data={z['source']:{'records':{r['target']:r for r in z['exact_all_shared_clocks_target_domain']},'fixed':z['fixed_clocks']} for z in old}
 for z in hs:
  s=z['source'];m=z['selected'];need(z['label']==[X[s]%128]+[Y[j][s]%2 for j in range(1,5)],'original native potentials authenticate helper label')
  data[s]={'records':{U:m['before'],V:m['after']},'fixed':m['common_clocks']}
 reservoir=set();cost={}
 for s,d in sorted(data.items()):
  can=int(nodes[s]['canon_a']);A=[[0]*8 for _ in range(8)]
  for k,(i,j) in enumerate(PAIRS):A[i][j]=(can>>(27-k))&1;A[j][i]=1-A[i][j]
  for target in (U,V):
   r=d['records'][target];path=r['path']
   need(sorted(path)==list(range(8)) and all(A[i][j] for i,j in zip(path,path[1:])),'actual Hamiltonian source path')
   m=sum(A[path[i]][path[j]]<<k for k,(i,j) in enumerate(CHORDS))
   need(m==r['mask'] and merged[owner[m]]==s and merged[owner[FULL^m]]==target,'actual source and target record')
   need(all(path[t]==o for t,o in d['fixed']),'all fixed common owner clocks')
   need((X[s]-X[target])%8==r['phase'],'chosen original X phase')
   words=orbit(m);need(len(words)==4 and not reservoir&words,'capacity-one disjoint normalized alternative reservoir');reservoir|=words
   rev=path[::-1];need(all(A[j][i] for i,j in zip(rev,rev[1:])),'actual converse Hamiltonian path')
   rm=sum(A[rev[j]][rev[i]]<<k for k,(i,j) in enumerate(CHORDS));need(rm==reflect(m),'converse reflected record')
   need(all(rev[7-t]==o for t,o in d['fixed']),'converse retained clocks')
  cost[s]=sum(a!=b for a,b in zip(d['records'][U]['path'],d['records'][V]['path']))
 need(len(data)==81 and len(reservoir)==648,'complete28+53 native reservation')
 need(label(U)==label(V) and X[U]==X[V]==146,'equal paired target charge and chosen integer X')
 def fibers(ids,r):
  out=[[] for _ in range(2048)]
  for t in combinations(ids,r):out[charge(t)].append(frozenset(t))
  return out
 def intersect(f,universe):
  out=set(universe)
  for t in f:out.intersection_update(t)
  return out
 banks={'40':supply['banks'][supply['main_bank_index']]['source_ids'],'46':json.loads((HERE/'robust_bank46.json').read_text())['source_ids'],'48':json.loads((HERE/'robust_bank48.json').read_text())['source_ids'],'53':sources}
 bank_info={};allfib={};minimal=[]
 for name,ids in banks.items():
  need(len(ids)==int(name) and len(set(ids))==len(ids) and set(ids)<=hset,'bank exact actual source set')
  fs=fibers(ids,3);allfib[name]=fs;need(all(fs),'every charge covered before restrictions')
  intersections=[intersect(f,ids) for f in fs]
  if name=='40':need(any(intersections),'compact40 not single-U resilient')
  else:need(not any(intersections),'every one-U deletion keeps every charge')
  bank_info[name]={'triples':sum(map(len,fs)),'minimum_multiplicity':min(map(len,fs)),'maximum_multiplicity':max(map(len,fs)),'one_U_resilient':not any(intersections)}
  if name=='46':
   for removed in ids:
    hostile=None
    for q,f in enumerate(fs):
     remaining=[t for t in f if removed not in t];inter=intersect(remaining,set(ids)-{removed})
     if inter:
      forced=min(inter);need(not any(removed not in t and forced not in t for t in f),'inclusion-minimal46 exact two-U witness');hostile={'removed':removed,'forced_U':forced,'missing_charge':q};break
    need(hostile is not None,'every46 deletion destroys1U resilience');minimal.append(hostile)
 q=305;expected=[{776,1182,2638},{1391,1948,3139},{1400,3111,3139},{2439,2728,3139}]
 need(set(allfib['53'][q])==set(map(frozenset,expected)),'complete exact four-triple hostile at305')
 need(all(t&{776,3139} for t in allfib['53'][q]),'two-U failure for full53 and every subbank')
 mixed=plus(label(776),q);need(mixed==1368,'shifted mixed1V1U hostile charge')
 need(not any(776 not in t and 3139 not in t for t in allfib['53'][q]),'mixed1V1U occupancy4 impossible')
 # All53 occupancy4 witness family. Certificates are checked using original charges.
 packing=json.loads((HERE/'quad_packings.json').read_text());need(packing['source_ids']==sources and len(packing['packings'])==2048,'exact packing universe and all charges')
 counts=[]
 for q,ms in enumerate(packing['packings']):
  seen=0
  for m in ms:
   need(isinstance(m,int) and 0<=m<(1<<53) and m.bit_count()==4,'quadruple is four actual distinct slots')
   need(not seen&m,'representations in this charge fibre are pairwise disjoint');seen|=m
   need(charge(sources[i] for i in range(53) if m>>i&1)==q,'quadruple exact original charge')
  need(len(ms)>=12,'twelve disjoint witnesses survive any eleven U restrictions');counts.append(len(ms))
 # One forcedV at occupancy4: residual triple avoids its specified vertex.
 for q,f in enumerate(allfib['53']):need(not intersect(f,sources),'forced-oneV padding of resilient triples')
 # Universal oneV at occupancy3 fails even granting the complete72 new labels.
 allnew=supply['all_new_source_labels'];need(len(allnew)==72,'all72 new-source alphabet')
 pairs={charge((a['source'],b['source'])) for a,b in combinations(allnew,2)}
 need(len(pairs)==1468 and 2 not in pairs,'universal forcedV requires unavailable residual pair charge2')
 # Causal routes: every total charge, adaptive destination restrictions, exact helper return.
 ids=banks['46'];fs=allfib['46'];paths=0;events=0;max_support=0;max_clock=0;evacuations=0
 owner_ids=sorted(data)+[U,V]
 def rows(active):
  out={s:[0]*6 for s in owner_ids}
  for s in data:
   t=V if s in active else U
   for j,potj in enumerate(pot):w=2*(potj[s]-potj[t]);out[s][j]+=w;out[t][j]-=w
  return out
 a0=frozenset(core[:14]);b0=frozenset(core[14:])
 for q in range(2048):
  h0=fs[minus(q,charge(a0))][0];cur=set(a0);prev=(a0,h0);baseline=rows(a0|h0)
  route=[]
  # Forward and backward swaps return the exact original complete records.
  pairs2=list(zip(sorted(a0-b0),sorted(b0-a0)))
  moves=pairs2+[(b,a) for a,b in reversed(pairs2)]
  for step,(a,b) in enumerate(moves):
   cur.remove(a);cur.add(b);debt=minus(q,charge(cur))
   if step+1==len(moves):
    h=h0;forbidden=next(s for s in ids if s not in h0)
   else:
    forbidden=min(prev[1]);h=next(t for t in fs[debt] if forbidden not in t);evacuations+=int(forbidden in prev[1])
   need(forbidden not in h and len(h)==3 and len(cur)==14,'arrival-indexed restriction and separate exact margins')
   need(plus(charge(cur),charge(h))==q,'every arrival retains exact marked total charge')
   now=rows(set(cur)|h)
   need(all((now[s][j]-baseline[s][j])%mods[j]==0 for s in owner_ids for j in range(6)),'all actual physical-row dyadic currents retain diagonals and factor2')
   changed=(prev[0]^frozenset(cur))|(prev[1]^h);support=len(changed);work=sum(cost[s] for s in changed)
   need(support<=8 and work<=50,'eight-source and50-clock primary budget')
   max_support=max(max_support,support);max_clock=max(max_clock,work);events+=1;prev=(frozenset(cur),h)
  need(prev==(a0,h0),'literal full helper and core endpoint return');paths+=1
 # Strong occupancy4 replay: eleven restrictions hit eleven disjoint options,
 # leaving the twelfth, followed by an admissible exact full-record return.
 strong_routes=0;strong_peak=0;strong_work=0
 quadsets=[[frozenset(sources[i] for i in range(53) if m>>i&1) for m in ms] for ms in packing['packings']]
 a4=frozenset(core[:14]);b4=(a4-{core[0]})|{core[14]}
 for q in range(2048):
  h0=quadsets[minus(q,charge(a4))][0]
  opts=quadsets[minus(q,charge(b4))];h1=opts[11];D={min(t) for t in opts[:11]}
  need(len(D)==11 and not D&h1,'eleven adaptive U restrictions leave twelfth quadruple')
  need(plus(charge(b4),charge(h1))==q,'occupancy4 exact debt')
  oldrows=rows(a4|h0);newrows=rows(b4|h1)
  need(all((newrows[t][j]-oldrows[t][j])%mods[j]==0 for t in owner_ids for j in range(6)),'occupancy4 original physical-row dyadic guards')
  changed=(a4^b4)|(h0^h1);sup=len(changed);work=sum(cost[t] for t in changed)
  need(sup<=10 and work<=64,'occupancy4 ten-source64-clock bound')
  Dreturn=set(t for t in sources if t not in h0);Dreturn=set(sorted(Dreturn)[:11])
  need(len(Dreturn)==11 and not h0&Dreturn,'exact return is permitted by final restriction')
  strong_routes+=1;strong_peak=max(strong_peak,sup);strong_work=max(strong_work,work)
 # Full mixed-pair census. A failing residual pair yields both directed1V/1U controls.
 mixed_controls=[];two_u_controls=[]
 for a,b in combinations(sources,2):
  missing=[q for q,f in enumerate(allfib['53']) if not any(a not in t and b not in t for t in f)]
  if missing:
   two_u_controls.append({'forced_U':[a,b],'missing_charges':missing})
   for q in missing:
    mixed_controls.extend([{'forced_V':a,'forced_U':b,'charge':plus(label(a),q)}, {'forced_V':b,'forced_U':a,'charge':plus(label(b),q)}])
 need(len(two_u_controls)==16 and all(len(z['missing_charges'])==1 for z in two_u_controls),'all16 failing unordered two-U pairs, each exactly one missing charge')
 need(len(mixed_controls)==32,'exact32 failing directed mixed pairs out of2756')
 # Sharpness of the uniform disjoint-witness certificate, not deletion resilience.
 total=charge(sources);inverse13=pow(13,-1,128)
 eligible={(((((total>>4)-(label(v)>>4))*inverse13)%128)<<4)|((total^label(v))&15) for v in sources}
 need(len(eligible)==53 and len(set(range(2048))-eligible)==1995,'thirteen disjoint quadruples are impossible at1995 charges')
 # Small coverage-vs-reachability hostile at lower budgets.
 labs=[0,0,0,1,2,3];small=[[set(t) for t in combinations(range(6),2) if labs[t[0]]^labs[t[1]]==q] for q in range(4)]
 need(all(not intersect(f,range(6)) for f in small),'small abstract bank is single-U resilient')
 need(small[1]==[{0,3},{1,3},{2,3},{4,5}],'exact hostile static charge fibre')
 need(all(len(t^{4,5})==4 for t in small[1][:-1]),'helper2 graph cannot reach isolated destination pair')
 out={'status':'FINITE-EXACT','checks':COUNT,'banks':bank_info,'bank46_inclusion_minimality':minimal,'full53_two_U_hostile':{'charge':305,'forced_U':[776,3139],'all_triples':[sorted(t) for t in expected]},'occupancy4_mixed_hostile':{'forced_V':[776],'forced_U':[3139],'charge':mixed},'quad_packing_histogram':dict(sorted(Counter(counts).items())),'guaranteed_occupancy4_U_restrictions':min(counts)-1,'sharp_uniform_disjoint_packing':12,'charges_excluding13packing':1995,'mixed_restriction_census':{'total_directed_pairs':2756,'full_coverage_pairs':2724,'failure_controls':mixed_controls,'unordered_two_U_controls':two_u_controls},'occupancy4_route_replays':{'routes':strong_routes,'maximum_support':strong_peak,'maximum_clock_work':strong_work,'forced_U_per_arrival':11,'exact_endpoint_return':True},'native':{'slots':81,'reservation_closed_words':648,'arrival_closed_words':324,'helper_clock_costs':dict(sorted(Counter(cost[s] for s in sources).items())),'occupancy3_support_bound':8,'occupancy3_clock_bound':50,'occupancy4_support_bound':10,'occupancy4_clock_bound':64},'causal_route_replays':{'total_charges':paths,'events':events,'actual_forced_evacuations':evacuations,'maximum_support':max_support,'maximum_clock_work':max_clock,'exact_endpoint_return':True},'source_sha256':{str(z.relative_to(R)):digest(z) for z in [supply_path,electrical_path,old_path,R/'compatible_charged_helpers_20261023/supply/check.py',R/'injective_native_blocks_20261023/electrical/ELECTRICAL_CODES.md']},'certificate_sha256':{n:digest(HERE/n) for n in ['robust_bank46.json','robust_bank48.json','quad_packings.json']},'scope':'Selected capacity-one order-eight packets; all original dyadic aggregate currents by complete paired-code theorem, not exact integer/odd-primary currents or every ordinary per-port phase; no all-order native supply claim.'}
 out=json.loads(json.dumps(out));dest=HERE/'receipt.json'
 if args.write:dest.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 else:need(out==json.loads(dest.read_text()),'frozen receipt exact match')
 print(json.dumps({'status':'PASS','checks':COUNT,'banks':bank_info,'quad_minimum':min(counts),'routes':paths,'events':events},sort_keys=True))
if __name__=='__main__':main()
