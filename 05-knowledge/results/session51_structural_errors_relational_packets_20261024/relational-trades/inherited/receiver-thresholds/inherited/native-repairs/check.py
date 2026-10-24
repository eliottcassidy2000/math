#!/usr/bin/env python3
"""Exact linked-native receiver repair and full-record route checker, stdlib only."""
from pathlib import Path
from itertools import combinations
from collections import defaultdict, Counter
from hashlib import sha256
from array import array
import argparse, csv, gzip, json, sys
P=Path(__file__).resolve().parent; I=P/'inherited/inputs'; CHECKS=0
U,V=1965,2234; FROZEN={176,1425,2071}; OLD_PAIR={1457,1841}
def need(ok,why):
 global CHECKS
 CHECKS+=1
 if not ok: raise RuntimeError(why)
def read(name):return json.loads((I/name).read_text())
def add(a,b):return (((a//16+b//16)%128)<<4)|((a^b)&15)
def sub(a,b):return (((a//16-b//16)%128)<<4)|((a^b)&15)
def total(xs):
 out=0
 for x in xs:out=add(out,x)
 return out
def source_pins():
 pins=json.loads((P/'SOURCEPINS.json').read_text())
 need(set(pins)=={p.relative_to(P).as_posix() for p in (P/'inherited').rglob('*') if p.is_file() and '__pycache__' not in p.parts},'complete inherited source closure')
 for name,row in pins.items():
  data=(P/name).read_bytes();need(len(data)==row['bytes'] and sha256(data).hexdigest()==row['sha256'],'source pin '+name)
def native():
 old=read('native_packet_decomposition_20261023/native/domains.json')['sources'];supply=read('compatible_charged_helpers_20261023/supply/selected_supply.json');repair=read('guarded_component_completion_20261024/sections/native_2125_receipt.json');host=read('guarded_component_completion_20261024/native/xor_tree_receipt.json')
 electric=read('metagraph_order8_electrical_20261023/receipt.json')['order8'];X=electric['higher_layers'][-1]['lifts'][0];Y=electric['radical_lifts_mod4']
 labels={s:((X[s]%128)<<4)|sum((Y[j][s]%2)<<(j-1) for j in range(1,5)) for s in range(len(X))}
 data={r['source']:{'records':sorted(r['exact_all_shared_clocks_target_domain'],key=lambda x:x['target']),'fixed':r['fixed_clocks'],'blocks':[]} for r in old};cores=sorted(data)
 for r in supply['helpers']:
  d=r['selected'];data[r['source']]={'records':[d['before'],d['after']],'fixed':d['common_clocks'],'blocks':[]}
 data[2125]={'records':repair['records'],'fixed':[],'blocks':repair['retained_blocks']}
 costs={s:sum(a!=b for a,b in zip(*(r['path'] for r in d['records']))) for s,d in data.items()}
 cheap=sorted(s for s in cores if costs[s]==2);outside=sorted(set(cores)-set(cheap));base=sorted(set(cheap)-FROZEN-{1457});free=sorted(set(data)-set(cores)-{368,776,3139})
 need((len(data),len(cheap),len(base),len(outside),len(free))==(82,20,16,8,51),'unchanged native inventories')
 prior=json.loads((P/'inherited/native_certificate.json').read_text());need(all(prior['labels'][str(s)]==labels[s] and prior['costs'][str(s)]==costs[s] for s in data),'join original electrical labels and actual clocks to full native authentication')
 nodes=list(csv.DictReader((I/'metagraph_order8_20261023/merged_nodes.tsv').open(),delimiter='\t'));owners=array('H');owners.frombytes(gzip.decompress((I/'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
 if sys.byteorder!='little':owners.byteswap()
 merged={int(r[k]):int(r['merged_id']) for r in nodes for k in ('class_a','class_b')};pairs=list(combinations(range(8),2));chords=[(i,j) for i in range(8) for j in range(i+2,8)];full=(1<<21)-1
 reflect=lambda m:sum(((m>>k)&1)<<chords.index((7-j,7-i)) for k,(i,j) in enumerate(chords))
 word=lambda A,p:sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(chords))
 reserved=set();literal={}
 for s,d in data.items():
  A=[[0]*8 for _ in range(8)];code=int(nodes[s]['canon_a'])
  for k,(i,j) in enumerate(pairs):A[i][j]=(code>>(27-k))&1;A[j][i]=1-A[i][j]
  literal[s]=[];frame=host['source_frames'][str(s)]
  if s in cores:need(all(not set(frame)&set(host['source_frames'][str(t)]) for t in data if t!=s),'complete private core frame')
  for bit,r in enumerate(d['records']):
   p=r['path'];m=r['mask'];idx={v:i for i,v in enumerate(p)};K=[[0 if i==j else A[i][j] if abs(idx[i]-idx[j])==1 else A[j][i] for j in range(8)] for i in range(8)]
   need(r['target']==(U,V)[bit] and all(A[a][b] for a,b in zip(p,p[1:])) and word(A,p)==m,'actual native path and target bit')
   need(merged[owners[m]]==s and merged[owners[full^m]]==(U,V)[bit],'source and target classification, not projected label')
   need(word(K,p)==full^m and word(list(map(list,zip(*A))),p[::-1])==reflect(m),'literal K lift and converse')
   need(all(p[t]==v for t,v in d['fixed']) and all(any(p[i:i+len(b)]==b for i in range(9-len(b))) for b in d['blocks']),'original full-record pins')
   orbit={m,reflect(m),full^m,full^reflect(m)};need(len(orbit)==4 and not orbit&reserved,'disjoint complete alternative reservation');reserved|=orbit
   literal[s].append({'words':orbit,'global_path':[frame[v] for v in p]})
 need(len(reserved)==656 and host['whole_source_host_order']==652,'unchanged reserved units and host')
 need([literal[368][0]['global_path'][1],literal[776][1]['global_path'][4],literal[3139][1]['global_path'][7]]==[4,4,4],'literal retained owner4 tree')
 return dict(data=data,cores=cores,cheap=cheap,outside=outside,base=base,free=free,labels=labels,costs=costs,literal=literal,potentials=[X]+Y,moduli=[256]+[4]*len(Y))
def compute(n):
 L=n['labels'];base=n['base'];Q=set(range(2048));z2=total(L[s] for s in OLD_PAIR)
 # Two independent support calculations: grade DP, then literal full subsets.
 dp=[{0}]
 for s in base:
  dp.append(set())
  for k in range(len(dp)-1,0,-1):dp[k]|={add(q,L[s]) for q in dp[k-1]}
 by=[defaultdict(list) for _ in range(17)];qs=[0]*(1<<16)
 for mask in range(1<<16):
  if mask:qs[mask]=add(qs[mask&(mask-1)],L[base[(mask&-mask).bit_length()-1]])
  by[mask.bit_count()][qs[mask]].append(mask)
 for k in range(17):need(set(by[k])==dp[k],'independent literal65536 subsets match gradeDP')
 def table(p,allowed):
  out=[{} for _ in range(20)]
  for a,b in allowed:
   shift=add(z2 if a else 0,L[p] if b else 0);extra=(OLD_PAIR if a else set())|({p} if b else set())
   for k,row in enumerate(by):
    for q,masks in row.items():
     qq=add(q,shift);out[k+2*a+b].setdefault(qq,frozenset(extra|{s for j,s in enumerate(base) if masks[0]>>j&1}))
  return out
 def support(p,allowed):
  out=[set() for _ in range(20)]
  for a,b in allowed:
   shift=add(z2 if a else 0,L[p] if b else 0)
   for k,row in enumerate(dp):out[k+2*a+b]|={add(q,shift) for q in row}
  return out
 old=[set() for _ in range(19)]
 for a in (0,1):
  for k,row in enumerate(dp):old[k+2*a]|={add(q,z2 if a else 0) for q in row}
 need({str(k):sorted(Q-old[k]) for k in (8,9,10)}=={'8':[477],'9':[],'10':[1098]},'old genuine receiver-reaching holes')
 rows=[];tables={};complement_checks=0
 for p in sorted(set(n['outside'])-{1841}):
  complete_charge=total(L[s] for s in set(base)|OLD_PAIR|{p})
  need(all(not x for x in support(p,[])),'empty relational table supplies no native record')
  for mask in range(1,16):
   allowed=[(a,b) for a in (0,1) for b in (0,1) if mask>>(2*a+b)&1];ds=support(p,allowed);literal=table(p,allowed)
   for k in range(20):need(ds[k]==set(literal[k]),'all literal allowed full-record branches agree, correct grade2a+b')
   complement=support(p,[(1-a,1-b) for a,b in allowed]);complement_mask=sum(1<<(3-i) for i in range(4) if mask>>i&1)
   need(all({sub(complete_charge,q) for q in ds[k]}==complement[19-k] for k in range(20)),'complement transforms full relation and grade, not a physical move');complement_checks+=1
   full=[k for k,row in enumerate(ds) if row==Q];bitmaps=b''.join(sum(1<<q for q in row).to_bytes(256,'little') for row in ds)
   rows.append({'source':p,'cost':n['costs'][p],'allowed_table_mask':mask,'complement_table_mask':complement_mask,'allowed_pair_new_bits':allowed,'coverage':[len(x) for x in ds],'complete_grades':full,'first_missing':[min(Q-x) if x!=Q else None for x in ds],'support_bitmap_sha256':sha256(bitmaps).hexdigest()})
   if mask==9:tables[p]=literal
   if mask==15:need(full==list(range(6,14)) if p==1923 else full==list(range(7,13)),'exact independent extra-slot band')
   if mask==9:need(full==({1532:[8,9,10,11],1923:[9,10],2087:[],2553:[],2835:[9,10],3015:[8,9,10,11],3360:[]}[p]),'exact native triple-equality band')
  # Fixed inventory adding a relation only removes assignments.
  free=support(p,[(0,0),(0,1),(1,0),(1,1)])
  for row in [r for r in rows if r['source']==p]:
   restricted=support(p,row['allowed_pair_new_bits']);need(all(x<=y for x,y in zip(restricted,free)),'relation monotonicity on common enlarged inventory')
 # Exact cost of changing the triple branch at fixed grade and same residual charge.
 p=3015;triple=OLD_PAIR|{p};z3=total(L[s] for s in triple);need(z3==1263,'independent triple label')
 boundaries=[];sharp=None
 for k in (8,9,10,11):
  hist=Counter();mins=[]
  for q in range(2048):
   left=by[k].get(q,[]);right=by[k-3].get(sub(q,z3),[])
   if not left or not right:continue
   d,a,b=min(((a^b).bit_count(),a,b) for a in left for b in right);cost=2*d+sum(n['costs'][s] for s in triple);hist[cost]+=1
   mins.append({'charge':q,'clock_minimum':cost,'slot_minimum':d+3,'U_branch_base_mask':a,'V_branch_base_mask':b,'branch_populations':[len(left),len(right)]})
   if k==8 and q==53:
    need((len(left),len(right),d,cost)==(1,1,13,34),'sharp two-state fixed-exterior native fibre')
    sharp={'receiver_grade':8,'receiver_charge':53,'triple_sources':sorted(triple),'states_V':[[s for j,s in enumerate(base) if a>>j&1],sorted(triple|{s for j,s in enumerate(base) if b>>j&1})],'exact_slots':16,'exact_clocks':34,'complete_fibre_size':2}
  boundaries.append({'grade':k,'both_branch_charges':len(mins),'exact_clock_histogram':dict(sorted(hist.items())),'minimizing_witnesses':mins})
 return {'frozen_cheap_sources':sorted(FROZEN),'base_sources':base,'base_coverage':[len(x) for x in dp],'old_pair_sources':sorted(OLD_PAIR),'old_pair_coverage':[len(x) for x in old],'old_pair_holes':{str(k):sorted(Q-old[k]) for k in (8,9,10)},'table_census':rows,'complement_table_checks':complement_checks,'triple3015_exact_boundaries':boundaries,'sharp_fixed_exterior_fibre':sharp},tables

def budgets(n):
 rows=[]
 for p in (1532,3015):
  R=set(n['base'])|OLD_PAIR|{p};rc=sorted((n['costs'][s] for s in R),reverse=True);W=lambda a:sum(rc[:min(a,19)]);c=n['costs'][p]
  for e in range(7):
   mx={name:0 for name in ('entry','adjacent','outside','helper')};support=dict(mx)
   for k in range(8,12+e):
    r0=max(8,k-e)
    for r in range(20):
     if not 0<=k-r<=e:continue
     if r<8:s=8
     elif r>11:s=11
     else:continue
     d=min(r+s,38-r-s);de=abs(r-s);mx['entry']=max(mx['entry'],W(d)+4*de);support['entry']=max(support['entry'],d+de)
    for r in range(r0+1,12):
     if 0<=k-r<e:
      d=min(2*r-1,39-2*r);mx['adjacent']=max(mx['adjacent'],W(d)+4);support['adjacent']=max(support['adjacent'],d+1)
    d=2*min(r0,19-r0);mx['helper']=max(mx['helper'],W(d)+15);support['helper']=max(support['helper'],d+2)
    if 0<k-r0<e:mx['outside']=max(mx['outside'],W(d)+8);support['outside']=max(support['outside'],d+2)
   need(max(mx.values())<=50+c and max(support.values())<=20,'every feasible exact-grade route bound')
   rows.append({'repair_source':p,'free_outside':e,'interior':[8,11+e],'clock_maxima':mx,'slot_maxima':support,'clock_bound':50+c,'slot_bound':20})
 need(sum(sorted((n['costs'][s] for s in n['cores']),reverse=True)[:14])==39,'actual top14 core boundary cost39')
 need(sum(sorted((n['costs'][s] for s in n['free']),reverse=True)[:4])==29,'actual top4 helper boundary cost29')
 freecore=set(n['cores'])-FROZEN
 need(len(freecore)==25 and sum(sorted((n['costs'][s] for s in freecore),reverse=True)[:24])==59,'single-event no-return core bound24slots59clocks')
 return rows

def routes(n,tables,sharp):
 L=n['labels'];C=n['costs'];charge=lambda ss:total(L[s] for s in ss);stats={'routes':0,'arrivals':0,'events':0,'maximum_slots':0,'maximum_clocks':0};literal=None;returning=None
 def physical(core,hp):
  active=set(core)|set(hp)|{776,3139};rows={s:[0]*len(n['potentials']) for s in list(n['data'])+[U,V]};used=set()
  for s in n['data']:
   bit=int(s in active);record=n['data'][s]['records'][bit];words=n['literal'][s][bit]['words'];need(not used&words,'arrival capacityone without extra kernel units');used|=words;t=record['target']
   for j,f in enumerate(n['potentials']):d=2*(f[s]-f[t]);rows[s][j]+=d;rows[t][j]-=d
  need(len(used)==328 and 368 not in active,'exact82 native occupied units owner4')
  return rows
 for p in (1532,3015):
  triple=OLD_PAIR|{p};R=set(n['base'])|triple;tabs=tables[p]
  for e in (0,1,6):
   E=set(sorted(set(n['outside'])-{1841,p})[:e]);fixed_core=set(n['cores'])-R-E
   for fixedV in (set(),set(sorted(fixed_core)[::2])):
    for k in range(8,12+e):
     for r in range(20):
      if not 0<=k-r<=e or not tabs[r]:continue
      hp0=frozenset(n['free'][:2]);hp1=frozenset(n['free'][-2:]);ec0=set(sorted(E)[:k-r]);receiver=next(iter(tabs[r].values()));start=(frozenset(fixedV|ec0|set(receiver)),hp0);q=charge(set(start[0])-fixedV|set(hp0))
      def section(ec,hp,grade):return(frozenset(fixedV|set(ec)|set(tabs[grade][sub(q,charge(set(ec)|set(hp)))])),frozenset(hp))
      r0=max(8,k-e);j=k-r0;ec1=set(sorted(E)[-j:]) if j else set();end=section(ec1,hp1,r0)
      cc,hp=start;ec=set(cc)&E;rr=len(set(cc)&R);path=[start]
      if rr<8:ec-=set(sorted(ec)[:8-rr]);rr=8;path.append(section(ec,hp,rr))
      elif rr>11:ec|=set(sorted(E-ec)[:rr-11]);rr=11;path.append(section(ec,hp,rr))
      while rr>r0:ec.add(min(E-ec));rr-=1;path.append(section(ec,hp,rr))
      for a,b in zip(sorted(ec-ec1),sorted(ec1-ec)):ec.remove(a);ec.add(b);path.append(section(ec,hp,r0))
      hp=set(hp)
      for a,b in zip(sorted(hp-set(hp1)),sorted(set(hp1)-hp)):hp.remove(a);hp.add(b);path.append(section(ec,hp,r0))
      path.append(end);path=[x for j,x in enumerate(path) if not j or x!=path[j-1]];need(path[0]==start and path[-1]==end,'literal endpoints and actual sections retained')
      first=physical(*path[0]);ev=[];changes=Counter()
      for core,hp in path:
       need(set(core)&fixed_core==fixedV and len(set(core)-fixedV)==k and len(hp)==2,'frozen full records and separate exact populations')
       actual_targets={n['data'][s]['records'][int(s in core)]['target'] for s in triple};need(len(actual_targets)==1,'native target-class equality on actual full records')
       need(charge((set(core)-fixedV)|set(hp))==q,'entire original paired charge retained')
       rows=physical(core,hp);need(all((rows[s][j]-first[s][j])%n['moduli'][j]==0 for s in rows for j in range(len(n['moduli']))),'every original generating physical source and target row, diagonal and paired factor2')
       stats['arrivals']+=1
      for a,b in zip(path,path[1:]):
       delta=(set(a[0])^set(b[0]))|(set(a[1])^set(b[1]));cost=sum(C[s] for s in delta);need(len(delta)<=20 and cost<=50+C[p],'native interior20slot50+cclock bound');changes.update(delta)
       ev.append({'changed_sources':sorted(delta),'slots':len(delta),'clocks':cost});stats['events']+=1;stats['maximum_slots']=max(stats['maximum_slots'],len(delta));stats['maximum_clocks']=max(stats['maximum_clocks'],cost)
      stats['routes']+=1
      record={'triple_sources':sorted(triple),'free_outside_sources':sorted(E),'fixed_core_V':sorted(fixedV),'free_core_grade':k,'residual_charge':q,'states':[{'core_V':sorted(c),'free_helper_V':sorted(h)} for c,h in path],'events':ev}
      if literal is None and p==3015 and e==6 and rr==r0 and len(path)>4:literal=record
      if returning is None and max(changes.values(),default=0)>1:returning={**record,'sources_returning_to_departed_option':sorted(s for s,t in changes.items() if t>1)}
 # Exact sharp two-state receiver fibre with every exterior record fixed.
 hp=frozenset(n['free'][:2]);x,y=[frozenset(v) for v in sharp['states_V']];rx,ry=physical(x,hp),physical(y,hp)
 need(all((rx[s][j]-ry[s][j])%n['moduli'][j]==0 for s in rx for j in range(len(n['moduli']))),'sharp native boundary all original generator rows')
 delta=x^y;need(len(delta)==16 and sum(C[s] for s in delta)==34,'sharp boundary literal per-slot clocks')
 sharp['full_native_states']=[{'core_V':sorted(v),'free_helper_V':sorted(hp),'helper_V':sorted(set(hp)|{776,3139})} for v in (x,y)]
 need(literal is not None and returning is not None,'literal successful route and temporal-return warning')
 return {'summary':stats,'literal_route':literal,'literal_returning_route':returning}

def main():
 source_pins();n=native();data,tables=compute(n);data['budgets']=budgets(n);data['routes']=routes(n,tables,data['sharp_fixed_exterior_fibre']);data['status']='FINITE-EXACT native receiver-spanning equality repair; conditional all-size mechanism not all-order supply';data['checks']=CHECKS
 return json.loads(json.dumps(data,sort_keys=True))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();out=main();target=P/'certificate.json'
 if args.write:target.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 else:need(out==json.loads(target.read_text()),'frozen independent finite certificate')
 print(json.dumps({'status':'PASS','checks':CHECKS,'table_cases':len(out['table_census']),'sharp_boundary':[16,34],'routes':out['routes']['summary']},sort_keys=True))
