#!/usr/bin/env python3
"""Independent support-bitset census, full path enumeration and route controls.
No imports from producer checkers. All universal claims remain proof obligations.
"""
from pathlib import Path
from itertools import combinations, permutations, product
from collections import Counter
from hashlib import sha256
import argparse, csv, gzip, json, random, struct

COUNT=0
def check(value, why):
 global COUNT
 COUNT+=1
 if not value: raise RuntimeError(why)
def main():
 global COUNT
 ap=argparse.ArgumentParser();ap.add_argument('packet',type=Path);ap.add_argument('--out',type=Path);a=ap.parse_args();P=a.packet;I=P/'inputs'
 read=lambda name:json.loads((I/name).read_text())
 manifest=json.loads((P/'MANIFEST.json').read_text())
 actual={str(x.relative_to(P)) for x in P.rglob('*') if x.is_file() and '__pycache__' not in x.parts and x!=P/'MANIFEST.json'}
 check(actual==set(manifest),'complete packet including any nested manifests')
 for name,row in manifest.items():
  data=(P/name).read_bytes();check(len(data)==row['bytes'] and sha256(data).hexdigest()==row['sha256'],'manifest hash '+name)
 sourcepins=json.loads((P/'SOURCEPINS.json').read_text())
 check(set(sourcepins)=={str(x.relative_to(P)) for x in I.rglob('*') if x.is_file()},'exact input closure')
 for name,row in sourcepins.items():check(sha256((P/name).read_bytes()).hexdigest()==row['sha256'],'source hash '+name)
 old=read('domains.json')['sources'];supply=read('supply.json')['helpers'];r2125=read('source2125.json');host=read('host.json')
 elect=read('electrical.json')['order8'];potentials=[elect['higher_layers'][-1]['lifts'][0]]+elect['radical_lifts_mod4'][:5];moduli=[256,4,4,4,4,4]
 label=lambda s:(potentials[0][s]%128,sum((potentials[j+2][s]%2)<<j for j in range(4)))
 plus=lambda x,y:((x[0]+y[0])%128,x[1]^y[1]);minus=lambda x,y:((x[0]-y[0])%128,x[1]^y[1])
 encode=lambda q:16*q[0]+q[1]
 def charge(ss):
  q=(0,0)
  for s in ss:q=plus(q,label(s))
  return q
 data={}
 for row in old:
  data[row['source']]={'rec':{x['target']:x for x in row['exact_all_shared_clocks_target_domain']},'pins':row['fixed_clocks'],'blocks':[],'raw':len(row['raw_all_shared_clocks_domain'])}
  check(label(row['source'])==(row['code'][0],sum(v<<j for j,v in enumerate(row['code'][1:]))),'original physical label')
 core=set(data)
 for row in supply:
  r=row['selected'];data[row['source']]={'rec':{x['target']:x for x in [r['before'],r['after']]},'pins':r['common_clocks'],'blocks':[],'raw':r['complete_pinned_count']}
 data[2125]={'rec':{x['target']:x for x in r2125['records']},'pins':r2125['fixed_clocks'],'blocks':r2125['retained_blocks'],'raw':r2125['raw_block_guard_count']}
 costs={s:sum(x!=y for x,y in zip(d['rec'][1965]['path'],d['rec'][2234]['path'])) for s,d in data.items()}
 cheap=sorted(s for s in core if costs[s]==2);outside=sorted(core-set(cheap));helpers=sorted(set(data)-core-{368,776,3139})
 check((len(data),len(core),len(cheap),len(outside),len(helpers))==(82,28,20,8,51),'inventory')
 check(sum(sorted([costs[s] for s in core],reverse=True)[:18])==47,'actual top18 cost')
 check(sum(sorted([costs[s] for s in helpers],reverse=True)[:4])==29,'actual helper top4 cost')
 check(max(costs[x]+costs[y] for x,y in combinations(helpers,2))==15,'actual helper edge max')
 # Reconstruct every source and exhaust its guarded permutations, not just selected paths.
 nodes=list(csv.DictReader((I/'merged_nodes.tsv').open(),delimiter='\t'))
 cl={int(r[k]):int(r['merged_id']) for r in nodes for k in ['class_a','class_b']}
 owner=struct.unpack('<'+'H'*(len(gzip.decompress((I/'owner_map.u16.gz').read_bytes()))//2),gzip.decompress((I/'owner_map.u16.gz').read_bytes()))
 chords=[(i,j) for i in range(8) for j in range(i+2,8)];allpairs=list(combinations(range(8),2));full=(1<<21)-1
 reflidx=[chords.index((7-j,7-i)) for i,j in chords]
 def reflection(w):return sum(((w>>j)&1)<<k for j,k in enumerate(reflidx))
 def word(A,p):return sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(chords))
 resources={};edges={};partner_options={};fullmenu=[];perms=list(permutations(range(8)));allwords=set()
 for s,d in sorted(data.items()):
  A=[[0]*8 for _ in range(8)];canon=int(nodes[s]['canon_a'])
  for k,(i,j) in enumerate(allpairs):A[i][j]=(canon>>(27-k))&1;A[j][i]=1-A[i][j]
  frame=host['source_frames'][str(s)]
  check(len(frame)==len(set(frame))==8,'injective source frame')
  for i,j in allpairs:
   u,v=frame[i],frame[j];edge=(min(u,v),max(u,v));bit=A[i][j] if u<v else A[j][i]
   check(edge not in edges or edges[edge]==bit,'coherent global source edge prescriptions');edges[edge]=bit
  found=[];raw=0
  for p in perms:
   if any(p[t]!=v for t,v in d['pins']):continue
   if any(not any(tuple(p[t:t+len(b)])==tuple(b) for t in range(9-len(b))) for b in d['blocks']):continue
   if any(not A[x][y] for x,y in zip(p,p[1:])):continue
   raw+=1;w=word(A,p);target=cl[owner[full^w]]
   if target in (1965,2234):found.append((target,p,w))
  check(raw==d['raw'],'all guarded paths count '+str(s))
  expected=sorted((t,tuple(r['path']),r['mask']) for t,r in d['rec'].items())
  check(sorted(found)==expected,'complete original U/V full-path menu '+str(s));fullmenu.append({'source':s,'raw_guarded_paths':raw,'retained_paths':len(found)})
  for t,p,w in found:
   check(cl[owner[w]]==s and cl[owner[full^w]]==t,'native source/target owner')
   check(word(list(map(list,zip(*A))),p[::-1])==reflection(w),'literal converse')
   order={v:i for i,v in enumerate(p)}
   K=[[0 if i==j else A[i][j] if abs(order[i]-order[j])==1 else A[j][i] for j in range(8)] for i in range(8)]
   check(word(K,p)==full^w,'literal K-partner')
   for i,j in allpairs:
    u,v=frame[i],frame[j];edge=(min(u,v),max(u,v));bit=K[i][j] if u<v else K[j][i]
    partner_options.setdefault(edge,{}).setdefault(s,set()).add(bit)
   for block in d['blocks']:check(any(tuple(p[::-1][i:i+len(block)])==tuple(block[::-1]) for i in range(9-len(block))),'converse moving blocks')
   orbit={w,reflection(w),full^w,full^reflection(w)}
   check(len(orbit)==4 and not orbit&allwords,'all reserved alternatives are disjoint');allwords|=orbit;resources[s,t]=orbit
   check((potentials[0][s]-potentials[0][t])%8==d['rec'][t]['phase'],'phase match')
  if s in core:check(all(not set(frame)&set(host['source_frames'][str(o)]) for o in data if o!=s),'every core source frame independent')
 for edge,by_source in partner_options.items():
  if len(by_source)>1:check(len(set.union(*by_source.values()))==1,'all alternative local K-partner edges glue')
 check(all(v[1965]==v[2234] for v in potentials),'literal common target potentials')
 check(len(allwords)==656,'164 reserved closed-orbit units')
 for s,t,clock in [(368,1965,1),(776,2234,4),(3139,2234,7)]:check(host['source_frames'][str(s)][data[s]['rec'][t]['path'][clock]]==4,'external owner4 fixed')
 # 16 cyclic bitsets of length128, indexed by the binary code.
 M=(1<<128)-1
 def shift(row,q):
  c,b=q
  return [((row[j^b]<<c)|(row[j^b]>>(128-c)))&M for j in range(16)]
 def support(ss):
  tab=[[0]*16 for _ in range(len(ss)+1)];tab[0][0]=1
  for n,s in enumerate(ss):
   for k in range(n+1,0,-1):tab[k]=[x|y for x,y in zip(tab[k],shift(tab[k-1],label(s)))]
  return tab
 def holes(row):return [16*c+b for c in range(128) for b in range(16) if not row[b]>>c&1]
 exceptions={p:[] for p in outside};mincov={p:{k:2048 for k in range(6,13)} for p in outside};grade6={};sixfreeze=[];sevenfreeze=[]
 for B in combinations(cheap,3):
  base=support([s for s in cheap if s not in B]);bad=[]
  for p in outside:
   allholes={}
   for k in range(6,13):
    row=[x|y for x,y in zip(base[k],shift(base[k-1],label(p)))];hh=holes(row);allholes[k]=hh
    mincov[p][k]=min(mincov[p][k],2048-len(hh))
    if k in (8,9,10):check(not hh,'independent direct grades8..10')
   check(bool(allholes[7])==bool(allholes[11]),'wide band complementary failure status')
   if allholes[7]:exceptions[p].append({'boundary':list(B),'grade7_missing':allholes[7],'grade11_missing':allholes[11]});bad.append(p)
   if allholes[6] and p not in grade6:grade6[p]={'boundary':list(B),'grade6_missing':allholes[6][0]}
  check(len(bad)<=2,'at least six wide repairs per triple')
  sixfreeze += [{'boundary':list(B),'surviving_pair':list(pair)} for pair in combinations(outside,2) if set(pair)<=set(bad)]
  sevenfreeze += [{'boundary':list(B),'survivor':p} for p in bad]
 native=json.loads((P/'native_receipt.json').read_text())['receiver_census']
 check({str(p):r for p,r in exceptions.items()}==native['exceptions'],'all independent exception tables exactly match')
 check([p for p in outside if not exceptions[p]]==[1532,1923,2087,2835,3015],'five universal repairs')
 check({str(p):r for p,r in grade6.items()}==native['grade6_failure_controls'],'independent grade6 failure controls')
 check(sixfreeze==[{'boundary':[1638,2071,2837],'surviving_pair':[1841,2553]}],'unique six frozen survivor failure')
 check(len(sevenfreeze)==4,'seven frozen failures')
 linked=plus(label(1841),label(1457));base=support([s for s in cheap if s not in {176,1425,2071,1457}]);internal={}
 for k in (8,9,10):internal[k]=holes([x|y for x,y in zip(base[k],shift(base[k-2],linked))])
 check(internal=={8:[477],9:[],10:[1098]},'receiver-internal actual equality holes')
 print('PASS independent full native menus and all1140x8 support census',flush=True)
 # Independent exhaustive grade / boundary arithmetic, all tables and all movable e.
 masks=list(range(1,256));gradecases=0;bridgecases=0;boundarycases=0
 for mask in masks:
  G=[b for b in range(8) if mask>>b&1];lo=7+max(b.bit_count() for b in G)
  for e in range(8):
   hi=11+e+min(b.bit_count() for b in G)
   for k in range(22+e):
    gradecases+=1
    if not lo<=k<=hi:
     check(2*min(k,21+e-k)<=18,'every out-of-band grade admits direct core bound');boundarycases+=1;continue
    intervals={b:[j for j in range(e+1) if 7<=k-b.bit_count()-j<=11] for b in G}
    check(all(intervals.values()),'every full boundary tuple at same total grade')
    for b,bp in combinations(G,2):
     common=set(intervals[b])&set(intervals[bp]);check(bool(common),'direct allowed tuple jump with common outside set')
     for j in common:
      r=k-b.bit_count()-j;s=k-bp.bit_count()-j;m=(b^bp).bit_count();d=min(r+s,36-r-s)
      check((r+s-m)%2==0 and d+m<=20,'parity-saving support')
      for c in (3,4):check(sum(sorted([2]*17+[c],reverse=True)[:d])+2*m<=38+c,'independent literal ternary bridge cost')
      bridgecases+=1
 # Generic clipped-interval identity by endpoint formula, larger range than producer.
 gaps=0
 for h in range(1,11):
  for l in range(h+1):
   for u in range(l,h+1):
    for e in range(8):
     for t,tp in combinations(range(8),2):
      for k in range(h+e+8):
       A=[j for j in range(e+1) if l<=k-t-j<=u];B=[j for j in range(e+1) if l<=k-tp-j<=u]
       if not A or not B:continue
       check(min(abs(x-y) for x in A for y in B)==max(0,tp-t-(u-l)),'exact clipped interval distance');gaps+=1
 # Genuine route controls from arbitrary noncanonical endpoints, including frozen V.
 rng=random.Random(481003);route_stats=Counter();route_examples=[];normalization=Counter();physical_checks=0
 def witnesses(ss):
  tabs=[{} for _ in range(len(ss)+1)];tabs[0][(0,0)]=frozenset()
  for i,s in enumerate(ss):
   for k in range(i+1,0,-1):
    for q,w in list(tabs[k-1].items()):tabs[k][plus(q,label(s))]=w|{s}
  return tabs
 def physical(C,H):
  nonlocal physical_checks
  active=C|H|{776,3139};rows={s:[0]*6 for s in set(data)|{1965,2234}};used=set()
  for s in data:
   t=2234 if s in active else 1965;w=resources[s,t]
   check(not used&w,'arrival normalized capacity');used|=w
   for j,v in enumerate(potentials):rows[s][j]+=2*(v[s]-v[t]);rows[t][j]-=2*(v[s]-v[t])
  check(len(used)==328,'arrival82 units');physical_checks+=1
  return rows
 Bs=[(176,1425,2071),(1638,2071,2837),(1425,2071,2837)]
 for Bt in Bs:
  B=set(Bt)
  for p in outside:
   if any(x['boundary']==list(Bt) for x in exceptions[p]):continue
   R=set(cheap)-B|{p};tab=witnesses(sorted(R));check(all(len(tab[k])==2048 for k in range(7,12)),'independent complete routing sections')
   rem=[s for s in outside if s!=p]
   for e in range(8):
    E=set(rem[:e]);frozen=set(rem[e:]);FV=frozenset(s for i,s in enumerate(sorted(frozen)) if i%2==0)
    for mask in (1,128,129,3,24,90,105,150,255):
     G=[frozenset(s for j,s in enumerate(Bt) if b>>j&1) for b in range(8) if mask>>b&1]
     lo=7+max(map(len,G));hi=11+e+min(map(len,G))
     for k in sorted({lo,hi,(lo+hi)//2}):
      b0=G[0];b1=G[-1];n=k-len(b0);r=rng.choice([r for r in range(19) if 0<=n-r<=e]);ec=set(rng.sample(sorted(E),n-r));rv=set(rng.sample(sorted(R),r));H=frozenset(rng.sample(helpers,2));start=(frozenset(FV|b0|ec|rv),H);q=charge(start[0]|H)
      def section(b,ec,h,r):
       fixed=FV|b|set(ec);return frozenset(fixed|tab[r][minus(q,charge(fixed|set(h)))]),frozenset(h)
      # Arbitrary actual receiver endpoint whenever possible, otherwise complete band.
      found=None
      for _ in range(20):
       h=frozenset(rng.sample(helpers,2));nn=k-len(b1);rr=rng.choice([r for r in range(19) if 0<=nn-r<=e]);ee=set(rng.sample(sorted(E),nn-rr));debt=minus(q,charge(FV|b1|ee|h))
       if debt in tab[rr]:found=section(b1,ee,h,rr);break
      if found is None:
       rr=max(7,k-len(b1)-e);ee=set(sorted(E)[:k-len(b1)-rr]);found=section(b1,ee,h,rr)
      end=found
      def migrate(x,y,b):
       target=k-len(b);rstar=max(7,target-e)
       def norm(z):
        C,H=z;ee=set(C)&E;rr=len(C&R);out=[z]
        if rr<7:
         normalization['below']+=1;ee-=set(sorted(ee)[:7-rr]);rr=7;out.append(section(b,ee,H,rr))
        if rr>11:
         normalization['above']+=1;ee|=set(sorted(E-ee)[:rr-11]);rr=11;out.append(section(b,ee,H,rr))
        while rr>rstar:
         ee.add(max(E-ee));rr-=1;out.append(section(b,ee,H,rr));normalization['consecutive']+=1
        return out
       left=norm(x);right=norm(y);out=left[:];ec=set(left[-1][0])&E;dest=set(right[-1][0])&E;hp=set(x[1])
       for v,w in zip(sorted(ec-dest),sorted(dest-ec)):
        ec=(ec-{v})|{w};out.append(section(b,ec,hp,rstar));normalization['outside']+=1
       for v,w in zip(sorted(hp-y[1]),sorted(y[1]-hp)):
        hp=(hp-{v})|{w};out.append(section(b,ec,hp,rstar));normalization['helper']+=1
       out+=list(reversed(right));return out
      j=max(0,k-len(b0)-11,k-len(b1)-11);ec=set(sorted(E)[:j]);a0=section(b0,ec,H,k-len(b0)-j);a1=section(b1,ec,H,k-len(b1)-j)
      route=migrate(start,a0,b0)+[a1]+migrate(a1,end,b1)[1:];route=[z for i,z in enumerate(route) if not i or z!=route[i-1]]
      check(route[0]==start and route[-1]==end,'restore exact noncanonical endpoints')
      baseline=physical(*start)
      for C,Hh in route:
       check(C&B in G and C&frozen==FV,'full ternary relation and frozen U/V records')
       check(len(C)-len(FV)==k and len(Hh)==2 and charge(C|Hh)==q,'residual grade and exact paired charge')
       rows=physical(C,Hh);check(all((rows[s][j]-baseline[s][j])%moduli[j]==0 for s in rows for j in range(6)),'all physical source and target rows including diagonals')
      for x,y in zip(route,route[1:]):
       diff=(x[0]^y[0])|(x[1]^y[1]);work=sum(costs[s] for s in diff);check(len(diff)<=20 and work<=49+costs[p],'actual interior support-clock bound')
       route_stats['events']+=1;route_stats['maximum_slots']=max(route_stats['maximum_slots'],len(diff));route_stats['maximum_clocks']=max(route_stats['maximum_clocks'],work)
      route_stats['routes']+=1;route_stats['arrivals']+=len(route);route_stats['frozen_V_routes']+=bool(FV)
      if Bt==(1638,2071,2837) and p==2087 and e==0 and mask==129 and k==10:
       route_examples.append({'boundary':list(Bt),'repair':p,'e':e,'frozen_V':sorted(FV),'k':k,'allowed_patterns':[0,7],'states':[{'core_V':sorted(C),'free_helper_V':sorted(Hh)} for C,Hh in route]})
 check(normalization['below'] and normalization['above'] and normalization['helper'] and route_stats['frozen_V_routes'],'all intended routing branches exercised')
 # Independent direct-boundary replacements, including arbitrary frozen V values.
 boundary_stats=Counter();pair_charge={}
 for hp in combinations(helpers,2):pair_charge.setdefault(charge(hp),[]).append(frozenset(hp))
 Bt=(1638,2071,2837);B=set(Bt);p=2087;remaining=[s for s in outside if s!=p]
 for e in range(8):
  E=set(remaining[:e]);frozen=set(remaining[e:]);FV=frozenset(s for i,s in enumerate(sorted(frozen)) if i%2==0);freecore=core-frozen;other=sorted(freecore-B)
  for mask in (1,128,129,3,24,90,105,150,255):
   G=[frozenset(s for j,s in enumerate(Bt) if b>>j&1) for b in range(8) if mask>>b&1]
   lo=7+max(map(len,G));hi=11+e+min(map(len,G))
   for k in range(22+e):
    if lo<=k<=hi:continue
    feasible=[b for b in G if 0<=k-len(b)<=len(other)]
    if not feasible:continue
    for repeat in range(2):
     b=rng.choice(feasible);C=FV|b|set(rng.sample(other,k-len(b)));H=frozenset(rng.sample(helpers,2));q=charge(C|H);D=Hh=None
     for attempt in range(100):
      bb=rng.choice(feasible);candidate=FV|bb|set(rng.sample(other,k-len(bb)));debt=minus(q,charge(candidate))
      if debt in pair_charge:D=candidate;Hh=rng.choice(pair_charge[debt]);break
     if D is None:D=C;Hh=H
     check(C&B in G and D&B in G and C&frozen==D&frozen==FV,'boundary literal relation and all frozen bits')
     check(len(C)==len(D)==k+len(FV) and charge(D|Hh)==q,'boundary residual grade and charge')
     lhs=physical(C,H);rhs=physical(D,Hh)
     check(all((lhs[s][j]-rhs[s][j])%moduli[j]==0 for s in lhs for j in range(6)),'direct boundary full physical rows')
     changed=(C^D)|(H^Hh);work=sum(costs[s] for s in changed)
     check(len(changed)<=22 and work<=76,'independent direct22/76 event bound')
     boundary_stats['events']+=1;boundary_stats['nonidentity']+=bool(changed);boundary_stats['frozen_V_events']+=bool(FV);boundary_stats['maximum_slots']=max(boundary_stats['maximum_slots'],len(changed));boundary_stats['maximum_clocks']=max(boundary_stats['maximum_clocks'],work)
 # Authenticate every producer literal route arrival and event independently.
 literal=json.loads((P/'native_receipt.json').read_text())['routing']['literal_route'];lbase=None
 for st in literal['states']:
  C=frozenset(st['core_V']);H=frozenset(st['free_helper_V']);check(len({2234 if s in C else 1965 for s in literal['boundary_sources']})==1,'actual native equality of targets')
  check(set(st['helper_V'])==H|{776,3139},'fixed external helper tree')
  rr=physical(C,H)
  if lbase is None:lbase=rr
  check(all((rr[s][j]-lbase[s][j])%moduli[j]==0 for s in rr for j in range(6)),'producer literal full physical rows')
 for x,y,event in zip(literal['states'],literal['states'][1:],literal['events']):
  diff=(set(x['core_V'])^set(y['core_V']))|(set(x['free_helper_V'])^set(y['free_helper_V']))
  check(sorted(diff)==event['changed_sources'] and len(diff)==event['support'] and sum(costs[s] for s in diff)==event['clock_occurrences'],'producer literal event independently measured')
 out={'status':'PASS independent audit controls; universal theorem checked separately in AUDIT.md','explicit_checks':COUNT,'manifest_files':len(manifest),'pinned_inputs':len(sourcepins),'full_path_menus':fullmenu,'all_native_costs':costs,'census':{'triples':1140,'repairs':outside,'minimum_coverages_grades6_to12':mincov,'exceptions':exceptions,'grade6_failure_controls':grade6,'six_frozen_certificate_failures':sixfreeze,'seven_frozen_certificate_failures':sevenfreeze},'receiver_internal_equality_holes':internal,'arithmetic':{'all_nonempty_tables':255,'grade_cases':gradecases,'direct_boundary_grade_cases':boundarycases,'bridge_common_count_cases':bridgecases,'generic_clipped_gap_cases':gaps},'independent_routes':dict(route_stats),'independent_direct_boundary_events':dict(boundary_stats),'normalization_branches':dict(normalization),'physical_arrivals':physical_checks,'frozen_V_equal_target_route':route_examples}
 if a.out:a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'status':'PASS','checks':COUNT,'routes':dict(route_stats),'normalization':dict(normalization)},sort_keys=True))
if __name__=='__main__':main()
