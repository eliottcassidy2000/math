#!/usr/bin/env python3
"""Portable exact native receiver census and relational full-record route replay."""
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import argparse,csv,gzip,json,sys,random
HERE=Path(__file__).resolve().parent;INPUT=HERE/'inputs';COUNT=0;U,V=1965,2234
def need(ok,why):
 global COUNT
 COUNT+=1
 if not ok:raise RuntimeError(why)
def add(a,b):return ((((a>>4)+(b>>4))%128)<<4)|((a^b)&15)
def sub(a,b):return ((((a>>4)-(b>>4))%128)<<4)|((a^b)&15)
def total(values):
 z=0
 for x in values:z=add(z,x)
 return z
def read(name):return json.loads((INPUT/name).read_text())
def verify_pins():
 pins=json.loads((HERE/'SOURCEPINS.json').read_text())
 need(set(pins)=={p.relative_to(HERE).as_posix() for p in INPUT.iterdir() if p.is_file()},'source pins equal complete input closure')
 for name,row in pins.items():need(sha256((HERE/name).read_bytes()).hexdigest()==row['sha256'],'source '+name)
def load_native():
 old=read('domains.json')['sources'];supply=read('supply.json');repair=read('source2125.json');host=read('host.json')
 nodes=list(csv.DictReader((INPUT/'merged_nodes.tsv').open(),delimiter='\t'))
 merged={int(z[k]):int(z['merged_id']) for z in nodes for k in ('class_a','class_b')}
 owners=array('H');owners.frombytes(gzip.decompress((INPUT/'owner_map.u16.gz').read_bytes()))
 if sys.byteorder!='little':owners.byteswap()
 electric=read('electrical.json')['order8'];X=electric['higher_layers'][-1]['lifts'][0];Y=electric['radical_lifts_mod4']
 labels={s:((X[s]%128)<<4)|sum((Y[j][s]%2)<<(j-1) for j in range(1,5)) for s in range(len(nodes))}
 data={}
 for row in old:
  rr={r['target']:r for r in row['exact_all_shared_clocks_target_domain']}
  data[row['source']]={'records':[rr[U],rr[V]],'fixed':row['fixed_clocks'],'blocks':[]}
  code=(row['code'][0]<<4)|sum(v<<j for j,v in enumerate(row['code'][1:]))
  need(labels[row['source']]==code,'old label independently rebuilt from original potentials')
 cores=sorted(data)
 for row in supply['helpers']:
  s=row['selected'];data[row['source']]={'records':[s['before'],s['after']],'fixed':s['common_clocks'],'blocks':[]}
 rr={r['target']:r for r in repair['records']};data[2125]={'records':[rr[U],rr[V]],'fixed':[],'blocks':repair['retained_blocks']}
 costs={s:sum(a!=b for a,b in zip(d['records'][0]['path'],d['records'][1]['path'])) for s,d in data.items()}
 cheap=sorted(s for s in cores if costs[s]==2);outside=sorted(set(cores)-set(cheap));free=sorted(set(data)-set(cores)-{368,776,3139})
 need((len(data),len(cheap),len(outside),len(free))==(82,20,8,51),'exact native inventory')
 for source in cores:
  frame=set(host['source_frames'][str(source)])
  need(all(not frame & set(host['source_frames'][str(other)]) for other in data if other!=source),'every original core frame is private from every other selected source')
 pairs=list(combinations(range(8),2));chords=[(i,j) for i in range(8) for j in range(i+2,8)];full=(1<<21)-1
 reflect=lambda m:sum(((m>>k)&1)<<chords.index((7-j,7-i)) for k,(i,j) in enumerate(chords))
 word=lambda A,p:sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(chords))
 literal={};all_words=set();prescribed={}
 for s,d in data.items():
  A=[[0]*8 for _ in range(8)];code=int(nodes[s]['canon_a'])
  for k,(i,j) in enumerate(pairs):A[i][j]=(code>>(27-k))&1;A[j][i]=1-A[i][j]
  em=host['source_frames'][str(s)]
  for i,j in pairs:
   a,b=em[i],em[j];edge=min(a,b),max(a,b);value=A[i][j] if a<b else A[j][i]
   need(edge not in prescribed or prescribed[edge]==value,'coherent source host');prescribed[edge]=value
  literal[s]=[]
  for bit,r in enumerate(d['records']):
   p,m=r['path'],r['mask'];at={v:i for i,v in enumerate(p)}
   partner=[[0 if i==j else A[i][j] if abs(at[i]-at[j])==1 else A[j][i] for j in range(8)] for i in range(8)]
   need(all(A[a][b] for a,b in zip(p,p[1:])),'literal directed native path')
   need(word(A,p)==m and merged[owners[m]]==s and merged[owners[full^m]]==(U,V)[bit],'literal source target and normalized word')
   need(word(partner,p)==full^m,'literal K partner')
   need(word([list(z) for z in zip(*A)],p[::-1])==reflect(m),'literal converse')
   need(all(p[t]==v for t,v in d['fixed']),'all inherited indexed pins')
   for block in d['blocks']:need(any(p[i:i+len(block)]==block for i in range(9-len(block))),'2125 block retained')
   words={m,reflect(m),full^m,full^reflect(m)}
   need(len(words)==4 and not words&all_words,'distinct normalized resources, not kernel multiplicities')
   all_words|=words;literal[s].append({'words':words,'global_path':[em[v] for v in p]})
   need((X[s]-X[(U,V)[bit]])%8==r['phase'],'original phase')
 need(len(all_words)==656 and host['whole_source_host_order']==652,'fixed reservoir and inherited host')
 need([literal[368][0]['global_path'][1],literal[776][1]['global_path'][4],literal[3139][1]['global_path'][7]]==[4,4,4],'immutable external owner4')
 need(sum(sorted((costs[s] for s in cores),reverse=True)[:18])==47,'top18 core clocks47')
 need(sum(sorted((costs[s] for s in free),reverse=True)[:4])==29,'top4 free helper clocks29')
 return dict(data=data,cores=cores,cheap=cheap,outside=outside,free=free,costs=costs,labels=labels,literal=literal,potentials=[X]+Y[:5],moduli=[256]+[4]*5)
def census(n):
 cheap,outside,labels=n['cheap'],n['outside'],n['labels']
 d=[[0]*2048 for _ in range(21)];d[0][0]=1
 for i,s in enumerate(cheap):
  for k in range(i+1,0,-1):
   for q,v in enumerate(d[k-1]):
    if v:d[k][add(q,labels[s])]+=v
 perms={s:[sub(q,labels[s]) for q in range(2048)] for s in cheap+outside}
 def remove(table,s):
  out=[table[0][:]];p=perms[s]
  for row in table[1:-1]:
   previous=out[-1];out.append([v-previous[p[q]] for q,v in enumerate(row)])
  return out
 exceptions={p:[] for p in outside};grade6_failure={};minimum={p:{r:2048 for r in (6,7,8,9)} for p in outside};triples=0
 for B in combinations(cheap,3):
  table=d
  for s in B:table=remove(table,s)
  need(len(table)==18 and all(sum(row)==comb(17,k) and min(row)>=0 for k,row in enumerate(table)),'exact nonnegative whole-slot deletion counts')
  remaining=[s for s in cheap if s not in B]
  for p in outside:
   coverage={};holes={};pp=perms[p]
   for k in (6,7,8,9):
    row,prev=table[k],table[k-1]
    missing=[q for q in range(2048) if row[q]+prev[pp[q]]==0]
    holes[k]=missing;coverage[k]=2048-len(missing);minimum[p][k]=min(minimum[p][k],coverage[k])
   need(coverage[8]==coverage[9]==2048,'all8 repairs preserve8..10')
   if holes[6] and p not in grade6_failure:grade6_failure[p]={'boundary':list(B),'grade6_missing':holes[6][0]}
   if holes[7]:
    S=total(labels[s] for s in remaining+[p])
    exceptions[p].append({'boundary':list(B),'grade7_missing':holes[7],'grade11_missing':sorted(sub(S,q) for q in holes[7])})
   need(sum(table[8][q]+table[7][pp[q]] for q in range(2048))==comb(18,8),'repaired exact population count')
  triples+=1
 need(triples==1140,'complete1140 triples')
 need(set(grade6_failure)==set(outside),'each fixed repair has a witnessed grade6 failure')
 universal=[p for p in outside if not exceptions[p]]
 need(universal==[1532,1923,2087,2835,3015],'five universal7..11 repairs')
 need([len(exceptions[p]) for p in outside]==[0,1,0,0,2,0,0,1],'all four narrow-band exceptions')
 all_bad=Counter(tuple(row['boundary']) for rows in exceptions.values() for row in rows)
 need(max(all_bad.values())==2,'at most two bad repairs per boundary')
 fresh_checks=[]
 for B,p in [(tuple(r['boundary']),p) for p,rows in exceptions.items() for r in rows]+[((1638,2071,2837),2087)]:
  ss=[s for s in cheap if s not in B]+[p];tabs=[set() for _ in range(19)];tabs[0]={0}
  for i,s in enumerate(ss):
   for k in range(i+1,0,-1):tabs[k]|={add(q,labels[s]) for q in tabs[k-1]}
  bad=next((r for r in exceptions[p] if tuple(r['boundary'])==B),None)
  need([len(tabs[k]) for k in (8,9,10)]==[2048]*3,'fresh enumeration of old band')
  if bad:need(sorted(set(range(2048))-tabs[7])==bad['grade7_missing'] and sorted(set(range(2048))-tabs[11])==bad['grade11_missing'],'fresh actual complementary holes')
  else:need(all(len(tabs[k])==2048 for k in range(7,12)),'fresh universal-source control')
  fresh_checks.append({'boundary':list(B),'repair':p,'coverage':[len(t) for t in tabs]})
 return {'triple_count':triples,'repair_sources':outside,'universal_wide_repairs':universal,'exceptions':exceptions,'minimum_coverage_grades6_to9':minimum,'grade6_failure_controls':grade6_failure,'maximum_bad_repairs_at_one_boundary':2,'fresh_enumeration_controls':fresh_checks}
def make_sections(n,B,p):
 R=sorted((set(n['cheap'])-set(B))|{p});tabs=[{} for _ in range(19)];tabs[0][0]=frozenset()
 for i,s in enumerate(R):
  for k in range(i+1,0,-1):
   for q,w in list(tabs[k-1].items()):tabs[k].setdefault(add(q,n['labels'][s]),w|{s})
 need(all(len(tabs[k])==2048 for k in range(7,12)),'actual wide receiver witnesses')
 return set(R),tabs
def route_tests(n):
 labels,costs=n['labels'],n['costs'];charge=lambda ss:total(labels[s] for s in ss)
 helper_pairs=list(combinations(n['free'],2));stats={'routes':0,'events':0,'arrivals':0,'maximum_slots':0,'maximum_clocks':0};literal_example=None
 def physical(core,hp):
  active=set(core)|set(hp)|{776,3139};rows={s:[0]*6 for s in list(n['data'])+[U,V]};used=set()
  for s in n['data']:
   bit=int(s in active);words=n['literal'][s][bit]['words']
   need(not used&words,'arrival capacity one');used|=words;t=(U,V)[bit]
   for j,f in enumerate(n['potentials']):
    v=2*(f[s]-f[t]);rows[s][j]+=v;rows[t][j]-=v
  need(len(used)==328 and 368 not in active,'82 occupied units and external tree')
  return rows
 for B in ((176,1425,2071),(1638,2071,2837)):
  p=2087;R,tabs=make_sections(n,B,p);B=set(B)
  for e in (0,7):
   E=set(sorted(set(n['outside'])-{p})[:e]);frozen=set(n['outside'])-{p}-E
   for pattern in ((0,7),(1,2,7),(0,3,5,6),tuple(range(8))):
    bitset=lambda b:frozenset(s for j,s in enumerate(sorted(B)) if (b>>j)&1)
    allowed=[bitset(b) for b in pattern];tmin=min(map(len,allowed));tmax=max(map(len,allowed));lo,hi=7+tmax,11+e+tmin
    for k in range(lo,hi+1):
     for q in (0,1,127,305,1024,1752,2047):
      def section(b,ec,hp,r):
       debt=sub(q,charge(set(b)|set(ec)|set(hp)))
       return frozenset(set(b)|set(ec)|set(tabs[r][debt])),frozenset(hp)
      def endpoint(b,hp,last=False):
       js=list(range(max(0,k-len(b)-11),min(e,k-len(b)-7)+1));j=js[-1 if last else 0]
       ec=set(sorted(E)[-j:]) if j and last else set(sorted(E)[:j])
       return section(b,ec,hp,k-len(b)-j)
      start=endpoint(allowed[0],helper_pairs[0]);end=endpoint(allowed[-1],helper_pairs[-1],True)
      def within(a,z,b):
       target=k-len(b);r0=max(7,target-e)
       def normalize(state):
        cc,hp=state;ec=set(cc)&E;r=len(set(cc)&R);out=[state]
        if r<7:ec-=set(sorted(ec)[:7-r]);r=7;out.append(section(b,ec,hp,r))
        elif r>11:ec|=set(sorted(E-ec)[:r-11]);r=11;out.append(section(b,ec,hp,r))
        while r>r0:
         ec.add(min(E-ec));r-=1;out.append(section(b,ec,hp,r))
        return out
       left,right=normalize(a),normalize(z);out=left[:];ec=set(left[-1][0])&E;dest=set(right[-1][0])&E;hp=set(a[1])
       for x,y in zip(sorted(ec-dest),sorted(dest-ec)):
        ec.remove(x);ec.add(y);out.append(section(b,ec,hp,r0))
       for x,y in zip(sorted(hp-set(z[1])),sorted(set(z[1])-hp)):
        hp.remove(x);hp.add(y);out.append(section(b,ec,hp,r0))
       out.append(right[-1]);out.extend(right[-2::-1]);return out
      b0,b1=allowed[0],allowed[-1];jlo=max(0,k-len(b0)-11,k-len(b1)-11);jhi=min(e,k-len(b0)-7,k-len(b1)-7)
      need(jlo<=jhi,'width4 admits common actual outside count')
      ec=set(sorted(E)[:jlo]);bridge0=section(b0,ec,start[1],k-len(b0)-jlo);bridge1=section(b1,ec,start[1],k-len(b1)-jlo)
      path=within(start,bridge0,b0)+[bridge1]+within(bridge1,end,b1)[1:]
      path=[x for i,x in enumerate(path) if i==0 or x!=path[i-1]]
      need(path[0]==start and path[-1]==end,'literal endpoint restoration')
      baseline=physical(*path[0]);events=[]
      for core,hp in path:
       need(set(core)&B in allowed and not set(core)&frozen,'actual ternary relation and frozen outside records')
       need(len(core)==k and len(hp)==2 and charge(set(core)|set(hp))==q,'both margins and full paired charge')
       current=physical(core,hp)
       need(all((current[s][j]-baseline[s][j])%n['moduli'][j]==0 for s in current for j in range(6)),'every physical source/target dyadic generating row')
       stats['arrivals']+=1
      for x,y in zip(path,path[1:]):
       changed=(set(x[0])^set(y[0]))|(set(x[1])^set(y[1]));clock=sum(costs[s] for s in changed)
       need(len(changed)<=20 and clock<=52,'native20slot52clock relational bound')
       events.append({'changed_sources':sorted(changed),'support':len(changed),'clock_occurrences':clock})
       stats['maximum_slots']=max(stats['maximum_slots'],len(changed));stats['maximum_clocks']=max(stats['maximum_clocks'],clock);stats['events']+=1
      stats['routes']+=1
      if literal_example is None and B=={1638,2071,2837} and e==0 and pattern==(0,7) and k==10 and q==0:
       literal_example={'boundary_sources':sorted(B),'allowed_bit_patterns':list(pattern),'repair_source':p,'frozen_U_sources':sorted(frozen),'core_grade':k,'residual_charge':q,'states':[{'core_V':sorted(c),'free_helper_V':sorted(h),'helper_V':sorted(set(h)|{776,3139})} for c,h in path],'events':events}
 need(literal_example is not None,'genuinely relational frozen-outside route retained')
 # Boundary-grade replay uses direct accepted endpoints. No complete section
 # outside the proved interior is assumed.
 rng=random.Random(48024);B={1638,2071,2837};other=sorted(set(n['cores'])-B)
 pair_by_q={}
 for hp in helper_pairs:pair_by_q.setdefault(charge(hp),[]).append(hp)
 boundary={'cases':0,'nonidentity':0,'maximum_slots':0,'maximum_clocks':0}
 for k in list(range(10))+list(range(19,29)):
  weights=[t for t in (0,3) if 0<=k-t<=25]
  for repeat in range(8):
   t=rng.choice(weights);ca=set(rng.sample(other,k-t))|(B if t else set());ha=set(rng.choice(helper_pairs));q=charge(ca|ha)
   for attempt in range(10000):
    tt=rng.choice(weights);cb=set(rng.sample(other,k-tt))|(B if tt else set());debt=sub(q,charge(cb))
    if debt in pair_by_q:hb=set(rng.choice(pair_by_q[debt]));break
   else:raise RuntimeError('no matched boundary endpoint found')
   need(len(ca)==len(cb)==k and ca&B in (set(),B) and cb&B in (set(),B),'boundary endpoints keep full ternary relation and grade')
   pa,pb=physical(ca,ha),physical(cb,hb)
   need(all((pa[z][j]-pb[z][j])%n['moduli'][j]==0 for z in pa for j in range(6)),'direct boundary all physical rows')
   changed=(ca^cb)|(ha^hb);work=sum(costs[z] for z in changed)
   need(len(changed)<=22 and work<=76,'native global22slot76clock endpoint event')
   boundary['cases']+=1;boundary['nonidentity']+=bool(changed);boundary['maximum_slots']=max(boundary['maximum_slots'],len(changed));boundary['maximum_clocks']=max(boundary['maximum_clocks'],work)
 return {'summary':stats,'boundary_controls':boundary,'literal_route':literal_example}
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
 verify_pins();n=load_native();tables=census(n);print('receiver census complete',flush=True);routes=route_tests(n)
 out={'status':'FINITE-EXACT native census and controls; universal quantifiers proved in THEOREM.md','receiver_census':tables,'routing':routes,'guard_scope':'All original aggregate nodewise dyadic currents by pinned paired-code equivalence; not individual marked-port phases, integer or odd-primary currents.','inheritance':'Full menu/kernel classification inherited from source-pinned domains, helper supply and XOR-tree certificate. Selected paths, converse, capacities, host restrictions and physical rows independently replayed.','checks':COUNT}
 out=json.loads(json.dumps(out,sort_keys=True));target=HERE/'native_receipt.json'
 if args.write:target.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 else:need(out==json.loads(target.read_text()),'frozen native receipt matches')
 print('PASS',COUNT,'checks;',routes['summary'])
if __name__=='__main__':main()
