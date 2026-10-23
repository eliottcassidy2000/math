#!/usr/bin/env python3
"""Construct and verify genuine nonuniqueness witnesses for four n8 bucket relaxations.
The counterfeits are relaxed transport tables, not realizable tournament wiring.
"""
import argparse,array,csv,gzip,hashlib,importlib.util,json,sys
from collections import Counter,defaultdict
from pathlib import Path
HERE=Path(__file__).resolve().parent
DEFAULT_INCIDENCE=HERE.parent/'session38-incidence'
DEFAULT_WIRING=HERE.parent/'session38-wiring'
N=8;FULL=(1<<21)-1;TILES=[(b,a) for b in range(6) for a in range(b+2,8)]
REF=[TILES.index((7-a,7-b)) for b,a in TILES]
STARTBITS=sum(1<<k for k,(b,a) in enumerate(TILES) if a==7)
ENDBITS=sum(1<<k for k,(b,a) in enumerate(TILES) if b==0)
checks=Counter()
def need(q,msg):
 checks[msg]+=1
 if not q:raise ValueError(msg)
def read_tsv(path):return list(csv.DictReader(gzip.open(path,'rt'),delimiter='\t'))
def clean(C):return Counter({k:v for k,v in C.items() if v})
def main():
 p=argparse.ArgumentParser();p.add_argument('--incidence',type=Path,default=DEFAULT_INCIDENCE);p.add_argument('--wiring',type=Path,default=DEFAULT_WIRING);p.add_argument('--write',action='store_true');args=p.parse_args()
 inc=args.incidence;wiring=args.wiring
 cls=read_tsv(wiring/'classes.tsv.gz');node={int(r['class_id']):int(r['node_id']) for r in cls};conv={int(r['class_id']):int(r['converse_class']) for r in cls};canon={int(r['class_id']):int(r['canon']) for r in cls};c3={int(r['class_id']):int(r['c3']) for r in cls}
 spec=importlib.util.spec_from_file_location('incidence_checker',inc/'check.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 sigcache={}
 def sig(ci,name):
  if ci not in sigcache:
   D,(c,E,v4,out,M)=mod.get_signatures(mod.matrix(canon[ci],8));D['E33']=mod.esig(c,E,out,3,3);sigcache[ci]=D
  return sigcache[ci][name]
 native=Counter()
 for r in read_tsv(wiring/'ports.tsv.gz'):
  key=tuple(int(r[k]) for k in ('source_class','target_class','start_score','end_score','blue'));native[key]=int(r['multiplicity'])
 need(sum(native.values())==1<<20,'native full transport mass')
 declared_open={}
 for rr in read_tsv(wiring/'endpoint_capacities.tsv.gz'):
  z=int(rr['open'])
  if z:declared_open[int(rr['class_id']),int(rr['start_score']),int(rr['end_score']),int(rr['blue'] in ('True','true','1'))]=z
 actual_open=Counter()
 for(i,j,a,b,blue),z in native.items():actual_open[j,8-a,6-b,blue]+=z
 need(dict(actual_open)==declared_open,'complete native open capacity ledger')
 def mirror(k):
  i,j,a,b,blue=k;return(conv[i],conv[j],7-b,7-a,blue)
 need(all(native[mirror(k)]==z for k,z in native.items()),'native full converse reflection')
 need(all(k[4] or k!=mirror(k) or z%2==0 for k,z in native.items()),'native reflection-fixed black cells even')
 merged=Counter();undirected=Counter()
 for(i,j,a,b,blue),z in native.items():merged[node[i],node[j],blue]+=z;undirected[min(node[i],node[j]),max(node[i],node[j]),blue]+=z
 need(all(blue or z%2==0 for(u,v,blue),z in merged.items()),'native literal black twins')
 def projections(delta):
  d=Counter();u=Counter()
  for(i,j,a,b,blue),z in delta.items():d[node[i],node[j],blue]+=z;u[min(node[i],node[j]),max(node[i],node[j]),blue]+=z
  return clean(d),clean(u)
 def constraints(delta,name,require_mirror):
  row=Counter();col=Counter()
  for(i,j,a,b,blue),z in delta.items():
   if native[i,j,a,b,blue]+z<0:return False
   row[i,a,b,blue,sig(j,name)]+=z;col[j,8-a,6-b,blue]+=z
  if any(row.values()) or any(col.values()):return False
  if require_mirror and any(delta[mirror(k)]!=v for k,v in delta.items()):return False
  d,u=projections(delta)
  if not u:return False
  if require_mirror and any(not blue and (merged[k]+z)%2 for k,z in d.items() for blue in [k[2]]):return False
  return True
 def serialize(delta,name,require_mirror,cycle):
  need(constraints(delta,name,require_mirror),'constructed counterfeit feasibility')
  d,u=projections(delta);need(bool(d) and bool(u),'nonzero directed and undirected projection')
  # Check each increment against an actual source row with the same target signature,
  # and an actual target endpoint capacity. This is not inferred from class counts.
  supportrow=set((i,a,b,blue,sig(j,name)) for(i,j,a,b,blue),z in native.items() if j in {k[1] for k in delta})
  targetcaps=Counter()
  for(i,j,a,b,blue),z in native.items():
   if j in {k[1] for k in delta}:targetcaps[j,8-a,6-b,blue]+=z
  for(i,j,a,b,blue),z in delta.items():
   need((i,a,b,blue,sig(j,name)) in supportrow,'literal source target-signature bucket exists')
   need(targetcaps[j,8-a,6-b,blue]>0,'literal target endpoint column exists')
   need(blue==0,'blue table unchanged')
  need(sum(delta.values())==0,'total line mass preserved')
  need(all(c3[j]-c3[i]==a-b-1 for i,j,a,b,blue in delta),'every changed cell obeys triangle law')
  dcurrent=Counter()
  for(i,j,a,b,blue),z in delta.items():
   dcurrent[node[i]]+=z*(c3[j]-c3[i]);dcurrent[node[j]]+=z*(c3[i]-c3[j])
  need(not any(dcurrent.values()),'merged triangle current unchanged')
  need(all(native[k]+v>=0 for k,v in delta.items()),'all changed cell bounds')
  if require_mirror:
   need(all(delta[mirror(k)]==v for k,v in delta.items()),'changed table full converse reflection')
   need(all(k[4] or k!=mirror(k) or (native[k]+z)%2==0 for k,z in delta.items()),'changed reflection-fixed black cells even')
   need(all(blue or (merged[key]+v)%2==0 for key,v in d.items() for blue in [key[2]]),'all changed black counts even')
  return {'model':'converse-symmetric and black-even' if require_mirror else 'ordinary endpoint/color/signature transport','cycle':cycle,'edits':[{'source_class':i,'target_class':j,'start_score':a,'end_score':b,'blue':blue,'before':native[i,j,a,b,blue],'delta':z,'after':native[i,j,a,b,blue]+z,'source_canon':canon[i],'target_canon':canon[j]} for(i,j,a,b,blue),z in sorted(delta.items())],'directed_merged_delta':[[*k,v] for k,v in sorted(d.items())],'undirected_merged_delta':[[*k,v] for k,v in sorted(u.items())],'all_source_signature_buckets_preserved':True,'all_target_endpoint_capacities_preserved':True,'blue_unchanged':True,'triangle_current_preserved':True,'converse_reflection_preserved':all(delta[mirror(k)]==v for k,v in delta.items()),'all_literal_black_multiplicities_even':all(blue or (merged[key]+v)%2==0 for key,v in d.items() for blue in [key[2]])}
 results={};deltas={};bases={};rejections=Counter();inputpaths=[wiring/'classes.tsv.gz',wiring/'ports.tsv.gz',wiring/'endpoint_capacities.tsv.gz',wiring/'class_of_mask.u16le.gz',inc/'check.py']
 for name in ('J4','E_edge2','E_neighbor2','E_neighbor3'):
  path=inc/f'residual_{name}.json';inputpaths.append(path);r=json.loads(path.read_text());rows={i:(key,z) for i,key,z in r['row_keys']};cols={j:(key,z) for j,key,z in r['column_keys']};cells={(i,j):z for i,j,z in r['allowed_entries']};adj=defaultdict(set)
  for(i,j),z in cells.items():adj[i].add(j)
  groups=defaultdict(list)
  for i,(key,z) in rows.items():groups[tuple(sorted(adj[i]))].append(i)
  # Independently authenticate every surviving residual cell and its capacity.
  expected={}
  for(i,j),z in cells.items():
   if z:expected[i]=sig(cols[j][0][0],name)
  for(i,j),z in cells.items():
   ci,a,b,blue,gid=rows[i][0];cj,ta,tb,tblue=cols[j][0]
   need((ta,tb,tblue)==(8-a,6-b,blue),'residual endpoint compatibility')
   need(sig(cj,name)==expected[i],'every residual cell target signature accepted')
   need(c3[cj]-c3[ci]==a-b-1,'every residual cell triangle law')
   need(native[ci,cj,a,b,blue]==z,'residual cell actual value')
  need(all(sum(z for(ii,j),z in cells.items() if ii==i)==cap for i,(key,cap) in rows.items()),'residual row capacities')
  need(all(sum(z for(i,jj),z in cells.items() if jj==j)==cap for j,(key,cap) in cols.items()),'residual column capacities')
  ordinary=None;symmetric=None;tested=0
  for coltuple,ris in sorted(groups.items()):
   if len(coltuple)<2:continue
   for c1 in coltuple:
    for c2 in coltuple:
     if c1>=c2:continue
     target1,target2=cols[c1][0][0],cols[c2][0][0]
     need(sig(target1,name)==sig(target2,name),'accepted target signatures truly equal')
     for r1 in ris:
      if not cells[r1,c1]:continue
      for r2 in ris:
       if r1==r2 or not cells[r2,c2]:continue
       source1,a,b,blue,_=rows[r1][0];source2,aa,bb,blue2,_=rows[r2][0]
       if blue or blue2:continue
       need((a,b,blue)==(aa,bb,blue2),'four-cycle endpoint and color agreement')
       delta=Counter({(source1,target1,a,b,blue):-1,(source1,target2,a,b,blue):1,(source2,target2,a,b,blue):-1,(source2,target1,a,b,blue):1});delta=clean(delta);tested+=1
       cycle=[['row',r1],['column',c2],['row',r2],['column',c1],['row',r1]]
       if ordinary is None and constraints(delta,name,False):ordinary=(delta,cycle)
       sym=delta.copy()
       for key,v in delta.items():sym[mirror(key)]+=v
       sym=clean(sym)
       if symmetric is None and constraints(sym,name,True):symmetric=(sym,cycle,delta)
       if ordinary is not None and symmetric is not None:break
      if ordinary is not None and symmetric is not None:break
     if ordinary is not None and symmetric is not None:break
    if ordinary is not None and symmetric is not None:break
   if ordinary is not None and symmetric is not None:break
  need(ordinary is not None,name+' ordinary counterfeit found')
  results[name]={'ordinary':serialize(ordinary[0],name,False,ordinary[1]),'candidate_four_cycles_tested':tested}
  deltas[name,'ordinary']=ordinary[0];bases[name,'ordinary']=ordinary[0]
  if symmetric is not None:
   results[name]['symmetric']=serialize(symmetric[0],name,True,symmetric[1]);deltas[name,'symmetric']=symmetric[0];bases[name,'symmetric']=symmetric[2]
   print(name,'nonunique in both models;',len(symmetric[0]),'symmetric changed cells; merged delta',results[name]['symmetric']['directed_merged_delta'],flush=True)
  else:
   need(name=='E_neighbor3','only intended symmetric uniqueness case')
   groups=defaultdict(list)
   for ci in canon:groups[sig(ci,name)].append(ci)
   bad=[ci for ci in canon if len(groups[sig(ci,name)])>1]
   need(all(len(groups[sig(conv[ci],name)])==1 for ci in bad),'every colliding class has singleton converse signature')
   need(all(len(groups[sig(ci,name)])==1 or len(groups[sig(conv[ci],name)])==1 for ci in canon),'singleton-or-converse theorem covers every target class')
   attempted=ordinary[0].copy()
   for k,z in ordinary[0].items():attempted[mirror(k)]+=z
   attempted=clean(attempted);rowchange=Counter()
   for(i,j,a,b,blue),z in attempted.items():rowchange[i,a,b,blue,sig(j,name)]+=z
   need(any(rowchange.values()),'naive reflected counterfeit violates retained source signature buckets')
   violations=[]
   for (i,a,b,blue,ss),z in rowchange.items():
    if z:violations.append({'source_class':i,'start_score':a,'end_score':b,'blue':blue,'target_classes_in_signature':groups[ss],'delta':z})
   results[name]['symmetric_uniqueness']={'status':'entire ordinary target-class transport table uniquely forced, hence merged wiring unique','criterion':'every class has a globally singleton S signature or a converse with globally singleton S signature','ordinary_classes':len(canon),'direct_singletons':len(canon)-len(bad),'additional_classes_forced_by_converse':len(bad),'colliding_classes':[{'class_id':i,'canon':canon[i],'converse_class':conv[i],'converse_canon':canon[conv[i]],'converse_signature_class_count':len(groups[sig(conv[i],name)])} for i in bad],'rejected_mirrored_switch_row_deltas':violations}
   print(name,'ordinary nonunique, but full symmetric table UNIQUE by singleton-or-converse certificate:',len(canon)-len(bad),'+',len(bad),'classes',flush=True)
 # Independent E_(3,3) recomputation: the same target partition as E_neighbor3.
 e33=defaultdict(list);exactedge=defaultdict(list)
 for ci in canon:e33[sig(ci,'E33')].append(ci);exactedge[sig(ci,'E_neighbor3')].append(ci)
 partition33=sorted(sorted(g) for g in e33.values());partition_exact=sorted(sorted(g) for g in exactedge.values())
 need(partition33==partition_exact,'E33 and E_neighbor3 class partitions identical')
 need(all(len(e33[sig(ci,'E33')])==1 or len(e33[sig(conv[ci],'E33')])==1 for ci in canon),'E33 singleton-or-converse theorem covers every class')
 need(constraints(deltas['E_neighbor3','ordinary'],'E33',False),'E33 ordinary counterfeit valid')
 e33certificate={'status':'ordinary table nonunique; converse-symmetric table unique','distinct_class_signatures':len(e33),'same_partition_as':'E_neighbor3','partition_sha256':hashlib.sha256(json.dumps(partition33,separators=(',',':')).encode()).hexdigest(),'singleton_or_converse_verified_on_all_classes':len(canon),'symbol_alphabet_size':9,'ordinary_counterfeit':'E_neighbor3 ordinary edits apply unchanged','minimality_scope':'within rectangular (c_target mod p, cyclic_edge_count mod q) alphabets with full converse constraints, p>=3 and q>=3 are necessary by the E_neighbor2 and E_edge2 symmetric counterfeits, and p=q=3 attains9 symbols; this is not a claim that all larger p,q succeed'}
 # Literal mask witnesses for every native changed cell, plus direct verification
 # that reflection realizes the claimed class/endpoint mirror.
 needed={k for name,ds in deltas.items() for k in ds if native[k]>0};witness={}
 mapping=array.array('H',gzip.decompress((wiring/'class_of_mask.u16le.gz').read_bytes()))
 if sys.byteorder!='little':mapping.byteswap()
 need(len(mapping)==1<<21,'literal map coverage')
 for mask in range(FULL+1):
  if not(mask>>5&1):continue
  a=7-(mask&STARTBITS).bit_count();b=(mask&ENDBITS).bit_count();mir=sum(((mask>>k)&1)<<REF[k] for k in range(21));key=(mapping[mask],mapping[mask^FULL],a,b,int(mask==mir))
  if key in needed:
   witness.setdefault(key,[])
   if len(witness[key])<max(1,-min(ds.get(key,0) for ds in deltas.values())):witness[key].append(mask)
  if all(len(witness.get(k,[]))>=max(1,-min(ds.get(k,0) for ds in deltas.values())) for k in needed):break
 need(set(witness)==needed,'literal witnesses for all nonzero native changed cells')
 for name,rs in results.items():
  for mode in ('ordinary','symmetric'):
   if mode not in rs:continue
   for e in rs[mode]['edits']:
    key=tuple(e[k] for k in ('source_class','target_class','start_score','end_score','blue'))
    e['native_masks']=witness.get(key,[])
    for mask in e['native_masks']:
     mir=sum(((mask>>k)&1)<<REF[k] for k in range(21));need(mapping[mir]==conv[key[0]] and mapping[mir^FULL]==conv[key[1]],'literal witness reflection classes')
    if e['delta']<0:need(len(e['native_masks'])>=-e['delta'],'distinct literal removed-port witnesses')
 # Exhibit the changes as literal alternate bijections of actual endpoint ports.
 def reflect(mask):return sum(((mask>>k)&1)<<REF[k] for k in range(21))
 for(name,mode),base in bases.items():
  removed=[k for k,z in base.items() if z<0]
  need(len(removed)==2 and all(base[k]==-1 for k in removed),'base cycle has two unit donors')
  m1,m2=(witness[k][0] for k in removed)
  pairing={m1:m2^FULL,m2:m1^FULL}
  if mode=='symmetric':
   for m,t in list(pairing.items()):
    need(reflect(m) not in pairing,'four distinct reflected source ports');pairing[reflect(m)]=reflect(t)
  need(set(pairing.values())=={m^FULL for m in pairing},'literal target ports permuted bijectively')
  observed=Counter();literal=[]
  for m,t in sorted(pairing.items()):
   i,jold,jnew=mapping[m],mapping[m^FULL],mapping[t];a=7-(m&STARTBITS).bit_count();b=(m&ENDBITS).bit_count()
   need(bool(m>>5&1) and not(t>>5&1),'literal closed to open matching')
   need(m!=reflect(m) and t!=reflect(t),'literal matching black color')
   need(sig(jold,name)==sig(jnew,name),'literal source target signature unchanged')
   need((7-(t&STARTBITS).bit_count(),(t&ENDBITS).bit_count())==(8-a,6-b),'literal target endpoint compatibility')
   need(t!=(m^FULL),'counterfeit differs from actual complement')
   observed[i,jold,a,b,0]-=1;observed[i,jnew,a,b,0]+=1
   literal.append({'source_closed_mask':m,'native_target_open_mask':m^FULL,'counterfeit_target_open_mask':t,'source_class':i,'native_target_class':jold,'counterfeit_target_class':jnew})
  need(clean(observed)==deltas[name,mode],'literal switch gives exact recorded table delta')
  if mode=='symmetric':need(all(pairing[reflect(m)]==reflect(t) for m,t in pairing.items()),'literal matching commutes with reflection')
  results[name][mode]['literal_port_matching']=literal
 out={'scope':'FINITE-EXACT n8: all four menus admit explicit alternate ordinary endpoint/color/signature transport. J4, E_edge2 and E_neighbor2 remain nonunique with full converse reflection, unchanged blue wiring and even black multiplicities. E_neighbor3 becomes unique under converse reflection by the singleton-or-converse certificate. No alternative tournament realizability is claimed.','results':results,'E33_certificate':e33certificate,'checks':dict(sorted(checks.items())),'input_sha256':{('wiring/' if p.parent==wiring else 'incidence/')+p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputpaths}}
 dest=HERE/'residual_counterfeits.json'
 if args.write:dest.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 else:need(json.loads(dest.read_text())==out,'frozen counterfeit receipt reproduction')
 print('PASS',sum(checks.values()),'checks;', 'receipt written' if args.write else 'receipt matched',flush=True)
if __name__=='__main__':main()
