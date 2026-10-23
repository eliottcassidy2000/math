#!/usr/bin/env python3
"""Literal converse-equivariant counterfeits for all255 binary E33 quotients."""
import argparse,array,csv,gzip,hashlib,itertools,json,sys
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;FULL=(1<<21)-1
TILES=[(b,a) for b in range(6) for a in range(b+2,8)];REF=[TILES.index((7-a,7-b)) for b,a in TILES]
START=sum(1<<k for k,(b,a) in enumerate(TILES) if a==7);END=sum(1<<k for k,(b,a) in enumerate(TILES) if b==0)
checks=Counter()
def need(x,msg):
 checks[msg]+=1
 if not x:raise ValueError(msg)
def clean(x):return Counter({k:v for k,v in x.items() if v})
def reflect(m):return sum(((m>>k)&1)<<REF[k] for k in range(21))
def main():
 p=argparse.ArgumentParser();p.add_argument('--wiring',type=Path,default=HERE.parent/'session38-wiring');p.add_argument('--write',action='store_true');a=p.parse_args();W=a.wiring
 classes=list(csv.DictReader(gzip.open(W/'classes.tsv.gz','rt'),delimiter='\t'));N=len(classes);conv=[int(r['converse_class']) for r in classes];node=[int(r['node_id']) for r in classes];c3=[int(r['c3']) for r in classes]
 need(N==6880,'complete ordinary class count');need(all(int(r['class_id'])==i for i,r in enumerate(classes)),'class order');need(all(conv[conv[i]]==i for i in range(N)),'converse involution')
 # Recompute every base symbol directly from literal class adjacency, no imported producer.
 D=[]
 for r in classes:
  code=int(r['canon']);A=[[0]*8 for _ in range(8)]
  for k,(i,j) in enumerate(itertools.combinations(range(8),2)):A[i][j]=(code>>(27-k))&1;A[j][i]=1-A[i][j]
  cv=[0]*8;ce=[[0]*8 for _ in range(8)]
  for i,j,k in itertools.combinations(range(8),3):
   if A[i][j]==A[j][k]==A[k][i]:
    for v in (i,j,k):cv[v]+=1
    for u,v in ((i,j),(i,k),(j,k)):ce[u][v]+=1;ce[v][u]+=1
  need(sum(cv)==3*int(r['c3']),'independent cyclic triple count')
  D.append([tuple(3*(cv[w]%3)+(ce[v][w]%3) for w in range(8) if A[v][w]) for v in range(8)])
 def sig(i,part):return tuple(sorted((len(row),sum(part[z]=='1' for z in row)) for row in D[i]))
 candidates=json.loads((HERE/'candidates.json').read_text());partitions={'0'+format(k,'08b') for k in range(1,256)}
 need(len(candidates)==255 and {r['partition'] for r in candidates}==partitions,'all binary unlabeled partitions exactly once')
 def mirror(k):i,j,x,y=k;return(conv[i],conv[j],7-y,7-x)
 deltas={};needed=set();donors=set();mergedkeys=set()
 for r in candidates:
  part=r['partition'];need(r['witness'] is not None,'candidate available');i,h,j,k,x,y=r['witness']
  need(len({node[i],node[h],node[j],node[k]})==4,'four distinct merged vertices');need(c3[j]==c3[k],'target triangle counts equal')
  need(sig(j,part)==sig(k,part),'target binary quotient equality');need(sig(conv[j],part)==sig(conv[k],part),'converse binary quotient equality')
  d=Counter({(i,j,x,y):-1,(i,k,x,y):1,(h,j,x,y):1,(h,k,x,y):-1})
  for cell,z in list(d.items()):d[mirror(cell)]+=z
  d=clean(d);deltas[part]=d;needed.update(d);donors.update(((i,j,x,y),(h,k,x,y)));mergedkeys.update((node[ii],node[jj]) for ii,jj,xx,yy in d)
 actual={k:0 for k in needed}
 for r in csv.DictReader(gzip.open(W/'ports.tsv.gz','rt'),delimiter='\t'):
  if int(r['blue']):continue
  key=tuple(int(r[k]) for k in ('source_class','target_class','start_score','end_score'))
  if key in actual:actual[key]=int(r['multiplicity'])
 baseline={k:0 for k in mergedkeys}
 for r in csv.DictReader(gzip.open(W/'edges.tsv.gz','rt'),delimiter='\t'):
  key=int(r['source']),int(r['target']);z=int(r['black']);need(z%2==0,'native merged black parity')
  if key in baseline:baseline[key]=z
 mapping=array.array('H',gzip.decompress((W/'class_of_mask.u16le.gz').read_bytes()))
 if sys.byteorder!='little':mapping.byteswap()
 need(len(mapping)==1<<21,'complete fixed-path map')
 masks={}
 for m in range(FULL+1):
  if not(m>>5&1):continue
  key=mapping[m],mapping[m^FULL],7-(m&START).bit_count(),(m&END).bit_count()
  if key in donors and key not in masks and m!=reflect(m):masks[key]=m
  if len(masks)==len(donors):break
 need(len(masks)==len(donors),'actual black donor mask for every witness')
 results=[]
 for r in candidates:
  part=r['partition'];i,h,j,k,x,y=r['witness'];m1=masks[i,j,x,y];m2=masks[h,k,x,y]
  matching={m1:m2^FULL,m2:m1^FULL}
  for m,t in list(matching.items()):
   need(reflect(m) not in matching,'four distinct reflected source masks');matching[reflect(m)]=reflect(t)
  need(set(matching.values())=={m^FULL for m in matching},'actual target port inventory preserved')
  d=Counter();literal=[]
  for m,t in sorted(matching.items()):
   ci,old,new=mapping[m],mapping[m^FULL],mapping[t];aa=7-(m&START).bit_count();bb=(m&END).bit_count()
   need(bool(m>>5&1) and not(t>>5&1),'literal closed to open')
   need(m!=reflect(m) and t!=reflect(t),'literal black color unchanged')
   need(matching[reflect(m)]==reflect(t),'literal matching commutes with reflection')
   need(sig(old,part)==sig(new,part),'literal source quotient label preserved')
   need((7-(t&START).bit_count(),(t&END).bit_count())==(8-aa,6-bb),'literal endpoint compatibility')
   need(c3[new]-c3[ci]==aa-bb-1,'literal pointwise triangle law')
   need(t!=(m^FULL),'literal alternate pairing genuinely differs')
   need(mapping[reflect(m)]==conv[ci] and mapping[reflect(t)]==conv[new],'literal class reflection')
   d[ci,old,aa,bb]-=1;d[ci,new,aa,bb]+=1
   literal.append({'closed_source_mask':m,'native_open_target_mask':m^FULL,'counterfeit_open_target_mask':t})
  d=clean(d);need(d==deltas[part],'literal table delta equals candidate')
  row=Counter();col=Counter();md=Counter();ud=Counter();current=Counter()
  for(ci,cj,aa,bb),z in d.items():
   need(actual[ci,cj,aa,bb]+z>=0,'all new ordinary cell counts nonnegative')
   need(d[mirror((ci,cj,aa,bb))]==z,'full class converse reflection')
   if (ci,cj,aa,bb)==mirror((ci,cj,aa,bb)):need((actual[ci,cj,aa,bb]+z)%2==0,'reflection-fixed black cells remain even')
   row[ci,aa,bb,sig(cj,part)]+=z;col[cj,8-aa,6-bb]+=z;md[node[ci],node[cj]]+=z;ud[min(node[ci],node[cj]),max(node[ci],node[cj])]+=z
   current[node[ci]]+=z*(c3[cj]-c3[ci]);current[node[cj]]+=z*(c3[ci]-c3[cj])
  md=clean(md);ud=clean(ud)
  need(not any(row.values()),'every source quotient bucket preserved');need(not any(col.values()),'every target endpoint capacity preserved');need(not any(current.values()),'every merged triangle current preserved')
  need(bool(md) and bool(ud),'nonzero directed and undirected wiring projection');need(all(z%2==0 and baseline[e]+z>=0 for e,z in md.items()),'changed literal black edges remain even and nonnegative')
  need(all(u!=v for u,v in ud),'all loop multiplicities unchanged')
  results.append({'partition':part,'source_classes':[i,h],'target_classes':[j,k],'start_score':x,'end_score':y,'literal_matching':literal,'ordinary_cell_delta':[[*key,actual[key],z,actual[key]+z] for key,z in sorted(d.items())],'directed_merged_delta':[[*key,z] for key,z in sorted(md.items())],'undirected_merged_delta':[[*key,z] for key,z in sorted(ud.items())]})
 out={'scope':'FINITE-EXACT n8: every binary quotient of the fixed9 E33 symbols admits a literal four-port counterfeit commuting with reflection and preserving endpoint/color/quotient signatures, triangle currents, black parity and all loops. No alternative complement-realizable tournament wiring is claimed.','base_symbol':'3*(cyclic_triangle_count_at_head mod3)+(cyclic_triangle_count_at_edge mod3)','partitions_tested':255,'verified_counterfeits':len(results),'distinct_literal_switches':len({tuple((d['closed_source_mask'],d['counterfeit_open_target_mask']) for d in r['literal_matching']) for r in results}),'unresolved':[],'results':results,'checks':dict(sorted(checks.items())),'input_sha256':{'candidates.json':hashlib.sha256((HERE/'candidates.json').read_bytes()).hexdigest(),**{name:hashlib.sha256((W/name).read_bytes()).hexdigest() for name in ('classes.tsv.gz','ports.tsv.gz','edges.tsv.gz','class_of_mask.u16le.gz')}}}
 dest=HERE/'receipt.json'
 if a.write:dest.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 else:need(json.loads(dest.read_text())==out,'frozen receipt reproduction')
 print(json.dumps({'status':'PASS','partitions':255,'literal_counterfeits':len(results),'unresolved':[],'checks':sum(checks.values()),'distinct_target_pairs':len({tuple(r['target_classes']) for r in results}),'distinct_donor_cells':len(donors)},indent=2))
if __name__=='__main__':main()
