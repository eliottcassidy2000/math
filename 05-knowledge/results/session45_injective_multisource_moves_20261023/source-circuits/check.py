#!/usr/bin/env python3
"""Independent endpoint-fibre census and literal three-source circuit audit.
No validation-critical assertions. Full black F/R resource keys are explicit.
"""
from array import array
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import argparse,csv,gzip,json,sys
HERE=Path(__file__).resolve().parent;CHECKS=0

def need(ok,why):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(why)

def mat(c):
 A=[[0]*8 for _ in range(8)]
 for b,(i,j) in enumerate(combinations(range(8),2)):A[i][j]=c>>(27-b)&1;A[j][i]=1-A[i][j]
 return tuple(map(tuple,A))

def native_paths(A):
 out=[];stack=[(i,) for i in range(8)]
 while stack:
  p=stack.pop()
  if len(p)==8:out.append(p);continue
  for v in range(8):
   if v not in p and A[p[-1]][v]:stack.append(p+(v,))
 return out

def literal_orbit(A,p):
 keep={frozenset((i,j)) for i,j in zip(p,p[1:])}
 B=tuple(tuple(0 if i==j else A[i][j] if frozenset((i,j)) in keep else 1-A[i][j] for j in range(8)) for i in range(8))
 return frozenset([(A,p),(B,p),(tuple(zip(*A)),p[::-1]),(tuple(zip(*B)),p[::-1])])

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 blobs={}
 for p in json.load(open(HERE/'SOURCE_PINS.json'))['files']:
  raw=(a.repo/p['path']).read_bytes();need(len(raw)==p['bytes'] and sha256(raw).hexdigest()==p['sha256'],'original source hash '+p['path']);blobs[Path(p['path']).name]=raw
 rows=list(csv.DictReader(blobs['merged_nodes.tsv'].decode().splitlines(),delimiter='\t'));need(len(rows)==3528,'full source universe')
 owner=array('H');owner.frombytes(gzip.decompress(blobs['owner_map.u16.gz']))
 if sys.byteorder!='little':owner.byteswap()
 mid={int(r[k]):int(r['merged_id']) for r in rows for k in ('class_a','class_b')};full=(1<<21)-1;chs=[(i,j) for i in range(8) for j in range(i+2,8)];rev=[chs.index((7-j,7-i)) for i,j in chs]
 def qorbit(q):
  r=sum((q>>k&1)<<rev[k] for k in range(21));return frozenset((q,r,q^full,r^full))
 D=json.loads(blobs['receipt.json'])['order8'];X=D['higher_layers'][-1]['lifts'][0];Y=D['radical_lifts_mod4'];mods=(2,2,2,2,128);labels=[tuple(y[i]%2 for y in Y[:4])+(X[i]%128,) for i in range(3528)]
 edges=[tuple(int(r[k]) for k in ('source','target','blue','black')) for r in csv.DictReader(gzip.decompress(blobs['merged_edges.tsv.gz']).decode().splitlines(),delimiter='\t')]
 need(len(edges)==438216,'full edge universe')
 def apply(x):
  out=[0]*3528
  for u,v,b,k in edges:z=(b+k)*(x[u]-x[v]);out[u]+=z;out[v]-=z
  return out
 for x,m in [(X,256)]+[(y,4) for y in Y]:need(all(z%m==0 for z in apply(x)),'original full harmonic '+str(m))
 need(any(z%512 for z in apply(X)),'fixed integer X is not harmonic mod512')
 need(all(z%512==0 for z in apply([2*v for v in X])),'actual scaled2X is harmonic mod512')
 # Independent complete census: group only by endpoint+target label, then inspect
 # shared adjacent clocks. Generator instead uses five indexed-arc fibres.
 actual={};source_counts={};total_paths=0
 for s,row in enumerate(rows):
  A=mat(int(row['canon_a']));paths=native_paths(A);source_counts[s]=len(paths);total_paths+=len(paths);fibres=defaultdict(list)
  for p in paths:
   q=sum(A[p[i]][p[j]]<<b for b,(i,j) in enumerate(chs));need(mid[owner[q]]==s,'complete source path authentication')
   t=mid[owner[q^full]];orb=qorbit(q)
   if t==s or len(orb)!=4:continue
   fibres[(p[0],p[-1],labels[t])].append((p,t,min(orb)))
  local=defaultdict(set)
  for rs in fibres.values():
   for p,t,q in rs:
    for r,u,z in rs:
     if t>=u:continue
     if any(p[k:k+2]==r[k:k+2] for k in range(1,6)):local[(t,u)].add((q,z))
  for targets,moves in local.items():actual[(s,*targets)]=moves
 frozen=json.load(open(HERE/'full-search/all_menus.json'));expected={(r['source'],*r['targets']):{(p['before']['resource'],p['after']['resource']) for p in r['moves']} for r in frozen}
 need(actual==expected,'independent complete endpoint census equals all frozen resource-edge menus')
 capacity2={k for k,moves in actual.items() if any(a[0]!=b[0] and a[1]!=b[1] for a,b in combinations(moves,2))}
 need(len(actual)==700 and len(capacity2)==8,'complete black source-targetpair menus and injective capacity2 census')
 groups=defaultdict(list)
 for s,t,u in actual:groups[(t,u)].append(s)
 twos=[];threes=[]
 for (t,u),ss in groups.items():
  ss=sorted(ss)
  for x,z in combinations(ss,2):
   if {x,z}&{t,u}:continue
   diff=[(c-d)%m for c,d,m in zip(labels[x],labels[z],mods)]
   if any(diff) and all(2*d%m==0 for d,m in zip(diff,mods)) and (x,t,u) in capacity2 and (z,t,u) in capacity2:twos.append([t,u,x,z])
   for y in ss:
    if y in (x,z,t,u) or (y,t,u) not in capacity2:continue
    if labels[x]!=labels[y] and all((c+d-2*b)%m==0 for c,d,b,m in zip(labels[x],labels[z],labels[y],mods)):threes.append([t,u,x,y,z])
 need(twos==[] and threes==[[2625,2628,1981,2458,2659]],'exact genuine four-slot two/three-source circuit census')
 packet=json.load(open(HERE/'three_source_packet.json'));states=[[],[]];resources=[[],[]];debts=[];allpins=[]
 for slot in packet['slots']:
  s=slot['source'];A=mat(int(rows[s]['canon_a']));need([list(r) for r in A]==slot['adjacency'],'explicit canonical physical source matrix')
  need(slot['local_to_physical_owners']==list(range(8)),'one exact identity owner frame')
  p0,p7,j,u,v=slot['pins'];records=[]
  for state,key in enumerate(('old','new')):
   rec=slot[key];p=tuple(rec['path']);need(sorted(p)==list(range(8)) and all(A[u][v] for u,v in zip(p,p[1:])),'actual Hamiltonian source path')
   q=sum(A[p[i]][p[j]]<<b for b,(i,j) in enumerate(chs));t=mid[owner[q^full]]
   need(q==rec['mask'] and t==rec['target'] and min(qorbit(q))==rec['resource'] and len(qorbit(q))==4,'original mask,target,black resource authentication')
   need(p[0]==p0 and p[-1]==p7 and p[j:j+2]==(u,v),'literal endpoint and indexed internal arc pins')
   need(all(p[k]==a for k,a in slot['preserved_owner_pins']),'every retained indexed owner pin')
   states[state].append(literal_orbit(A,p));resources[state].append(rec['resource']);records.append((s,t,(X[s]-X[t])%8))
  need(records[0][2]==records[1][2],'per-slot original mod8 phase retained')
  allpins.append(dict(source=s,old_phase=records[0][2],new_phase=records[1][2],pins=slot['preserved_owner_pins']))
 need(packet['slots'][2]['adjacency']==packet['slots'][3]['adjacency'] and packet['slots'][2]['local_to_physical_owners']==packet['slots'][3]['local_to_physical_owners'],'middle resource slots share ONE whole source and owner frame')
 for st in states:need(all(len(o)==4 for o in st) and len(set().union(*st))==16,'four distinct physical K/converse resources per endpoint')
 need(len(set().union(*(states[0]+states[1])))==32 and len(set(resources[0]+resources[1]))==8,'all eight black resources distinct across endpoint union')
 def debt(char,indices):
  d=Counter()
  for i in indices:
   slot=packet['slots'][i];s=slot['source']
   for key,sgn in [('old',-1),('new',1)]:t=slot[key]['target'];v=2*sgn*(char[s]-char[t]);d[s]+=v;d[t]-=v
  return dict(sorted((k,v) for k,v in d.items() if v))
 xdebt=debt(X,range(4));need(xdebt=={2625:256,2628:-256},'exact nonzero integral X debt is256 target difference')
 generators=[(X,256)]+[(y,4) for y in Y]
 for char,m in generators:need(all(v%m==0 for v in debt(char,range(4)).values()),'joint packet keeps every generating harmonic modulo its actual modulus')
 indiv=[]
 for inds in ([0,2],[1,3]):
  dx=debt(X,inds);ys=[debt(y,inds) for y in Y]
  need(any(v%256 for v in dx.values()),'each constituent rectangle fails original higher harmonic256')
  need(all(v%8==0 for v in dx.values()),'each constituent keeps chosen original mod8 character')
  indiv.append(dict(slots=inds,integer_X_debt=dx,mod4_Y_debts=ys))
 need(all(v%512==0 for v in debt([2*z for z in X],range(4)).values()),'actual scaled512 harmonic has zero512 debt')
 # Module relation itself is direct, not inferred from weighted total zero.
 t,u=packet['targets'];coeffs={1981:1,2458:-2,2659:1}
 need(labels[t]==labels[u] and all(sum(c*labels[s][j] for s,c in coeffs.items())%m==0 for j,m in enumerate(mods)),'exact mixed-module charge circuit')
 receipt=dict(status='FINITE-EXACT PASS',checks=CHECKS,universe=dict(canonical_sources=3528,actual_source_paths=total_paths,black_endpoint_arc_menus=len(actual),capacity2_menus=[list(k) for k in sorted(capacity2)],genuine_two_source_four_slot_repairs=twos,genuine_three_source_four_slot_repairs=threes),packet_sha256=sha256((HERE/'three_source_packet.json').read_bytes()).hexdigest(),source_path_counts={str(s):source_counts[s] for s in coeffs},source_labels={str(s):labels[s] for s in coeffs},target_labels={str(v):labels[v] for v in (t,u)},source_X={str(s):X[s] for s in coeffs},target_X=[X[t],X[u]],shared_pins=allpins,literal_resource_counts=dict(per_endpoint=4,union=8,ports_per_endpoint=16,union_ports=32),original_resource_keys=resources,individual_rectangles=indiv,joint_integer_X_debt=xdebt,mod256_preserved=True,raw_X_mod512_harmonic=False,scaled2X_mod512_harmonic_and_preserved=True,scope='Partial selected paired-black resource replacement. Finite minimum source support applies only to four-slot common equal-mixed-label-target-pair moves with endpoints and an indexed internal arc retained. All-dyadic full-bank interpretation uses the independently audited complete electrical filtration and actual harmonic scaling, not exact integer X equality; no global arbitrary remainder insertion and no full per-port bank preservation.')
 a.out.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(dict(status=receipt['status'],checks=CHECKS,receipt_sha256=sha256(a.out.read_bytes()).hexdigest())))
if __name__=='__main__':main()
