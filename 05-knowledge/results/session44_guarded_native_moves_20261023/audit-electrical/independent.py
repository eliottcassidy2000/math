#!/usr/bin/env python3
"""Independent arithmetic/atlas audit; no producer import and no floating point."""
from array import array
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse, csv, gzip, hashlib, json, sys

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--producer',type=Path,default=HERE.parent/'session44-electrical-duality')
parser.add_argument('--repo',type=Path,default=None)
args=parser.parse_args()
PACK=args.producer.resolve()
INPUT=PACK/'inputs'
CHECKS=0

def need(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(msg)

def get(name):
 p=INPUT/('05-knowledge__results__'+name.replace('/','__'))
 if not p.exists():
  bases=([args.repo] if args.repo else [])+list(HERE.parents)
  p=next((b/'05-knowledge/results'/name for b in bases if (b/'05-knowledge/results'/name).exists()),p)
 raw=p.read_bytes()
 pins=json.loads((PACK/'SOURCE_PINS.json').read_text())
 need(hashlib.sha256(raw).hexdigest()==pins[name], 'Authenticated inherited source bytes')
 return raw

def matvec(a,x,m):
 return tuple(sum(v*w for v,w in zip(row,x))%m for row in a)

def connected(n,edges,skip=None):
 seen={0};todo=[0]
 while todo:
  v=todo.pop()
  for i,j,w in edges:
   if skip==(i,j):continue
   u=j if i==v else i if j==v else None
   if u is not None and u not in seen:seen.add(u);todo.append(u)
 return len(seen)==n

def solve(a,b):
 a=[[F(x) for x in row]+[F(v)] for row,v in zip(a,b)]
 for j in range(len(a)):
  p=next(i for i in range(j,len(a)) if a[i][j])
  a[j],a[p]=a[p],a[j];v=a[j][j];a[j]=[x/v for x in a[j]]
  for i in range(len(a)):
   if i!=j:
    v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[j])]
 return [row[-1] for row in a]

def small_graphs():
 total=0;pairs=0;aliases=0
 for n in (2,3,4):
  ee=list(combinations(range(n),2))
  for weights in product(range(3),repeat=len(ee)):
   edges=[(a,b,w) for (a,b),w in zip(ee,weights) if w]
   if not connected(n,edges):continue
   total+=1
   L=[[0]*n for _ in range(n)]
   for a,b,w in edges:L[a][a]+=w;L[b][b]+=w;L[a][b]-=w;L[b][a]-=w
   bridge_edges=[(a,b) for a,b,w in edges if w==1 and not connected(n,edges,(a,b))]
   def bridge_connected(a,b):
    seen={a};todo=[a]
    while todo:
     v=todo.pop()
     for u,w in bridge_edges:
      z=w if u==v else u if w==v else None
      if z is not None and z not in seen:seen.add(z);todo.append(z)
    return b in seen
   A=[row[:-1] for row in L[:-1]]
   for a,b in combinations(range(n),2):
    charge=[int(i==a)-int(i==b) for i in range(n-1)]
    f=solve(A,charge)
    integral=all(v.denominator==1 for v in f)
    need(integral==bridge_connected(a,b),'Integer image iff literal bridge path')
    pairs+=1;aliases+=integral
 # Single modulus is not all moduli: a unit triangle is invisible modulo8.
 L=[[2,-1,-1],[-1,2,-1],[-1,-1,2]]
 H8=[(a,b,0) for a,b in product(range(8),repeat=2) if matvec(L,(a,b,0),8)==(0,0,0)]
 H9=[(a,b,0) for a,b in product(range(9),repeat=2) if matvec(L,(a,b,0),9)==(0,0,0)]
 need(H8==[(0,0,0)] and any(x[0]!=x[1] for x in H9),'Fixed-modulus alias hostile')
 # Rowwise criterion over Z8, including nonunit conductances and nonsymmetric delta.
 rowcases=0
 for w in range(1,9):
  L=[[w,-w],[-w,w]];vec=list(product(range(8),repeat=2))
  H=[x for x in vec if x[1]==0 and matvec(L,x,8)==(0,0)]
  image={matvec(L,x,8) for x in vec}
  for a,b in product(range(8),repeat=2):
   delta=[[a,-a],[b,-b]]
   annihilates=all(matvec(delta,x,8)==(0,0) for x in H)
   rows_in_image=tuple(v%8 for v in delta[0]) in image and tuple(v%8 for v in delta[1]) in image
   need(annihilates==rows_in_image,'Composite rowwise criterion')
   rowcases+=1
 need(matvec([[0,0],[0,0]],(1,0),8)==(0,0) and matvec([[1,-1],[-1,1]],(1,0),8)!=(0,0),
      'Delta=-L preserves old bank but can enlarge it')
 return {'connected_weighted_graphs':total,'vertex_pairs':pairs,'all_modulus_aliases':aliases,
         'rowwise_mod8_cases':rowcases,'single_modulus_cycle_hostile':True,'bank_inclusion_not_equality_hostile':True}

def atlas():
 N=3528;root=N-1
 nodes=list(csv.DictReader(get('metagraph_order8_20261023/merged_nodes.tsv').decode().splitlines(),delimiter='\t'))
 edges=[];adj=[set() for _ in range(N)];rows=[0]*N
 for row in csv.DictReader(gzip.decompress(get('metagraph_order8_20261023/merged_edges.tsv.gz')).decode().splitlines(),delimiter='\t'):
  i,j=int(row['source']),int(row['target']);w=int(row['blue'])+int(row['black'])
  edges.append((i,j,w));adj[i].add(j);adj[j].add(i)
  if i!=j and w%2:rows[i]^=(1<<i)|(1<<j);rows[j]^=(1<<i)|(1<<j)
 # Different computation from producer: reduce vertex unit charges modulo the
 # full binary ROW SPACE, rather than constructing a grounded kernel basis.
 echelon={}
 for r in reversed(rows):
  while r:
   p=r.bit_length()-1
   if p in echelon:r^=echelon[p]
   else:echelon[p]=r;break
 def residue(r):
  for p in sorted(echelon,reverse=True):
   if (r>>p)&1:r^=echelon[p]
  return r
 binary=[residue(1<<i) for i in range(N)]
 need(len(echelon)==164,'Independent full binary rank')
 d=json.loads(get('metagraph_order8_electrical_20261023/receipt.json'))['order8']
 Y=d['radical_lifts_mod4'];Z=d['eighth_radical_lifts_mod8'][0];X=d['higher_layers'][-1]['lifts'][0]
 need(d['root']==root,'Actual original grounding')
 def lap(x):
  out=[0]*N
  for a,b,w in edges:
   out[a]+=w*(x[a]-x[b]);out[b]+=w*(x[b]-x[a])
  return out
 for x,m in [(y,4) for y in Y]+[(Z,8),(X,256)]:
  need(x[root]==0 and len(x)==N and all(v%m==0 for v in lap(x)),'Independent full modular lift authentication')
 signatures={2:binary,4:[(binary[i],tuple(y[i]%4 for y in Y)) for i in range(N)],
             8:[(binary[i],tuple(y[i]%4 for y in Y),Z[i]%8) for i in range(N)]}
 partitions={}
 for m,sigs in signatures.items():
  buckets=defaultdict(list)
  for i,sig in enumerate(sigs):buckets[sig].append(i)
  partitions[m]={'classes':len(buckets),'nonsingletons':sorted(v for v in buckets.values() if len(v)>1)}
 pairs=[[0,3527],[53,3492],[109,3523]]
 need(partitions[2]['classes']==3521 and partitions[4]['classes']==3524 and partitions[8]['classes']==3525,
      'Independent dyadic vertex partition counts')
 need(partitions[8]['nonsingletons']==pairs,'Independent exact vertex aliases')
 for parent,leaf in pairs:
  x=[0]*N;x[leaf]=1;actual=lap(x);expected=[0]*N;expected[leaf]=1;expected[parent]=-1
  need(actual==expected,'Integer all-modulus alias certificate')
  need(adj[leaf]=={parent} and parent not in adj[parent],'No co-target leaf-parent pair')
 # Independent complete presentation replay, retain directed source target sets.
 raw=gzip.decompress(get('metagraph_order8_20261023/owner_map.u16.gz'))
 owners=array('H');owners.frombytes(raw)
 if sys.byteorder!='little':owners.byteswap()
 need(len(owners)==2**21 and max(owners)<6880,'Complete inherited owner-map universe')
 merged=[None]*6880
 for row in nodes:
  for key in ('class_a','class_b'):merged[int(row[key])]=int(row['merged_id'])
 need(all(x is not None for x in merged),'Every source class has merged owner')
 targetmasks=[0]*N;counts=Counter();last=len(owners)-1
 for q in range(len(owners)):
  a=merged[owners[q]];b=merged[owners[last-q]]
  targetmasks[a]|=1<<b
  counts[min(a,b)*N+max(a,b)]+=1
 need(len(counts)==len(edges)==438216,'All exported unordered pair rows reconstructed')
 need(all(counts[a*N+b]==2*w for a,b,w in edges),'Every native pair multiplicity matches')
 oldambiguous=0;paircollisions=0
 reports=[]
 for s,mask in enumerate(targetmasks):
  targets=[]
  while mask:
   bit=mask&-mask;mask-=bit;targets.append(bit.bit_length()-1)
  need(set(targets)==adj[s],'Exact source neighborhood')
  need(len({signatures[8][t] for t in targets})==len(targets),'Every source has injective universal target signature')
  h=Counter((X[s]-X[t])%8 for t in targets)
  lost=sum(v*(v-1)//2 for v in h.values())
  oldambiguous+=bool(lost)
  if s in (1238,2242,2246,2315):reports.append([s,len(targets),lost])
 nodehist=Counter(x%8 for x in X)
 paircollisions=sum(v*(v-1)//2 for v in nodehist.values())
 need(oldambiguous==3518 and paircollisions==777399,'Original single-character losses')
 # Read the final result solely to check reporting consistency, not as evidence.
 producer=json.loads((PACK/'receipt.json').read_text())
 need(producer['hostile']['source']==1238,'Hostile enclosing source must equal authenticated path sources')
 need(all(row['source']==1238 for row in producer['hostile']['paths']),'Hostile path source records')
 return {'vertices':N,'binary_full_rank':len(echelon),'partitions':partitions,
         'presentation_masks':len(owners),'pair_rows':len(edges),'sources_tested':len(targetmasks),
         'all_source_target_signatures_injective':True,'single_character_ambiguous_sources':oldambiguous,
         'single_character_vertex_false_equalities':paircollisions-3,'source_reports':reports,
         'mod8_charge_quotient_bits':3358+2*4+3,'minimum_module_generators':3358+4+1}

def main():
 out={'status':'INDEPENDENT PASS','small':small_graphs(),'atlas':atlas(),'checks':CHECKS}
 text=json.dumps(out,indent=2,sort_keys=True)+'\n'
 print(text,end='')

if __name__=='__main__':main()
