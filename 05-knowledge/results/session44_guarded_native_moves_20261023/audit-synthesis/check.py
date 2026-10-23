#!/usr/bin/env python3
"""Independent synthesis cross-check. No producer or prior-auditor imports.
Finite arithmetic controls do not establish all-size theorems.
"""
from pathlib import Path
from array import array
from collections import Counter,defaultdict
from functools import lru_cache
from fractions import Fraction
from itertools import product,permutations,combinations
import argparse,csv,gzip,hashlib,json,sys
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,default=HERE.parent);p.add_argument('--repo',type=Path);args=p.parse_args()
W=args.workspace.resolve()
candidates=([args.repo] if args.repo else [])+[W/'session35-publication/clean-repo']+list(HERE.parents)
repo=next((q for q in candidates if (q/'05-knowledge/results/metagraph_order8_20261023/owner_map.u16.gz').exists()),None)
if repo is None:raise RuntimeError('Supply --repo pointing to the original math repository')
repo=repo.resolve();R=repo/'05-knowledge/results';checks=0;hashes={}
def need(b,m):
 global checks
 checks+=1
 if not b:raise RuntimeError(m)
def read(path,key=None):
 raw=path.read_bytes()
 if key is not None:pass
 elif path.is_relative_to(repo):key='repo/'+str(path.relative_to(repo))
 elif path.is_relative_to(W):key=str(path.relative_to(W));key=key.removeprefix('session44-')
 else:key=str(path)
 hashes[key]=hashlib.sha256(raw).hexdigest();return raw
def packetpath(rel):
 q=W/rel
 if q.exists():return q
 logical=rel.removeprefix('session44-').removeprefix('root/')
 return W/logical
def load(rel):return json.loads(read(packetpath(rel),rel.removeprefix('session44-').removeprefix('root/')))
def harmonic(L,m):
 return [x for tail in product(range(m),repeat=len(L)-1) for x in [(0,)+tail] if all(sum(a*b for a,b in zip(row,x))%m==0 for row in L)]
def current(s,t,x,m):
 return tuple((x[s]-x[t])*(int(v==s)-int(v==t))%m for v in range(len(x)))
# The scalar-phase/current-vector implication really fails, even modulo8.
L=((2,-1,0,-1),(-1,3,-2,0),(0,-2,4,-2),(-1,0,-2,3))
H=harmonic(L,8)
need(H==[(0,0,0,0),(0,0,4,0),(0,4,2,4),(0,4,6,4)],'complete weighted-C4 bank')
need(all(x[1]==x[3] for x in H),'universal phase equality')
need(current(0,1,H[2],8)!=(current(0,3,H[2],8)),'target-owned current differs')
controls={'weighted_C4_mod8_harmonics':H,'distinct_currents':[current(0,j,H[2],8) for j in (1,3)]}
# Native cross-packet reconstruction from original immutable atlas bytes.
pins=load('session44-electrical-duality/SOURCE_PINS.json')
for name,h in pins.items():
 raw=read(R/name);need(hashlib.sha256(raw).hexdigest()==h,'original input SHA256 '+name)
f=R/'metagraph_order8_20261023'
nodes=list(csv.DictReader((f/'merged_nodes.tsv').open(),delimiter='\t'))
owner=array('H');owner.frombytes(gzip.decompress((f/'owner_map.u16.gz').read_bytes()))
if sys.byteorder!='little':owner.byteswap()
merged={int(row[k]):int(row['merged_id']) for row in nodes for k in ('class_a','class_b')}
x=json.loads((R/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']['higher_layers'][-1]['lifts'][0]
pairs=list(combinations(range(8),2));chords=[(i,j) for i,j in pairs if j>i+1]
def matrix(code):
 A=[[0]*8 for _ in range(8)]
 for b,(i,j) in enumerate(pairs):A[i][j]=(code>>(27-b))&1;A[j][i]=1-A[i][j]
 return A
def partner(A,P):
 keep={frozenset((a,b)) for a,b in zip(P,P[1:])}
 return [[0 if a==b else A[a][b] if frozenset((a,b)) in keep else A[b][a] for b in range(len(A))] for a in range(len(A))]
adaptive=load('session44-actual-guards/receipt.json');records=load('session44-actual-guards/path_records.json');endpoints=load('session44-actual-guards/endpoint.normal.json')
reports=[];native={};total_fibres=0
for source,want,ends in zip((1238,2242,2246,2315),adaptive['sources'],endpoints['sources']):
 A=matrix(int(nodes[source]['canon_a']));need(A==want['adjacency'],'source matrix')
 rows=[]
 for P in permutations(range(8)):
  if not all(A[a][b] for a,b in zip(P,P[1:])):continue
  word=sum(A[P[i]][P[j]]<<b for b,(i,j) in enumerate(chords));target=merged[owner[word^((1<<21)-1)]]
  need(merged[owner[word]]==source,'source ownership')
  rows.append((P,target,(x[source]-x[target])%8,word))
 need(len(rows)==want['H'],'literal path inventory')
 actual={(P,t,ph,w) for P,t,ph,w in rows}
 claimed={(tuple(r['path']),r['target'],r['phase'],r['mask']) for r in records[str(source)]}
 need(actual==claimed,'complete records cross-packet match')
 native[source]=(A,rows)
 outputs=[[r[2] for r in rows],[r[1] for r in rows],[r[0] for r in rows]]
 names=('phase','target','path')
 @lru_cache(None)
 def solve(S):
  homogeneous=[all(v[i]==v[S[0]] for i in S) for v in outputs]
  if all(homogeneous):return ((0,0),(0,0),(0,0))
  choices=[]
  for j in range(8):
   blocks=defaultdict(list)
   for i in S:blocks[rows[i][0][j]].append(i)
   if len(blocks)==1:continue
   child=[solve(tuple(b)) for b in blocks.values()]
   choices.append(tuple((1+max(c[q][0] for c in child),len(S)+sum(c[q][1] for c in child)) for q in range(3)))
  return tuple((0,0) if homogeneous[q] else (min(c[q][0] for c in choices),min(c[q][1] for c in choices)) for q in range(3))
 S=tuple(range(len(rows)));costs=solve(S);fixed=[]
 for q,key in enumerate(names):
  need(costs[q][0]==want['objectives'][key]['worst_depth'],'minimax optimum '+key)
  need(str(Fraction(costs[q][1],len(S)))==want['objectives'][key]['optimal_mean'],'minimum mean '+key)
  for k in range(8):
   found=False
   for clocks in combinations(range(8),k):
    seen={};good=True
    for i in S:
     observation=tuple(rows[i][0][j] for j in clocks)
     if observation in seen and seen[observation]!=outputs[q][i]:good=False;break
     seen[observation]=outputs[q][i]
    if good:found=True;break
   if found:fixed.append(k);break
  need(fixed[-1]==6,'exact nonadaptive footprint')
 byend=defaultdict(list)
 for i,(P,*_) in enumerate(rows):byend[P[0],P[-1]].append(i)
 distributions=[Counter() for _ in range(3)]
 for group in byend.values():
  vals=solve(tuple(group))
  for q in range(3):distributions[q][vals[q][0]]+=1
 for q,key in enumerate(names):
  need({str(k):v for k,v in sorted(distributions[q].items())}==ends['extra_adaptive_depth_distributions'][key],'endpoint histogram '+key)
  need(max(distributions[q])==3,'endpoint sharp max three '+key)
 total_fibres+=len(byend)
 reports.append({'source':source,'paths':len(rows),'depths':dict(zip(names,[c[0] for c in costs])),'minimum_mean':dict(zip(names,[str(Fraction(c[1],len(S))) for c in costs])),'fixed_depths':fixed,'endpoint_fibres':len(byend),'endpoint_depths':dict(zip(names,[dict(sorted(d.items())) for d in distributions]))})
need(total_fibres==210,'all210 endpoint fibres')
# Selected two-slot and three-slot packets: keep resource and operation types.
coupled=load('session44-coupled-moves/receipt.json')
A,rows=native[1238];slots=[]
for start,end,target in ((3,2,2397),(2,3,2790)):
 eligible=[P for P,t,ph,w in rows if (P[0],P[-1],t)==(start,end,target)]
 need(len(eligible)==2,'complete target+endpoint fibre');slots.append(eligible)
states=[s for s in product(*slots) if s[0][1:7]==s[1][1:7]]
need(len(states)==2 and all(states[0][i]!=states[1][i] for i in (0,1)),'necessary two-slot move')
need(len({P for slot in slots for P in slot})==4,'four distinct resources')
g=(4,0,2,3,1,5,6,7)
for i in (0,1):
 B,C=(partner(A,states[j][i]) for j in (0,1))
 need(all(B[a][b]==C[g[a]][g[b]] for a,b in pairs),'actual complement-target transporter')
A2,rows2=native[2242];tm=coupled['ternary_matching'];ps=[tuple(z['path']) for z in tm['paths']]
need({P for P,t,ph,w in rows2 if t==3442}==set(ps),'complete three-port target fibre')
eligible=[[i for i,P in enumerate(ps) if [P[j] for j in clocks]==pins_] for clocks,pins_ in zip(tm['clocksets'],tm['pin_values'])]
matchings=[s for s in product(*eligible) if len(set(s))==3]
need(matchings==[(0,1,2),(1,2,0)],'exact assignment cycle')
need(sum(a!=b for a,b in zip(*matchings))==3,'no pairwise assignment step')
# Full harmonic rectangle certificates checked against all original edges.
images=load('session44-rectangle-bank/image_witnesses.json')
rectangle=load('session44-rectangle-bank/receipt.json')
electrical=json.loads((R/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
ys=[[2*a for a in y] for y in electrical['radical_lifts_mod4'][:2]]
zs=[z['potential_mod8'] for z in images['solutions']]
vectors=ys+zs;lap=[[0]*len(nodes) for _ in vectors]
edge_count=0
with gzip.open(f/'merged_edges.tsv.gz','rt') as stream:
 for edge in csv.DictReader(stream,delimiter='\t'):
  a,b=int(edge['source']),int(edge['target']);w=int(edge['blue'])+int(edge['black']);edge_count+=1
  for v,Lv in zip(vectors,lap):
   term=w*(v[a]-v[b]);Lv[a]+=term;Lv[b]-=term
need(edge_count==438216,'all original rectangle edge rows')
for Lv in lap[:2]:need(all(a%8==0 for a in Lv),'genuine entire-network harmonic separator')
for i,(z,claim) in enumerate(zip(images['solutions'],rectangle['charge_orders'])):
 a,b=z['charge'];rhs=[4*(int(j==a)-int(j==b)) for j in range(len(nodes))]
 need(all((u-v)%8==0 for u,v in zip(lap[i+2],rhs)),'direct four-charge image witness')
 correction=[[j,(u-v)//8] for j,(u,v) in enumerate(zip(lap[i+2],rhs)) if u!=v]
 need(correction==z['integer_correction_div8_nonzero'],'entire integer correction')
 y=ys[claim['lower_separator_mod4_lift_index']]
 need((2*(y[a]-y[b]))%8==4,'twice charge excluded: exact order4')
_,base_rows=native[1238]
for target,count in ((881,2),(2762,1),(2650,2),(2843,1)):
 need(sum(t==target for P,t,ph,w in base_rows)==count,'complete rectangle source-port capacity')
# Spectral cross-check from incidence receipts, not a new native enumeration.
import sympy as sp
fractions=load('session44-fractions-ternary/receipt.json');zvar=sp.Symbol('z');degrees={}
for order in ('4','5'):
 for name,case in fractions['spectra'][order]['cases'].items():
  incidence=sp.Matrix(case['incidence']);nr,nc=incidence.shape
  rm=[sum(incidence[i,j] for j in range(nc)) for i in range(nr)]
  cm=[sum(incidence[i,j] for i in range(nr)) for j in range(nc)]
  K=sp.diag(*[1/v for v in rm])*incidence*sp.diag(*[1/v for v in cm])*incidence.T
  cp=K.charpoly(zvar).as_poly()
  need(cp.all_coeffs()==[sp.Rational(c) for c in case['charpoly_descending']],'independent exact incidence characteristic polynomial')
  if order=='5' and name!='refinement':
   factors=sp.factor_list(cp)[1];big=[q for q,multiplicity in factors if q.degree()>1]
   need(len(big)==1 and big[0].degree()==(5 if name=='paired_compressed' else 7),'irreducible residual degree')
   cert=fractions['spectra']['algebraic_certificates'][name]['root_isolation'];lo,hi=map(sp.Rational,cert['largest_root_interval'])
   need(big[0].count_roots(lo,hi)==1 and big[0].count_roots(hi,1)==0,'largest-root exact interval')
   degrees[name]=big[0].degree()
# Actual C9 lift: exact supplied moves, not every native metagraph edge.
V=[(a,b) for a,b in product(range(9),repeat=2) if (a-b)%3==0]
def comps(generators):
 unseen=set(V);out=[]
 while unseen:
  root=min(unseen);seen={root};todo=[root];unseen.remove(root)
  while todo:
   a,b=todo.pop()
   for c,d in generators:
    z=((a+c)%9,(b+d)%9)
    if z not in seen:seen.add(z);todo.append(z);unseen.discard(z)
  out.append(len(seen))
 return sorted(out)
need(comps(((1,1),(-1,-1)))==[9,9,9],'three native lift components')
need(comps(((1,1),(-1,-1),(3,0),(-3,0)))==[27],'kernel-direction repair')
need(sum(1 for a,b,c,d in product(range(9),repeat=4) if a%3==b%3==c%3==d%3==0 and b==c)==27,'joint gluing count at q0')
# Native continuation: count binomial residues independently, all r3..160.
from math import comb
for r in range(3,161):
 counts=lambda n:[sum(comb(n,k) for k in range(delta,n+1,3)) for delta in range(3)]
 den=counts(r);num=counts(r-1);J=(2**r-(-1)**r)//3
 for d in range(3):
  pair=sorted((Fraction(num[d],den[d]),Fraction(num[(d-1)%3],den[d])))
  need(pair==[Fraction(1,2)]*2 or pair==[Fraction(J-1,2*J),Fraction(J+1,2*J)],'Jacobsthal conditional law')
  need(Fraction(7,15) not in pair and Fraction(8,15) not in pair,'not order4 spectral pair')
  if pair[0]!=Fraction(1,2):
   need(all(abs(2*z.numerator-z.denominator)==1 for z in pair),'individual Farey adjacency')
   need(abs(pair[0].numerator*pair[1].denominator-pair[1].numerator*pair[0].denominator)==J>1,'mutual nonadjacency at r>=3')
# Nilpotency top ideal exact regular columns, using only cyclic shifts.
def shift(v,i,r):
 out=[0]*len(v)
 for index,value in enumerate(v):
  digit=(index//(3**i))%3;j=index+3**i if digit<2 else index-2*3**i
  out[j]=value
 return out
def compositions(n,r):
 if r==1:yield (n,);return
 for i in range(n+1):
  for tail in compositions(n-i,r-1):yield (i,)+tail
newton_count=0
for r in range(1,4):
 for k in range(1,5):
  N=2*(k+r)-1;m=3**k;top=[]
  for degree in (N-1,N):
   for js in compositions(degree,r):
    v=[1]+[0]*(3**r-1)
    for i,j in enumerate(js):
     for _ in range(j):v=[(a-b)%m for a,b in zip(shift(v,i,r),v)]
    expected=sum((j-1)//2 for j in js if j)>=(k)
    need((not any(v))==expected,'exact augmentation monomial valuation')
    if degree==N-1 and any(v):top.append(v)
    newton_count+=1
  z=[((-1)**(k-1)*3**(k-1))%m]*(3**r)
  need(top and all(v==z for v in top),'all top monomials equal signed norm')
# Root's stable-padding construction, audited directly before publication.
padded=[]
for k,suffix in list(product(range(7),range(9)))+[(32,0),(0,32),(16,16),(17,23)]:
 n=8+k+suffix
 B=[[0]*n for _ in range(n)]
 for a,b in combinations(range(n),2):
  B[a][b]=A[a-k][b-k] if k<=a<b<k+8 else 1;B[b][a]=1-B[a][b]
 gp=tuple(range(k))+tuple(k+v for v in g)+tuple(range(k+8,n))
 source_degrees=list(map(sum,B))
 need(source_degrees[:k]==list(range(n-1,n-k-1,-1)),'unique prefix source scores')
 need(source_degrees[k+8:]==list(range(suffix-1,-1,-1)),'unique suffix source scores')
 need(all(suffix<=d<=suffix+7 for d in source_degrees[k:k+8]),'separate core score interval')
 padded_slots=[]
 for slot in slots:
  paths=[tuple(range(k))+tuple(k+a for a in P)+tuple(range(k+8,n)) for P in slot]
  padded_slots.append(paths)
  targets=[partner(B,P) for P in paths]
  need(all(all(B[a][b] for a,b in zip(P,P[1:])) for P in paths),'padded paths Hamiltonian')
  need(all(targets[0][a][b]==targets[1][gp[a]][gp[b]] for a,b in combinations(range(n),2)),'padded literal target isomorphism')
  need(paths[0][k]==paths[1][k] and paths[0][k+7]==paths[1][k+7],'preserved core starts and ends')
  need(all(P[k+2:k+4]==(k+6,k+7) and P[k+5:k+7]==(k+5,k+4) for P in paths),'padded fixed shared arcs')
  need([i for i in range(n) if paths[0][i]!=paths[1][i]]==[k+1,k+4],'exact padded support')
  for P,C in zip(paths,targets):
   transpose=lambda M:list(map(list,zip(*M)))
   need(partner(transpose(B),P[::-1])==transpose(C),'padded converse square')
 compatible=[z for z in product(*padded_slots) if z[0][k+1:k+7]==z[1][k+1:k+7]]
 need(len(compatible)==2 and all(compatible[0][j]!=compatible[1][j] for j in (0,1)),'padded two-state coupling and one-slot obstruction')
 padded.append({'r':k,'s':suffix,'n':n})
need(all(not all(A[u][v]==A[p[u]][p[v]] for u,v in pairs) for p in permutations(range(8)) if p!=tuple(range(8))),'base source has trivial automorphisms')
result={'status':'PASS','checks':checks,'scope':'Independent native reconstructions and finite cross-checks, plus proof review in AUDIT.md; no new full-atlas classification or native wiring claim.','mod8_current_hostile':controls,'adaptive':reports,'endpoint_fibres':total_fibres,'selected_matching_assignments':matchings,'lift_components':[9,9,9],'repaired_components':[27],'spectral_residual_degrees':degrees,'rectangle_charge_orders':[4,4,4],'continuation_clocks':[3,160],'newton_boundary_monomials':newton_count,'stable_padding_hosts':padded,'input_sha256':hashes}
print(json.dumps(result,indent=2,sort_keys=True))
