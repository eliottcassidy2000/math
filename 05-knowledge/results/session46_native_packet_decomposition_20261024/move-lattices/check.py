"""Direct integer witnesses and complete low-cardinality native charge fibers."""
import argparse,collections,hashlib,itertools,json
from fractions import Fraction
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--supply',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();HERE=Path(__file__).resolve().parent;checks=0

def need(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)

raw=(a.repo/'05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json').read_bytes();D=json.loads(raw)['order8'];Y=D['radical_lifts_mod4'];X=D['higher_layers'][-1]['lifts'][0]
rows=json.loads(a.supply.read_text())['records'];ids=[r['source'] for r in rows];n=len(ids);need(n==28 and ids==sorted(set(ids)),'actual ordered28')
L=[]
for r in rows:
 s=r['source'];label=[Y[j][s]%2 for j in range(4)]+[X[s]%128];need(label==r['source_label'],'source charge authenticated');L.append((sum(label[j]<<j for j in range(4)),label[4]))
need(len(set(L))==n,'no two-slot relation between distinct labels')

def key(aa):
 q=0;x=0
 for i in aa:q^=L[i][0];x+=L[i][1]
 return q,x%128

def moves(k):
 groups=collections.defaultdict(list)
 for aa in itertools.combinations(range(n),k):groups[key(aa)].append(aa)
 out=[]
 for pp in groups.values():
  for aa,bb in itertools.combinations(pp,2):
   if set(aa).isdisjoint(bb):out.append([int(i in aa)-int(i in bb) for i in range(n)])
 return out

four=moves(2);six=moves(3);need(len(four)==27 and len(six)==1917,'complete move counts')
for v in four+six:
 need(sum(v)==0,'zero total margin')
 need(sum(v[i]*L[i][1] for i in range(n))%128==0 and all(sum(v[i]*((L[i][0]>>j)&1) for i in range(n))%2==0 for j in range(4)),'complete paired charge guard')

def rank(A):
 A=[[Fraction(x) for x in r] for r in A];r=0
 for j in range(len(A[0])):
  k=next((k for k in range(r,len(A)) if A[k][j]),None)
  if k is None:continue
  A[r],A[k]=A[k],A[r];d=A[r][j];A[r]=[x/d for x in A[r]]
  for k in range(r+1,len(A)):
   if A[k][j]:d=A[k][j];A[k]=[x-d*y for x,y in zip(A[k],A[r])]
  r+=1
  if r==len(A):break
 return r
need(rank(four)==23,'four-slot rational rank23')
pair=(ids.index(1638),ids.index(1974));need(all(v[pair[0]]+v[pair[1]]==0 for v in four),'critical integer invariant')
# The degree-zero charge map has image all2048 elements: finite BFS from differences.
gens=[(L[i][0]^L[-1][0],(L[i][1]-L[-1][1])%128) for i in range(n-1)];seen={(0,0)};todo=[(0,0)]
while todo:
 p0,x0=todo.pop()
 for p1,x1 in gens:
  z=(p0^p1,(x0+x1)%128)
  if z not in seen:seen.add(z);todo.append(z)
need(len(seen)==2048,'full zero-margin group image')
W=json.loads((HERE/'integer_basis_witness.json').read_text());V=(four+six)[:W['prefix_moves']];B=W['basis_columns'];C=W['integer_combinations'];need(len(V)==48 and len(B)==len(C)==27,'certificate dimensions')
for b,c in zip(B,C):
 need(len(c)==48 and len(b)==27 and all(type(x)is int for x in b+c),'exact integer certificate')
 need([sum(c[k]*V[k][i] for k in range(48)) for i in range(27)]==b,'direct integer linear combination')

def det(A):
 A=[r[:] for r in A];last=1;sign=1
 for k in range(len(A)-1):
  j=next((j for j in range(k,len(A)) if A[j][k]),None)
  if j is None:return 0
  if j!=k:A[k],A[j]=A[j],A[k];sign=-sign
  pivot=A[k][k]
  for i in range(k+1,len(A)):
   for j in range(k+1,len(A)):
    num=A[i][j]*pivot-A[i][k]*A[k][j];need(num%last==0,'Bareiss exact division');A[i][j]=num//last
   A[i][k]=0
  last=pivot
 return sign*A[-1][-1]
need(abs(det(list(map(list,zip(*B)))))==2048,'generated sublattice determinant equals full charge-kernel index')
# All fibers at cardinality2,3,4; legal <=6 moves are exactly pair distances<=6.
reports=[];bad4=[]
for k in (2,3,4):
 groups=collections.defaultdict(list)
 for aa in itertools.combinations(range(n),k):groups[key(aa)].append(sum(1<<i for i in aa))
 dist=collections.Counter()
 for charge,states in groups.items():
  remaining=set(states);components=[]
  while remaining:
   seed=remaining.pop();got={seed};stack=[seed]
   while stack:
    x=stack.pop();adj={y for y in remaining if (x^y).bit_count()<=6};remaining.difference_update(adj);got.update(adj);stack.extend(adj)
   components.append(sorted(got))
  dist[len(components)]+=1
  if k==4 and len(components)>1:bad4.append({'charge':list(charge),'components':components})
 reports.append({'selected':k,'states':sum(map(len,groups.values())),'fibers':len(groups),'component_distribution':dict(sorted(dist.items()))})
need(reports[0]['component_distribution']=={1:351},'cardinality2 fibers')
need(reports[1]['component_distribution']=={1:1646},'cardinality3 fibers')
need(reports[2]['component_distribution']=={1:1851,2:175,3:21,4:1},'full cardinality4 components')
h=next(z for z in bad4 if z['charge']==[7,15]);need(sorted(h['components'])==[[32901],[136318978]],'complete two-state capacity hostile')
u,v=32901,136318978;need((u^v).bit_count()==8 and (u&v)==0,'sharp eight-slot move')
result={'status':'PASS','checks':checks,'source_order':ids,'four_slot_moves':27,'six_slot_moves':1917,'four_move_rank':23,'four_move_integer_invariant_sources':[1638,1974],'zero_margin_group_image':2048,'integer_move_lattice_index':2048,'integer_move_lattice_equals_full_charge_kernel':True,'low_cardinality_fibers':reports,'binary_hostile':{'charge':[7,15],'masks':[u,v],'selected_sources':[[ids[i] for i in range(n) if mask>>i&1] for mask in(u,v)],'sharp_move_support':8},'pins':{'original_electrical_receipt':hashlib.sha256(raw).hexdigest(),'inherited_supply':hashlib.sha256(a.supply.read_bytes()).hexdigest(),'integer_basis_witness':hashlib.sha256((HERE/'integer_basis_witness.json').read_bytes()).hexdigest()},'scope':'Complete28-source target-margin binary menu, with native realization supplied by separately audited pin-domain certificate. Integer lattice equality does not imply binary-box connectivity; no all-fiber6 connectivity claim.'}
a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
