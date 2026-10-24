#!/usr/bin/env python3
"""Independent exact audit: no imports of producer or original checker code."""
import argparse,csv,gzip,json,hashlib,random
from pathlib import Path
from collections import Counter,defaultdict
from itertools import product,combinations

COUNT=0
def require(ok,why):
 global COUNT
 COUNT+=1
 if not ok: raise RuntimeError(why)

def bits(x):
 while x:
  y=x&-x;yield y.bit_length()-1;x-=y

def binary_kernel(rows,n):
 """High pivots, independent of the source checker's low-pivot routine."""
 piv={}
 for z in rows:
  while z:
   p=z.bit_length()-1
   if p in piv:z^=piv[p]
   else:piv[p]=z;break
 free=[j for j in range(n)if j not in piv];basis=[]
 for j in free:
  x=1<<j
  for p,row in sorted(piv.items()):
   if (row&x).bit_count()%2:x^=1<<p
  basis.append(x)
 return basis,free,len(piv)

def rank(rows):return binary_kernel(rows,max((x.bit_length() for x in rows),default=0))[2]

def independent_filtration(repo,packet):
 base=repo/'05-knowledge/results';edgepath=base/'metagraph_order8_20261023/merged_edges.tsv.gz';receiptpath=base/'metagraph_order8_electrical_20261023/receipt.json'
 receipt=json.loads(receiptpath.read_text())['order8'];n=3528;root=n-1
 edges=[];parents=list(range(n));binary=[0]*root
 def find(i):
  while parents[i]!=i:parents[i]=parents[parents[i]];i=parents[i]
  return i
 for r in csv.DictReader(gzip.open(edgepath,'rt'),delimiter='\t'):
  i,j,b,k=(int(r[c])for c in['source','target','blue','black'])
  require(0<=i<=j<n and min(b,k)>=0 and k%2==0,'authentic edge schema/paired black')
  w=b+k;edges.append((i,j,w))
  if i==j or not w:continue
  parents[find(i)]=find(j)
  if w%2:
   if i<root:binary[i]^=1<<i
   if j<root:binary[j]^=1<<j
   if j<root:binary[i]^=1<<j;binary[j]^=1<<i
 require(len({find(i)for i in range(n)})==1,'full graph connected')
 basis,free,r0=binary_kernel(binary,root);d=len(basis)
 require(r0==164 and d==3363,'independent full binary matrix rank and nullity')
 for x in basis:require(all((row&x).bit_count()%2==0 for row in binary),'independent full binary kernel vector')
 labels=[0]*n
 for k,x in enumerate(basis):
  for i in bits(x):labels[i]|=1<<k
 # Q is B^T L B / 2 mod 2. Split every actual weight into 2 floor(w/2)+(w mod 2).
 # The even part is rank-one over F2. Accumulate the odd part over Z before dividing.
 Q=[0]*d;odd_positions=sorted({k for i,j,w in edges if w%2 for k in bits(labels[i]|labels[j])})
 loc={k:p for p,k in enumerate(odd_positions)};O=[[0]*len(loc)for _ in loc]
 for i,j,w in edges:
  z=labels[i]^labels[j]
  if (w//2)%2:
   for k in bits(z):Q[k]^=z
  if w%2:
   diff=[(loc[k],1 if labels[i]>>k&1 else -1)for k in bits(z)]
   for k,s in diff:
    for l,t in diff:O[k][l]+=s*t
 for k,row in enumerate(O):
  for l,value in enumerate(row):
   require(value%2==0,'odd-edge half energy integral before reduction')
   if value//2%2:Q[odd_positions[k]]^=1<<odd_positions[l]
 qbasis,qfree,r1=binary_kernel(Q,d)
 require(r1==3358 and len(qbasis)==5,'independent half-energy full bilinear rank')
 require(all(((Q[i]>>j)&1)==((Q[j]>>i)&1)for i in range(d)for j in range(i)),'full independent half form symmetric')
 def apply(x):
  out=[0]*n
  for i,j,w in edges:
   v=w*(x[i]-x[j]);out[i]+=v;out[j]-=v
  return out
 Y=receipt['radical_lifts_mod4'];Z=receipt['eighth_radical_lifts_mod8'][0];X=receipt['higher_layers'][-1]['lifts'][0]
 require(len(Y)==5 and all(len(x)==n and x[root]==0 for x in Y+[Z,X]),'rooted literal lift schema')
 AY=[apply(y)for y in Y];AZ=apply(Z);AX=apply(X)
 for y,ay in zip(Y,AY):
  require(all(t%4==0 for t in ay),'literal original L Y = 0 mod4')
  code=sum((y[i]%2)<<j for j,i in enumerate(free));rebuilt=0
  for k in bits(code):rebuilt^=basis[k]
  require(rebuilt==sum((v%2)<<i for i,v in enumerate(y[:root])),'Y belongs to independently recomputed binary basis')
  require(all((row&code).bit_count()%2==0 for row in Q),'Y parity lies in independent full half-form radical')
 require(rank([sum((v%2)<<i for i,v in enumerate(y))for y in Y])==5,'five independent Y parities span radical')
 require(all(v%8==0 for v in AZ)and all(v%256==0 for v in AX),'literal original Z and X harmonic equations')
 R=[]
 for y in Y:
  row=0
  for j,ay in enumerate(AY):
   e=sum(a*b for a,b in zip(y,ay));require(e%4==0,'literal quarter energy integral');row|=(e//4%2)<<j
  R.append(row)
 require(rank(R)==4,'independent quarter form rank four')
 require(all((X[j]-Y[0][j]-Y[4][j])%2==0 for j in range(n)),'X=Y0+Y4 mod2')
 require(all((X[j]-Z[j]-2*Y[0][j]-2*Y[3][j])%4==0 for j in range(n)),'literal X/Z mod4 relation')
 require(any(x%2 for x in X)and all((row&17).bit_count()%2==0 for row in R),'primitive X spans the quarter radical')
 energy=sum(x*v for x,v in zip(X,AX));require(energy%256==0 and energy//256%2==1,'independent final divided energy terminates at exponent256')
 ranks=[r1,rank(R)]+[energy//(1<<k)%2 for k in range(3,9)]
 require(ranks==[3358,4,0,0,0,0,0,1],'full independently certified dyadic valuation ranks')
 # Enumerate the marked harmonic generator modules as actual full vectors, not abstract coordinate names.
 tower=[]
 for k in range(2,9):
  m=1<<(k-1);half=m//2
  generated={bytes((half*sum(Y[h][j]for h in range(4)if mask>>h&1)+t*X[j])%m for j in range(n))for mask in range(16)for t in range(m)}
  target_log=sum(r*min(k,s)for s,r in enumerate(ranks,1))-d
  require(len(generated)==1<<target_log==1<<(k+3),'complete reduced harmonic tower generator cardinality')
  tower.append({'k':k,'modulus':1<<k,'reduced_image_order':len(generated),'predicted_log2_order':target_log})
 lab=[tuple([*(Y[h][j]%2 for h in range(4)),X[j]%128])for j in range(n)]
 saved=json.loads((packet/'node_labels.json').read_text());require(saved['labels']==[list(x)for x in lab]and saved['moduli']==[2,2,2,2,128],'all marked node coordinates')
 fibres=defaultdict(list)
 for j,l in enumerate(lab):fibres[l].append(j)
 dist=dict(sorted(Counter(map(len,fibres.values())).items()));require(len(fibres)==1695 and dist=={1:642,2:545,3:310,4:139,5:44,6:15},'independent 1695 occupied label fibres')
 require(saved['fibres']==[{'label':list(l),'vertices':v}for l,v in sorted(fibres.items())],'exact full fibre membership')
 return {'vertices':n,'edge_rows':len(edges),'binary_rank':r0,'binary_nullity':d,'odd_edge_kernel_coordinates':len(loc),'half_energy_rank':r1,'quarter_form_rows':R,'valuation_ranks':ranks,'higher_primitive_integer_energy':energy,'higher_primitive_energy_div_256':energy//256,'critical_two_primary':'(Z/2)^3358 + (Z/4)^4 + Z/256','paired_charge_type':'(Z/2)^4 + Z/128','paired_charge_order':2048,'tower':tower,'node_label_fibres':len(fibres),'fibre_size_distribution':dist,'independent_half_form_sha256':hashlib.sha256(json.dumps(Q).encode()).hexdigest()},lab

def lap(n,edges):
 A=[[0]*n for _ in range(n)]
 for i,j,w in edges:A[i][i]+=w;A[j][j]+=w;A[i][j]-=w;A[j][i]-=w
 return A

def matvec(A,x):return [sum(a*b for a,b in zip(row,x))for row in A]

def rectangle(n,s,t,u,v,a):
 U=[int(i==s)-int(i==t)for i in range(n)];V=[int(i==u)-int(i==v)for i in range(n)]
 return [[a*(U[i]*V[j]+V[i]*U[j])for j in range(n)]for i in range(n)]

def carry(A,X,width):
 c=[0]*len(A)
 for k in range(width):
  z=matvec(A,[(x>>k)&1 for x in X]);z=[a+b for a,b in zip(c,z)]
  if any(v%2 for v in z):return False,k
  c=[v//2 for v in z]
 return True,None

def finite_controls():
 # Exhaustive row-image vs full harmonic-bank tests, including odd-primary factors and overlapping roles.
 families=[(3,[(0,1,1),(1,2,1),(2,0,1)]),(3,[(0,1,3),(1,2,2)]),(4,[(0,1,1),(1,2,2),(2,3,2),(3,0,1)]),(4,[(0,1,2),(1,2,2),(2,3,2),(3,0,2)])]
 ntests=0;cuttests=0;carrytests=0
 for n,edges in families:
  A=lap(n,edges)
  for m in [2,4,8]:
   allvec=list(product(range(m),repeat=n));im={tuple(v%m for v in matvec(A,x))for x in allvec};H=[x for x in allvec if x[0]==0 and not any(v%m for v in matvec(A,x))]
   for (s,t),(u,v) in product(combinations(range(n),2),repeat=2):
    for a in [-2,-1,1,2]:
     F=rectangle(n,s,t,u,v,a)
     for b in range(4):
      direct=all(not any((1<<b)*z%m for z in matvec(F,x))for x in H)
      criterion=all(tuple((1<<b)*z%m for z in row)in im for row in F)
      require(direct==criterion,'all-size b row criterion finite exhaustive control');ntests+=1
     if m!=2:continue
     for y in product([0,1],repeat=n):
      degree=[sum(F[i][j]for j in range(n)if y[i]!=y[j])%2 for i in range(n)]
      require(degree==[v%2 for v in matvec(F,y)],'physical-row cut parity including loops/overlap');cuttests+=1
     for X in product(range(4),repeat=n):
      require(carry(F,X,2)[0]==(not any(v%4 for v in matvec(F,X))),'signed exact carry equals divisibility');carrytests+=1
 return {'row_image_full_harmonic_cases':ntests,'cut_cases':cuttests,'carry_cases':carrytests}

def original_menu_tests(lab):
 mod=[2,2,2,2,128];n=len(lab)
 def defects(rects,k=8):
  touched=sorted({x for r in rects for x in r[:4]});ix={x:i for i,x in enumerate(touched)};A=[[0]*len(touched)for _ in touched]
  for s,t,u,v,a in rects:
   B=rectangle(len(touched),ix[s],ix[t],ix[u],ix[v],a)
   A=[[x+y for x,y in zip(ar,br)]for ar,br in zip(A,B)]
  rows=[[sum(A[i][j]*lab[v][h]for j,v in enumerate(touched))%(2 if h<4 else 1<<min(k-1,7))for h in range(5)]for i in range(len(touched))]
  return rows,A,touched
 safe=[(913,2701,398,528,1)];old=[(1238,2246,881,2762,1)];shared=[(0,1,398,528,1),(55,46,398,528,1)]
 require(not any(any(r)for r in defects(safe)[0]),'claimed safe rectangle full bank')
 require(any(any(r)for r in defects(old)[0]),'old primitive-mod8 false-positive rejected')
 require(not any(any(r)for r in defects(shared)[0]),'two source-pair shared-target circuit')
 require(all(any(any(row)for row in defects([r])[0])for r in shared),'each circuit component individually unsafe')
 hostile=[(0,28,398,528,1)];rows,A,touch=defects(hostile)
 require(not any(any(r[:4])for r in rows)and carry(A,[lab[j][-1]for j in touch],7)==(False,1),'higher carry hostile')
 rng=random.Random(451024);trials=0
 for _ in range(120):
  vertices=rng.sample(range(n),rng.randrange(2,9));rects=[]
  for __ in range(rng.randrange(1,8)):
   s,t=rng.sample(vertices,2);u,v=rng.sample(vertices,2);rects.append((s,t,u,v,rng.choice([-3,-2,-1,0,1,2,3])))
  for k in range(2,9):
   rows,A,touch=defects(rects,k)
   c=carry(A,[lab[j][-1]for j in touch],min(k-1,7))
   require((not any(any(r)for r in rows))==(not any(any(r[:4])for r in rows)and c[0]),'tower label debt equals signed shared-row carry');trials+=1
 return {'safe_rectangle':safe,'old_rectangle_rejected':old,'shared_circuit':shared,'higher_carry_failure_bit':1,'random_multirectangle_tower_cases':trials}

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--packet',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 pins=json.loads((a.packet/'SOURCE_PINS.json').read_text());sources=[]
 for entry in pins['files']:
  data=(a.repo/entry['path']).read_bytes();require(hashlib.sha256(data).hexdigest()==entry['sha256']and len(data)==entry['bytes'],'pinned original source bytes');sources.append(entry)
 for path in sorted(a.packet.iterdir()):
  if path.is_file():sources.append({'path':'producer/'+path.name,'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 filtration,lab=independent_filtration(a.repo,a.packet);controls=finite_controls();menus=original_menu_tests(lab)
 result={'status':'PASS: independently certified original filtration, complete marked dyadic tower and shared-row cut/carry laws','checks':COUNT,'filtration':filtration,'finite_controls':controls,'original_menu_controls':menus,'source_hashes':sources,'scope':'All-dyadic inclusion of the original harmonic banks. No native resource authentication, odd-primary preservation, kernel equality, or unconditional all-order native supply claim.'}
 a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'status':'PASS','checks':COUNT,'output':str(a.out)}))
if __name__=='__main__':main()
