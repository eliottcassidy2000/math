#!/usr/bin/env python3
"""Complete marked paired-dyadic bank certificate; no optimization-sensitive asserts."""
import argparse,csv,gzip,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
HERE=Path(__file__).resolve().parent
CHECKS=0

def need(ok,why):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(why)

def rank_binary(vectors):
 piv={}
 for row in vectors:
  while row:
   p=row.bit_length()-1
   if p in piv:row^=piv[p]
   else:piv[p]=row;break
 return len(piv)

def compile_rectangles(rectangles):
 """A single integer row per physical node, even if roles overlap."""
 F=defaultdict(lambda:defaultdict(int))
 for r in rectangles:
  s,t=r['sources'];u,v=r['targets'];a=r.get('amplitude',1)
  if type(a)is not int or s==t or u==v:raise ValueError('invalid oriented rectangle')
  U={s:1,t:-1};V={u:1,v:-1}
  for i,x in U.items():
   for j,y in V.items():F[i][j]+=a*x*y;F[j][i]+=a*x*y
 F={i:{j:v for j,v in row.items() if v} for i,row in F.items() if any(row.values())}
 for i,row in F.items():
  need(sum(row.values())==0,'zero-row-sum compile')
  for j,v in row.items():need(F.get(j,{}).get(i,0)==v,'symmetric shared-row compile')
 return F

def apply_sparse(F,x):return {i:sum(c*x[j] for j,c in row.items()) for i,row in sorted(F.items())}

def label_defects(F,labels,moduli):
 return {i:[sum(c*labels[j][k] for j,c in row.items())%m for k,m in enumerate(moduli)]for i,row in sorted(F.items())}

def carry_certificate(F,X,bits):
 c={i:0 for i in F};layers=[]
 for t in range(bits):
  numerator={i:c[i]+sum(w*((X[j]>>t)&1) for j,w in row.items())for i,row in sorted(F.items())}
  bad=[[i,z%2]for i,z in numerator.items()if z%2]
  layers.append({'bit':t,'incoming_carry_nonzero':[[i,v]for i,v in sorted(c.items())if v],'numerator_nonzero':[[i,v]for i,v in numerator.items()if v],'odd_rows':bad})
  if bad:return {'passes':False,'first_failed_bit':t,'layers':layers}
  c={i:z//2 for i,z in numerator.items()}
 return {'passes':True,'layers':layers,'final_carry_nonzero':[[i,v]for i,v in sorted(c.items())if v]}

def cut_ledger(F,y):
 cross=defaultdict(int)
 for i,row in F.items():
  for j,w in row.items():
   if i<j and w%2 and (y[i]-y[j])%2:cross[i]^=1;cross[j]^=1
 literal={i:v%2 for i,v in apply_sparse(F,y).items()}
 need(all(literal[i]==cross[i]for i in F),'binary cut degree agrees with owned current')
 return sorted(i for i in F if cross[i])

def dense_laplacian(n,edges):
 L=[[0]*n for _ in range(n)]
 for i,j,w in edges:L[i][i]+=w;L[j][j]+=w;L[i][j]-=w;L[j][i]-=w
 return L

def small_controls():
 cases=0
 for n,edges in [(3,[(0,1,1),(1,2,2)]),(4,[(0,1,1),(0,3,1),(1,2,2),(2,3,2)]),(4,[(0,1,2),(1,2,2),(2,3,2),(3,0,2)])]:
  L=dense_laplacian(n,edges)
  for m in [4,8]:
   vectors=list(product(range(m),repeat=n))
   image={tuple(sum(a*b for a,b in zip(row,x))%m for row in L)for x in vectors}
   H=[x for x in vectors if x[0]==0 and all(sum(a*b for a,b in zip(row,x))%m==0 for row in L)]
   for a,b in combinations(list(combinations(range(n),2)),2):
    for amp in [-1,1,2]:
     F=compile_rectangles([{'sources':a,'targets':b,'amplitude':amp}])
     literal=all(all(2*z%m==0 for z in apply_sparse(F,x).values())for x in H)
     rowcrit=all(tuple(2*row.get(j,0)%m for j in range(n))in image for row in F.values())
     need(literal==rowcrit,'finite full-bank image criterion including overlapping roles')
     for y in product(range(2),repeat=n):cut_ledger(F,y)
     for X in product(range(4),repeat=n):
      carry=carry_certificate(F,X,2)
      need(carry['passes']==all(z%4==0 for z in apply_sparse(F,X).values()),'integer carry exactness including signed coefficients')
     cases+=1
 # Weighted-C4 hostile: scalar equality is not owned current equality.
 L=dense_laplacian(4,[(0,1,1),(0,3,1),(1,2,2),(2,3,2)])
 H=[(0,)+x for x in product(range(8),repeat=3)if all(sum(a*b for a,b in zip(row,(0,)+x))%8==0 for row in L)]
 need(H==[(0,0,0,0),(0,0,4,0),(0,4,2,4),(0,4,6,4)],'weighted C4 entire harmonic bank')
 need(all(x[1]==x[3]for x in H),'weighted C4 scalar alias')
 def current(t,x):
  out=[0]*4;out[0]=(x[0]-x[t])%8;out[t]=(-out[0])%8;return out
 need(current(1,H[2])==[4,4,0,0]and current(3,H[2])==[4,0,0,4],'weighted C4 target-owned current hostile')
 F=compile_rectangles([{'sources':[0,1],'targets':[0,1]}])
 need(all(z%8==0 for z in apply_sparse(F,[0,4]).values()),'overlap two terms cancel')
 return {'finite_ring_cases':cases,'weighted_C4_bank':H,'weighted_C4_currents':[current(1,H[2]),current(3,H[2])]}

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path);a=p.parse_args()
 pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
 for entry in pins['files']:need(sha256((a.repo/entry['path']).read_bytes()).hexdigest()==entry['sha256'],'source hash '+entry['path'])
 base=a.repo/'05-knowledge/results';E=base/'metagraph_order8_20261023/merged_edges.tsv.gz';D=json.loads((base/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
 N=D['vertices'];root=D['root'];edges=[tuple(int(r[k])for k in['source','target','blue','black'])for r in csv.DictReader(gzip.decompress(E.read_bytes()).decode().splitlines(),delimiter='\t')]
 def apply(x):
  z=[0]*N
  for i,j,b,k in edges:w=(b+k)*(x[i]-x[j]);z[i]+=w;z[j]-=w
  return z
 need(N==3528 and root==3527 and len(edges)==438216,'exact original graph universe')
 Y=D['radical_lifts_mod4'];Z=D['eighth_radical_lifts_mod8'][0];X=D['higher_layers'][-1]['lifts'][0]
 need(len(Y)==5,'five saved radical lifts')
 for i,y in enumerate(Y):need(len(y)==N and y[root]==0 and all(v%4==0 for v in apply(y)),'literal Y harmonic '+str(i))
 need(Z[root]==X[root]==0 and all(z%8==0 for z in apply(Z))and all(z%256==0 for z in apply(X)),'literal Z and X harmonic')
 rank=rank_binary([sum((v%2)<<j for j,v in enumerate(y))for y in Y]);need(rank==5,'five independent mod2 lifts')
 for j in range(N):
  need((X[j]-Y[0][j]-Y[4][j])%2==0,'X parity relation')
  need((X[j]-Z[j]-2*Y[0][j]-2*Y[3][j])%4==0,'X/Z integral reduction relation')
 need(any(v%2 for v in X),'primitive X has odd coordinate')
 need(D['complete_two_primary']and D['critical_two_rank']==3363 and D['exact_valuation_one']==3358 and D['exact_valuation_two']==4 and D['exact_valuation_three']==0,'pinned complete filtration counts')
 need([(r['level'],r['dimension'],r['rank'],r['next_dimension'])for r in D['higher_layers']]==[(4,1,0,1),(5,1,0,1),(6,1,0,1),(7,1,0,1),(8,1,1,0)],'filtration terminates at exponent256')
 reduction_log_order=3358+2*4+8-3363;need(reduction_log_order==11,'complete H256 reduction cardinality')
 # Obtain size2048 without storing huge numbers of full harmonic vectors.
 even_group={bytes(2*sum(y[j]for i,y in enumerate(Y)if b>>i&1)%4 for j in range(N))for b in range(32)}
 need(len(even_group)==32 and bytes(2*v%4 for v in X)in even_group,'cyclic/even generator intersection order2')
 reduced8={bytes((e[j]+t*Z[j])%4 for j in range(N))for e in even_group for t in range(4)}
 need(len(reduced8)==64,'actual complete H8 reduced bank cardinality64')
 labels=[[*(Y[i][j]%2 for i in range(4)),X[j]%128]for j in range(N)];moduli=[2,2,2,2,128]
 saved=json.loads((HERE/'node_labels.json').read_text());need(saved['labels']==labels and saved['moduli']==moduli,'full node coordinate map matches direct original lifts')
 fibres=defaultdict(list)
 for j,l in enumerate(labels):fibres[tuple(l)].append(j)
 distribution=dict(sorted(Counter(map(len,fibres.values())).items()))
 need(len(fibres)==1695 and distribution=={1:642,2:545,3:310,4:139,5:44,6:15},'complete paired-charge label census')
 need(saved['fibres']==[{'label':list(k),'vertices':v}for k,v in sorted(fibres.items())],'full fibres equal recomputed map')
 controls=small_controls()
 def certificate(rectangles):
  F=compile_rectangles(rectangles);defects=label_defects(F,labels,moduli);cuts=[cut_ledger(F,[y[j]%2 for j in range(N)])for y in Y[:4]]
  carry=carry_certificate(F,X,7)
  good=not any(any(d)for d in defects.values())
  need(good==(not any(cuts)and carry['passes']),'marked module equals cut and carry compiler')
  # Direct generated reduced H256 verification at all physical touched rows.
  all_good=True;generators_checked=0
  for b in range(16):
   for t in range(128):
    for i,row in F.items():
     debt=sum(w*(64*sum(Y[h][j]for h in range(4)if b>>h&1)+t*X[j])for j,w in row.items())
     if 2*debt%256:all_good=False
    generators_checked+=1
  need(all_good==good,'direct entire reduced H256 bank equals compiler')
  return {'rectangles':rectangles,'integer_F_rows':[[i,sorted(row.items())]for i,row in sorted(F.items())],'passes_all_original_dyadic_banks':good,'row_label_defects':defects,'four_binary_cut_odd_rows':cuts,'carry':carry,'direct_reduced_characters_checked':generators_checked}
 old=certificate([{'sources':[1238,2246],'targets':[881,2762]}]);need(not old['passes_all_original_dyadic_banks'],'old unit rectangle fails full bank')
 F=compile_rectangles(old['rectangles']);need(all(2*z%8==0 for z in apply_sparse(F,X).values()),'old selected primitive modulus8 still passes')
 safe=certificate([{'sources':[913,2701],'targets':[398,528]}]);need(safe['passes_all_original_dyadic_banks'],'new candidate is all-dyadic electrically safe')
 # Pick65 distinct source labels, avoiding the safe target pair; pair-sum collision.
 reps=[];seen=set()
 for j,label in enumerate(labels):
  if j in [398,528]or tuple(label)in seen:continue
  seen.add(tuple(label));reps.append(j)
  if len(reps)==65:break
 sums={};witness=None
 for p,q in combinations(reps,2):
  label=tuple((labels[p][i]+labels[q][i])%moduli[i]for i in range(5))
  if label in sums:witness=(*sums[label],p,q);break
  sums[label]=(p,q)
 need(witness is not None and len(set(witness))==4,'finite module pair-sum collision has disjoint source pairs')
 A,B,C,D=witness;rects=[{'sources':[A,C],'targets':[398,528]},{'sources':[B,D],'targets':[398,528]}]
 circuit=certificate(rects);individual=[certificate([r])for r in rects]
 need(circuit['passes_all_original_dyadic_banks']and not any(c['passes_all_original_dyadic_banks']for c in individual),'two individually unsafe source-pair rectangles cancel at shared target rows')
 carry_hostile=certificate([{'sources':[0,28],'targets':[398,528]}]);need(not carry_hostile['passes_all_original_dyadic_banks']and not any(carry_hostile['four_binary_cut_odd_rows'])and carry_hostile['carry']['first_failed_bit']==1,'higher carry hostile passes all lower parity cuts then fails bit1')
 # Direct fixed-source tariff on the old two disjoint target rectangles.
 tariff=0
 for a1,a2 in product(range(-2,3),repeat=2):
  F=compile_rectangles([{'sources':[1238,2246],'targets':[881,2762],'amplitude':a1},{'sources':[1238,2246],'targets':[2650,2843],'amplitude':a2}])
  bad=any(any(x)for x in label_defects(F,labels,[2,2,2,2,4]).values())
  need((not bad)==(a1%2==0 and a2%2==0),'exact original two-rectangle fullmod8 parity tariff');tariff+=1
 result={'status':'PROVED compressed dyadic module; FINITE-EXACT original atlas controls','checks':CHECKS,'source_pins_sha256':sha256((HERE/'SOURCE_PINS.json').read_bytes()).hexdigest(),'node_labels_sha256':sha256((HERE/'node_labels.json').read_bytes()).hexdigest(),'original_edge_rows':len(edges),'paired_charge_group':{'type':'(Z/2)^4 + Z/128','order':2048,'full_H256_log2_order':3374,'reduction_kernel_log2_order':3363,'binary_Y_rank':rank,'full_H8_reduced_order':64,'native_label_fibres':len(fibres),'fibre_size_distribution':distribution},'small_controls':controls,'old_rectangle':old,'safe_rectangle_electrical_only':safe,'shared_target_source_circuit_algebraic_only':{'source_partition':[[A,B],[C,D]],'combined':circuit,'individual':individual},'higher_carry_hostile':carry_hostile,'old_two_rectangle_tariff_cases':tariff,'scope':'All dyadic original-bank inclusion. No odd-primary claim. Native authentication and injective resource availability are separate certificates.'}
 raw=json.dumps(result,indent=2,sort_keys=True)+'\n'
 if a.out:a.out.write_text(raw)
 else:print(raw,end='')
if __name__=='__main__':main()
