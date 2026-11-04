"""Exact physical small-model extraction and complete modulo-eight type controls."""
from pathlib import Path
from itertools import product,combinations
from functools import lru_cache
import hashlib,importlib.util,json,random,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('u',ROOT.parent/'uniform_block_starts'/'verify.py');u=importlib.util.module_from_spec(s);s.loader.exec_module(u)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def physical(mask,blocks,offs):return sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,offs)) if mask>>i&1)
def run_count(R,k):
 counts=tuple(k)+(1,1)
 @lru_cache(None)
 def go(left,last):
  if not any(left):return 1
  return sum(go(left[:v]+(left[v]-1,)+left[v+1:],v) for v,t in enumerate(left) if t and (last<0 or R[last]>>v&1))
 return go(counts,-1)
def run_words(R,k):
 left=list(k)+[1,1];out=[]
 def go(word):
  if not any(left):out.append(tuple(word));return
  for v,t in enumerate(left):
   if t and (not word or R[word[-1]]>>v&1):
    left[v]-=1;go(word+[v]);left[v]+=1
 go([]);return out
def evaluate(coeff,types):
 states={'S':(1,0,0),'P':(0,0,0),'E':(0,1,0),'O':(0,1,1),'C':(0,0,1)}
 ss=[states[t] for t in types];v=coeff['K0']
 for i,(s,b,c) in enumerate(ss):v+=(2+4*b-2*s)*coeff['K2'][i]+2*c*coeff['K3'][i]
 for i,j,k in coeff['K22']:v+=4*(1-ss[i][0])*(1-ss[j][0])*k
 return v%8
def rank2(vectors):
 basis={}
 for row in vectors:
  x=sum((z&1)<<i for i,z in enumerate(row))
  while x:
   p=x.bit_length()-1
   if p not in basis:basis[p]=x;break
   x^=basis[p]
 return len(basis)
def main():
 rng=random.Random(188335);small=[]
 for n in range(1,4):
  for code in range(1<<(n*(n-1)//2)):small.append((f'all{n}-{code}',u.rows(n,code),list(product(range(1<<n),repeat=2))))
 extras=[('prime5a',[20,1,2,7,14]),('prime5b',[18,8,3,5,14]),('prime5c',[24,17,3,6,12]),('prime8_endpoint_blind',[238,4,240,118,195,146,34,74])]
 extras += [(f'random{n}',u.rows(n,rng.randrange(1<<(n*(n-1)//2)))) for n in (4,6,7)]
 extras += [('transitive8',u.tt(8))]
 for name,Q in extras:
  m=len(Q);M=(1<<m)-1;qs=[(0,0),(M,M),(1,0),(0,1),(1,3),(M^1,M^3)]
  while len(qs)<12:
   q=(rng.randrange(M+1),rng.randrange(M+1))
   if q not in qs:qs.append(q)
  small.append((name,Q,qs))
 cases=[];model_records=[]
 for name,Q,queries in small:
  m=len(Q);models=[('base',[1]*m)]
  for i in range(m):
   for r in (2,3):sizes=[1]*m;sizes[i]=r;models.append((f'{i}:{r}',sizes))
  for i,j in combinations(range(m),2):sizes=[1]*m;sizes[i]=sizes[j]=2;models.append((f'{i},{j}:2',sizes))
  ids=[]
  for label,sizes in models:
   blocks=[u.tt(n) for n in sizes];T,offs=u.blow(Q,blocks);U,_=u.blow(u.op(Q),blocks);pq=[(physical(a,blocks,offs),physical(c,blocks,offs)) for a,c in queries];ids.append(len(cases));cases.append((T,U,pq))
  model_records.append({'name':name,'quotient_rows':Q,'block_queries':queries,'models':[{'label':label,'sizes':sizes} for label,sizes in models],'case_ids':ids})
 frozen=json.loads((ROOT.parent/'five_profile_types'/'normal.json').read_text())
 # Actual noncanonical children, with odd Hamilton units and singleton states.
 variants=[]
 variant_inputs=[]
 for r in frozen['actual_profile_transfer_controls']:variant_inputs.append(('cyclic_larger', [2,4,1], r['child_rows'],list(product(range(8),repeat=2))))
 # Quotients not covered by the cyclic classification; keep physical order modest.
 for name,Q in extras:
  m=len(Q);blocks=[u.tt(1) for _ in Q]
  blocks[0]=[2,4,1];blocks[-1]=[8,1,3,6]
  if m<=5:blocks[1]=u.tt(2)
  queries=next(z[2] for z in small if z[0]==name)
  variant_inputs.append((name,Q,blocks,queries))
 for name,Q,blocks,queries in variant_inputs:
  T,offs=u.blow(Q,blocks);U,_=u.blow(u.op(Q),blocks);pq=[(physical(a,blocks,offs),physical(c,blocks,offs)) for a,c in queries];qs=[u.profile(B) for B in blocks];types=[];unit=1
  for B,q in zip(blocks,qs):
   unit=unit*q[1]%8;r=[(q[k] if k<len(q) else 0)*pow(q[1],-1,8)%8 for k in (1,2,3)];typ=next(t for t,v in zip('SPEOC',[[1,0,0],[1,2,0],[1,6,4],[1,6,6],[1,2,2]]) if r==v);types.append(typ)
  variants.append({'name':name,'quotient_rows':Q,'child_rows':blocks,'profiles':qs,'types':types,'odd_H_product_mod8':unit,'block_queries':queries,'physical_queries':pq,'case_id':len(cases)});cases.append((T,U,pq))
 lines=[str(len(cases))]
 for T,U,qs in cases:lines += [f'{len(T)} {len(qs)}',' '.join(map(str,T)),' '.join(map(str,U))]+[f'{a} {c}' for a,c in qs]
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'literal.cpp'),'-o',str(exe)],check=True);values=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 at=0;vals=[]
 for T,U,qs in cases:
  vals.append([values[at+2*i:at+2*i+2] for i in range(len(qs))]);at+=2*len(qs)
 need(at==len(values),'every physical count consumed');lookup={};word_coefficient_checks=0
 for r in model_records:
  vv=[vals[i] for i in r.pop('case_ids')]
  for md,v in zip(r['models'],vv):md['physical_values']=v
  Q=r['quotient_rows'];m=len(Q);out=[]
  for qi,(a,c) in enumerate(r['block_queries']):
   pair=[]
   for side,P in enumerate((Q,u.op(Q))):
    f0=vv[0][qi][side];K2=[];K3=[]
    for i in range(m):
     f2=vv[1+2*i][qi][side];f3=vv[2+2*i][qi][side];need((f2-f0)%2==0,'two-run exact extraction');need((f3-3*f2+2*f0)%12==0,'three-run even extraction');K2.append((f2-f0)//2);K3.append((f3-3*f2+2*f0)//6)
    K22=[]
    for j,(i,k) in enumerate(combinations(range(m),2)):
     f22=vv[1+2*m+j][qi][side];num=f22-vv[1+2*i][qi][side]-vv[1+2*k][qi][side]+f0;need(num%4==0,'two-child exact extraction');K22.append([i,k,num//4])
    coeff={'K0':f0,'K2':K2,'K3':K3,'K22':K22};need(min([f0]+K2+K3+[z[2] for z in K22])>=0,'nonnegative literal run counts');pair.append(coeff)
    if m<=3 or qi<2:
     R=u.extend(P,a,c);patterns=[([1]*m,f0)]
     for i in range(m):
      for visits,cnt in [(2,K2[i]),(3,K3[i])]:ks=[1]*m;ks[i]=visits;patterns.append((ks,cnt))
     for i,j,cnt in K22:ks=[1]*m;ks[i]=ks[j]=2;patterns.append((ks,cnt))
     for ks,cnt in patterns:need(run_count(R,ks)==cnt,'independent run-state DP equals physical extraction');word_coefficient_checks+=1
   out.append({'block_query':[a,c],'coefficients':pair});lookup[(tuple(Q),a,c)]=pair
  r['extracted_coefficients']=out
 # Full canonical banks retain both integer source tables and residue predictions.
 canonical=[];normalized={}
 for record in frozen['canonical_records']:
  types=record['types'];units=pow(3,types.count('C'),8);rows=[];nvals=[]
  for (a,c),vv in zip(product(range(8),repeat=2),record['values']):
   cs=lookup[(tuple([2,4,1]),a,c)];pred=[units*evaluate(z,types)%8 for z in cs];need(pred==[v%8 for v in vv],'all16000 canonical integer counts reconstructed modulo8');rows.append([a,c,*vv,*pred]);nvals.append([(v*pow(units,-1,8))%8 for v in vv])
  canonical.append({'types':types,'odd_H_product_mod8':units,'records':rows});normalized[tuple(types)]=nvals
 for r in variants:
  r['physical_values']=vals[r.pop('case_id')];pred=[]
  for (a,c),vv in zip(r['block_queries'],r['physical_values']):
   cs=lookup[(tuple(r['quotient_rows']),a,c)];w=[r['odd_H_product_mod8']*evaluate(z,r['types'])%8 for z in cs];need(w==[z%8 for z in vv],'noncanonical actual profile physical restoration');pred.append(w)
  r['predictions_mod8']=pred
 # Independently enumerate, transform and certify every retained return word.
 inv=[];witnesses={};invwords=0
 for name,Q,queries in small:
  m=len(Q)
  if m>3:continue
  for a,c in queries:
   R=u.extend(Q,a,c)
   for i in range(m):
    patterns=[('i3',None,[3 if v==i else 1 for v in range(m)])]
    for j in range(m):
     if i==j:continue
     for times in (2,3):ks=[1]*m;ks[i]=3;ks[j]=times;patterns.append((f'i3j{times}',j,ks))
    for kind,j,ks in patterns:
     words=run_words(R,ks);W=set(words);need(len(W)==len(words) and len(W)%2==0,'word multiplicity and parity');mapped=[]
     for w in words:
      p=[k for k,x in enumerate(w) if x==i];l0,l1,l2=p;z=w[:l0+1]+w[l1+1:l2]+(i,)+w[l0+1:l1]+w[l2:];need(z in W and z!=w,'actual fixed-point-free return involution');pp=[k for k,x in enumerate(z) if x==i];r0,r1,r2=pp;zz=z[:r0+1]+z[r1+1:r2]+(i,)+z[r0+1:r1]+z[r2:];need(zz==w,'involution squares to identity');mapped.append(z)
     need(len(set(mapped))==len(W),'return map bijective');invwords+=len(W)
     digest=hashlib.sha256(json.dumps(words,separators=(',',':')).encode()).hexdigest();inv.append([name,a,c,i,j,kind,len(W),digest])
     if words and kind not in witnesses:witnesses[kind]={'quotient_rows':Q,'block_query':[a,c],'visits':ks,'words':words,'images':mapped}
 need(set(witnesses)=={'i3','i3j2','i3j3'},'positive controls for each involution clause')
 # The complete cyclic contrast has six independent nontrivial type derivatives.
 base=normalized[('P',)*3];basecontrast=[(a-b)%8 for a,b in base];derivatives=[]
 for i in range(3):
  for typ in ('E','C'):
   ts=['P']*3;ts[i]=typ;dv=[((a-b)-basecontrast[j])%8 for j,(a,b) in enumerate(normalized[tuple(ts)])];need(set(dv)<=set([0,4]),'binary contrast derivatives');derivatives.append([v//4 for v in dv])
 rank=rank2(derivatives);need(rank==6,'cyclic six-bit contrast rank');bankset=set();blind=[]
 for ts in product('PEOC',repeat=3):
  contrast=tuple((a-b)%8 for a,b in normalized[ts]);bankset.add(contrast)
  if not any(contrast):blind.append(ts)
 need(len(bankset)==64 and blind==[('E','E','E')],'all nontrivial normalized cyclic banks distinct and unique blind')
 dep=['five_profile_types/PROOF.md','five_profile_types/normal.json','uniform_block_starts/verify.py']
 return {'status':'PASS','source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','literal.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in dep},'physical_H_counts':len(values),'small_model_controls':model_records,'independent_run_coefficient_checks':word_coefficient_checks,'canonical_full_bank_columns':['a','c','H_source','H_target','predicted_source_mod8','predicted_target_mod8'],'canonical_full_banks':canonical,'actual_profile_variants':variants,'involution_columns':['quotient_name','a','c','i','j','pattern','words','sha256'],'involution_instances':inv,'total_mapped_return_words':invwords,'positive_involution_witnesses':witnesses,'cyclic_nontrivial_contrast_derivatives':derivatives,'cyclic_nontrivial_contrast_rank':rank,'distinct_cyclic_nontrivial_banks':len(bankset),'unique_cyclic_blind_type':blind,'scope':'Complete block-uniform response reduction for arbitrary quotients and nonempty children modulo8. Exact small-model extraction and scalar restoration, not a separation or arbitrary-ambient nonannihilation theorem.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
