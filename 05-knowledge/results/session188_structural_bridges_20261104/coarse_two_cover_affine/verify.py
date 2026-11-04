"""Complete physical two-cover type compression with exact path-swap controls."""
from pathlib import Path
from itertools import product,combinations
from functools import lru_cache
import hashlib,importlib.util,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('u',ROOT.parent/'uniform_block_starts'/'verify.py');u=importlib.util.module_from_spec(s);s.loader.exec_module(u)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def physical(mask,blocks,offs):return sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,offs)) if mask>>i&1)
def coeff_count(R,k,h=2):
 @lru_cache(None)
 def go(left,last,cuts):
  if not any(left):return int(cuts==h-1)
  z=0
  for v,t in enumerate(left):
   if not t:continue
   ll=left[:v]+(t-1,)+left[v+1:]
   if last<0 or R[last]>>v&1:z+=go(ll,v,cuts)
   if last>=0 and cuts<h-1:z+=go(ll,v,cuts+1)
  return z
 return go(tuple(k)+(1,1),-1,0)
def word_covers(R,k,h=2):
 left=list(k)+[1,1];out=[]
 def go(words):
  if not any(left):
   if len(words)==h:out.append(tuple(tuple(w) for w in words))
   return
  for v,t in enumerate(left):
   if not t:continue
   left[v]-=1
   if not words:go([[v]])
   else:
    if R[words[-1][-1]]>>v&1:go(words[:-1]+[words[-1]+[v]])
    if len(words)<h:go(words+[[v]])
   left[v]+=1
 go([]);return out
def evaluate(coeff,types):
 cbit={'S':0,'P':0,'E':0,'O':1,'C':1}
 return (coeff['B0']+2*sum((t!='S')*b+cbit[t]*c for t,b,c in zip(types,coeff['B2'],coeff['B3'])))%8
def main():
 old=json.loads((ROOT.parent/'coarse_mod8_affine'/'normal.json').read_text());frozen=json.loads((ROOT.parent/'five_profile_types'/'normal.json').read_text());cases=[];small=[]
 for r in old['small_model_controls']:
  Q=r['quotient_rows'];m=len(Q);qs=r['block_queries'];models=[md for md in r['models'] if sum(n>1 for n in md['sizes'])<=1];ids=[]
  for md in models:
   blocks=[u.tt(n) for n in md['sizes']];T,offs=u.blow(Q,blocks);U,_=u.blow(u.op(Q),blocks);pq=[(physical(a,blocks,offs),physical(c,blocks,offs)) for a,c in qs];ids.append(len(cases));cases.append((T,U,pq))
  small.append({'name':r['name'],'quotient_rows':Q,'block_queries':qs,'models':[{'sizes':md['sizes'],'label':md['label']} for md in models],'ids':ids})
 canonical=[]
 for r in frozen['canonical_records']:
  canonical.append({'types':r['types'],'child_rows':r['child_rows'],'physical_queries':r['physical_queries'],'id':len(cases),'expected_H_values':r['values']});cases.append((r['source_old_rows'],r['target_old_rows'],r['physical_queries']))
 variants=[]
 for r in old['actual_profile_variants']:
  blocks=r['child_rows'];Q=r['quotient_rows'];T,offs=u.blow(Q,blocks);U,_=u.blow(u.op(Q),blocks);variants.append({k:r[k] for k in ['name','quotient_rows','child_rows','profiles','types','odd_H_product_mod8','block_queries','physical_queries']});variants[-1]['id']=len(cases);cases.append((T,U,r['physical_queries']))
 lines=[str(len(cases))]
 for T,U,qs in cases:lines += [f'{len(T)} {len(qs)}',' '.join(map(str,T)),' '.join(map(str,U))]+[f'{a} {c}' for a,c in qs]
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);values=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 vals=[];at=0
 for T,U,qs in cases:vals.append([values[at+4*i:at+4*i+4] for i in range(len(qs))]);at+=4*len(qs)
 need(at==len(values),'all actual physical H and q2 counts consumed');lookup={};independent_checks=0
 for r in small:
  vv=[vals[i] for i in r.pop('ids')]
  for md,v in zip(r['models'],vv):md['values_H_q2_source_target']=v
  Q=r['quotient_rows'];m=len(Q);records=[]
  for qi,(a,c) in enumerate(r['block_queries']):
   cs=[]
   for side,P in enumerate((Q,u.op(Q))):
    pos=2*side+1;v0=vv[0][qi][pos];b2=[];b3=[];need(v0%2==0,'base coefficient even')
    for i in range(m):
     v2=vv[1+2*i][qi][pos];v3=vv[2+2*i][qi][pos];need((v2-v0)%4==0,'two-run numerator divisible4');need((v3-3*v2+2*v0)%12==0,'three-run numerator divisible12');b2.append((v2-v0)//2);b3.append((v3-3*v2+2*v0)//6)
    cc={'B0':v0,'B2':b2,'B3':b3};need(min([v0]+b2+b3)>=0,'nonnegative extracted coefficients');cs.append(cc)
    if m<=3 or qi<2:
     R=u.extend(P,a,c);need(coeff_count(R,[1]*m)==v0,'independent cut-state base');independent_checks+=1
     for i in range(m):
      for k,v in [(2,b2[i]),(3,b3[i])]:ks=[1]*m;ks[i]=k;need(coeff_count(R,ks)==v,'independent cut-state extraction');independent_checks+=1
   lookup[(tuple(Q),a,c)]=cs;records.append({'block_query':[a,c],'coefficients':cs})
  r['extracted_coefficients']=records
 for r in canonical:
  vv=vals[r.pop('id')];expected=r.pop('expected_H_values');need([[v[0],v[2]] for v in vv]==expected,'fresh complete canonical Hamilton controls');unit=pow(3,r['types'].count('C'),8);pred=[]
  for (a,c),v in zip(product(range(8),repeat=2),vv):
   w=[unit*evaluate(z,r['types'])%8 for z in lookup[(tuple([2,4,1]),a,c)]];need(w==[v[1]%8,v[3]%8],'complete16000 canonical two-cover values');pred.append(w)
  r['values_H_q2_source_target']=vv;r['predicted_q2_source_target_mod8']=pred;r['odd_H_product_mod8']=unit
 for r in variants:
  vv=vals[r.pop('id')];pred=[]
  for (a,c),v in zip(r['block_queries'],vv):
   w=[r['odd_H_product_mod8']*evaluate(z,r['types'])%8 for z in lookup[(tuple(r['quotient_rows']),a,c)]];need(w==[v[1]%8,v[3]%8],'larger actual child physical restoration');pred.append(w)
  r['values_H_q2_source_target']=vv;r['predicted_q2_mod8']=pred
 # Explicit global-path exchange orbits, including repeated labels at a cut.
 instances=[];total=0;witnesses={}
 for Q in [[0],[2,0],[2,4,1],u.tt(3)]:
  m=len(Q)
  for a,c in product(range(1<<m),repeat=2):
   R=u.extend(Q,a,c)
   for ks in product(range(1,4),repeat=m):
    if sum(k>1 for k in ks)>2:continue
    words=word_covers(R,ks);W=set(words);need(len(W)==len(words) and len(W)%2==0,'even distinct ordered word pairs')
    for w in words:need((w[1],w[0]) in W and w[0]!=w[1],'global-path exchange has no fixed point')
    need(coeff_count(R,ks)==len(W),'explicit words match cut-state DP');total+=len(W)
    tag='same_label_at_cut' if any(w[0][-1]==w[1][0] for w in words) else 'ordinary'
    if words and tag not in witnesses:witnesses[tag]={'quotient_rows':Q,'block_query':[a,c],'visits':ks,'words':words}
    instances.append({'quotient_rows':Q,'block_query':[a,c],'visits':ks,'count':len(W),'sha256':hashlib.sha256(json.dumps(words,separators=(',',':')).encode()).hexdigest()})
 need('same_label_at_cut' in witnesses,'actual same-child boundary control')
 # Three global paths have non-free permutation action: exact coefficient9.
 models=[u.profile(u.extend(u.tt(n),0,0)) for n in (1,2,3)];g=[q[3] for q in models];num=g[2]-3*g[1]+2*g[0];need(num==54,'odd three-cover coefficient numerator');R=u.extend([0],0,0);three=word_covers(R,[3],3);need(len(three)==9 and ((1,2,0),(0,),(0,)) in three,'literal fixed three-cover pattern');need(coeff_count(R,[3],3)==9,'three-cover cut-state boundary')
 deps=['five_profile_types/PROOF.md','five_profile_types/normal.json','coarse_mod8_affine/PROOF.md','coarse_mod8_affine/normal.json','uniform_block_starts/verify.py']
 return {'status':'PASS','source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','count.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in deps},'physical_graph_occurrences':len(values)//2,'physical_H_counts':len(values)//2,'physical_q2_counts':len(values)//2,'independent_cut_state_coefficient_checks':independent_checks,'small_model_controls':small,'canonical_controls':canonical,'larger_child_controls':variants,'path_swap_instances':instances,'total_swapped_word_pairs':total,'positive_word_witnesses':witnesses,'three_cover_boundary':{'actual_model_profiles':models,'actual_three_cover_values':g,'extracted_coefficient':num//6,'all_nine_abstract_covers':three,'fixed_under_last_two_path_swap':[[1,2,0],[0],[0]]},'scope':'Every full two-cover response modulo8 is affine only in singleton and odd-cycle type bits after restoring the actual odd Hamilton scalar. Arbitrary quotient and nonempty child sizes; no non-annihilation or arbitrary-ambient separation is asserted.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
