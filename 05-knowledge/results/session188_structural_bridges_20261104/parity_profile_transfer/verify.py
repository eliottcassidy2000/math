"""Universal endpoint-parity transfer, with complete named cubic/quartic control."""
from pathlib import Path
from itertools import product,combinations
from functools import lru_cache
import importlib.util,json,hashlib,random,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('starts',ROOT.parent/'uniform_block_starts'/'verify.py');u=importlib.util.module_from_spec(s);s.loader.exec_module(u)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def pbit(Q):
 S=u.paths(Q);E=u.paths(u.op(Q));n=len(Q)
 need(all(sum(row[i] for row in S)%2==E[-1][i]%2 and sum(row[i] for row in E)%2==S[-1][i]%2 for i in range(n)),'all-path and Hamilton endpoint parity bridge')
 return sum(((a+b)%2)<<i for i,(a,b) in enumerate(zip(S[-1],E[-1])))
def prime(Q):
 n=len(Q)
 return all(any(Q[v]&M not in(0,M) for v in range(n) if not(M>>v&1)) for M in range(1,(1<<n)-1) if M.bit_count()>1)
def strong(Q):
 n=len(Q);R=[x|(1<<i) for i,x in enumerate(Q)]
 for k in range(n):
  for i in range(n):
   if R[i]>>k&1:R[i]|=R[k]
 return all(x==(1<<n)-1 for x in R)
def phys(mask,blocks,offs):return sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,offs)) if mask>>i&1)
def repeats(R,z):
 n=len(R);full=((1<<n)-1)^(1<<z)
 @lru_cache(None)
 def f(S,k,last):
  out=[int(S==full and k==2),int(S==full and k==3)];rem=R[last]
  while rem:
   bit=rem&-rem;rem-=bit;v=bit.bit_length()-1
   if v==z:
    if k==3:continue
    t=f(S,k+1,v)
   else:
    if S&bit:continue
    t=f(S|bit,k,v)
   out[0]^=t[0];out[1]^=t[1]
  return tuple(out)
 out=[0,0]
 for v in range(n):
  r=f(0,1,v) if v==z else f(1<<v,0,v);out[0]^=r[0];out[1]^=r[1]
 return out

def main():
 profiles=[];quotients=[];counts={};rng=random.Random(188163)
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   T=u.rows(n,code);q=u.profile(T);need(q[1]%2==1,'Redei profile parity')
   need(all(q[k]%4==(2 if k==2 and n>=2 else 2 if k==3 and n%2 else 0) for k in range(2,n+1)),'universal higher profile residues')
   profiles.append([n,code,q]);p=pbit(T);pr=prime(T) and strong(T);quotients.append([n,code,p,pr]);counts[n]=counts.get(n,0)+bool(p)
 cases=[];labels=[];involution=0
 def add_case(name,Q,blocks,queries,law=True,ambient=False):
  nonlocal involution
  T,o=u.blow(Q,blocks);U,_=u.blow(u.op(Q),blocks);P=pbit(Q);physical=[(phys(a,blocks,o),phys(c,blocks,o)) for a,c in queries]
  if ambient:
   # Put the core at vertex0 of an outer C3, with two singleton exterior names.
   n=len(T);T,_=u.blow([2,4,1],[T,[0],[0]]);U,_=u.blow([2,4,1],[U,[0],[0]]);E=1<<(n+1);physical=[(a|E,c|E) for a,c in physical]
  cases.append((T,U,physical));labels.append({'name':name,'quotient':Q,'child_rows':blocks,'block_queries':queries,'physical_queries':physical,'endpoint_parity_mask':P,'linear_law':law,'ambient':ambient})
  if not ambient and len(Q)<=4:
   for a,c in queries[:4]:
    for A in(Q,u.op(Q)):
     R=u.extend(A,a,c)
     for z in range(len(Q)):need(repeats(R,z)[1]==0,'three-return skeleton involution parity');involution+=1
 for n in range(1,4):
  for code in range(1<<(n*(n-1)//2)):
   Q=u.rows(n,code);add_case('all_small_pairs',Q,[u.tt(2)]*n,list(product(range(1<<n),repeat=2)))
 for code in range(64):
  Q=u.rows(4,code);queries=[(rng.randrange(16),rng.randrange(16)) for _ in range(4)];add_case('four_quotient_pairs',Q,[u.tt(2)]*4,queries)
 for blocks in([u.tt(2),[2,4,1],[8,1,3,6]],[u.tt(3),[2,4,1],u.tt(2)],[[2,4,1]]*3):add_case('mixed_nontrivial_C3',[2,4,1],blocks,list(product(range(8),repeat=2)))
 Q=[18,8,3,5,14];need(prime(Q) and strong(Q) and pbit(Q)==9,'prime-five endpoint witness')
 blocks=[u.tt(2)]*5;T,o=u.blow(Q,blocks);U,_=u.blow(u.op(Q),blocks);n=len(T);half=[m for m in range(1<<n) if m.bit_count()<=2];bank_queries=list(product(half,repeat=2));coarse=list(product(range(32),repeat=2));actualcoarse=[(phys(a,blocks,o),phys(c,blocks,o)) for a,c in coarse];bank_index=len(cases);cases.append((T,U,bank_queries+actualcoarse));labels.append({'name':'complete_named_ten_bank'})
 for blocks in([u.tt(2),[2,4,1],u.tt(2),u.tt(2),u.tt(2)],[[0],u.tt(2),[2,4,1],[0],u.tt(2)]):
  allnontrivial=all(len(b)>1 for b in blocks);queries=[(1,0),(8,0),(9,0),(0,0)]
  add_case('arbitrary_prime_children',Q,blocks,queries,allnontrivial)
  add_case('ambient_prime_children',Q,blocks,queries,allnontrivial,True)
 # Disjoint highest prime changes: endpoint-visible five-quotient core plus another C3.
 M,_=u.blow(Q,[u.tt(2)]*5);Mr,_=u.blow(u.op(Q),[u.tt(2)]*5);Top,_=u.blow([2,4,1],[M,[2,4,1],[0]]);Mid,_=u.blow([2,4,1],[Mr,[2,4,1],[0]]);End,_=u.blow([2,4,1],[Mr,[4,1,2],[0]]);E=1<<13;highquery=(3|E,E);high_index=len(cases);cases.append((Top,Mid,[highquery]));labels.append({'name':'highest_core_only'});cases.append((Mid,End,[highquery]));labels.append({'name':'other_highest_restoration'})
 lines=[str(len(cases))]
 for T,U,qs in cases:lines.extend([f'{len(T)} {len(qs)}',' '.join(map(str,T)),' '.join(map(str,U))]+[f'{a} {c}' for a,c in qs])
 with tempfile.TemporaryDirectory() as d:
  exe=Path(d)/'literal';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'literal.cpp'),'-o',str(exe)],check=True);raw=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 need(len(raw)==2*sum(len(qs) for T,U,qs in cases),'literal output count');cursor=0;case_values=[]
 for T,U,qs in cases:
  size=len(qs);case_values.append([raw[cursor+2*i:cursor+2*i+2] for i in range(size)]);cursor+=2*size
 controls=[]
 for j,(meta,vals) in enumerate(zip(labels,case_values)):
  if 'block_queries' not in meta:continue
  p=meta['endpoint_parity_mask']
  for (a,c),(f,g) in zip(meta['block_queries'],vals):
   if meta['linear_law']:need((f-g)%4==2*(((a^c)&p).bit_count()%2),'complete coarse linear mod4 law')
   elif c==0 and a.bit_count()==1:need((f-g)%4==2*bool(a&p),'singleton-allowed endpoint query')
  controls.append({**meta,'query_values':vals})
 values=case_values[bank_index];lookup={key:pair for key,pair in zip(bank_queries,values)};records=[];hist={}
 for a,c in bank_queries:
  co=[0,0];aa=a
  while True:
   cc=c
   while True:
    sign=(-1)**((a^aa).bit_count()+(c^cc).bit_count());pair=lookup[aa,cc]
    for k in range(2):co[k]+=sign*pair[k]
    if cc==0:break
    cc=(cc-1)&c
   if aa==0:break
   aa=(aa-1)&a
  delta=co[0]-co[1];degree=a.bit_count()+c.bit_count()
  if delta:hist[degree]=hist.get(degree,0)+1
  records.append([a,c,*lookup[a,c],*co])
 need(hist.get(4,0)>0 and hist.get(3,0)>0,'full named cubic and quartic retained')
 for (a,c),(aa,cc),(f,g) in zip(coarse,actualcoarse,values[len(bank_queries):]):
  got=[sum(r[4+k] for r in records if r[0]&aa==r[0] and r[1]&cc==r[1]) for k in range(2)]
  need(got==[f,g],'full named interpolation at every coarse query');need((f-g)%4==2*(((a^c)&9).bit_count()%2),'prime-five complete coarse mod4 law')
 high=case_values[high_index:high_index+2];need((high[0][0][0]-high[0][0][1])%4==2 and high[1][0][0]==high[1][0][1],'multiple highest prime restoration')
 return {'status':'PASS','scope':'Universal modulo4 transfer and one-query endpoint criterion; balanced endpoint fingerprints do not imply exact or full named collisions.','profile_controls':profiles,'quotient_endpoint_records':quotients,'quotients_with_nonzero_endpoint_mask_by_order':counts,'three_return_involution_checks':involution,'literal_count_total':len(raw),'literal_controls':controls,'full_named_ten':{'quotient':Q,'T':cases[bank_index][0],'U':cases[bank_index][1],'basis_size':len(records),'basis_records_a_c_FT_FU_coefT_coefU':records,'omission_convention':'All Boolean coefficients with either port degree above2 are zero. Every remaining named coordinate, including zero and all cubic/quartic entries, is listed.','nonzero_contrast_counts_by_degree':hist,'all1024_coarse_values':values[len(bank_queries):]},'multiple_highest_control':{'T':Top,'middle':Mid,'U':End,'query':highquery,'values':high},'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT/'literal.cpp',ROOT.parent/'uniform_block_starts'/'verify.py']}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
