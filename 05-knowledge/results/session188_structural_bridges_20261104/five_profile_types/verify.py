"""Complete canonical residue certificate and actual named response controls."""
from pathlib import Path
from itertools import product
from functools import lru_cache
import hashlib,importlib.util,json,random,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('u',ROOT.parent/'uniform_block_starts'/'verify.py');u=importlib.util.module_from_spec(s);s.loader.exec_module(u)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def subsets(S):
 R=S
 while True:
  yield R
  if not R:break
  R=(R-1)&S
@lru_cache(None)
def profile(T):return tuple(u.profile(list(T)))
def cycles(T):
 n=len(T)
 if n==1:return 0
 d=[[0]*n for _ in range(1<<n)];d[1][0]=1
 for S in range(1,1<<n,2):
  for v in range(n):
   if not d[S][v]:continue
   R=T[v]&~S
   while R:
    b=R&-R;R-=b;d[S|b][b.bit_length()-1]+=d[S][v]
 return sum(d[-1][v] for v in range(1,n) if T[v]&1)
def typeid(T,z):
 n=len(T)
 return 0 if n==1 else 1 if n==2 else 2 if n%2==0 else 4 if z%2 else 3
def phys(mask,blocks,offs):return sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,offs)) if mask>>i&1)
def lifted(T,U,queries,Q):
 n=len(T);A,_=u.blow(Q,[T]+[[0] for _ in Q[1:]]);B,_=u.blow(Q,[U]+[[0] for _ in Q[1:]]);outside=sum(1<<(n+j-1) for j in range(1,len(Q)) if Q[j]&1)
 return A,B,[(a|outside,c|outside) for a,c in queries]

def main():
 can=[u.tt(1),u.tt(2),u.tt(4),u.tt(3),[2,4,1]];names=['S','P','E','O','C'];cp=[profile(tuple(T)) for T in can];normalized=[[1,0,0],[1,2,0],[1,6,4],[1,6,6],[1,2,2]];rng=random.Random(188805);profile_records=[]
 tournaments=[u.rows(n,code) for n in range(1,6) for code in range(1<<(n*(n-1)//2))]
 for n in range(6,10):
  tournaments += [u.rows(n,rng.randrange(1<<(n*(n-1)//2))) for _ in range(12)]
 for T in tournaments:
  q=profile(tuple(T));z=cycles(T) if len(T)%2 else 0;t=typeid(T,z);H=q[1];unit=H*pow(cp[t][1],-1,8)%8;need(H%2 and all((q[k]-unit*(cp[t][k] if k<len(cp[t]) else 0))%8==0 for k in range(1,len(q))),'complete actual profile canonical reduction');need([(q[k] if k<len(q) else 0)*pow(H,-1,8)%8 for k in range(1,4)]==normalized[t],'correct normalized exceptional types');profile_records.append({'rows':T,'q':q,'odd_Hamilton_cycles':z,'type':names[t],'unit':unit})
 Q=[2,4,1];O=u.op(Q);blockqueries=list(product(range(8),repeat=2));cases=[];canonical=[]
 for ts in product(range(5),repeat=3):
  blocks=[can[t] for t in ts];T,offs=u.blow(Q,blocks);U,_=u.blow(O,blocks);queries=[(phys(a,blocks,offs),phys(c,blocks,offs)) for a,c in blockqueries];idx=len(cases);cases.append((T,U,queries));canonical.append({'types':[names[t] for t in ts],'child_rows':blocks,'source_old_rows':T,'target_old_rows':U,'physical_queries':queries,'case':idx})
 # Complete named physical bank for a usual-probe-blind but fully visible type.
 blocks=[u.tt(3)]*3;T,offs=u.blow(Q,blocks);U,_=u.blow(O,blocks);n=len(T);half=[a for a in range(1<<n) if a.bit_count()<=2];bank=list(product(half,repeat=2));extra=[(rng.randrange(1<<n),rng.randrange(1<<n)) for _ in range(32)];bankcase=len(cases);cases.append((T,U,bank+extra));bankT,bankU=T,U
 # Arbitrary actual profiles with scalar H residues and both odd-cycle types.
 odd5=next(r['rows'] for r in profile_records if len(r['rows'])==5 and r['odd_Hamilton_cycles']%2 and r['q'][1]%8!=3)
 variants=[]
 for blocks in [[u.tt(5),u.tt(2),[2,4,1]],[[8,1,3,6],u.tt(3),u.tt(2)],[odd5,u.tt(2),u.tt(3)]]:
  T,offs=u.blow(Q,blocks);U,_=u.blow(O,blocks);qs=[profile(tuple(B)) for B in blocks];zs=[cycles(B) if len(B)%2 else 0 for B in blocks];ts=[typeid(B,z) for B,z in zip(blocks,zs)];unit=1
  for q,t in zip(qs,ts):unit=unit*q[1]*pow(cp[t][1],-1,8)%8
  queries=[(phys(a,blocks,offs),phys(c,blocks,offs)) for a,c in blockqueries];idx=len(cases);cases.append((T,U,queries));variants.append({'child_rows':blocks,'profiles':qs,'odd_cycles':zs,'types':[names[t] for t in ts],'unit':unit,'case':idx})
 # Even-total core: use an exact known canonical query selected below after counts.
 ambient=[]
 for a,c in [(1,3),(1,2),(3,1)]:
  blocks=[u.tt(5),u.tt(2),[2,4,1]];T,offs=u.blow(Q,blocks);U,_=u.blow(O,blocks);query=(phys(a,blocks,offs),phys(c,blocks,offs))
  for outer in [[0],[2,4,1],[6,12,24,17,3]]:
   A,B,queries=lifted(T,U,[query],outer);idx=len(cases);cases.append((A,B,queries));ambient.append({'block_query':[a,c],'outer_rows':outer,'source_old_rows':A,'target_old_rows':B,'physical_query':queries[0],'case':idx})
 # The explicit odd-total modular loss remains a boundary, not a collision.
 hostile=[];blocks=[u.tt(2),u.tt(2),u.tt(3)];T,offs=u.blow(Q,blocks);U,_=u.blow(O,blocks);query=(3,12)
 for outer in [[0],[2,4,1]]:
  A,B,queries=lifted(T,U,[query],outer);idx=len(cases);cases.append((A,B,queries));hostile.append({'outer_rows':outer,'physical_query':queries[0],'case':idx})
 lines=[str(len(cases))]
 for A,B,queries in cases:lines += [f'{len(A)} {len(queries)}',' '.join(map(str,A)),' '.join(map(str,B))]+[f'{a} {c}' for a,c in queries]
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'literal.cpp'),'-o',str(exe)],check=True);nums=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 vals=[];at=0
 for A,B,queries in cases:
  count=len(queries);vals.append([nums[at+2*j:at+2*j+2] for j in range(count)]);at+=2*count
 need(at==len(nums),'all physical counts consumed');blind=[];hist={};lookup={};nontrivial=0
 for r in canonical:
  vv=vals[r.pop('case')];contrasts=[(a-b)%8 for a,b in vv];visible=sum(bool(z) for z in contrasts);hist[visible]=hist.get(visible,0)+1;ts=tuple(r['types']);lookup[ts]=vv;r['values']=vv;r['contrasts_mod8']=contrasts;r['nonzero_query_count']=visible
  if not visible:blind.append(ts)
  else:r['first_separator_block_query']=list(blockqueries[next(i for i,z in enumerate(contrasts) if z)])
  if 'S' not in ts:
   need(set(contrasts)<=set([0,4]),'nontrivial-child residue-four law');nontrivial+=bool(visible)
 need(blind==[('E','E','E')] and nontrivial==63,'complete125-type classification and63 nontrivial separators')
 allO=next(r for r in canonical if r['types']==['O']*3);need(allO['nonzero_query_count']==12 and allO['contrasts_mod8'][blockqueries.index((1,2))]==0 and allO['contrasts_mod8'][blockqueries.index((1,3))]==4,'usual-probe blindness is not full coarse blindness')
 for r in variants:
  vv=vals[r.pop('case')];cc=lookup[tuple(r['types'])];need(all((x-r['unit']*y)%8==0 for a,b in zip(vv,cc) for x,y in zip(a,b)),'literal arbitrary-profile full-bank transfer');r['physical_values']=vv
 for r in ambient:
  r['values']=vals[r.pop('case')][0];r['contrast_mod8']=(r['values'][0]-r['values'][1])%8
 for i in range(0,len(ambient),3):
  local=ambient[i]['contrast_mod8'];need(all((r['contrast_mod8']!=0)==(local!=0) for r in ambient[i:i+3]),'even-total ambient nonzero transfer')
 for r in hostile:r['values']=vals[r.pop('case')][0];r['integer_contrast']=r['values'][0]-r['values'][1]
 need([r['integer_contrast'] for r in hostile]==[-68,-600],'odd-total explicit modular-loss boundary')
 values={q:vals[bankcase][i] for i,q in enumerate(bank)};coeff={};bankrecords=[];degreehist={k:0 for k in range(5)}
 for a,c in bank:
  z=[0,0]
  for x in subsets(a):
   for y in subsets(c):
    sign=(-1)**(a.bit_count()+c.bit_count()-x.bit_count()-y.bit_count())
    for j,v in enumerate(values[x,y]):z[j]+=sign*v
  coeff[a,c]=z;bankrecords.append([a,c,*values[a,c],*z]);degreehist[a.bit_count()+c.bit_count()]+=z[0]!=z[1]
 full=[]
 for i,(a,c) in enumerate(extra):
  z=[sum(v[j] for (x,y),v in coeff.items() if x&a==x and y&c==y) for j in range(2)];need(z==vals[bankcase][len(bank)+i],'complete named polynomial reconstruction');full.append([a,c,*z])
 need(len(bankrecords)==2116 and degreehist[3]>0,'complete named degree-four bank retained')
 deps=['mod8_transfer/PROOF.md','uniform_block_starts/verify.py']
 return {'status':'PASS','canonical_type_order':names,'canonical_child_rows':can,'canonical_profiles':cp,'normalized_profiles_q1_q2_q3':normalized,'actual_profile_controls':profile_records,'canonical_type_triples':125,'complete_coarse_queries_per_triple':64,'canonical_records':canonical,'only_blind_triples':blind,'nonzero_count_histogram':hist,'literal_counts':len(nums),'actual_profile_transfer_controls':variants,'even_total_ambient_controls':ambient,'odd_total_boundary':hostile,'complete_named_TT3_bank':{'source_old_rows':bankT,'target_old_rows':bankU,'columns':['a','c','F_source','F_target','coefficient_source','coefficient_target'],'records':bankrecords,'nonzero_difference_histogram':degreehist,'full_mask_reconstructions':full},'source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','literal.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in deps},'scope':'Complete five-type modulo8 reduction for block-uniform responses; exact125-type cyclic classification; arbitrary ambient lift only for even total core order. Not an exact coarse or full named collision.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
