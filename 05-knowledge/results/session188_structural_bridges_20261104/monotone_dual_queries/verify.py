"""Exact dual-query monotonicity, with complete named cubic/quartic records."""
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
def modules(T):
 return [M for M in range(1,1<<len(T)) if M.bit_count()>1 and all(T[v]&M in(0,M) for v in range(len(T)) if not(M>>v&1))]
def reverse(T,selected):
 R=list(T)
 for M in selected:
  for i in range(len(T)):
   if M>>i&1:R[i]^=M^(1<<i)
 return R

def main():
 Q=[4,13,24,17,3];anchor=2;S=u.paths(Q)[-1];E=u.paths(u.op(Q))[-1];need((S[anchor]+E[anchor])%2==1,'specified endpoint anchor')
 allcases=[];records=[];baseline_data=[]
 specs=[[u.tt(3),[0],[0],[0],[0]],[u.tt(3),u.tt(2),u.tt(2),[0],[0]],[u.tt(3),u.tt(3),[0],[0],[0]]]
 for sid,blocks in enumerate(specs):
  T,offs=u.blow(Q,blocks);V,_=u.blow(u.op(Q),blocks);n=len(T);M=(1<<n)-1;I=((1<<len(blocks[anchor]))-1)<<offs[anchor];corequeries=[(I,0),(M,M^I)]
  childmods=[modules(B) for B in blocks];choices=[list(range(1<<len(mm))) for mm in childmods]
  # Check the complete augmented-profile duality on the smallest source.
  if sid==0:
   qs=[profile(tuple(u.extend(A,a,c))) for A in(T,V) for a,c in corequeries];need(qs[0]==qs[3] and qs[1]==qs[2],'all-coordinate converse-complement duality');baseline_data.append({'source_old_rows':T,'baseline_old_rows':V,'core_queries':corequeries,'complete_augmented_profiles':qs})
  for outer in([0],[2,4,1]):
   source,_=u.blow(outer,[T]+[[0] for _ in outer[1:]]);base,_=u.blow(outer,[V]+[[0] for _ in outer[1:]]);ext=sum(1<<(n+j-1) for j in range(1,len(outer)) if outer[j]&1);queries=[(a|ext,c|ext) for a,c in corequeries];bi=len(allcases);allcases.append((source,base,queries));families=[]
   for masks in product(*choices):
    targetblocks=[];strict=False;selected=[]
    for B,mm,bits in zip(blocks,childmods,masks):
     chosen=[m for j,m in enumerate(mm) if bits>>j&1];C=u.op(reverse(B,chosen));qa,qb=profile(tuple(B)),profile(tuple(C));need(all(b>=a for a,b in zip(qa,qb)),'hereditary child profile monotonicity');strict|=qb!=qa;targetblocks.append(C);selected.append(chosen)
    U,_=u.blow(u.op(Q),targetblocks);target,_=u.blow(outer,[U]+[[0] for _ in outer[1:]]);ci=len(allcases);allcases.append((source,target,queries));families.append({'selected_original_child_modules':selected,'target_child_rows':targetblocks,'profile_strict':strict,'case_index':ci})
   records.append({'source_spec':sid,'quotient':Q,'anchor':anchor,'child_rows':blocks,'outer_rows':outer,'source_old_rows':source,'baseline_old_rows':base,'physical_queries':queries,'baseline_case':bi,'families':families})
 # Complete seven-name bank, retaining scalar-strict crossing changes inside TT3.
 blocks=specs[0];T,offs=u.blow(Q,blocks);V,_=u.blow(u.op(Q),blocks);cross=u.op(reverse(blocks[0],[3,6]));need(profile(tuple(cross))[1]>profile(tuple(blocks[0]))[1],'crossing child creates strict Hamilton increase');W,_=u.blow(u.op(Q),[cross]+blocks[1:]);n=len(T);half=[a for a in range(1<<n) if a.bit_count()<=2];bank=list(product(half,repeat=2));rng=random.Random(188740);extra=[(rng.randrange(1<<n),rng.randrange(1<<n)) for _ in range(40)];b0=len(allcases);allcases.append((T,V,bank+extra));b1=len(allcases);allcases.append((T,W,bank+extra))
 # Independently check the concrete augmented-core lift on every family.
 for r in records:
  n=sum(map(len,r['child_rows']));oldn=len(r['source_old_rows']);coremask=((1<<n)-1)|(1<<oldn)|(1<<(oldn+1))
  for A in (r['source_old_rows'],r['baseline_old_rows']):
   for a,c in r['physical_queries']:
    R=u.extend(A,a,c);need(all(R[v]&coremask in(0,coremask) for v in range(n,oldn)),'common augmented module under physical lift')
 lines=[str(len(allcases))]
 for A,B,queries in allcases:lines += [f'{len(A)} {len(queries)}',' '.join(map(str,A)),' '.join(map(str,B))]+[f'{a} {c}' for a,c in queries]
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'literal.cpp'),'-o',str(exe)],check=True);nums=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 out=[];at=0
 for A,B,queries in allcases:
  count=len(queries);out.append([nums[at+2*i:at+2*i+2] for i in range(count)]);at+=2*count
 need(at==len(nums),'all literal outputs consumed');famcount=strictcount=0
 for r in records:
  baseline=out[r.pop('baseline_case')];d=[a-b for a,b in baseline];need(d[0]==-d[1] and d[0]%4==2,'opposite nonzero ambient baseline contrasts');chosen=0 if d[0]<0 else 1;r['baseline_values']=baseline;r['fixed_negative_query_index']=chosen
  for f in r['families']:
   vals=out[f.pop('case_index')];need(all(v[1]>=b[1] for v,b in zip(vals,baseline)),'literal target profile monotonicity');need(vals[chosen][0]-vals[chosen][1]<=d[chosen]<0,'one selected query handles every family');f['physical_values']=vals;famcount+=1;strictcount+=f['profile_strict']
 # Boolean coefficients retain every coordinate, including all zeros.
 values={key:[out[b0][i][0],out[b0][i][1],out[b1][i][1]] for i,key in enumerate(bank)};coeff={}
 for a,c in bank:
  z=[0,0,0]
  for x in subsets(a):
   for y in subsets(c):
    sign=(-1)**((a.bit_count()+c.bit_count())-(x.bit_count()+y.bit_count()))
    for j,v in enumerate(values[x,y]):z[j]+=sign*v
  coeff[a,c]=z
 histogram={d:0 for d in range(5)};bankrecords=[]
 for a,c in bank:
  z=coeff[a,c];d=a.bit_count()+c.bit_count();histogram[d]+=z[0]!=z[2];bankrecords.append([a,c,*values[a,c],*z])
 reconstruction=[]
 for i,(a,c) in enumerate(extra):
  z=[sum(v[j] for (x,y),v in coeff.items() if x&a==x and y&c==y) for j in range(3)];actual=[out[b0][len(bank)+i][0],out[b0][len(bank)+i][1],out[b1][len(bank)+i][1]];need(z==actual,'full cubic/quartic Boolean reconstruction');reconstruction.append([a,c,*actual])
 need(len(bank)==841 and histogram[3]>0 and histogram[4]>0,'complete named bank and genuine cubic/quartic changes')
 dependencies=['parity_profile_transfer/PROOF.md','module_decomposition/PROOF.md','uniform_block_starts/verify.py']
 return {'status':'PASS','scope':'Two fixed physical queries, one baseline-sign-selected member for all profile-increasing quotient-converse targets; fixed ambient and original-module root-family corollary. No universal reconstruction or exact collision.','quotient_start_counts':S,'quotient_end_counts':E,'family_controls':records,'families':famcount,'profile_strict_families':strictcount,'all_coordinate_duality':baseline_data,'complete_named_bank':{'source_old_rows':T,'baseline_old_rows':V,'strict_crossing_target_old_rows':W,'columns':['a','c','F_source','F_baseline','F_target','coefficient_source','coefficient_baseline','coefficient_target'],'records':bankrecords,'source_target_nonzero_histogram':histogram,'full_mask_reconstructions':reconstruction},'literal_counts':len(nums),'source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','literal.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in dependencies}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
