"""Complete even-child mod16 upper-domain certificate and literal controls."""
from pathlib import Path
from itertools import product
from functools import lru_cache
import hashlib,importlib.util,json,random,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
r=load('rank',ROOT.parent/'coarse_response_rank'/'verify.py');u=load('u',ROOT.parent/'uniform_block_starts'/'verify.py')
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
@lru_cache(None)
def profile(T):return tuple(u.profile(list(T)))
def subsets(S):
 R=S
 while True:
  yield R
  if not R:break
  R=(R-1)&S

def main():
 formal=[[0,1,b,c,8,d] for b,c,d in product([6,14],[4,12],[0,8])];queries=list(product(range(8),repeat=2));Gt=list(map(list,zip(*r.G)));numerators={q:r.add(r.numerator(r.G,*q),r.scale(r.numerator(Gt,*q),-1)) for q in queries};records=[];hist={};common=set(range(64))
 for ts in product(range(8),repeat=3):
  values=[r.functional(p,[formal[t] for t in ts]) for p in numerators.values()];res=[v%16 for v in values];nonzero=[i for i,v in enumerate(res) if v];need(nonzero and set(res)<=set([0,8]),'complete formal residue non-annihilation');need(res[queries.index((1,3))]==8,'one fixed universal physical query');hist[len(nonzero)]=hist.get(len(nonzero),0)+1;common &= set(nonzero);records.append({'profile_type_indices':ts,'integer_functionals':values,'residues_mod16':res,'first_separator':queries[nonzero[0]]})
 need(hist=={24:256,28:192,12:64},'exact full512-class histogram');need([queries[i] for i in sorted(common)]==[(1,3),(1,5),(2,3),(2,6),(4,5),(4,6)],'six universal queries')
 rng=random.Random(1881616);actual=[u.rows(4,c) for c in range(64)]
 for n,count in[(6,24),(8,24),(10,8)]:actual += [u.rows(n,rng.randrange(1<<(n*(n-1)//2))) for _ in range(count)]
 actualprofiles=[]
 for T in actual:
  q=profile(tuple(T));inv=pow(q[1],-1,16);norm=[0]+[(q[k] if k<len(q) else 0)*inv%16 for k in range(1,6)];need(norm in formal and all(v%16==0 for v in q[6:]),'complete actual normalized even profile');actualprofiles.append({'rows':T,'q':q,'normalized_first_five':norm,'type':formal.index(norm)})
 cases=[];controls=[]
 oddjoin=u.blow([2,0],[[2,4,1],[0]])[0];strong4=[8,1,3,6]
 variants=[[u.tt(4)]*3,[strong4,oddjoin,u.tt(4)],[oddjoin]*3,[u.tt(4),u.tt(4),u.tt(6)],[u.tt(4),oddjoin,u.tt(6)]]
 for j,blocks in enumerate(variants):
  T,o=u.blow([2,4,1],blocks);U,_=u.blow([4,1,2],blocks);qq=queries if j<3 else[(1,3),(1,2),(4,6)];physical=[(r.physical_mask(a,blocks,o),r.physical_mask(c,blocks,o)) for a,c in qq];idx=len(cases);cases.append((T,U,physical));qs=[profile(tuple(B)) for B in blocks];ts=[formal.index([0]+[(q[k] if k<len(q) else 0)*pow(q[1],-1,16)%16 for k in range(1,6)]) for q in qs];controls.append({'child_rows':blocks,'profiles':qs,'type_indices':ts,'block_queries':qq,'physical_queries':physical,'source_old_rows':T,'target_old_rows':U,'case':idx})
 # Complete named bank on the all-even residual's smallest physical core.
 blocks=variants[0];T,o=u.blow([2,4,1],blocks);U,_=u.blow([4,1,2],blocks);n=len(T);half=[m for m in range(1<<n) if m.bit_count()<=2];bank=list(product(half,repeat=2));extra=[(rng.randrange(1<<n),rng.randrange(1<<n)) for _ in range(32)];bankindex=len(cases);cases.append((T,U,bank+extra));bankT,bankU=T,U
 # Every original-module family inside the first TT4 child, root also selected.
 B=blocks[0];mods=[M for M in range(1,16) if M.bit_count()>1 and all(B[v]&M in(0,M) for v in range(4) if not M>>v&1)];need(len(mods)==6,'complete original TT4 interval modules');familyrecords=[];pair=[(1,3),(4,6)];pairphys=[(r.physical_mask(a,blocks,o),r.physical_mask(c,blocks,o)) for a,c in pair];baseline=len(cases);cases.append((T,U,pairphys))
 for bits in range(64):
  child=B[:];chosen=[M for i,M in enumerate(mods) if bits>>i&1]
  for M in chosen:
   for i in range(4):
    if M>>i&1:child[i]^=M^(1<<i)
  targetblocks=[u.op(child),u.op(B),u.op(B)];targetprofiles=[profile(tuple(C)) for C in targetblocks];need(all(b>=a for p in targetprofiles for a,b in zip(profile(tuple(B)),p)),'hereditary target profile monotonicity');V,_=u.blow([4,1,2],targetblocks);idx=len(cases);cases.append((T,V,pairphys));familyrecords.append({'selected_original_first_child_modules':chosen,'target_child_rows':targetblocks,'target_profiles':targetprofiles,'profile_strict':targetprofiles[0]!=profile(tuple(B)),'case':idx})
 lines=[str(len(cases))]
 for A,B,qq in cases:lines += [f'{len(A)} {len(qq)}',' '.join(map(str,A)),' '.join(map(str,B))]+[f'{a} {c}' for a,c in qq]
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'literal.cpp'),'-o',str(exe)],check=True);nums=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 values=[];at=0
 for A,B,qq in cases:values.append([nums[at+2*j:at+2*j+2] for j in range(len(qq))]);at+=2*len(qq)
 need(at==len(nums),'consume all actual counts')
 for z in controls:
  vv=values[z.pop('case')];expected=[r.functional(numerators[tuple(q)],z['profiles']) for q in z['block_queries']];need([a-b for a,b in vv]==expected,'full literal physical versus exact profile functional');need(all((a-b)%16==8 for q,(a,b) in zip(z['block_queries'],vv) if tuple(q) in[(1,3),(4,6)]),'fixed actual query residue');z['physical_values']=vv
 d=[a-b for a,b in values[baseline]];need(d[0]==-d[1] and d[0]%16==8,'opposite baseline exact signs');chosen=0 if d[0]<0 else 1
 for z in familyrecords:
  vv=values[z.pop('case')];need(vv[chosen][0]-vv[chosen][1]<=d[chosen]<0,'same dual member detects every root family');z['physical_values']=vv;z['negative_query_index']=chosen
 bv={q:values[bankindex][i] for i,q in enumerate(bank)};coeff={};bankrecords=[];degreehist={k:0 for k in range(5)}
 for a,c in bank:
  z=[0,0]
  for x in subsets(a):
   for y in subsets(c):
    sign=(-1)**(a.bit_count()+c.bit_count()-x.bit_count()-y.bit_count())
    for j,v in enumerate(bv[x,y]):z[j]+=sign*v
  coeff[a,c]=z;bankrecords.append([a,c,*bv[a,c],*z]);degreehist[a.bit_count()+c.bit_count()]+=z[0]!=z[1]
 full=[]
 for i,(a,c) in enumerate(extra):
  z=[sum(v[j] for(x,y),v in coeff.items() if x&a==x and y&c==y) for j in range(2)];need(z==values[bankindex][len(bank)+i],'complete named polynomial reconstruction');full.append([a,c,*z])
 need(len(bankrecords)==6241 and degreehist[3]>0 and degreehist[4]>0,'complete residual named cubic/quartic data')
 deps=['coarse_response_rank/PROOF.md','coarse_response_rank/verify.py','five_profile_types/PROOF.md','mod8_transfer/PROOF.md','monotone_dual_queries/PROOF.md','module_decomposition/PROOF.md','uniform_block_starts/verify.py']
 return {'status':'PASS','formal_profile_upper_domain':formal,'block_query_order':queries,'complete_numerator_polynomials':[[list(q),[[*e,c] for e,c in sorted(p.items())]] for q,p in numerators.items()],'formal_triples':512,'complete_formal_records':records,'nonzero_query_histogram':hist,'universal_queries':[queries[i] for i in sorted(common)],'actual_profile_controls':actualprofiles,'actual_physical_controls':controls,'root_family_records':familyrecords,'root_family_baseline_values':values[baseline],'root_family_block_queries':pair,'complete_named_TT4_bank':{'source_old_rows':bankT,'target_old_rows':bankU,'columns':['a','c','F_source','F_target','coefficient_source','coefficient_target'],'records':bankrecords,'nonzero_difference_histogram':degreehist,'full_mask_reconstructions':full},'literal_counts':len(nums),'source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in['PROOF.md','verify.py','literal.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in deps},'scope':'Universal fixed coarse modulo16 separator for arbitrary even children>=4, completing local C3 quotient-converse separation for all nonempty children. No arbitrary mod16 ambient lift in this packet.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
