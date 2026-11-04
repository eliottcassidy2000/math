"""Complete singleton residue bridge and unified cyclic ambient controls."""
from pathlib import Path
from itertools import product
from functools import lru_cache
import hashlib,importlib.util,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('u',ROOT.parent/'uniform_block_starts'/'verify.py');u=importlib.util.module_from_spec(s);s.loader.exec_module(u)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
@lru_cache(None)
def profile(T):return tuple(u.profile(list(T)))
def modules(T,proper=False):
 top=(1<<len(T))-1
 return[M for M in range(1,top+1) if not(proper and M==top) and all(T[v]&M in(0,M) for v in range(len(T)) if not M>>v&1)]
def flip(T,family):
 R=list(T)
 for M in family:
  for v in range(len(T)):
   if M>>v&1:R[v]^=M^(1<<v)
 return R
def build(blocks,outer,query,targetblocks=None):
 T,offs=u.blow([2,4,1],blocks);U,_=u.blow([4,1,2],blocks if targetblocks is None else targetblocks);n=len(T);a,c=[sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,offs)) if m>>i&1)for m in query];outside=sum(1<<(n+j-1)for j in range(1,len(outer))if outer[j]&1);T,_=u.blow(outer,[T]+[[0]for _ in outer[1:]]);U,_=u.blow(outer,[U]+[[0]for _ in outer[1:]]);a|=outside;c|=outside
 for V in(T,U):
  G=u.extend(V,a,c);M=((1<<n)-1)|(1<<len(V))|(1<<(len(V)+1));need(all(G[v]&M in(0,M)for v in range(n,len(V))),'actual common augmented-core lift')
 return T,U,a,c

def main():
 data=json.loads((ROOT.parent/'five_profile_types'/'normal.json').read_text());rows=[r for r in data['canonical_records']if'S'in r['types']];need(len(rows)==61,'all singleton-containing types');queries=list(product(range(8),repeat=2));jobs=[];single=[];chosen={}
 def add(T,U,a,c,K):
  i=len(jobs);jobs.extend([(u.extend(T,a,c),K),(u.extend(U,a,c),K)]);return i
 for r in rows:
  residues=[v%4 for v in r['contrasts_mod8']];j=next(i for i,v in enumerate(residues)if v);need(residues[j]==2,'every singleton type has residue2');query=queries[j];chosen[tuple(r['types'])]=query;a,c=r['physical_queries'][j];idx=add(r['source_old_rows'],r['target_old_rows'],a,c,3);single.append({'types':r['types'],'complete_mod4_bank':residues,'block_query':query,'physical_query':[a,c],'source_old_rows':r['source_old_rows'],'target_old_rows':r['target_old_rows'],'frozen_H_values':r['values'][j],'jobs':idx})
 larger=[]
 variants=[(['S','O','E'],[[0],u.tt(5),[8,1,3,6]]),(['C','S','E'],[[2,4,1],[0],u.tt(6)]),(['E','C','S'],[u.tt(4),[2,4,1],[0]]),(['S','S','E'],[[0],[0],u.tt(8)])]
 for types,blocks in variants:
  query=chosen[tuple(types)]
  for outer in[[0],[2,4,1],[6,12,24,17,3]]:
   T,U,a,c=build(blocks,outer,query);idx=add(T,U,a,c,1);larger.append({'types':types,'child_rows':blocks,'outer_rows':outer,'block_query':query,'physical_query':[a,c],'source_old_rows':T,'target_old_rows':U,'jobs':idx})
 # Every original child-module family, including crossing TT3 intervals.
 blocks=[[0],u.tt(3),u.tt(3)];mods=[[m for m in modules(B)if m.bit_count()>1]for B in blocks];sourceprofiles=[profile(tuple(B))for B in blocks];z=chosen['S','O','O'];dual=(7^z[1],7^z[0]);families=[]
 for outer in[[0],[2,4,1]]:
  baseline=[]
  for query in[z,dual]:T,U,a,c=build(blocks,outer,query);baseline.append(add(T,U,a,c,1))
  for masks in product(*[range(1<<len(mm))for mm in mods]):
   selected=[[m for j,m in enumerate(mm)if bits>>j&1]for mm,bits in zip(mods,masks)];target=[u.op(flip(B,ff))for B,ff in zip(blocks,selected)];tp=[profile(tuple(B))for B in target];need(all(b>=a for x,y in zip(sourceprofiles,tp)for a,b in zip(x,y)),'hereditary singleton-family profile monotonicity');ids=[]
   for query in[z,dual]:T,U,a,c=build(blocks,outer,query,target);ids.append(add(T,U,a,c,1))
   families.append({'outer_rows':outer,'selected_original_child_modules':selected,'target_child_rows':target,'target_profiles':tp,'profile_strict':tp!=sourceprofiles,'baseline_jobs':baseline,'target_jobs':ids})
 # Two highest changed cyclic quotients, with exact exterior restoration.
 blocks=[u.tt(2),[2,4,1],[0]];M,offs=u.blow([2,4,1],blocks);V,_=u.blow([4,1,2],blocks);T,_=u.blow([2,4,1],[M,[2,4,1],[0]]);B,_=u.blow([2,4,1],[V,[2,4,1],[0]]);U=flip(T,[63,448]);query=chosen['P','C','S'];a,c=[sum(((1<<len(C))-1)<<o for i,(C,o)in enumerate(zip(blocks,offs))if bits>>i&1)|(1<<9)for bits in query];hi=len(jobs);jobs.extend([(u.extend(A,a,c),1)for A in(T,B,U)]);allmods=modules(T);strong=[m for m in allmods if not any((m&x)and(m&~x)and(x&~m)for x in allmods)];nodes=[];changed=[]
 for m in strong:
  children=[x for x in strong if x!=m and x&m==x and not any(x!=y and y!=m and x&y==x and y&m==y for y in strong)]
  if len(children)<2:continue
  reps=[(x&-x).bit_length()-1 for x in children];flag=any(bool(T[i]>>j&1)!=bool(U[i]>>j&1)for i in reps for j in reps if i!=j)
  if flag:changed.append(m)
  nodes.append({'mask':m,'children':children,'changed_quotient':flag})
 highest=[m for m in changed if not any(m!=x and m&x==m for x in changed)];need(sorted(highest)==[63,448],'two genuine highest original strong quotients');need(u.hp(T)==u.hp(U),'multiple-highest scalar equality')
 # There is no prime strongly connected quotient of order four.
 primes4=[c for c in range(64)if not any(1<m.bit_count()<4 for m in modules(u.rows(4,c),True))];need(not primes4,'no four-name prime quotient')
 payload=[str(len(jobs))]+[f'{len(T)} {K} '+' '.join(map(str,T))for T,K in jobs]
 with tempfile.TemporaryDirectory()as td:
  exe=Path(td)/'physical';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'physical.cpp'),'-o',str(exe)],check=True);values=[list(map(int,line.split()))for line in subprocess.check_output([str(exe)],input=('\n'.join(payload)+'\n').encode()).decode().splitlines()]
 need(len(values)==len(jobs),'all literal physical counts')
 for r in single:
  i=r.pop('jobs');x,y=values[i:i+2];need([x[0],y[0]]==r['frozen_H_values'],'canonical Hamilton replay');need([(a-b)%4 for a,b in zip(x,y)]==[2,0,0],'first-three same-order cover residues');r['source_first_three']=x;r['target_first_three']=y
 for r in larger:
  i=r.pop('jobs');r['values']=[values[i][0],values[i+1][0]];need((r['values'][0]-r['values'][1])%4==2,'actual larger singleton-child arbitrary ambient')
 for r in families:
  base=[[values[i][0],values[i+1][0]]for i in r.pop('baseline_jobs')];vv=[[values[i][0],values[i+1][0]]for i in r.pop('target_jobs')];diff=[a-b for a,b in base];need(diff[0]==-diff[1] and diff[0]%4==2,'singleton-family exact baseline duality');k=0 if diff[0]<0 else 1;need(vv[k][0]-vv[k][1]<=diff[k]<0,'one fixed member handles every crossing family');r['baseline_values']=base;r['target_values']=vv;r['negative_query_index']=k
 hv=[values[hi+j][0]for j in range(3)];need(hv[1]==hv[2]and(hv[0]-hv[1])%4==2,'highest-node restoration and physical separator')
 deps=['five_profile_types/PROOF.md','five_profile_types/normal.json','parity_profile_transfer/PROOF.md','cyclic_odd_profile_lift/PROOF.md','cyclic_odd_profile_lift/normal.json','even_cyclic_mod16/PROOF.md','even_cyclic_mod16_ambient/PROOF.md','even_cyclic_mod16_ambient/normal.json','monotone_dual_queries/PROOF.md','module_decomposition/PROOF.md','highest_quotient_reduction/PROOF.md','uniform_block_starts/verify.py']
 return{'status':'PASS','singleton_canonical_triples':61,'singleton_records':single,'larger_child_ambient_controls':larger,'root_family_contexts':len(families),'profile_strict_family_contexts':sum(r['profile_strict']for r in families),'root_family_block_queries':[z,dual],'root_family_records':families,'multiple_highest_control':{'source_rows':T,'baseline_rows':B,'target_rows':U,'original_strong_nodes':nodes,'highest_changed_masks':highest,'physical_query':[a,c],'values':hv},'prime_four_codes':primes4,'physical_graph_counts':len(jobs),'physical_profile_coordinates':sum(K for T,K in jobs),'source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in['PROOF.md','verify.py','physical.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest()for f in deps},'scope':'Every nonempty cyclic child triple has a block-uniform quotient-converse separator in every fixed ambient; negative-dual monotonicity detects every original-module root family. Any remaining full-response collision has only highest prime quotients of order>=5.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
