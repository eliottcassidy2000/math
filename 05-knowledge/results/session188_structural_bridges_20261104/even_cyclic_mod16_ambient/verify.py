"""Exact special cyclic all-cover residue certificate and ambient controls."""
from pathlib import Path
from itertools import product
import hashlib,importlib.util,json,subprocess,tempfile
from cover_words import tensor
ROOT=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
u=load('u',ROOT.parent/'uniform_block_starts'/'verify.py')
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def build(blocks,outer,query=(1,3),targetblocks=None):
 T,offs=u.blow([2,4,1],blocks);U,_=u.blow([4,1,2],blocks if targetblocks is None else targetblocks);n=len(T);a,c=[sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,offs)) if mask>>i&1) for mask in query];outside=sum(1<<(n+j-1) for j in range(1,len(outer)) if outer[j]&1);T,_=u.blow(outer,[T]+[[0] for _ in outer[1:]]);U,_=u.blow(outer,[U]+[[0] for _ in outer[1:]]);a|=outside;c|=outside
 for V in(T,U):
  X=u.extend(V,a,c);M=((1<<n)-1)|(1<<len(V))|(1<<(len(V)+1));need(all(X[v]&M in(0,M) for v in range(n,len(V))),'same actual augmented core module')
 return T,U,a,c

def main():
 local=json.loads((ROOT.parent/'even_cyclic_mod16'/'normal.json').read_text());formal=local['formal_profile_upper_domain'];A=tensor(1,3);B=tensor(1,3,True);records=[];j=local['block_query_order'].index([1,3])
 def evaluate(K,ts):return[sum(v*formal[ts[0]][i]*formal[ts[1]][j]*formal[ts[2]][k] for(h,i,j,k),v in K.items() if h==hh) for hh in range(1,6)]
 for old,ts in zip(local['complete_formal_records'],product(range(8),repeat=3)):
  x,y=evaluate(A,ts),evaluate(B,ts);diff=[a-b for a,b in zip(x,y)];need(diff[0]==old['integer_functionals'][j],'exact Hamilton agreement with independent rational run formula');need([v%16 for v in diff]==[8,0,0,0,0],'every512 full first-five contrast');records.append({'profile_type_indices':ts,'source_first_five':x,'target_first_five':y,'integer_contrasts':diff})
 jobs=[];controls=[];core_records=[];oddjoin=u.blow([2,0],[[2,4,1],[0]])[0];variants=[[u.tt(4)]*3,[[8,1,3,6],oddjoin,u.tt(4)],[u.tt(4),u.tt(4),u.tt(6)],[oddjoin]*3]
 for blocks in variants:
  profiles=[u.profile(B) for B in blocks];normalized=[[0]+[(q[k] if k<len(q) else 0)*pow(q[1],-1,16)%16 for k in range(1,6)]for q in profiles];types=[formal.index(q) for q in normalized];unit=1
  for q in profiles:unit=unit*q[1]%16
  outers=[[0],[2,0],[2,4,1]]+([[6,12,24,17,3]] if sum(map(len,blocks))==12 else [])
  for outer in outers:
   T,U,a,c=build(blocks,outer);K=5 if len(outer)==1 else 1;idx=len(jobs);jobs.extend([(u.extend(T,a,c),K),(u.extend(U,a,c),K)]);controls.append({'child_rows':blocks,'profiles':profiles,'types':types,'unit':unit,'outer_rows':outer,'source_old_rows':T,'target_old_rows':U,'physical_query':[a,c],'coordinates_counted':K,'jobs':idx})
 # A crossing original-module target, with the same dual query pair in a strong exterior.
 blocks=variants[0];child=blocks[0][:]
 for M in[3,6]:
  for v in range(4):
   if M>>v&1:child[v]^=M^(1<<v)
 targetblocks=[u.op(child),u.op(blocks[1]),u.op(blocks[2])];qp=[u.profile(B)for B in targetblocks];need(all(b>=a for q in qp for a,b in zip(u.profile(blocks[0]),q)) and qp[0][1]>1,'strict crossing child profile growth');family=[]
 for query in[(1,3),(4,6)]:
  for target in[None,targetblocks]:
   T,U,a,c=build(blocks,[2,4,1],query,target);idx=len(jobs);jobs.extend([(u.extend(T,a,c),1),(u.extend(U,a,c),1)]);family.append({'block_query':query,'profile_increased':target is not None,'physical_query':[a,c],'source_old_rows':T,'target_old_rows':U,'jobs':idx})
 lines=[str(len(jobs))]+[f'{len(T)} {K} '+' '.join(map(str,T))for T,K in jobs]
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'physical';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'physical.cpp'),'-o',str(exe)],check=True);out=subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).decode().splitlines();values=[list(map(int,line.split()))for line in out]
 need(len(values)==len(jobs),'all physical output rows')
 for r in controls:
  i=r.pop('jobs');x,y=values[i:i+2];need(len(x)==len(y)==r['coordinates_counted'],'physical coordinate width');r['source_values']=x;r['target_values']=y;r['integer_contrasts']=[a-b for a,b in zip(x,y)];need(r['integer_contrasts'][0]%16==8,'every actual ambient gap8')
  if r['coordinates_counted']==5:
   need([v%16 for v in r['integer_contrasts']]==[8,0,0,0,0],'actual complete higher-coordinate vanishing');X,Y=evaluate(A,r['types']),evaluate(B,r['types']);need(all((v-r['unit']*w)%16==0 for zz,ww in[(x,X),(y,Y)]for v,w in zip(zz,ww)),'exact actual-to-formal all-cover odd unit')
 for r in family:
  i=r.pop('jobs');r['values']=[values[i][0],values[i+1][0]];r['contrast']=r['values'][0]-r['values'][1]
 need(family[0]['contrast']==-family[2]['contrast'] and family[0]['contrast']%16==8,'opposite actual ambient baseline signs');chosen=0 if family[0]['contrast']<0 else 2;need(family[chosen+1]['contrast']<=family[chosen]['contrast']<0,'one baseline-negative member survives crossing growth')
 deps=['even_cyclic_mod16/PROOF.md','even_cyclic_mod16/normal.json','monotone_dual_queries/PROOF.md','uniform_block_starts/verify.py']
 return {'status':'PASS','source_cover_tensor':[[*k,v]for k,v in sorted(A.items())],'target_cover_tensor':[[*k,v]for k,v in sorted(B.items())],'formal_profiles':formal,'all512_profile_evaluations':records,'physical_controls':controls,'crossing_root_family_control':family,'negative_family_query_index':chosen//2,'physical_graph_counts':len(jobs),'physical_profile_coordinates':sum(k for T,k in jobs),'source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in['PROOF.md','verify.py','cover_words.py','physical.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest()for f in deps},'scope':'Special fixed query on arbitrary even>=4 cyclic children has augmented profile contrast(8,0,...)mod16 and hence arbitrary fixed ambient gap8. Not a general mod16 transfer theorem.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
