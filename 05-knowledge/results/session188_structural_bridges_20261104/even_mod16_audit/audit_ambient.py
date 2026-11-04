#!/usr/bin/env python3
"""Independent exact cover-cut automaton and physical ambient certificate audit."""
from pathlib import Path
from itertools import product
from functools import lru_cache
import hashlib,json,subprocess,tempfile
from cover_automaton import cover_kernel
from cover_physical import main as small_physical
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'even_cyclic_mod16_ambient';BASE=ROOT.parent
PINS={'PROOF.md':'afa692b80873b3f49e05767c8ff316405ea3df0bf2e2738ff99ae491d4118b83','verify.py':'97b1acb85e6dab0c593809362a4007a9fceee5974dff9901f6f9d9f0b2ec66a4','cover_words.py':'9055ff08ea35ac726cd6bb04a24b3a3822c81a9fdfbb379af80c024e5d932a21','physical.cpp':'4d980994dfafa599ddb065768ce34fdf54abacc936a9fb89e177834e38b3152c','normal.json':'ba3c4ac5b4a1e2d8faf485244a60a31193938bce5f30b1fc421458e6d882b62b','optimized.json':'ba3c4ac5b4a1e2d8faf485244a60a31193938bce5f30b1fc421458e6d882b62b'}
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def reverse(T,mm):
 R=list(T)
 for M in mm:
  for i in range(len(T)):
   if M>>i&1:R[i]^=M^(1<<i)
 return R
def blow(Q,B):
 off=[sum(map(len,B[:i])) for i in range(len(B))];R=[0]*sum(map(len,B))
 for i,C in enumerate(B):
  for v,r in enumerate(C):
   R[off[i]+v]=r<<off[i]
   for j,D in enumerate(B):
    if Q[i]>>j&1:R[off[i]+v]|=((1<<len(D))-1)<<off[j]
 return R,off
def augment(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for i in range(n):
  if a>>i&1:R[i]|=1<<n
  else:R[n]|=1<<i
  if c>>i&1:R[i]|=1<<(n+1)
  else:R[n+1]|=1<<i
 return tuple(R)
def endpoints(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for M in range(1,1<<n):
  for i in range(n):
   if M>>i&1 and M!=(1<<i):D[M][i]=sum(D[M^(1<<i)][j] for j in range(n) if T[j]>>i&1)
 return D
@lru_cache(None)
def profile(T):
 n=len(T);N=1<<n;H=list(map(sum,endpoints(T)));q=[[0]*N for _ in range(n+1)];q[0][0]=1
 for h in range(1,n+1):
  for M in range(1,N):
   part=M
   while part:q[h][M]+=H[part]*q[h-1][M^part];part=(part-1)&M
 return [q[h][-1] for h in range(n+1)]
def modules(T):return [M for M in range(1,1<<len(T)) if M.bit_count()>1 and all(T[i]&M in(0,M) for i in range(len(T)) if not M>>i&1)]
def subsets(M):
 S=M
 while True:
  yield S
  if not S:break
  S=(S-1)&M

def build(blocks,outer,query=(1,3),target=None):
 T,off=blow([2,4,1],blocks);U,_=blow([4,1,2],blocks if target is None else target);n=len(T);a,c=[sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,off)) if mask>>i&1) for mask in query];ext=sum(1<<(n+j-1) for j in range(1,len(outer)) if outer[j]&1);T,_=blow(outer,[T]+[[0]]*(len(outer)-1));U,_=blow(outer,[U]+[[0]]*(len(outer)-1));return T,U,a|ext,c|ext,n

def main():
 for f,h in PINS.items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'final source '+f)
 d=json.loads((SRC/'normal.json').read_text());need(d['status']=='PASS','source status')
 for f,h in d['dependency_pins'].items():need(hashlib.sha256((BASE/f).read_bytes()).hexdigest()==h,'dependency '+f)
 types=[(0,1,b,c,8,e) for b,c,e in product((6,14),(4,12),(0,8))];need(types==list(map(tuple,d['formal_profiles'])),'complete eight-profile upper domain')
 K=[cover_kernel(1,3,r) for r in(False,True)];tensors=[]
 for kk,key in zip(K,('source_cover_tensor','target_cover_tensor')):
  tensor={(h,*e):v for h in range(1,6) for e,v in kk[h].items() if min(e)>0};need(tensor=={tuple(z[:4]):z[4] for z in d[key]},'full exact integer cover tensor, not only residues');tensors.append(tensor)
 local=json.loads((BASE/'even_cyclic_mod16'/'normal.json').read_text());at=local['block_query_order'].index([1,3]);need(len(d['all512_profile_evaluations'])==512,'complete all-cover domain');evaluations=[]
 for ids,r,old in zip(product(range(8),repeat=3),d['all512_profile_evaluations'],local['complete_formal_records']):
  need(tuple(r['profile_type_indices'])==ids==tuple(old['profile_type_indices']),'every ordered profile triple');P=[types[i] for i in ids];vals=[[sum(v*P[0][i]*P[1][j]*P[2][k] for (h,i,j,k),v in t.items() if h==hh) for hh in range(1,6)] for t in tensors];need(vals==[r['source_first_five'],r['target_first_five']],'all exact first-five source/target profile vectors');gap=[a-b for a,b in zip(*vals)];need(gap==r['integer_contrasts'] and gap[0]==old['integer_functionals'][at] and [z%16 for z in gap]==[8,0,0,0,0],'exact Hamilton link and entire residue identity');evaluations.append(dict(types=ids,source=vals[0],target=vals[1],contrast=gap))
 # Separate actual named cover convolution checks the cut interpretation.
 small=small_physical();raw=(json.dumps(small,sort_keys=True,indent=2)+'\n').encode();need(raw==(ROOT/'cover_physical.normal.json').read_bytes()==(ROOT/'cover_physical.optimized.json').read_bytes(),'independent small physical cover replay')
 jobs=[];expected=[]
 def retain(T,a,c,N,values,ncore):
  G=augment(T,a,c);M=((1<<ncore)-1)|(1<<len(T))|(1<<(len(T)+1));need(all(G[v]&M in(0,M) for v in range(ncore,len(T))),'same actual augmented core module');jobs.append((G,N));expected.append(values)
 for r in d['physical_controls']:
  B=r['child_rows'];outer=r['outer_rows'];T,U,a,c,n=build(B,outer);need(T==r['source_old_rows'] and U==r['target_old_rows'] and [a,c]==r['physical_query'],'fixed physical masks and unchanged exterior');qs=[profile(tuple(b)) for b in B];need(qs==r['profiles'],'exact actual child profiles');ts=[];unit=1
  for q in qs:
   inv=pow(q[1],-1,16);norm=(0,)+tuple((q[k] if k<len(q) else 0)*inv%16 for k in range(1,6));need(norm in types and all(z%16==0 for z in q[6:]),'actual entire profile domain');ts.append(types.index(norm));unit=unit*q[1]%16
  need(ts==r['types'] and unit==r['unit'],'profile types and retained odd unit');vals=[r['source_values'],r['target_values']];gap=[s-t for s,t in zip(*vals)];need(gap==r['integer_contrasts'] and gap[0]%16==8,'actual ambient residue8');N=r['coordinates_counted'];need(N==(5 if len(outer)==1 else 1),'profile versus Hamilton job semantics')
  for V,v in zip((T,U),vals):retain(V,a,c,N,v,n)
  if N==5:
   rec=evaluations[ts[0]*64+ts[1]*8+ts[2]];need(all((v-unit*w)%16==0 for vv,ww in zip(vals,(rec['source'],rec['target'])) for v,w in zip(vv,ww)) and [z%16 for z in gap]==[8,0,0,0,0],'physical full-profile normalization')
 B=d['physical_controls'][0]['child_rows'];child=reverse(B[0],[3,6]);target=[op(child),op(B[1]),op(B[2])];need(all(all(y>=x for x,y in zip(profile(tuple(b)),profile(tuple(c)))) for b,c in zip(B,target)) and profile(tuple(target[0]))[1]>profile(tuple(B[0]))[1],'crossing profile monotonicity')
 family=d['crossing_root_family_control']
 for index,r in enumerate(family):
  query=(1,3) if index<2 else(4,6);growth=bool(index%2);need(tuple(r['block_query'])==query and r['profile_increased']==growth,'complete dual baseline/growth cases');T,U,a,c,n=build(B,[2,4,1],query,target if growth else None);need(T==r['source_old_rows'] and U==r['target_old_rows'] and [a,c]==r['physical_query'],'actual crossing original-module target');need(r['contrast']==r['values'][0]-r['values'][1],'literal sign accounting')
  for V,v in zip((T,U),r['values']):retain(V,a,c,1,[v],n)
 need(family[0]['contrast']==-family[2]['contrast'] and family[0]['contrast']%16==8,'opposite actual ambient baseline signs');chosen=0 if family[0]['contrast']<0 else 2;need(d['negative_family_query_index']==chosen//2 and family[chosen+1]['contrast']<=family[chosen]['contrast']<0,'fixed negative member handles crossing growth')
 need(len(jobs)==d['physical_graph_counts']==38 and sum(N for G,N in jobs)==d['physical_profile_coordinates']==70,'full physical occurrence/coordinate census')
 with tempfile.TemporaryDirectory(prefix='even-cover-ambient-audit-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'ambient_count.cpp'),'-o',str(exe)],check=True);wire=str(len(jobs))+'\n'+'\n'.join(str(len(G))+' '+str(N)+' '+' '.join(map(str,G)) for G,N in jobs)+'\n';rows=[list(map(int,line.split())) for line in subprocess.check_output([str(exe)],input=wire.encode()).decode().splitlines()];need(rows==expected,'all actual profile coordinates via predecessor HP and ordered subset convolution')
 return dict(status='PASS',source_sha256=PINS,dependencies=d['dependency_pins'],complete_source_tensor=[[*k,v] for k,v in sorted(tensors[0].items())],complete_target_tensor=[[*k,v] for k,v in sorted(tensors[1].items())],formal_profile_evaluations=evaluations,independent_small_physical_graphs=small['actual_augmented_graphs'],independent_small_cover_coordinates=small['integer_cover_coordinates'],small_physical_receipt_sha256=hashlib.sha256(raw).hexdigest(),actual_graph_occurrences=len(jobs),actual_profile_coordinates=70,complete_actual_values=rows)
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
