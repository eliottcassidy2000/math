from pathlib import Path
from itertools import combinations,permutations
from math import comb
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
DEPS=['prime_five_frontier/PROOF.md','prime_five_frontier/normal.json','prime_six_reconnaissance/PROOF.md','prime_six_reconnaissance/endpoint_table.tsv','prime_six_reconnaissance/normal.json','module_decomposition/PROOF.md','highest_quotient_reduction/PROOF.md','monotone_dual_queries/PROOF.md','cyclic_ambient_synthesis/PROOF.md','prime_five_all_ambient/PROOF.md']
EXPECTED={'prime_five_frontier/PROOF.md':'3ede97d419ae29f637d8b4ff338cf42a4fbe575f57ec44706828bdcc5ad9feb6','prime_six_reconnaissance/PROOF.md':'f270aec9732d87d2b7f9ec685b9655d7efeea00b89fbb052d5b73af6dfcfe894','prime_six_reconnaissance/endpoint_table.tsv':'902de128fc24d1460d611fdb1136bce64453a5d891d100f1c38ee46d399a4d56','prime_six_reconnaissance/normal.json':'83ad3788c8f58ce7633c764e9cc07a9ddf337d0878c00a819cc69a1d32a551f7'}
def need(c,m):
 if not c:raise RuntimeError(m)
def op(T):return [((1<<len(T))-1)^(1<<i)^v for i,v in enumerate(T)]
def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def blow(Q,B):
 off=[];n=0
 for b in B:off.append(n);n+=len(b)
 T=[]
 for i,b in enumerate(B):T.extend((v<<off[i])|sum(((1<<len(z))-1)<<off[j] for j,z in enumerate(B) if Q[i]>>j&1) for v in b)
 return T,off
def restrict(T,M):
 vs=[i for i in range(len(T)) if M>>i&1];return [sum(((T[i]>>j)&1)<<k for k,j in enumerate(vs)) for i in vs]
def mods(T):return [M for M in range(1,1<<len(T)) if all(T[i]&M in(0,M) for i in range(len(T)) if not M>>i&1)]
def overlaps(A,B):return bool(A&B) and bool(A&~B) and bool(B&~A)
def quotient(T,children):
 reps=[(M&-M).bit_length()-1 for M in children];return [sum(((T[i]>>j)&1)<<k for k,j in enumerate(reps)) for i in reps]
def structure(T,U):
 ms=mods(T);target=set(mods(U));strong=[M for M in ms if not any(overlaps(M,N) for N in ms)];nodes=[];changed=[]
 for M in strong:
  if M.bit_count()<2:continue
  proper=[N for N in strong if N!=M and N&M==N];children=[N for N in proper if not any(N!=K and N&K==N for K in proper)];children.sort(key=lambda N:(N&-N).bit_length());need(all(N in target for N in children),'persistent original child modules');a=quotient(T,children);b=quotient(U,children)
  if a!=b:changed.append(M)
  nodes.append(dict(mask=M,children=children,source_quotient=a,target_quotient=b))
 highest=[M for M in changed if not any(M!=N and M&N==M for N in changed)]
 return ms,strong,nodes,highest

def profile(T):
 n=len(T);q=[0]*(n+1)
 for w in permutations(range(n)):
  bad=sum(not(T[i]>>j&1) for i,j in zip(w,w[1:]));good=n-1-bad
  for h in range(bad+1,n+1):q[h]+=comb(good,h-1-bad)
 return q
def aug(T,a,c):
 n=len(T);U=T[:]+[1<<(n+1),0]
 for i in range(n):
  if a>>i&1:U[i]|=1<<n
  else:U[n]|=1<<i
  if c>>i&1:U[i]|=1<<(n+1)
  else:U[n+1]|=1<<i
 return U
def xor(T,family):
 U=T[:]
 for i,j in combinations(range(len(T)),2):
  if sum(bool(M>>i&1 and M>>j&1) for M in family)%2:U[i]^=1<<j;U[j]^=1<<i
 return U
def replace(T,M,B):
 vs=[i for i in range(len(T)) if M>>i&1];U=T[:]
 for i,v in enumerate(vs):U[v]=(U[v]&~M)|sum(((B[i]>>j)&1)<<w for j,w in enumerate(vs))
 return U
def main():
 deps={f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in DEPS}
 for f,h in EXPECTED.items():need(deps[f]==h,'frozen dependency '+f)
 for p in ['prime_five_frontier/normal.json','prime_six_reconnaissance/normal.json']:need(json.loads((ROOT.parent/p).read_text())['status']=='PASS','dependency receipt status')
 Q6=[34,8,3,5,15,30];Q5=[24,17,3,6,12];Q8=[238,4,240,118,195,146,34,74];C3=[2,4,1];records=[];jobs=[]
 boundary=json.loads((ROOT.parent/'prime_six_reconnaissance/normal.json').read_text())['prime8_boundary'];need(boundary['rows']==Q8 and boundary['p']==0 and all((a+b)%2==0 for a,b in zip(boundary['S'],boundary['E'])),'pinned endpoint-balanced other prime node')
 def job(T):jobs.append(T);return len(jobs)-1
 for title,children,other,descendant in [('six_with_descendant_and_five',[tt(3)]+[tt(1)]*5,Q5,True),('six_and_endpoint_balanced_eight',[tt(1)]*6,Q8,False)]:
  Mcore,offs=blow(Q6,children);T,outeroffs=blow(C3,[Mcore,other,[0]]);m=len(Mcore);M=(1<<m)-1;N=((1<<len(other))-1)<<m;family=[M,N]+([3] if descendant else[]);U=xor(T,family);ms,strong,nodes,highest=structure(T,U);need(all(F in ms for F in family),'all selected sets original modules');need(set(highest)=={M,N},'exact two highest nodes');need(M in strong and N in strong,'both highest nodes originally strong')
  chosen=next(r for r in nodes if r['mask']==M);need(chosen['source_quotient']==Q6 and chosen['target_quotient']==op(Q6),'six-child quotient converse on fixed partition')
  p=0;S=[0]*6;E=S[:]
  for w in permutations(range(6)):
   if all(Q6[i]>>j&1 for i,j in zip(w,w[1:])):S[w[0]]+=1;E[w[-1]]+=1
  for i in range(6):p|=((S[i]+E[i])%2)<<i
  need(p==63,'literal quotient endpoint vector');anchor=0;A=chosen['children'][anchor];a=A|(1<<(len(T)-1));c=1<<(len(T)-1)
  # V has the actual target M interior only; B restores all chosen children to source interiors.
  V=replace(T,M,restrict(U,M));baselinecore,_=blow(op(Q6),children);B=replace(T,M,baselinecore);need(replace(V,N,restrict(U,N))==U,'exact restoration reaches actual target')
  profile_records=[]
  for label,K in [('other_highest',N)]+[(f'child{i}',K) for i,K in enumerate(chosen['children'])]+[('chosen_whole',M)]:
   qs=[profile(restrict(X,K)) for X in (T,U)];need(qs[0]==qs[1],'required whole/child complete profile equality');profile_records.append([label,K,qs[0]])
  # Ports see the other highest node uniformly in every physical graph.
  need(a&N in(0,N) and c&N in(0,N),'uniform physical other-node ports')
  ids=[job(aug(X,a,c)) for X in (T,B,V,U)];oldids=[job(X) for X in (T,B,V,U)];records.append(dict(name=title,source=T,restored_baseline=B,chosen_target_only=V,target=U,original_family=family,original_modules=ms,original_strong_nodes=strong,original_tree_nodes=nodes,highest=highest,chosen=M,other=N,chosen_anchor=anchor,query=[a,c],profile_equalities=profile_records,query_ids=ids,old_ids=oldids))
 with tempfile.TemporaryDirectory(prefix='frontier-seven-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);wire=str(len(jobs))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in jobs);vv=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()))
 need(len(vv)==len(jobs),'all literal physical counts')
 for r in records:
  r['physical_query_counts']=[vv[i] for i in r.pop('query_ids')];r['old_H_counts']=[vv[i] for i in r.pop('old_ids')];vals=r['physical_query_counts'];need(len(set(r['old_H_counts']))==1,'scalar equality');need(vals[1]==vals[2]==vals[3] and (vals[0]-vals[3])%4==2,'actual descendant and disjoint-node restoration preserves separator')
 pins={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','count.cpp']}
 return dict(status='PASS',source_pins=pins,dependency_pins=deps,multiple_highest_node_controls=records,physical_graphs=len(jobs),maximum_augmented_order=max(map(len,jobs)),scope='Necessary highest changed prime quotient order at least seven for nontrivial original-module full-response collisions. Source-relative family restrictions and fixed exterior are explicit; no general reconstruction claim.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
