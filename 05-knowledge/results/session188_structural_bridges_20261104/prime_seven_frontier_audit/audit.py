"""Independent original-tree, subset-cover and forward physical audit."""
from pathlib import Path
from itertools import combinations
import hashlib,json,subprocess,tempfile
HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'prime_seven_frontier'
def need(c,m):
 if not c:raise RuntimeError(m)
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def inflate(Q,blocks):
 labels=[(i,j) for i,B in enumerate(blocks) for j in range(len(B))];R=[]
 for i,u in labels:
  R.append(sum(1<<w for w,(j,v) in enumerate(labels) if (blocks[i][u]>>v&1 if i==j else Q[i]>>j&1)))
 return R

def modules(T):
 full=(1<<len(T))-1;answer=[]
 for S in range(1,full+1):
  base=(S&-S).bit_length()-1;outside=full^S
  if all(not((T[base]^T[v])&outside) for v in range(len(T)) if S>>v&1):answer.append(S)
 return answer

def tree(T,U):
 ms=modules(T);strong=[S for S in ms if all(not(S&K and S&~K and K&~S) for K in ms)];nodes={};changed=[]
 for S in strong:
  if S.bit_count()==1:continue
  cs=[]
  for K in strong:
   if K==S or K&~S:continue
   if not any(K!=L and L!=S and not(K&~L) and not(L&~S) for L in strong):cs.append(K)
  cs.sort(key=lambda s:(s&-s).bit_length());reps=[(s&-s).bit_length()-1 for s in cs]
  qr=lambda R:[sum(((R[u]>>v)&1)<<j for j,v in enumerate(reps)) for u in reps]
  a,b=qr(T),qr(U);nodes[S]={'mask':S,'children':cs,'source_quotient':a,'target_quotient':b}
  if a!=b:changed.append(S)
 highest=[S for S in changed if not any(S!=K and not(S&~K) for K in changed)]
 return ms,strong,nodes,highest

def restrict(T,S):
 vs=[i for i in range(len(T)) if S>>i&1]
 return [sum(((T[u]>>v)&1)<<j for j,v in enumerate(vs)) for u in vs]
def profiles(T):
 n=len(T);N=1<<n;d=[[0]*n for _ in range(N)]
 for v in range(n):d[1<<v][v]=1
 for S in range(1,N):
  for v,z in enumerate(d[S]):
   if z:
    for w in range(n):
     if T[v]>>w&1 and not S>>w&1:d[S|1<<w][w]+=z
 h=[sum(row) for row in d];q=[[0]*(n+1) for _ in range(N)];q[0][0]=1
 # Unordered partitions, choosing the block of the least remaining physical name.
 for S in range(1,N):
  v=S&-S;A=S
  while A:
   if A&v:
    for k in range(1,n+1):q[S][k]+=h[A]*q[S^A][k-1]
   A=(A-1)&S
 fac=1
 for k in range(1,n+1):fac*=k;q[-1][k]*=fac
 return q[-1]
def xor(T,family):
 U=T[:]
 for i,j in combinations(range(len(T)),2):
  if sum(((S>>i)&1)*((S>>j)&1) for S in family)%2:U[i]^=1<<j;U[j]^=1<<i
 return U
def replace(T,S,B):
 vs=[i for i in range(len(T)) if S>>i&1];R=T[:]
 for i,u in enumerate(vs):R[u]=(T[u]&~S)|sum(((B[i]>>j)&1)<<v for j,v in enumerate(vs))
 return R
def aug(T,a,c):
 n=len(T);R=T[:]+[1<<(n+1),0]
 for v in range(n):
  for bits,w in [(a,n),(c,n+1)]:
   if bits>>v&1:R[v]|=1<<w
   else:R[w]|=1<<v
 return R

def main():
 d=json.loads((SRC/'normal.json').read_text());need(d['status']=='PASS','source receipt')
 for f,h in d['dependency_pins'].items():need(hashlib.sha256((HERE.parent/f).read_bytes()).hexdigest()==h,'frozen dependency '+f)
 for name in ['prime_five_frontier_audit','prime_six_endpoint_audit']:
  need(json.loads((HERE.parent/name/'normal.json').read_text())['status']=='PASS','separate reviewed input '+name)
 Q6=[34,8,3,5,15,30];Q5=[24,17,3,6,12];Q8=[238,4,240,118,195,146,34,74];C3=[2,4,1];jobs=[];expect=[];records=[]
 for case,r in enumerate(d['multiple_highest_node_controls']):
  blocks=([[6,4,0]]+[[0]]*5) if case==0 else [[0]]*6;other=Q5 if case==0 else Q8;core=inflate(Q6,blocks);T=inflate(C3,[core,other,[0]]);m=len(core);M=(1<<m)-1;N=((1<<len(other))-1)<<m;family=[M,N]+([3] if case==0 else []);U=xor(T,family);need(T==r['source'] and U==r['target'] and family==r['original_family'],'independent physical source and original-family XOR')
  ms,strong,nodes,highest=tree(T,U);need(ms==r['original_modules'] and strong==r['original_strong_nodes'],'all modules and original strong nodes');need(nodes=={x['mask']:x for x in r['original_tree_nodes']},'entire original tree and quotient action');need(set(highest)==set(r['highest'])=={M,N},'exact highest-node antichain')
  need(all(S in ms for S in family),'selected sets are original modules');need(all(S in modules(U) for S in strong),'common original recursive nodes');need(nodes[M]['source_quotient']==Q6 and nodes[M]['target_quotient']==op(Q6),'six quotient converse');B=replace(T,M,inflate(op(Q6),blocks));V=replace(T,M,restrict(U,M));need(B==r['restored_baseline'] and V==r['chosen_target_only'] and replace(V,N,restrict(U,N))==U,'both restoration stages')
  coverchecks=[]
  for label,S,q in r['profile_equalities']:
   found=[profiles(restrict(W,S)) for W in [T,U]];need(found[0]==found[1]==q,'independent unordered-partition complete profile '+label);coverchecks.append([label,S,found[0]])
  a,c=r['query'];anchor=nodes[M]['children'][0];need(a==anchor|(1<<(len(T)-1)) and c==1<<(len(T)-1),'same physical root-incidence endpoint query');need((a&N) in(0,N) and (c&N) in(0,N),'uniform other highest-node ports')
  for W,val in zip([T,B,V,U],r['physical_query_counts']):jobs.append(aug(W,a,c));expect.append(val)
  for W,val in zip([T,B,V,U],r['old_H_counts']):jobs.append(W);expect.append(val)
  need(len(set(r['old_H_counts']))==1 and r['physical_query_counts'][1:]==[r['physical_query_counts'][1]]*3 and (r['physical_query_counts'][0]-r['physical_query_counts'][3])%4==2,'source-relative exact restoration and gap')
  records.append({'name':r['name'],'highest':highest,'strong_nodes':strong,'profiles':coverchecks,'query':[a,c],'counts':r['physical_query_counts'],'old_counts':r['old_H_counts']})
 wire=''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in jobs)
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'f';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'forward.cpp'),'-o',str(exe)],check=True);actual=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()))
 need(actual==expect,'all16 actual counts by independent forward extension');need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired frozen source receipts')
 print(json.dumps({'status':'PASS','source_pins':{f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','count.cpp','normal.json','optimized.json']},'dependency_pins':d['dependency_pins'],'audit_pins':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['audit.py','forward.cpp']},'controls':records,'physical_graphs':len(jobs),'maximum_order':max(map(len,jobs)),'physical_values':actual},sort_keys=True,indent=2))
if __name__=='__main__':main()
