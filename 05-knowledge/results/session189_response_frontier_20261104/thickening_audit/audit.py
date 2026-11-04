"""Independent pairwise module tests and complete strong-tree inflation audit."""
from pathlib import Path
import hashlib,json,math
D=Path(__file__).resolve().parent;P=D.parent/'thickening'
if not P.exists():P=D.parent/'session189-thickening'
def need(x,s):
 if not x:raise RuntimeError(s)
def modules(T):
 n=len(T);full=(1<<n)-1;out=[]
 for M in range(1,full+1):
  u=(M&-M).bit_length()-1;left=M^(1<<u);external=full^M;valid=True
  while left:
   b=left&-left;left-=b;v=b.bit_length()-1
   if (T[u]^T[v])&external:valid=False;break
  if valid:out.append(M)
 return out
def structure(T):
 ms=modules(T);ss=[A for A in ms if all(not(A&B) or A&B==A or A&B==B for B in ms)];ps=[];ls=[]
 for M in ss:
  if M.bit_count()<2:continue
  proper=[A for A in ss if A!=M and A&M==A];children=[A for A in proper if not any(A!=B and A&B==A for B in proper)];children.sort(key=lambda A:A&-A);need(sum(children)==M,'maximal proper strong partition')
  v=[(A&-A).bit_length()-1 for A in children];Q=[sum(1<<j for j,w in enumerate(v) if T[u]>>w&1) for u in v];tri=any((Q[i]>>j&1) and (Q[j]>>k&1) and (Q[k]>>i&1) for i in range(len(Q)) for j in range(len(Q)) for k in range(len(Q)))
  rec=dict(node=M,children=children,rows=Q)
  if tri:
   need(set(modules(Q))=={1<<i for i in range(len(Q))}|{(1<<len(Q))-1},'nonlinear strong quotient prime');ps.append(rec)
  else:ls.append(rec)
 return ms,ss,ps,ls

def inflate(T,N):
 n=len(T)*N;R=[]
 for u in range(n):
  a,i=divmod(u,N);row=0
  for v in range(n):
   b,j=divmod(v,N)
   if (a==b and i<j) or (a!=b and T[a]>>b&1):row|=1<<v
  R.append(row)
 return R

def profiles_mod16(T):
 n=len(T);dp=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n);H[0]=1
 for i in range(n):dp[1<<i][i]=1
 for S in range(1,1<<n):
  v=S
  while v:
   b=v&-v;v-=b;i=b.bit_length()-1;R=S^b
   if R:dp[S][i]=sum(dp[R][j] for j in range(n) if R>>j&1 and T[j]>>i&1)%16
  H[S]=sum(dp[S])%16
 C=[[1,0,0,0,0,0]]+[None]*((1<<n)-1)
 for S in range(1,1<<n):
  low=S&-S;out=[0]*6;A=S
  while A:
   if A&low:
    for k in range(5):out[k+1]=(out[k+1]+H[A]*C[S^A][k])%16
   A=(A-1)&S
  C[S]=out
 return [C[-1][k]*math.factorial(k)%16 for k in range(1,6)]
def tournament(n,code):
 R=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>k&1:R[i]|=1<<j
   else:R[j]|=1<<i
   k+=1
 return R

def main():
 data=json.loads((P/'structural_controls.json').read_text());records=[];subset_checks=prime_nodes=0
 for rec in data:
  T=rec['old_rows'];N=rec['inflation'];U=inflate(T,N);need(U==rec['physical_rows'],'all direct physical arcs in inflation');cells=[list(range(N*i,N*(i+1))) for i in range(len(T))];need(cells==rec['physical_cells'],'all fixed physical fiber names')
  old,os,op,ol=structure(T);new,ns,np,nl=structure(U);need(old==rec['old_modules'] and new==rec['inflated_modules'],'all modules by independent pairwise exterior-difference test');need(os==rec['old_strong'] and ns==rec['inflated_strong'],'all strong modules independently')
  def key(p):return p['node'],tuple(p['children']),tuple(p['rows'])
  need(sorted(map(key,op))==sorted(map(key,rec['old_prime_nodes'])) and sorted(map(key,np))==sorted(map(key,rec['inflated_prime_nodes'])),'all prime quotient records')
  def lift(M):return sum(((1<<N)-1)<<(i*N) for i in range(len(T)) if M>>i&1)
  expected=[dict(node=lift(p['node']),children=[lift(c) for c in p['children']],rows=p['rows']) for p in op];need(sorted(map(key,expected))==sorted(map(key,np)),'exact prime quotient and child-set persistence');prime_nodes+=len(np);subset_checks+=(1<<len(T))+(1<<len(U))-2
  records.append(dict(label=rec['label'],N=N,old_rows=T,physical_rows=U,physical_cells=cells,old_modules=old,new_modules=new,old_strong=os,new_strong=ns,old_prime_nodes=op,new_prime_nodes=np,old_linear_nodes=ol,new_linear_nodes=nl))
 # Independent actual profile closure checks. TT4 itself has g profile mod16.
 closure=[]
 for n in range(1,4):
  for code in range(1<<(n*(n-1)//2)):
   T=tournament(n,code);U=inflate(T,4);q=profiles_mod16(U);need(q==[1,14,4,8,0],'actual substituted TT4 complete g profile mod16');closure.append(dict(base_rows=T,N=4,actual_rows=U,q_1_through_5_mod16=q,h_ge_6_zero_by_factorial=True))
 outputs={'structural_recount.json':records,'profile_closure_controls.json':closure}
 for name,obj in outputs.items():
  txt=json.dumps(obj,separators=(',',':'))+'\n';p=D/name
  if p.exists():need(p.read_text()==txt,'independent artifact stable')
  else:p.write_text(txt)
 return {'status':'PASS','literal_inflations':len(records),'all_nonempty_subsets_checked':subset_checks,'prime_node_child_comparisons':prime_nodes,'actual_complete_profile_closure_cases':len(closure),'maximum_structural_order':max(len(r['physical_rows']) for r in records),'proof_sha256':hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest(),'producer_data_sha256':hashlib.sha256((P/'structural_controls.json').read_bytes()).hexdigest(),'audit_source_sha256':hashlib.sha256((D/'audit.py').read_bytes()).hexdigest(),'output_pins':{name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs},'scope':'All-size proof reviewed separately. Small structural inflations are not claimed to meet every Seidel exponent; actual TT4 closure controls verify the profile lemma modulo16, not an arbitrary-size response census.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
