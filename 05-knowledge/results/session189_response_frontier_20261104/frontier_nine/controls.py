from pathlib import Path
import json,math
D=Path(__file__).parent

def op(T):return tuple(((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T))
def blow(Q,As):
 B=[];n=0
 for A in As:B.append(list(range(n,n+len(A))));n+=len(A)
 T=[0]*n
 for i,A in enumerate(As):
  for j,r in enumerate(A):
   for k in range(len(A)):
    if r>>k&1:T[B[i][j]]|=1<<B[i][k]
  for k,C in enumerate(As):
   if Q[i]>>k&1:
    for u in B[i]:
     for v in B[k]:T[u]|=1<<v
 return T,B

def modules(T):
 n=len(T);full=(1<<n)-1
 return [M for M in range(1,full+1) if all((T[v]&M) in (0,M) for v in range(n) if not M>>v&1)]
def strong(T):
 ms=modules(T);return [A for A in ms if all(not(A&B) or A&B==A or A&B==B for B in ms)]
def quotient(T,children):
 v=[(C&-C).bit_length()-1 for C in children]
 return [sum(1<<j for j,w in enumerate(v) if T[u]>>w&1) for u in v]
def highest(T,U):
 ss=strong(T);changed=[];records=[]
 for M in ss:
  if M.bit_count()<2:continue
  sub=[A for A in ss if A!=M and A&M==A];ch=sorted([A for A in sub if not any(A!=B and A&B==A for B in sub)])
  a=quotient(T,ch);b=quotient(U,ch);records.append(dict(node=M,children=ch,source_quotient=a,target_quotient=b,changed=a!=b))
  if a!=b:changed.append(M)
 return [A for A in changed if not any(A!=B and A&B==A for B in changed)],records

def xor(T,family):
 U=T.copy()
 for i in range(len(T)):
  for j in range(i+1,len(T)):
   if sum(bool(M>>i&1 and M>>j&1) for M in family)%2:U[i]^=1<<j;U[j]^=1<<i
 return U

def Hsub(T):
 n=len(T);dp=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n);H[0]=1
 for i in range(n):dp[1<<i][i]=1
 for M in range(1,1<<n):
  for i in range(n):
   if M>>i&1:
    v=dp[M][i];todo=T[i]&~M
    while todo:
     b=todo&-todo;todo-=b;j=b.bit_length()-1;dp[M|b][j]+=v
  H[M]=sum(dp[M])
 return H

def profile(T):
 n=len(T);H=Hsub(T);P=[None]*(1<<n);P[0]=[1]
 for M in range(1,1<<n):
  v=[0]*(M.bit_count()+1);p=M&-M;A=M
  while A:
   if A&p:
    for k,c in enumerate(P[M^A]):v[k+1]+=H[A]*c
   A=(A-1)&M
  P[M]=v
 return [[v*math.factorial(k) for k,v in enumerate(row)] for row in P]
def induce(T,M):
 vs=[i for i in range(len(T)) if M>>i&1]
 return [sum(1<<j for j,w in enumerate(vs) if T[v]>>w&1) for v in vs]
def ext(T,a,c):
 n=len(T);R=[r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(T)]
 return R+[(1<<(n+1))|(((1<<n)-1)^a),((1<<n)-1)^c]
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def build():
 Q0=[4,17,42,131,45,75,31,119];Q1=[16,5,41,67,46,139,55,95];C3=[2,4,1];cases=[]
 T,blocks=blow([2,0],[Q0,C3]);M=sum(1<<v for v in blocks[0]);N=sum(1<<v for v in blocks[1]);family=[M,N];U=xor(T,family);V=xor(T,[M]);a=0;c=128
 cases.append(dict(name='ordered_endpoint_zero_prime8_and_triangle',source=T,intermediate=V,target=U,family=family,core=M,other=N,query=[a,c],core_child_masks=[1<<i for i in range(8)]))
 core,child=blow(Q1,[[2,0]]+[[0] for _ in range(7)]);T,blocks=blow(C3,[core,C3,[0]]);M=sum(1<<v for v in blocks[0]);N=sum(1<<v for v in blocks[1]);family=[M,N];U=xor(T,family);V=xor(T,[M]);outside=sum(1<<v for v in blocks[2]);a=outside;c=sum(1<<v for v in child[0])|outside
 cases.append(dict(name='cyclic_endpoint_zero_prime8_changed_descendant_and_other_highest',source=T,intermediate=V,target=U,family=family,core=M,other=N,query=[a,c],core_child_masks=[sum(1<<v for v in b) for b in child]))
 for r in cases:
  T,U,V=r['source'],r['target'],r['intermediate'];n=len(T);ms=modules(T);check(all(M in ms for M in r['family']),'every selected set is an original module');check(xor(T,r['family'])==U,'simultaneous edge XOR')
  hh,nodes=highest(T,U);check(set(hh)=={r['core'],r['other']},'exact highest-node antichain');r['original_modules']=ms;r['strong_node_records']=nodes;r['highest_changed_nodes']=hh
  pT,pU=profile(T),profile(U);check(pT==pU,'all fixed induced subsets preserve complete profiles');r['hereditary_profiles_source']=pT;r['hereditary_profiles_target']=pU
  a,c=r['query'];hs=[Hsub(ext(R,a,c))[-1] for R in [T,V,U]];check(hs[1]==hs[2],'exact other-node restoration');check((hs[0]-hs[1])%4==2,'actual ambient mod4 separator');r['actual_response_values']=hs
  z=r['core'];outside=((1<<n)-1)^z
  for v in range(n):
   if outside>>v&1:
    sample=(z&-z).bit_length()-1;want=bool(T[v]>>sample&1);check(bool(a>>v&1)==want and bool(c>>v&1)==want,'both physical ports inherit core exterior incidences')
  N=r['other'];check((a&N) in (0,N) and (c&N) in (0,N),'same physical query sees other node uniformly');check(induce(T,r['core'])!=induce(U,r['core']),'changed core retained')
  if 'descendant' in r['name']:check(induce(T,r['core_child_masks'][0])!=induce(U,r['core_child_masks'][0]),'selected root leaves a changed profile-equal pair child')
 json.dump(cases,open(D/'controls.json','w'),separators=(',',':'))
 models=[dict(case=r['name'],role=role,rows=r[role]) for r in cases for role in ['source','intermediate','target']]
 json.dump(models,open(D/'models.json','w'))
 with open(D/'models.in','w') as f:
  print(len(models),file=f)
  for r in models:print(len(r['rows']),*r['rows'],file=f)
 print(json.dumps([{'name':r['name'],'highest':r['highest_changed_nodes'],'responses':r['actual_response_values']} for r in cases],indent=2))
if __name__=='__main__':build()
