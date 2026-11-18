from itertools import combinations
import json

def require(x):
 if not x:raise AssertionError

def tables(out):
 n=len(out);D=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n);H[0]=1
 for i in range(n):D[1<<i][i]=1
 for m in range(1,1<<n):
  for v in range(n):
   if m>>v&1:
    for u in range(n):
     if m>>u&1 and out[u]>>v&1:D[m][v]+=D[m^(1<<v)][u]
  H[m]=sum(D[m])
 return H,D

def conv(f,g):
 z=[0]*len(f)
 for m in range(len(f)):
  s=m
  while True:
   z[m]+=f[s]*g[m^s]
   if not s:break
   s=(s-1)&m
 return z

def slots(H,qmax):
 ans=[[1]+[0]*(len(H)-1)]
 for q in range(qmax):ans.append(conv(ans[-1],H))
 return [a[-1] for a in ans]

def converse(out):return [sum(((out[j]>>i)&1)<<j for j in range(len(out))) for i in range(len(out))]

def shell(T,v=None):
 n=len(T);a=n;b=n+1;out=T.copy()+[0,1<<a]
 for i in range(n):out[i]|=1<<b
 out[a]=(1<<n)-1
 if v is not None:out[a]^=1<<v;out[v]|=1<<a
 return out

def reflection(n,S):return [sum(1<<j for j in range(n) if j!=i and ((i<j)!=(i+j+2 in S))) for i in range(n)]

def nested_poly(m):
 f=[0]*(m%2)+[1]
 from math import comb
 for _ in range(m//2):
  z=[0]*(len(f)+2)
  for j,c in enumerate(f):
   z[j+2]+=c
   for k in range(1,j+1,2):z[j-k+1]+=2*c*comb(j,k)
  f=z
 return f

def profile(T,v):
 H,E=tables(T);_,S=tables(converse(T));F=len(H)-1
 ks=sum(S[m][v]*H[F^m] for m in range(F+1) if m>>v&1)
 ke=sum(E[m][v]*H[F^m] for m in range(F+1) if m>>v&1)
 return H[-1],conv(H,H)[-1],S[-1][v],E[-1][v],ks,ke
shellchecks=brokenchecks=0
for n in range(1,6):
 pairs=list(combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  T=[0]*n
  for k,(i,j) in enumerate(pairs):
   if bits>>k&1:T[i]|=1<<j
   else:T[j]|=1<<i
  H,_=tables(T)
  if n<=4:
   P=slots(H,4);Q=slots(tables(shell(T))[0],3)
   for q in range(1,4):require(Q[q]==q*q*P[q]+q*(P[q+1]-P[q-1]));shellchecks+=1
  for v in range(n):
   h,p,s,e,ks,ke=profile(T,v);actual=tables(shell(T,v))[0][-1]
   require(actual==h+p+ke-ks-e-s);brokenchecks+=1
single=[]
for n in range(1,13):
 for s in range(2,2*n+1):
  m=min(s-1,2*n-s+1);H=tables(reflection(n,{s}))[0][-1]
  require(H==sum(nested_poly(m)));single.append((n,s,H))
chain=[]
for n in range(2,13):
 T=reflection(n,{n,n+1});h=tables(T)[0][-1]
 if n>=3:
  old=reflection(n-2,{n-2,n-1});U=shell(old,n-3)
  # shell order [old,a,b] corresponds [2,...,n-1,1,n]
  order=list(range(1,n-1))+[0,n-1]
  remap=[sum(((T[order[i]]>>order[j])&1)<<j for j in range(n)) for i in range(n)]
  require(U==remap)
 chain.append((n,h))
print(json.dumps({'shell_identity_checks':shellchecks,'broken_shell_checks':brokenchecks,'single_reflection_checks':len(single),'two_reflection_counts':chain,'nested_counts':[(n,sum(nested_poly(n))) for n in range(13)],'status':'PASS'},indent=2))
