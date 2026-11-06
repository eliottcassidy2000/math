"""Exact all-word physical four-register certificate. Python 3 + SymPy."""
import itertools,json
from collections import Counter
from pathlib import Path
import sympy as S
from literal_forest import forest,read

def require(x,msg):
 if not x:raise RuntimeError(msg)
def patterns(labels,b):
 if b<0:return set()
 if not labels:return {tuple([(-1,-1)]*q)for q in range(b+1)}
 u,*rest=labels;out=set()
 for p in patterns(tuple(rest),b):
  for q in [p,p+((u,u),),p+((u,-1),)]:out.add(tuple(sorted(q)))
 for p in patterns(tuple(rest),b-1):out.add(tuple(sorted(p+((-1,u),))))
 for v in rest:
  for p in patterns(tuple(x for x in rest if x!=v),b):
   out.add(tuple(sorted(p+((u,v),))));out.add(tuple(sorted(p+((v,u),))))
 return out
def intro(parts,mask,new):
 out=Counter();left=[-1]+[i for i,(_,e)in enumerate(parts)if e==-1 or not(mask>>e&1)];right=[-1]+[i for i,(s,_)in enumerate(parts)if s>=0 and mask>>s&1]
 for i in left:
  for j in right:
   if i==j and i!=-1:continue
   q=[p for k,p in enumerate(parts)if k not in (i,j)];q.append((parts[i][0]if i!=-1 else new,parts[j][1]if j!=-1 else new));out[tuple(sorted(q))]+=1
 return out
P=sorted(patterns((0,1,2),2));I={p:i for i,p in enumerate(P)}
# Published operations 90 and 106, written without index dependence.
ops=[(3,3,(0,1,2)),(4,1,(3,0,2))]
M=[]
for mask,drop,remaining in ops:
 ren={x:i for i,x in enumerate(remaining)};ren[drop]=-1;ren[-1]=-1;m=S.zeros(len(P))
 for col,p in enumerate(P):
  for q,n in intro(p,mask,3).items():
   q=tuple(sorted((ren[a],ren[b])for a,b in q))
   if sum(a==-1 for a,b in q)<=2:m[I[q],col]+=n
 M.append(m)
C=S.zeros(2,len(P))
for col,p in enumerate(P):
 for q,a in intro(p,3,3).items():
  for z,b in intro(q,3,4).items():
   ss=int(len(z)==1 and z[0][0]==2);B=int(len(z)==1)-ss+int(len(z)in(1,2)and any(e==2 for _,e in z));C[0,col]+=a*b*B;C[1,col]+=a*b*ss
D=3*C[0,:]+2*C[1,:];s=C[1,:]
seed=(4,5,48,7,11,27);f,_=forest(seed,(0,1,2));x=S.Matrix([f.get(p,0)for p in P])
R=S.Matrix.vstack(D,D*M[0]-2*D,D*M[1]-2*D,s)
A=[S.Matrix([[2,1,0,0],[0,1,0,0],[S.Rational(1,6),S.Rational(7,6),0,-36],[0,0,0,1]]),S.Matrix([[2,0,1,0],[0,2,-1,-36],[S.Rational(1,6),S.Rational(13,6),-1,-78],[0,0,0,1]])]

def span(seed,mats):
 vectors=[];bas=[];piv=[];todo=[seed]
 while todo:
  v=todo.pop(0);r=list(v)
  for b,p in zip(bas,piv):
   t=r[p];r=[a-t*c for a,c in zip(r,b)]
  if not any(r):continue
  p=next(i for i,z in enumerate(r)if z);t=r[p];bas.append([z/t for z in r]);piv.append(p);vectors.append(v)
  todo.extend(m*v for m in mats)
 return S.Matrix.hstack(*vectors)
V=span(x,M);require(V.cols==8,'full reachable dimension');require((R*V).rank()==4,'observable dimension')
for i in range(2):require((R*M[i]-A[i]*R)*V==S.zeros(4,V.cols),'all-word intertwining')
require(R*x==S.Matrix([563,-23,-17,1]),'initial register')
# Cross-run identities proved by annihilation of an entire invariant span.
for i,j,start in [(0,1,1),(1,0,0)]:
 W=span(M[i]**start*x,[M[i]])
 require(D*M[j]*(M[i]**2-3*M[i]+2*S.eye(len(P)))*W==S.zeros(1,W.cols),'cross-run recurrence')
for n in range(1,3):require((D*M[1]*M[0]**n*x)[0]==1125*2**n-13,'zero-run initial')
for n in range(2):require((D*M[0]*M[1]**n*x)[0]==1050*2**n+53,'one-run initial')
# Literal graph controls use subset paths and unordered set partitions.
def append(rows,ports,mask):
 n=len(rows);rows=list(rows)+[0]
 for i in range(n):
  if i in ports and (mask>>ports.index(i)&1):rows[n]|=1<<i
  else:rows[i]|=1<<n
 return rows
checks=0
for k in range(3):
 for word in itertools.product(range(2),repeat=k):
  rows=list(seed);ports=[0,1,2];state=x
  for bit in word:
   mask,drop,rem=ops[bit];n=len(rows);rows=append(rows,ports,mask);ports=[(ports+[n])[a]for a in rem];state=M[bit]*state
  root=ports[2]
  for mask in [3,3]:n=len(rows);rows=append(rows,ports,mask);ports=ports+[n]
  ff,_=forest(rows,tuple(ports));h,ss,L,B=read(ff,2)
  require((B,ss)==tuple(C*state),'literal closing');checks+=1
print(json.dumps(dict(status='PASS',formal_dimension=len(P),reachable_dimension=V.cols,observable_dimension=4,literal_checks=checks,root_start=1)))
