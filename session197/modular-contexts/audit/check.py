#!/usr/bin/env python3
"""Independent: literal odd-cycle collections, exact basis recovery, finite machines."""
import math,json,hashlib
from functools import lru_cache
from itertools import combinations,product
from collections import deque

def ck(x,msg):
 if not x:raise RuntimeError(msg)
def trim(p):
 p=list(p)
 while len(p)>1 and p[-1]==0:p.pop()
 return p
def plus(a,b):
 p=[0]*max(len(a),len(b))
 for i,v in enumerate(a):p[i]+=v
 for i,v in enumerate(b):p[i]+=v
 return trim(p)
def times(a,b):
 p=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):p[i+j]+=x*y
 return trim(p)
def scale(a,b):return trim([x*b for x in a])
def ev(p,x):
 z=0
 for c in p[::-1]:z=z*x+c
 return z
def tour(n,c):
 A=[0]*n
 for i,(u,v)in enumerate(combinations(range(n),2)):
  if c>>i&1:A[u]|=1<<v
  else:A[v]|=1<<u
 return tuple(A)
def tr(n):return tuple(sum(1<<j for j in range(i+1,n))for i in range(n))
def subst(Q,r,T):
 labels=[(i,0)for i in range(len(Q))if i!=r]+[(r,j)for j in range(len(T))]
 return tuple(sum(1<<j for j,(v,b)in enumerate(labels)if i!=j and ((T[a]>>b&1)if u==v else(Q[u]>>v&1)))for i,(u,a)in enumerate(labels))
def cy(A,B,C):
 sizes=[len(A),len(B),len(C)];blocks=[A,B,C];offset=[0,sizes[0],sizes[0]+sizes[1]];R=[]
 for k in range(3):
  nxt=(k+1)%3;mask=((1<<sizes[nxt])-1)<<offset[nxt]
  R.extend((z<<offset[k])|mask for z in blocks[k])
 return tuple(R)

@lru_cache(None)
def prof(A):
 # Every cyclic word is rooted at its least vertex, so rotations occur exactly once.
 n=len(A);cycles={}
 for first in range(n):
  allowed=((1<<n)-1)^((1<<(first+1))-1)
  def dfs(last,mask,length):
   if length>=3 and length%2 and A[last]>>first&1:cycles[mask]=cycles.get(mask,0)+1
   options=A[last]&allowed&~mask
   while options:
    bit=options&-options;options-=bit;dfs(bit.bit_length()-1,mask|bit,length+1)
  dfs(first,1<<first,1)
 byvertex=[[]for _ in A]
 for mask,c in cycles.items():
  for v in range(n):
   if mask>>v&1:byvertex[v].append((mask,c))
 @lru_cache(None)
 def collect(S):
  if not S:return (1,)
  b=S&-S;v=b.bit_length()-1;p=[0]+list(collect(S^b))
  for C,c in byvertex[v]:
   if C&S==C:p=plus(p,[0]+scale(collect(S^C),2*c))
  return tuple(trim(p))
 return collect((1<<n)-1)

def basis(k,odd):
 p=[0,1]if odd else[0,0,1]
 for j in range(1,k+1):p=times(p,[-j*j,0,1])
 return p
def weight(k,odd):return math.factorial(k)*math.factorial(2*k+1)*(1 if odd else k+1)
def central(p,odd):
 p=list(p);d=(len(p)-1-(1 if odd else 2))//2+1;v=[0]*max(0,d)
 for k in range(d-1,-1,-1):
  b=basis(k,odd);v[k]=p[len(b)-1]if len(p)>=len(b)else 0;p=plus(p,scale(b,-v[k]))
 ck(trim(p)==[0],'central remainder');return v
def wvec(p,L):
 values=[ev(p,j)for j in range(L)];out=[]
 for r in range(L):
  if r:out.append(math.factorial(r-1)*values[0])
  values=[values[j+1]-values[j]for j in range(len(values)-1)]
 return out
def cutoff(M,central_odd=None):
 k=0 if central_odd is not None else 1
 while (weight(k,central_odd)if central_odd is not None else math.factorial(k-1)*math.factorial(k))%M:k+=1
 return k

def contextmatrix(Q,r,L):
 images={s:prof(subst(Q,r,tr(s)))for s in range(1,L)};A=[[0]*(L-1)for _ in range(L-1)];fall=[1]
 for k in range(1,L):
  fall=times(fall,[-(k-1),1]);p=[0]
  for s,c in enumerate(fall):
   if c:p=plus(p,scale(images[s],c))
  w=wvec(p,L);den=math.factorial(k-1)*math.factorial(k)
  for j,z in enumerate(w):ck(z%den==0,'weighted Newton column integrality');A[j][k-1]=z//den
 return A

def wrapmatrix(M,odd,c):
 K=cutoff(M,odd);outodd=odd if c=='D'else not odd;J=cutoff(M,outodd);A=[[0]*K for _ in range(J)]
 for k in range(K):
  p=[0]
  for s,v in enumerate(basis(k,odd)):
   if v:p=plus(p,scale(prof(cy(tr(s),tr(1 if c=='D'else 2),tr(1))),v))
  u=central(p,outodd)
  for j in range(J):
   num=(u[j]if j<len(u)else 0)*weight(j,outodd);den=weight(k,odd);ck(num%den==0,'physical central column integrality');A[j][k]=num//den%M
 return A

def finite(M):
 K=cutoff(M,True);ck(K==cutoff(M,False),'equal dimensions');mat={(o,c):wrapmatrix(M,o,c)for o in (False,True)for c in'DE'}
 states=[(o,v)for o in(False,True)for v in product(range(M),repeat=K)]
 def step(s,c):
  o,v=s;return(o if c=='D'else not o,tuple(sum(a*b for a,b in zip(row,v))%M for row in mat[o,c]))
 start=(True,(1,)+(0,)*(K-1));dist={start:0};todo=deque([start])
 while todo:
  s=todo.popleft()
  for c in'DE':
   z=step(s,c)
   if z not in dist:dist[z]=dist[s]+1;todo.append(z)
 rev={s:[]for s in states}
 for s in states:
  for c in'DE':rev[step(s,c)].append(s)
 back={start};todo=deque([start])
 while todo:
  for s in rev[todo.popleft()]:
   if s not in back:back.add(s);todo.append(s)
 ck(len(dist)==(48 if M==5 else 686),'all physical orbit states')
 ck(max(dist.values())==(7 if M==5 else 14),'maximum shortest word')
 ck(back=={s for s in states if any(s[1])},'full nonzero communicating class')
 if M==5:
  for c in'DE':ck(len({step(s,c)for s in states})==len(states),'both letters permutations')
 J,V=(4,10)if M==5 else(5,13);words=[]
 for l in range(1,J+1):
  for w in product('DE',repeat=l):
   if sum(2 if c=='D'else 3 for c in w)<=V:words.append(w)
 steering=0
 for s in states:
  if not any(s[1]):continue
  got=set()
  for w in words:
   z=s
   for c in w:z=step(z,c)
   got.add(z[1][0])
  ck(got==set(range(M)),'joint positive word/vertex steering');steering+=M
 fs=[]
 for odd in(False,True):
  found=set()
  for v in product(range(M),repeat=K):
   p=[0]
   for k,a in enumerate(v):p=plus(p,scale(basis(k,odd),a*pow(weight(k,odd),-1,M)))
   found.add(tuple(ev(p,x)%M for x in range(M)))
  ck(len(found)==M**K,'same parity injectivity');fs.append(found)
 ck(fs[0]&fs[1]=={(0,)*M},'cross-parity overlap exactly zero')
 return {'modulus':M,'matrices':{str(o)+c:A for(o,c),A in mat.items()},'orbit':len(dist),'maximum_shortest_word':max(dist.values()),'communicating_nonzero':len(back),'steering_targets':steering,'global_classes':len(fs[0]|fs[1])}

def main():
 moduli=[2,3,4,5,7,9,16,25,49];L=max(cutoff(M)for M in moduli);inputs=[tour(n,c)for n in range(1,5)for c in range(1<<(n*(n-1)//2))]
 counts=checks=integercols=0
 for n in range(1,4):
  for code in range(1<<(n*(n-1)//2)):
   Q=tour(n,code)
   for r in range(n):
    A=contextmatrix(Q,r,L);integercols+=L-1
    for T in inputs:
     p=prof(T);q=prof(subst(Q,r,T));counts+=1
     for M in moduli:
      R=cutoff(M);w=wvec(p,R);out=wvec(q,R);pred=[sum(A[i][j]*w[j]for j in range(R-1))%M for i in range(R-1)]
      ck(pred==[z%M for z in out],'actual full modular context transition');checks+=1
 # Polynomial symmetry/function identification and zero representative from literal cycles.
 C5=tuple(sum(1<<((i+j)%5)for j in(1,2))for i in range(5));ck(all(z%5==0 for z in wvec(prof(C5),5)),'literal zero class')
 quotients=[finite(5),finite(7)]
 # Independent integer D matrices recovered from actual substituted profiles through degree8.
 recurrence=[]
 for odd in(False,True):
  A=[[0]*3 for _ in range(3)]
  for k in range(3):
   p=[0]
   for s,c in enumerate(basis(k,odd)):
    if c:p=plus(p,scale(prof(cy(tr(s),tr(1),tr(1))),c))
   u=central(p,odd)
   for j in range(3):
    num=(u[j]if j<len(u)else 0)*weight(j,odd);ck(num%weight(k,odd)==0,'exact Jacobi integer');A[j][k]=num//weight(k,odd)
  def mm(X,Y):return[[sum(X[i][k]*Y[k][j]for k in range(3))for j in range(3)]for i in range(3)]
  I=[[int(i==j)for j in range(3)]for i in range(3)];A2=mm(A,A);A3=mm(A2,A);co=(336,-231,32)if odd else(708,-345,38)
  ck(A3==[[co[0]*I[i][j]+co[1]*A[i][j]+co[2]*A2[i][j]for j in range(3)]for i in range(3)],'characteristic cubic')
  if odd:ck(all((A3[i][j]-4*A2[i][j])%7==0 for i in range(3)for j in range(3)),'odd recurrence')
  else:
   P=I
   for _ in range(16):P=[[z%7 for z in row]for row in mm(A,P)]
   ck(P==I,'even period16')
  recurrence.append(A)
 for p in[3,5,7,11,13]:
  for odd in(False,True):
   ck(cutoff(p,odd)==(p-1)//2 and cutoff(p*p,odd)==(p if odd else p-1),'prime/square exact cutoff')
   for e in range(1,8):ck(cutoff(p**e,odd)<=p*((e+2)//3),'prime power finite bound')
 print(json.dumps({'status':'PASS','literal_context_profiles':counts,'modular_transitions':checks,'integer_columns':integercols,'cached_actual_graphs':prof.cache_info().currsize,'full_quotients':quotients,'exact_jacobi_matrices':recurrence,'zero_C5_profile':prof(C5)},indent=2,sort_keys=True))
if __name__=='__main__':main()
