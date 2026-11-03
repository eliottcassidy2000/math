def require(x):
 if not x: raise AssertionError("audit failed")
from collections import deque
Q=(11,4,18,7);B=(5,2,8,3);z=(25,34);I=(1,0,0,1)
def mul(A,B,m):return tuple(sum(A[2*i+k]*B[2*k+j] for k in range(2))%m for i in range(2) for j in range(2))
def power(A,k,m):
 R=I
 while k:
  if k&1:R=mul(R,A,m)
  A=mul(A,A,m);k//=2
 return R
def act(A,v,m):return ((A[0]*v[0]+A[1]*v[1])%m,(A[2]*v[0]+A[3]*v[1])%m)
def inv(A,m):
 a,b,c,d=A; delta=pow((a*d-b*c)%m,-1,m);return tuple(x*delta%m for x in (d,-b,-c,a))
P=7*13*19;T=19152
require(power(Q,T,P)==I)
m=31; A=power(Q,T,m); C=mul(mul(B,A,m),inv(B,m),m)
print('base',P,'Q^T identity; A,C mod31',A,C)
start=tuple(x%m for x in z); vis={start}; todo=deque([start])
while todo:
 v=todo.popleft()
 for g in (A,C):
  w=act(g,v,m)
  if w not in vis:vis.add(w);todo.append(w)
require(vis=={(a,b) for a in range(m) for b in range(m) if a or b})
print('31-kernel orbit',len(vis),'full set of nonzero vectors')
for p in (7,13,19,31):
 m=p*p; L=power(B,576,m); M=mul(mul(Q,L,m),inv(Q,m),m)
 require(tuple(x%p for x in L)==tuple(x%p for x in M)==I)
 v=[tuple(((a-b)%m)//p for a,b in zip(act(R,z,m),z)) for R in (L,M)]
 d=(v[0][0]*v[1][1]-v[0][1]*v[1][0])%p;require(d)
 print('p',p,'vectors',v,'det',d)
