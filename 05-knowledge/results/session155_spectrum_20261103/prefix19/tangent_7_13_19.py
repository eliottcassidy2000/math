from math import prod
Q=(11,4,18,7); B=(5,2,8,3);z=(25,34); I=(1,0,0,1)
def mul(A,B,m):return tuple(sum(A[2*i+k]*B[2*k+j] for k in range(2))%m for i in range(2) for j in range(2))
def power(A,k,m):
 R=I
 while k:
  if k&1:R=mul(R,A,m)
  A=mul(A,A,m);k//=2
 return R
def act(A,v,m):return ((A[0]*v[0]+A[1]*v[1])%m,(A[2]*v[0]+A[3]*v[1])%m)
for p in (7,13,19):
 m=p*p; L=power(B,144,m); inv=pow(5,-1,m);Qi=tuple(t*inv%m for t in (7,-4,-18,11));M=mul(mul(Q,L,m),Qi,m)
 if tuple(t%p for t in L)!=I or tuple(t%p for t in M)!=I: raise AssertionError("identity")
 v=[tuple(((x-y)%m)//p for x,y in zip(act(A,z,m),z)) for A in (L,M)]
 det=(v[0][0]*v[1][1]-v[0][1]*v[1][0])%p
 if not det: raise AssertionError("determinant")
 print('p',p,'B_order',next(k for k in range(1,1000) if power(B,k,p)==I),'L',L,'M',M,'vectors',v,'det',det)
