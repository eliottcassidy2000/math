from math import isqrt, comb
ATOMS={19:(38,4,8,19,39,14),217:(70,60,24,17,33,77,158,63),
       223:(230,28,8,17,37,78,158,62),599:(178,100,137,19,38,76,157,58)}
SHIFTS={19:8,217:12,223:12,599:13}
SHELLS=[
 ((8,1,3,6,15,31,63),3,(128,1)),
 ((2,132,9,35,15,23,63,125),38,(896,7)),
 ((32,5,81,7,75,30,171,63),100,(1024,-1)),
 ((8,1,3,166,15,23,63,119),80,(224,-7)),
]
def dp(A):
 n=len(A); B=[[0]*n for _ in range(1<<n)]
 for i in range(n):B[1<<i][i]=1
 for s in range(1,len(B)):
  for i,c in enumerate(B[s]):
   if not c:continue
   m=A[i]&~s
   while m:
    j=m&-m;m-=j;B[s|j][j.bit_length()-1]+=c
 return B
def read(A):
 H=[1]+list(map(sum,dp(A)[1:]));F=len(H)-1
 return H[F],sum(H[s]*H[F^s] for s in range(F+1))
def flip(A,i,j):
 A=list(A);A[i]^=1<<j;A[j]^=1<<i;return tuple(A)
def shellRead(A,I):
 n=len(A);F=(1<<n)-1;V=dp(A);W=dp(tuple(r^F^(1<<i) for i,r in enumerate(A)))
 G=Q=0
 for s in range(1,F+1):
  d=sum(x for i,x in enumerate(V[s]) if I>>i&1);t=F^s
  G+=d*(sum(W[t]) if t else 1)
  Q+=d*(sum(x for i,x in enumerate(W[t]) if not I>>i&1) if t else 1)
 return G,sum(V[F])+Q-G
def shell(A,I):
 n=len(A)
 return tuple(r|(((I>>i)&1)<<(n+1)) for i,r in enumerate(A))+((1<<n)-1,(((1<<n)-1)^I)|(1<<n))
def strong(A):
 F=(1<<len(A))-1
 for i in range(len(A)):
  r=1<<i
  while True:
   t=r
   for j,v in enumerate(A):
    if r>>j&1:t|=v
   if r==t:break
   r=t
  if r!=F:return False
 return True
def prime(n):return n>1 and all(n%d for d in range(2,isqrt(n)+1))
def inverse(p):
 for sign,c in ((1,7),(-1,10)):
  K=p.bit_length()-int(sign==1);D=sign*(p-(1<<K))
  if D<1:continue
  A=[];U=D
  for d in ATOMS:
   a=0
   while U%d==0:a+=1;U//=d
   A.append(a)
  e=K-c-sum(SHIFTS[d]*a for d,a in zip(ATOMS,A))
  if U==1 and e>=0 and any(A+[e]):yield sign,K,D,tuple(A+[e])
def evaluate(F,t):return sum(v*t**i for i,v in enumerate(F))
def wrap(F):
 G=[0]*(len(F)+2)
 for j,a in enumerate(F):
  G[j+2]+=a
  for i in range((j-1)%2,j,2):G[i+1]+=2*a*comb(j,i)
 return tuple(G)
def inverse_pw(F,target):
 h,q=evaluate(F,1),evaluate(F,2)
 if h==target:return ()
 if h>target:return None
 a=0
 while h+q*4**a<=target:
  seq=inverse_pw(wrap((0,)*(2*a)+tuple(F)),target)
  if seq is not None:return (a,)+seq
  a+=1
if __name__=='__main__':
 for A,I,x in SHELLS:
  assert shellRead(A,I)==x and strong(shell(A,I))
 for d,A in ATOMS.items():assert read(A)==(d,1<<SHIFTS[d])
 for p in (1073741843,524269,1048793,1048799,134218327,268434857,549755809651,268439579):
  opts=list(inverse(p));assert len(opts)==1 and prime(p)
  print('exact prime cofactor',p,'sevenfold',7*p,'inverse',opts[0])
 M=(244,81,210,39,8,150,56,90)
 assert read(M)==(383,5960) and read(flip(M,0,5))==(389,5960)
 values=tuple(128*(5960<<3)+d for d in (383,389));assert values==(6103423,6103429)
 assert all(map(prime,values))
 print('literal prime bridge',values,tuple(7*p for p in values))
 f=(0,0,1)
 for a in (2,1,3):f=wrap((0,)*(2*a)+f)
 assert inverse_pw((0,0,1),evaluate(f,1)) is not None
 print('exact physical-word inverse control PASS')
