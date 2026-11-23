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

DIGITS={
 1:((2,12,25,17,3),16,2,(20,11)),
 3:((2,12,25,17,35,15),34,0,(80,13)),
 5:((2,4,41,51,39,3),32,1,(40,15)),
 7:((2,4,9,19,39,15),32,2,(20,17)),
 9:((2,12,41,17,39,11),12,0,(80,19)),
}
SEAMS=[
 ((720,733,153,129,776,31,188,304,623,236),(2,6),(6173,6175)),
 ((148,181,176,391,136,89,415,32,183),(0,5),(493,495)),
 ((194,180,9,19,37,77,158,60),(1,6),(495,497)),
]
def build(A,I,X):
 n=len(A);m=len(X); b=n+m
 return tuple(r|((I>>i&1)<<b) for i,r in enumerate(A))+tuple((r<<n)|((1<<n)-1) for r in X)+(((1<<n)-1)^I|(((1<<m)-1)<<n),)
def t(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
def exact_inverse(N):
 h,d=divmod(N,10)
 x=h-1
 if d in DIGITS and x>=16 and x.bit_count()==1:return x.bit_length()-4,d
def decode_smooth(N):
 m=N.bit_length()-7
 if m<1:return
 R=N-80*(1<<m)
 if R<=0:return
 for d in DIGITS:
  D,r=divmod(R,10+d)
  if r:continue
  n=D;ex=[]
  for k in ATOMS:
   a=0
   while n%k==0:a+=1;n//=k
   ex.append(a)
  e=m-sum(SHIFTS[k]*x for k,x in zip(ATOMS,ex))
  if n==1 and e>=0 and any(ex+[e]):yield d,m,D,tuple(ex+[e])
if __name__=='__main__':
 for d,(A,I,s,w) in DIGITS.items():
  assert shellRead(A,I)==w and strong(shell(A,I))
  assert (w[0]<<s,w[1])==(80,10+d)
  for m in (1,2,3):
   N=80*2**m+10+d
   assert read(build(A,I,t(m+s)))[0]==N
   assert exact_inverse(N)==(m,d)
   assert list(decode_smooth(N))==[(d,m,1,(0,0,0,0,m))]
  print('digit',d,'coefficient',w,'pad',s)
 for d,A in ATOMS.items():assert read(A)==(d,1<<SHIFTS[d])
 for A,e,(a,b) in SEAMS:
  B=flip(A,*e)
  assert strong(A) and strong(B) and (read(A)[0],read(B)[0])==(a,b)
  print('literal strong seam',len(A),a,b,'edge',e)
