from math import gcd,lcm

def need(c,w):
 if not c:raise AssertionError(w)

def factor(n):
 d=2;f={}
 while d*d<=n:
  while n%d==0:f[d]=f.get(d,0)+1;n//=d
  d+=1
 if n>1:f[n]=f.get(n,0)+1
 return f

def crt(a,A,b,B):return a if B==1 else a+A*((b-a)*pow(A,-1,B)%B)

def hmod(a,m):return ((10*pow(4,a,9*m)+6*a-1)%(9*m))//9 if m>1 else 0

def Hmod(k,a,n,m):
 if m==1:return 0
 q=1<<k;v=pow(q-1,-1,m)
 b=(pow(q,n+1,m)-1)*v % m
 S=(q*(pow(q,n,m)-1)*v*v-n*v)%m
 C=2*(4**a-1)//3;D=(4*4**a+6*a-13)//9
 return (2+n+(q//2-1)*S+C*b+D)%m

def construct(M,r):
 f=factor(M);e=f.get(2,0);m=M//(2**e)
 need(not(e and r%2==0),('ineligible',M,r))
 if e:
  a=next(i for i in range(1,2**e+e+3) if hmod(i,2**e)==r%(2**e))
 else:a=1
 ps=factor(m);L=1
 for p in ps:
  d=1
  while pow(2,d,p)!=1:d+=1
  L=lcm(L,d)
 k=1+L*max(1,(max(2,e+1)-1+L-1)//L)
 need(k>=max(2,e+1),('k',k))
 need(gcd((1<<k)-1,M)==1,('q-1',M,k))
 coordinates=factor(L)
 for p,v in ps.items():coordinates[p]=max(coordinates.get(p,0),v)
 R,Q=0,1
 for p,d in sorted(coordinates.items()):
  v=ps.get(p,0);W=p**v
  if v:
   choices=[x for x in range(W) if Hmod(k,a,crt(R,Q,x,W)+Q*W,W)==r%W]
   need(choices,('local',M,r,p,R,Q))
   x=choices[0]
  else:x=0
  Z=p**d;R=crt(R,Q,x,Z);Q*=Z
 n=R+Q
 need(Hmod(k,a,n,M)==r%M,('output',M,r,k,a,n,Hmod(k,a,n,M)))
 return k,a,n

def rows(k,a,n):
 w=(1,)+((0,)*(k-1)+(1,))*n+(1,0)*a
 N=len(w)+1;A=[0]*N
 for i in range(1,N):
  for j in range(i+1,N):A[i]|=1<<j
 for i,x in enumerate(w,1):
  if x:A[0]|=1<<i
  else:A[i]|=1
 return A

def strong(A):
 N=len(A)
 for v in range(N):
  R=1<<v
  while True:
   old=R
   for i in range(N):
    if R>>i&1:R|=A[i]
   if R==old:break
  if R!=(1<<N)-1:return False
 return True

def dp(A):
 N=len(A);D=[[0]*N for _ in range(1<<N)]
 for i in range(N):D[1<<i][i]=1
 for S in range(1,1<<N):
  for i in range(N):
   if D[S][i]:
    for j in range(N):
     if not(S>>j&1) and A[i]>>j&1:D[S|1<<j][j]+=D[S][i]
 return sum(D[-1])

for k,a,n in ((2,1,1),(3,1,1),(3,2,1),(2,1,2),(4,1,1),(3,1,2)):
 A=rows(k,a,n);actual=dp(A)
 need(strong(A),('strong',k,a,n))
 need(Hmod(k,a,n,1000003)==actual%1000003,('DP',k,a,n,actual))
print('Literal strong graphs and independent subset DP: six small parameter triples')
for M in list(range(1,251))+[256,343,385,512,625,675,1001,1024,1575,2*3*5*7*11]:
 for r in range(M):
  if M%2==0 and r%2==0:continue
  construct(M,r)
print('Constructive ordered CRT: all eligible residues mod M<=250 plus mixed and prime-power controls')
