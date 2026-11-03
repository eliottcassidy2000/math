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

def fmod(n,m,s):
 if m==1:return 0
 coeff=80 if s==1 else 192
 c=67 if s==1 else 53
 return ((coeff*pow(8,n,49*m)+28*n+c)%(49*m))//49

def construct(m,r):
 need(gcd(m,7)==1 and not(m%2==0 and r%2==0),('input',m,r))
 s=2 if m%4==0 and r%4==1 else 1
 fs=factor(m);L=1
 for p in fs:
  if p==2:continue
  k=1
  while pow(8,k,p)!=1:k+=1
  L=lcm(L,k)
 ps=factor(L)
 for p,e in fs.items():ps[p]=max(ps.get(p,0),max(0,e-2) if p==2 else e)
 R,Q=0,1
 for p,a in sorted(ps.items()):
  if not a:continue
  e=fs.get(p,0);precision=max(0,e-2) if p==2 else e
  W=p**precision
  if precision:
   choices=[]
   for x in range(W):
    n=crt(R,Q,x,W)+Q*W*(m.bit_length()+2)
    if fmod(n,p**e,s)==r%(p**e):choices.append(x)
   need(choices,('local',m,r,s,p,R,Q))
   x=choices[0]
  else:x=0
  Z=p**a;R=crt(R,Q,x,Z);Q*=Z
 n=R+Q*(m.bit_length()+3)
 need(fmod(n,m,s)==r%m,('output',m,r,s,n,fmod(n,m,s)))
 return n,s

def rows(n):
 w=(1,)+(0,0,1)*n+(0,0);N=len(w)+1;A=[0]*N
 for i in range(1,N):
  for j in range(i+1,N):A[i]|=1<<j
 for i,b in enumerate(w,1):
  if b:A[0]|=1<<i
  else:A[i]|=1
 return A

def dp(A):
 N=len(A);d=[[0]*N for _ in range(1<<N)]
 for i in range(N):d[1<<i][i]=1
 for S in range(1,1<<N):
  for i in range(N):
   if d[S][i]:
    for j in range(N):
     if not(S>>j&1) and A[i]>>j&1:d[S|1<<j][j]+=d[S][i]
 return sum(d[-1])
for n in range(1,4):
 h=(192*8**n+28*n+53)//49
 need(dp(rows(n))==h,('companion DP',n,h))
print('Companion independent graph DP n=1..3')
for m in list(range(1,251))+[256,625,675,1001,1024,1575]:
 if m%7==0:continue
 for r in range(m):
  if m%2==0 and r%2==0:continue
  construct(m,r)
print('All eligible classes verified by constructive CRT for M<=250, 7 not dividing M, plus extra moduli')
