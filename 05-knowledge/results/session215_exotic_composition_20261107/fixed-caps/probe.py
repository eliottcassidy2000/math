from math import factorial,comb
from itertools import product
import json

def tournament(n,bits):
 a=[[0]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):
   a[i][j]=(bits>>k)&1;a[j][i]=1-a[i][j];k+=1
 return a

def profile(a):
 n=len(a);N=1<<n
 dp=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for s in range(1,N):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    dp[s][v]=1 if t==0 else sum(dp[t][u] for u in range(n) if t>>u&1 and a[u][v])
  h[s]=sum(dp[s])
 pc=[[0]*(n+1) for _ in range(N)];pc[0][0]=1
 for s in range(1,N):
  bit=s&-s;t=s
  while t:
   if t&bit:
    for k in range(1,n+1):pc[s][k]+=h[t]*pc[s^t][k-1]
   t=(t-1)&s
 return [pc[-1][k]*factorial(k) for k in range(n+1)]

def q(f,j):return sum(comb(j,k)*v for k,v in enumerate(f) if k<=j)
def cap(f,m):
 a=[0]+[factorial(k)*sum((-1)**(k-j)*comb(k,j)*j**m for j in range(k+1))//factorial(k) for k in range(1,m+1)]
 # f(TT_m,k)=k! S(m,k), using Q(t)=t^m
 ans=0
 for x in range(1,m+1):
  for y in range(1,len(f)):
   for z in range(1,m+1):
    if max(x,y,z)-min(x,y,z)<=1:ans+=(3 if x==y==z else 1)*a[x]*f[y]*a[z]
 return ans
if __name__=='__main__':
 for n in range(3,7):
  seen={};hit=None
  for bits in range(1<<(n*(n-1)//2)):
   f=profile(tournament(n,bits));key=(f[1],q(f,2))
   if key in seen and seen[key][1]!=f:
    hit=(seen[key],(bits,f));break
   seen[key]=(bits,f)
  print(n,'collision',hit,flush=True)
  if hit:break
