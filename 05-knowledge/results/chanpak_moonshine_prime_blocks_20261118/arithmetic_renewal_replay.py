from itertools import permutations
from collections import Counter

def require(b,msg):
 if not b: raise RuntimeError(msg)
def rows(n,c):
 r=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if (c>>k)&1:r[i]|=1<<j
   else:r[j]|=1<<i
   k+=1
 return r
def canon(r):
 n=len(r);out=[]
 for p in permutations(range(n)):
  c=0;k=0
  for i in range(n):
   for j in range(i+1,n):
    c|=((r[p[i]]>>p[j])&1)<<k;k+=1
  out.append(c)
 return min(out)
def converse(r):
 n=len(r);return [((1<<n)-1)^r[i]^(1<<i) for i in range(n)]
def count(r):
 n=len(r)
 if not n:return 1
 d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m in range(1,1<<n):
  for i in range(n):
   if not d[m][i]:continue
   a=r[i]&~m
   while a:
    q=a&-a;j=q.bit_length()-1;a-=q;d[m|q][j]+=d[m][i]
 return sum(d[-1])
def strong(r):
 n=len(r)
 for s in range(n):
  reach=1<<s
  while True:
   old=reach
   for i in range(n):
    if reach>>i&1:reach|=r[i]
   if reach==old:break
  if reach!=(1<<n)-1:return False
 return n>0
B=5
A=[Counter() for _ in range(B+1)];S=[Counter() for _ in range(B+1)];C=[Counter() for _ in range(B+1)];U=[Counter() for _ in range(B+1)]
A[0][1]=C[0][1]=1
for n in range(1,B+1):
 seen=set()
 for c in range(1<<(n*(n-1)//2)):
  r=rows(n,c);v=canon(r)
  if v in seen:continue
  seen.add(v);h=count(r);sc=canon(converse(r))==v
  A[n][h]+=1;C[n][h]+=sc
  if strong(r):
   S[n][h]+=1;U[n][h]+=sc
   require(n==1 or h>=n,'Camion bound')
for n in range(1,B+1):
 a=Counter();c=Counter()
 for k in range(1,n+1):
  for h,x in S[k].items():
   for g,y in A[n-k].items():a[h*g]+=x*y
 require(a==A[n],('Dirichlet convolution',n,a,A[n]))
 if n%2==0:
  for h,x in A[n//2].items():c[h*h]+=x
 for k in range(1,n+1):
  if (n-k)%2:continue
  for h,x in U[k].items():
   for g,y in A[(n-k)//2].items():c[h*g*g]+=x*y
 require(c==C[n],('Converse square divisor law',n,c,C[n]))
 for h in [3,5,7,11,13,21,35]:
  expected=sum((n-k+1)*S[k][h] for k in range(1,n+1))
  require(A[n][h]==expected,('prime-seven lane',n,h))
print('PASS literal isomorphism/strong/converse census through order5')
print('class totals',[sum(a.values()) for a in A])
print('strong inventories', [dict(s) for s in S])
print('PASS Dirichlet renewal, converse square-divisor rule, prime/sevenfold exact linear kernels')
