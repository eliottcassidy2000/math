from math import isqrt
from itertools import product

def check(condition):
 if not condition: raise AssertionError("check failed")

def make(bits):
 u=3;k0=0;rev=[(3,1)];ks=[]; H=3
 for b in bits:
  s=isqrt(2*u+k0+1)+1+b
  v=s*s-u;k=v-u-1
  check(k>k0>=0 and u+v==s*s)
  rev.append((v,u));ks.append(k);H*=1+2**k
  u,k0=v,k
 return u,rev,ks,H

def decode(H):
 H//=3;ks=[]
 while H>1:
  x=H-1;k=(x&-x).bit_length()-1
  f=1+2**k
  check(H%f==0)
  H//=f;ks.append(k)
 return ks

def dp(n,rev):
 edges=set(rev);out=[0]*n
 for i in range(1,n+1):
  for j in range(1,n+1):
   if i!=j and ((i,j) in edges or i<j and (j,i) not in edges):out[i-1]|=1<<(j-1)
 d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,1<<n):
  for i in range(n):
   z=d[S][i]
   if not z: continue
   avail=out[i]&~S
   while avail:
    t=avail&-avail;j=t.bit_length()-1;d[S|t][j]+=z;avail-=t
 return sum(d[-1])

for bits in [(),(0,),(1,),(0,0)]:
 n,e,k,H=make(bits); d=dp(n,e); print(bits,n,k,H,d);check(H==d)
for m in range(9):
 for bits in product((0,1),repeat=m):
  n,e,k,H=make(bits)
  check(decode(H)==k)
  check(all(a*a==u+v for v,u in e for a in [isqrt(u+v)]))
print('all good')
