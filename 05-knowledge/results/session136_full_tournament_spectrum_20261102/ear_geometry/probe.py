from itertools import permutations
import json

def adj(n,code):
 a=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
   k+=1
 return a

def h_all(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)];hs=[0]*(1<<n);hs[0]=1
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for j in range(n):
   z=dp[s][j]
   if z:
    bs=a[j]&~s
    while bs:
     bit=bs&-bs;bs-=bit;dp[s|bit][bit.bit_length()-1]+=z
  hs[s]=sum(dp[s])
 return hs

def kernel(a,hs):
 n=len(a);V=(1<<n)-1;r=[[0]*n for _ in range(n)]
 def rec(start,end,mask):
  if mask.bit_count()%2==0:r[start][end]+=hs[V^mask]
  bs=a[end]&~mask
  while bs:
   bit=bs&-bs;bs-=bit;rec(start,bit.bit_length()-1,mask|bit)
 for i in range(n):rec(i,i,1<<i)
 return r

def cut(r,s):
 n=len(r)
 return sum(r[i][j] for i in range(n) if s>>i&1 for j in range(n) if not(s>>j&1))

def strong(a):
 n=len(a)
 for i in range(n):
  reach=1<<i
  while True:
   z=reach
   for j in range(n):
    if reach>>j&1:z|=a[j]
   if z==reach:break
   reach=z
  if reach!=(1<<n)-1:return False
 return True

def run(n):
 count=0;miss=[];misscount=0;maxmult=(0,None);minsizes={};hist={}
 for code in range(1<<(n*(n-1)//2)):
  a=adj(n,code)
  if not strong(a):continue
  count+=1;hs=h_all(a);r=kernel(a,hs);vals=[cut(r,s) for s in range(1,2**n-1)];mn=min(vals);masks=[s for s in range(1,2**n-1) if vals[s-1]==mn]
  if len(masks)>maxmult[0]:maxmult=(len(masks),code)
  if not any(s.bit_count() in (1,n-1) for s in masks):
   misscount+=1
   if len(miss)<3:miss.append((code,hs[-1],mn,masks,r))
  sizes=tuple(sorted(set(min(s.bit_count(),n-s.bit_count()) for s in masks)))
  minsizes[sizes]=minsizes.get(sizes,0)+1
  hist[len(masks)]=hist.get(len(masks),0)+1
 print(json.dumps({'n':n,'strong':count,'no_singleton_min_count':misscount,'no_singleton_min':miss,'max_multiplicity':maxmult,'size_sets':str(minsizes),'multiplicity_histogram':hist}),flush=True)

if __name__=='__main__':
 for n in range(3,7):run(n)
