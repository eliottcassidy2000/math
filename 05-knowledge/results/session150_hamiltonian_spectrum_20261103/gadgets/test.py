#!/usr/bin/env python3
"""Exhaustive small-tournament checks of rooted one-sided multi-insertion."""

def endpoints(a):
 n=len(a)
 if not n:return 1,[]
 d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if d[mask][i]:
    for j in range(n):
     if not(mask>>j&1) and a[i][j]:d[mask|1<<j][j]+=d[mask][i]
 return sum(d[-1]),d[-1]
def induced(a,S):return [[a[i][j] for j in S] for i in S]
def reverse(a):return [list(row) for row in zip(*a)]
def add(a,v):
 n=len(a);return [row+[int(i==v)] for i,row in enumerate(a)]+[[int(i!=v) for i in range(n)]+[0]]
def graph(n,bits):
 a=[[0]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):
   a[i][j]=(bits>>k)&1;a[j][i]=1-a[i][j];k+=1
 return a
def profile(a,v):
 n=len(a);u=[i for i in range(n) if i!=v]; h,_=endpoints(a);_,starts=endpoints(reverse(a));s=starts[v];hu,_=endpoints(induced(a,u));f=0
 for mask in range(1<<len(u)):
  left=[x for i,x in enumerate(u) if mask>>i&1];right=[x for i,x in enumerate(u) if not(mask>>i&1)]
  if not left: l=1
  else:
   _,ep=endpoints(induced(a,left));l=sum(ep[i] for i,x in enumerate(left) if a[x][v])
  hr,_=endpoints(induced(a,right));f+=l*hr
 return h,s,hu,f

def run():
 total=0; witnessed=None; rays=set()
 for n in range(2,6):
  for code in range(1<<(n*(n-1)//2)):
   a=graph(n,code)
   for v in range(n):
    h,s,u,f=profile(a,v); b=a
    for k in range(1,4):
     b=add(b,v);actual,_=endpoints(b);_,starts=endpoints(reverse(b));pred=h-s+u+(2**k-1)*(f-u)
     if (actual,starts[v])!=(pred,u): raise AssertionError((n,code,v,k,(actual,starts[v]),(pred,u),(h,s,u,f)))
     total+=1
    if n==2 and code==0 and v==0:witnessed=(h,s,u,f)
    rays.add((h-s+2*u-f, f-u))
 print('rooted gadget checks:',total)
 print('counterexample n=2, bitmask=0, v=0, (H,s,H(T-v),F):',witnessed)
 print('distinct affine rays through order 5:',len(rays))
if __name__=='__main__':run()
