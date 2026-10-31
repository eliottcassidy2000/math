#!/usr/bin/env python3
import itertools,json
import sympy as S
D=36
idx=lambda o,p,q,r:18*o+6*p+2*q+r

def matrix(c):
 M=S.zeros(D); bits=[c&1,(c>>1)&1]
 for o in range(2):
  a,b=(0,1) if o==0 else (1,0)
  for p,q,r in itertools.product(range(3),range(3),range(2)):
   j=idx(o,p,q,r)
   M[idx(o,2 if bits[a] else 1,q,r),j]+=1
   if q or bits[a]: M[idx(o,p,2 if bits[b] else 1,r),j]+=1
   if r or bits[b]: M[idx(o,p,q,1),j]+=1
 return M
MS=[matrix(c) for c in range(4)]
INIT=S.zeros(D,1)
for o in range(2):INIT[idx(o,0,0,0)]=1
OUT=S.zeros(1,D)
for o in range(2):
 for p,q,r in itertools.product((0,1),((0,1) if o==0 else (1,)),range(2)): OUT[0,idx(o,p,q,r)]=1

def state(w):
 v=INIT
 for c in w:v=MS[c]*v
 return v

def direct(w):
 n=len(w)+2; adj=[0]*n;adj[0]|=2
 for i,c in enumerate(w,2):
  for a in range(2):
   if c>>a&1:adj[a]|=1<<i
   else:adj[i]|=1<<a
  for j in range(2,i):adj[j]|=1<<i
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if not dp[mask][i]:continue
   nxt=adj[i]&~mask
   while nxt:
    bit=nxt&-nxt;nxt-=bit;j=bit.bit_length()-1
    dp[mask|bit][j]+=dp[mask][i]
 return sum(dp[-1])

def main():
 sparse=[[(i,j,int(M[i,j])) for i in range(D) for j in range(D) if M[i,j]] for M in MS]
 frontier=[((),list(map(int,INIT)))];tested=0
 for m in range(7):
  nxt=[]
  for w,v in frontier:
   h=sum(int(OUT[0,i])*v[i] for i in range(D))
   hd=direct(w)
   if h!=hd:raise RuntimeError((w,h,hd))
   tested+=1
   if m<6:
    for c in range(4):
     nv=[0]*D
     for i,j,z in sparse[c]:nv[i]+=z*v[j]
     nxt.append((w+(c,),nv))
  frontier=nxt
 print(json.dumps(dict(exhaustive_words_through_length=6,direct_checks=tested)),flush=True)
if __name__=='__main__':main()
