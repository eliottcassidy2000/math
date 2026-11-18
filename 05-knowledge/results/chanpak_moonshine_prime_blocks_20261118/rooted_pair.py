from sys import path
from pathlib import Path
path.insert(0,str(Path(__file__).resolve().parent.parent/'moonshine_prime_blocks'))
from replay import dp,transpose
from itertools import combinations
from functools import lru_cache
from collections import deque
from math import isqrt
A=(8,5,9,18,7,31,63,127);r=4
C=(8,17,11,50,69,23,47,127);t=5

def rooted(T,b):
   E=dp(T); S=dp(transpose(T)); f=len(E)-1; H=[1]+[sum(Z) for Z in E[1:]]
   R=sum(E[s][b]*H[f^s] for s in range(1,f+1))
   W=sum(S[s][b]*H[f^s] for s in range(1,f+1))
   M=f^(1<<b);a=M;v=0
   while 1:
      v+=H[a]*H[a^M]
      if not a:break
      a=(a-1)&M
   return H[f],E[f][b],W,S[f][b],R,v
def append(A,a,B,b):
   m=len(A);n=len(B);high=((1<<n)-1)<<m
   T=[r|high for r in A]+[r<<m for r in B]
   T[a]^=1<<(m+b)
   T[m+b]|=1<<a
   return tuple(T),m+b
def transform(x,w,n):
   h,e,a=x;g,f,s,S,R,V=w
   return h*g-e*S+a*R,h*f-(e if n==1 else 0),h*s+a*V
def TT(n):
   return tuple((1<<n)-(1<<(j+1)) for j in range(n))
def decode(i,n):
   B=[0]*n
   for j,(a,b) in enumerate(combinations(range(n),2)):
      if i>>j&1:B[a]|=1<<b
      else:B[b]|=1<<a
   return tuple(B)

def run():
   assert rooted(A,r)[:3]==(11,1,27) and rooted(C,t)[:3]==(77,7,189)
   for plan in [((TT(3),2),(TT(1),0)),((TT(2),1),((2,4,1),2),(TT(2),1)),((TT(3),0),(TT(2),1),(TT(2),0))]:
      P,R=A,r;Q,S=C,t;state=rooted(P,R)[:3]
      for B,k in plan:
         state=transform(state,rooted(B,k),len(B))
         P,R=append(P,R,B,k);Q,S=append(Q,S,B,k)
         assert rooted(P,R)[:3]==state
         assert rooted(Q,S)[:3]==tuple(7*j for j in state)
      print('literal layers', [(len(B),k) for B,k in plan],state[0],7*state[0])

   # Choose the prime first. Exact search uses a finite, stated palette.
   P=sorted({(T,j) for n in range(1,5) for i in range(1<<(n*(n-1)//2)) for T in [decode(i,n)] for j in range(n)})
   U={(*rooted(T,j),len(T)):(T,j) for T,j in P}
   rr=(11,1,27)
   q=deque([rr]);old={rr:None}
   targets={30011,30013,10007,10009,100003,100019};M=max(targets)
   byH={}
   while q:
      w=q.popleft(); h=w[0]
      if h in targets and h not in byH:byH[h]=w
      if len(byH)==len(targets):break
      for (g,e,s,S,R,V,n),B in U.items():
         z=transform(w,(g,e,s,S,R,V),n)
         if z in old or z[0]>M or (z[2]>M and z[0] not in targets):continue
         old[z]=w,B
         if z[2]<=M:q.append(z)
         if z[0] in targets:byH.setdefault(z[0],z)
   assert set(byH) >= {100003,100019}
   for h in sorted(byH):
      plan=[];w=byH[h]
      while old[w]:
         w,B=old[w];plan.append(B)
      assert all(h%d for d in range(2,isqrt(h)+1))
      seq=list(reversed(plan))
      P,R=A,r;Q,S=C,t
      z=rooted(P,R)[:3];w=rooted(Q,S)[:3]
      for B,k in seq:
         z=transform(z,rooted(B,k),len(B))
         w=transform(w,rooted(B,k),len(B))
         P,R=append(P,R,B,k);Q,S=append(Q,S,B,k)
      assert (z[0],w[0])==(h,7*h)
      print(h,7*h,seq,'orders',len(P),len(Q))
def mixed(fixtures=False,other=False):
    seeds=((4,9,18,21,35,15,63,127),(4,9,66,53,39,135,59,95))
    path=((TT(3),2,True),((2,4,1),0,False),((0,),0,False),((12,9,2,4),1,False))
    W=((29,11,31,13),(153,29,153,29),(889,29,1749,87),
       (2609,87,3411,87),(30013,87,108543,435))
    if other:
        path=(((0,),0,False),(TT(3),2,False),(TT(2),1,False),(TT(2),1,False),((2,4,1),0,False))
        W=((29,11,31,13),(49,13,49,13),(245,13,379,13),(1003,13,1503,13),
           (4009,13,5999,13),(30011,13,71949,39))
    for scale,A in zip((1,7),seeds):
        root=0;x=(*rooted(A,root)[:3],rooted(tuple(r>>1 for r in A[1:]),0)[0])
        assert x==tuple(scale*v for v in W[0])
        if fixtures:print(len(A),root,*A)
        for ii,(B,b,migrate) in enumerate(path,1):
            h,e,t,d=x; g,f,S1,S,R,V=rooted(B,b)
            D=dp(B);M=((1<<len(B))-1)^(1<<b)
            g0=sum(D[M]) if M else 1
            ab=[1]+[sum(y) for y in D[1:]]
            BB=sum(ab[k]*ab[len(ab)-1-k] for k in range(len(ab)))
            T,w=append(A,root,B,b)
            hh=g*h-S*e+R*t
            if migrate:x=(hh,f*h-(e if len(B)==1 else 0),S1*h+V*t,h*g0);root=w
            else:x=(hh,d*f,t*BB-d*S1,d*g)
            A=T
            assert x==tuple(scale*v for v in W[ii])
            if fixtures:print(len(A),root,*A)
    if not fixtures:
        assert all(W[-1][0]%d for d in range(2,isqrt(W[-1][0])+1))
        print('mixed keep/migrate',W[-1][0],7*W[-1][0])

if __name__=='__main__':
    import sys
    if '--mixed-fixtures' in sys.argv:mixed(True)
    elif '--fixtures30011' in sys.argv:mixed(True,True)
    else:
        mixed()
        run()
