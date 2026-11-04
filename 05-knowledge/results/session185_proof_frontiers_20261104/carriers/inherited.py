"""Copied DP helpers from Session184 spectrum/verify.py; inherited, not novel."""
import itertools, json, hashlib
from pathlib import Path

def need(b,m):
 if not b: raise RuntimeError(m)
def signature(word):
 out=None
 for kind,a in word:
  z=(1,2,a,a,0,a,a) if kind==1 else (3,12,3*a,3*a,a*(a-1),a,a)
  if out is None:out=z;continue
  h,p,f,g,j,e,s=out;H,P,F,G,J,E,S=z
  out=(h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H)
 return out

def words_with_F(F):
 """All strong-attachment component words at exact positive F, all orders."""
 k=0
 while F%3**k==0:
  u=F//3**k
  def rec(pos,left,w):
   if pos>=u.bit_length():
    if left==0 and w[-1][1]>0:yield tuple(w)
    return
   yield from rec(pos+1,left,w+[(1,(u>>pos)&1)])
   if left:yield from rec(pos+2,left-1,w+[(3,(u>>pos)&3)])
  yield from rec(0,k,[])
  k+=1

def normalized_candidates(F):
 while F%2==0:F//=2
 for G in (F,2*F):
  for w in words_with_F(G):
   if w[0][1]:yield w

def graph(word):
 rows=[];ports=set()
 for kind,a in word:
  n=len(rows);block=[0] if kind==1 else [2,4,1]
  rows=[x|(((1<<kind)-1)<<n) for x in rows]+[x<<n for x in block]
  ports.update(range(n,n+a))
 return rows,ports

def tables(rows):
 n=len(rows);N=1<<n;h=[0]*N;h[0]=1;e=[[0]*n for _ in range(N)]
 for v in range(n):e[1<<v][v]=1
 for S in range(1,N):
  for v in range(n):
   if S>>v&1:
    T=S^(1<<v)
    if T:e[S][v]=sum(e[T][u] for u in range(n) if rows[u]>>v&1)
  h[S]=sum(e[S])
 return h,e

def direct_signature(rows,ports):
 n=len(rows);h,e=tables(rows);reverse=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)];_,s=tables(reverse);N=(1<<n)-1
 ep=[sum(v[i] for i in ports) for v in e];sp=[sum(v[i] for i in ports) for v in s]
 return (h[N],sum(h[x]*h[N^x] for x in range(N+1)),sum(ep[x]*h[N^x] for x in range(N+1)),sum(h[x]*sp[N^x] for x in range(N+1)),sum(ep[x]*sp[N^x] for x in range(N+1)),ep[N],sp[N])

def attachment(word,t):
 rows,W=graph(word);b=len(rows);pref=[sum(1<<j for j in range(i+1,t))|(((1<<b)-1)<<t) for i in range(t)]
 rows=pref+[v<<t for v in rows];root=t+b;rows.append((1<<t)-1)
 for i in range(b):
  if i in W:rows[t+i]|=1<<root
  else:rows[root]|=1<<(t+i)
 return rows

def strong(rows):
 N=(1<<len(rows))-1
 for v in range(len(rows)):
  S=1<<v
  while True:
   T=S
   for i in range(len(rows)):
    if S>>i&1:T|=rows[i]
   if S==T:break
   S=T
  if S!=N:return False
 return True
