from fractions import Fraction
from itertools import permutations
import json

def coefficients(rows):
 n=len(rows); N=1<<n; dp=[[0]*n for _ in range(N)]; sp=[[0]*n for _ in range(N)]
 for v in range(n):dp[1<<v][v]=sp[1<<v][v]=1
 for S in range(1,N):
  for v in range(n):
   if not S>>v&1:continue
   R=S^(1<<v)
   if not R:continue
   dp[S][v]=sum(dp[R][u] for u in range(n) if rows[u]>>v&1)
   sp[S][v]=sum(sp[R][u] for u in range(n) if rows[v]>>u&1)
 h=[sum(a) for a in dp];h[0]=1
 u=[0]*n;v=[0]*n;d=[0]*n;e=[0]*n;b=[[0]*n for _ in range(n)]
 for S in range(1,N):
  C=(N-1)^S;sg=(-1)**S.bit_count()
  for i in range(n):
   u[i]+=sp[S][i]*h[C];v[i]+=dp[S][i]*h[C]
   d[i]+=sg*sp[S][i]*h[C];e[i]+=sg*dp[S][i]*h[C]
   for j in range(n):b[i][j]+=dp[S][i]*sp[C][j]
 return h[-1],sum(h[S]*h[(N-1)^S] for S in range(N)),u,v,d,e,b,sp[-1],dp[-1]

def packet(rows,x):
 h,Q,u,v,d,e,b,s,t=coefficients(rows);n=len(rows)
 H=sum(x[i]*s[i]+(1-x[i])*t[i] for i in range(n))+sum((1-x[i])*x[j]*b[i][j] for i in range(n) for j in range(n))
 return (H,h,Q,h+sum(x[i]*u[i] for i in range(n)),h+sum((1-x[i])*v[i] for i in range(n)),h+sum(x[i]*d[i] for i in range(n)),h+sum((1-x[i])*e[i] for i in range(n)))

def roots(rows):
 n=len(rows);cc=coefficients(rows); h,Q,u,v,d,e,b,s,t=cc
 for mask in range(1<<n):
  x=[mask>>i&1 for i in range(n)]
  H=sum(x[i]*s[i]+(1-x[i])*t[i] for i in range(n))+sum((1-x[i])*x[j]*b[i][j] for i in range(n) for j in range(n))
  yield mask,(h+sum(x[i]*u[i] for i in range(n)),h+sum((1-x[i])*v[i] for i in range(n)),h+sum(x[i]*d[i] for i in range(n)),h+sum((1-x[i])*e[i] for i in range(n))),H

def allrows(n):
 edges=[(i,j) for i in range(n) for j in range(i+1,n)]
 for code in range(1<<len(edges)):
  rows=[0]*n
  for k,(i,j) in enumerate(edges):
   if code>>k&1:rows[i]|=1<<j
   else:rows[j]|=1<<i
  yield code,rows

if __name__=='__main__':
 for n in range(1,7):
  nf=0;found=[]
  for code,rows in allrows(n):
   fibers={}
   for mask,response,H in roots(rows):
    if response in fibers:
     mask0,H0=fibers[response]
     if H!=H0:
      found.append(dict(n=n,code=code,rows=rows,cuts=[mask0,mask],H=[H0,H],response=response));break
    else:fibers[response]=(mask,H)
   if found:break
  print(json.dumps({'order':n,'first':found}),flush=True)
