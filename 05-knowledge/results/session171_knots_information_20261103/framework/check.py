#!/usr/bin/env python3
"""Exact square-zero response and contextual-observability controls."""
import itertools,json,random

def add(a,b):
 c=a.copy()
 for s,v in b.items():
  c[s]=c.get(s,0)+v
  if not c[s]: del c[s]
 return c

def mul(a,b):
 c={}
 for s,u in a.items():
  for t,v in b.items():
   if not s&t:c[s|t]=c.get(s|t,0)+u*v
 return {s:v for s,v in c.items() if v}

def zero(n):return [[{} for _ in range(n)] for _ in range(n)]
def identity(n):
 a=zero(n)
 for i in range(n):a[i][i]={0:1}
 return a

def mm(a,b):
 c=[[{} for _ in b[0]] for _ in a]
 for i in range(len(a)):
  for k in range(len(b)):
   if a[i][k]:
    for j in range(len(b[0])):
     if b[k][j]:c[i][j]=add(c[i][j],mul(a[i][k],b[k][j]))
 return c

def ma(a,b):return [[add(x,y) for x,y in zip(ar,br)] for ar,br in zip(a,b)]
def star(a,bound):
 r=p=identity(len(a))
 for _ in range(bound):
  p=mm(p,a);r=ma(r,p)
 if any(any(row) for row in mm(p,a)):raise RuntimeError('nilpotence bound failed')
 return r

def eliminate(a,indices,bound):
 keep=[i for i in range(len(a)) if i not in indices]
 if not indices:return a,keep
 b=[[a[i][j] for j in keep] for i in keep]
 if not keep:raise RuntimeError('empty boundary')
 bi=[[a[i][j] for j in indices] for i in keep]
 ii=[[a[i][j] for j in indices] for i in indices]
 ib=[[a[i][j] for j in keep] for i in indices]
 return ma(b,mm(mm(bi,star(ii,bound)),ib)),keep

def adjacency(n,mask):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if mask>>k&1:a[i][j]=1
  else:a[j][i]=1
 return a

def augmented(a):
 n=len(a);m=zero(n+2)
 for j in range(n):m[0][j+1]={1<<j:1}
 for i in range(n):
  m[i+1][n+1]={0:1}
  for j in range(n):
   if a[i][j]:m[i+1][j+1]={1<<j:a[i][j]}
 return m

def brute(a):
 n=len(a);r={}
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(a[p[i]][p[i+1]] for i in range(k-1)):
    s=sum(1<<i for i in p);r[s]=r.get(s,0)+1
 return r

def main():
 totals={};staged=0
 for n in range(1,6):
  count=0
  for mask in range(1<<(n*(n-1)//2)):
   a=adjacency(n,mask);m=augmented(a);want=brute(a)
   if star(m,n+1)[0][-1]!=want:raise RuntimeError(('direct',n,mask))
   for imask in range(1<<n):
    inds=[i+1 for i in range(n) if imask>>i&1]
    red,keep=eliminate(m,inds,n+1)
    if star(red,n+1)[0][-1]!=want:raise RuntimeError(('response',n,mask,imask))
    cut=len(inds)//2;first=inds[:cut];rest=inds[cut:]
    r1,k1=eliminate(m,first,n+1)
    r2,k2=eliminate(r1,[k1.index(i) for i in rest],n+1)
    if r2!=red:raise RuntimeError(('staged',n,mask,imask))
    count+=1;staged+=1
  totals[n]=count
 # ker C is not invariant: C=(1,0), swap sends blind e2 to visible e1.
 C=(1,0);e2=(0,1);swap_e2=(1,0)
 if sum(a*b for a,b in zip(C,e2))!=0 or sum(a*b for a,b in zip(C,swap_e2))!=1:raise RuntimeError('hostile')
 # Curvature: Ch=C but CAh e2 differs from CA e2 over Q.
 h_e2=(0,-1)
 A_h_e2=(h_e2[1],h_e2[0])
 if C[0]*A_h_e2[0]+C[1]*A_h_e2[1]!=-1:raise RuntimeError('curvature hostile')
 # Sharp horizon d-r for C selecting coordinate 0 and shift e_j -> e_(j-1).
 horizons={}
 for d in range(1,9):
  seen={0};steps=0
  while len(seen)<d:
   seen|={j+1 for j in seen if j+1<d};steps+=1
  if steps!=d-1:raise RuntimeError('horizon')
  horizons[d]=steps
 print(json.dumps({'status':'PASS','exhaustive_tournament_boundary_cases':totals,'staged_eliminations':staged,'sharp_observation_horizons':horizons,'hostile_currently_blind_then_visible':True},indent=2,sort_keys=True))
if __name__=='__main__':main()
