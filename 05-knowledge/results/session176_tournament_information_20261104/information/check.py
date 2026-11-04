#!/usr/bin/env python3
import itertools,json,math
from fractions import Fraction

def require(b,m):
 if not b:raise RuntimeError(m)

def tournament(n,c):
 g=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if c>>k&1:g[i]|=1<<j
   else:g[j]|=1<<i
   k+=1
 return g

def counts(g):
 n=len(g);d=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n)
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  h[s]=sum(d[s])
  for v,x in enumerate(d[s]):
   rem=g[v]&~s
   while x and rem:
    b=rem&-rem;rem-=b;d[s|b][b.bit_length()-1]+=x
 return h

def profile(g):
 h=counts(g);N=len(h);n=len(g);d=[[0]*(n+1) for _ in range(N)];d[0][0]=1
 for s in range(1,N):
  t=s
  while t:
   for k in range(1,n+1):d[s][k]+=h[t]*d[s^t][k-1]
   t=(t-1)&s
 return d[-1][1:]

def probe(g,r):
 n=len(g);G=g[:]+[0]*(2*r)
 for i in range(2*r):
  v=n+i
  for j in range(i+1,2*r):G[v]|=1<<(n+j)
  if i%2:G[v]|=(1<<n)-1
  else:
   for u in range(n):G[u]|=1<<v
 return G

def kernel(r):
 # DP counts exterior permutations interleaved with indistinguishable X tokens.
 m=2*r;cap=r+1;full=(1<<m)-1;D={(0,-1,0):1}
 for length in range(m+cap):
  for (s,last,k),val in list(D.items()):
   if s.bit_count()+k!=length:continue
   if k<cap and (last==-1 or last>=0 and last%2==1):
    key=(s,-2,k+1);D[key]=D.get(key,0)+val
   for j in range(m):
    if s>>j&1:continue
    if last==-1 or last==-2 and j%2==0 or last>=0 and last<j:
     key=(s|1<<j,j,k);D[key]=D.get(key,0)+val
 return [sum(v for (s,last,k),v in D.items() if s==full and k==j) for j in range(1,cap+1)]

def literal(g):
 return sum(all(g[p[i]]>>p[i+1]&1 for i in range(len(g)-1)) for p in itertools.permutations(range(len(g))))

def rank(rows):
 a=[[Fraction(x) for x in row] for row in rows];k=0
 for col in range(len(a[0])):
  pivot=next((i for i in range(k,len(a)) if a[i][col]),None)
  if pivot is None:continue
  a[k],a[pivot]=a[pivot],a[k];v=a[k][col];a[k]=[x/v for x in a[k]]
  for i in range(len(a)):
   if i!=k:
    v=a[i][col];a[i]=[x-v*y for x,y in zip(a[i],a[k])]
  k+=1
 return k

def main():
 K=[kernel(r) for r in range(7)]
 for r,row in enumerate(K):require(row[-1]==math.factorial(r),'diagonal')
 profiles={};checks=0;collision=None;ranks={}
 for n in range(1,6):
  P={}
  for c in range(1<<(n*(n-1)//2)):
   g=tournament(n,c);q=profile(g);P.setdefault(tuple(q),g)
   require(q[0]==literal(g),'literal H')
   for r in range(2):
    expected=sum(x*y for x,y in zip(K[r],q));require(counts(probe(g,r))[-1]==expected,'probe');checks+=1
  ranks[n]={'distinct_profiles':len(P),'linear_rank':rank(list(P))}
  profiles[n]=P
  byH={}
  for q,g in P.items():
   if q[0] in byH and collision is None:
    old,og=byH[q[0]];rr=next(i for i in range(n) if sum(x*y for x,y in zip(K[i],q))!=sum(x*y for x,y in zip(K[i],old)))
    collision={'n':n,'H':q[0],'A':og,'B':g,'qA':old,'qB':q,'probe_r':rr,'observedA':counts(probe(og,rr))[-1],'observedB':counts(probe(g,rr))[-1]}
   byH[q[0]]=(q,g)
 for n,P in profiles.items():
  for q,g in P.items():
   obs=[]
   for r in range(n):
    val=sum(x*y for x,y in zip(K[r],q));obs.append(val)
    if n+2*r<=11:require(counts(probe(g,r))[-1]==val,'deep probe');checks+=1
   recovered=[]
   for r,val in enumerate(obs):recovered.append((val-sum(K[r][i]*recovered[i] for i in range(r)))//math.factorial(r))
   require(recovered==list(q),'recovery')
 require(collision is not None,'scalar collision')
 print(json.dumps({'status':'PASS','probe_checks':checks,'kernels':K,'realized_profile_ranks_through_5':ranks,'same_H_distinguished_by_uniform_context':collision},indent=2))
if __name__=='__main__':main()
