#!/usr/bin/env python3
from verify import *
from math import factorial,gcd,comb
from random import Random
from functools import reduce

def induced(a,vs):return tuple(sum(((a[u]>>v)&1)<<j for j,v in enumerate(vs))for u in vs)
def substitute(q,r,t):
 vs=[v for v in range(len(q))if v!=r];N=len(vs);a=[0]*(N+len(t))
 for i,u in enumerate(vs):
  for j,v in enumerate(vs):
   if q[u]>>v&1:a[i]|=1<<j
  for j in range(len(t)):
   if q[u]>>r&1:a[i]|=1<<(N+j)
   else:a[N+j]|=1<<i
 for i,row in enumerate(t):a[N+i]|=row<<N
 validate(a);return tuple(a)
def exterior_runs(q,r):
 vs=[v for v in range(len(q))if v!=r];a=induced(q,vs);N=len(a);M=1<<N
 O=sum(1<<i for i,v in enumerate(vs)if q[r]>>v&1);I=(M-1)^O
 dp=[[0]*N for _ in range(M)];od=[[0]*N for _ in range(M)]
 for i in range(N):dp[1<<i][i]=1;od[1<<i][i]=int(bool(O>>i&1))
 end=[0]*M;start=[0]*M;inside=[0]*M;end[0]=start[0]=1
 for S in range(1,M):
  for v in range(N):
   if I>>v&1:end[S]+=dp[S][v];inside[S]+=od[S][v]
   start[S]+=od[S][v]
   todo=a[v]&~S
   while todo:
    b=todo&-todo;todo-=b;j=b.bit_length()-1;dp[S|b][j]+=dp[S][v];od[S|b][j]+=od[S][v]
 outer=[0]*M
 for S in range(M):
  A=S
  while True:
   outer[S]+=end[A]*start[S^A]
   if not A:break
   A=(A-1)&S
 coeff=[];F=[0]*M;F[0]=1
 for t in range(N//2+1):
  coeff.append(sum(F[S]*outer[(M-1)^S]for S in range(M)))
  nxt=[0]*M
  for S in range(1,M):
   A=S
   while A:
    nxt[S]+=inside[A]*F[S^A];A=(A-1)&S
  F=nxt
 return coeff

def bank(r):
 k=r-1;a=[0]*(2*k+1)
 for i in range(k):a[0]|=1<<(1+i);a[1+k+i]|=1
 for i,j in combinations(range(k),2):a[1+i]|=1<<(1+j);a[1+k+i]|=1<<(1+k+j)
 for i,j in product(range(k),repeat=2):
  if i==j:a[1+i]|=1<<(1+k+j)
  else:a[1+k+j]|=1<<(1+i)
 validate(a);return tuple(a)
def contents(c):
 f=[factorial(j)*c[j]for j in range(1,len(c))];w=[factorial(j)*f[j]for j in range(len(f))]
 return reduce(gcd,f),reduce(gcd,w),f,w

def newton(P):
 vals=[value(P,x)for x in range(len(P))];out=[]
 while vals:
  out.append(vals[0]);vals=[b-a for a,b in zip(vals,vals[1:])]
 return out

def central_checks():
 columns=formal=0;rng=Random(1971427)
 for parity in (0,1):
  for k in range(13):
   gamma=factorial(2*k+1)*(1 if parity else k+1);P=basis(k,parity);f=newton(P)
   for r,z in enumerate(f):
    choose=lambda n,j:comb(n,j)if 0<=j<=n else 0
    N=choose(k,2*k+1-r)if parity else choose(k+1,2*k+2-r)+choose(k,2*k+2-r)
    require(z==gamma*N,'normalized exact binomial/Newton column')
    require((r>=k+1 or z==0)and(r!=k+1 or z==gamma),'central triangular boundary')
   columns+=1
  for d in range(1,10):
   for _ in range(40):
    u=[rng.randrange(-10,11)for k in range(d)];u[-1]=rng.randrange(1,11)
    P=add(*(scale(basis(k,parity),v)for k,v in enumerate(u)));f=newton(P)[1:];w=[factorial(j)*v for j,v in enumerate(f)]
    gu=[factorial(2*k+1)*(1 if parity else k+1)*v for k,v in enumerate(u)]
    g=reduce(gcd,f);G=reduce(gcd,w)
    require(g==reduce(gcd,gu)and G==reduce(gcd,[factorial(k)*v for k,v in enumerate(gu)]),'signed formal central contents')
    require(G==reduce(gcd,w[:d])and g**3%(G**2)==0,'short bank and valuation bound');formal+=1
 return columns,formal

def run():
 columns,formal=central_checks()
 contexts=divchecks=clones=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   q=tournament(n,code)
   for r in range(n):
    a=exterior_runs(q,r);contexts+=1
    for j,co in enumerate(a,1):require(co%factorial(j-1)==0,'free internal fragment factorial');divchecks+=1
    for m in (1,2,3):
     c=cover_polynomial(transitive(m));expected=sum(co*factorial(j)*c[j]for j,co in enumerate(a,1)if j<len(c))
     require(h_subsets(substitute(q,r,transitive(m)))[-1]==expected,'actual transitive clone versus outside fragments');clones+=1
 B=[]
 for r in range(1,7):
  a=exterior_runs(bank(r),0);require(len(a)==r and a[-1]==factorial(r-1),'physical matching bank leading coefficient')
  B.append([a[j]//factorial(j)for j in range(r)])
 require(all(B[r][r]==1 for r in range(6)),'integer unit triangular diagonal')
 inputs=physical=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   t=tournament(n,code);c=cover_polynomial(t);g,G,f,w=contents(c)
   vals=[sum(B[r][j]*w[j]for j in range(r+1))for r in range(n)]
   require(reduce(gcd,vals)==G and reduce(gcd,vals[:(n+1)//2])==G and G%g==0,'shortened physical test gcd equality');inputs+=1
   P=bin_polynomial(t);u=coordinates(P,n%2);gu=[factorial(2*k+1)*(1 if n%2 else k+1)*v for k,v in enumerate(u)]
   require(g==reduce(gcd,gu)and G==reduce(gcd,[factorial(k)*v for k,v in enumerate(gu)])and g**3%(G**2)==0,'physical central-factorial contents and sharp valuation bound')
   for p in (2,3,5,7,11):
    require((G%p==0)==(g%p==0)and(G%(p*p)==0)==(g%(p*p)==0),'radical and square-divisibility scope')
   if n<=4:
    for r in range(1,n+1):require(h_subsets(substitute(bank(r),0,t))[-1]==vals[r-1],'actual finite-bank substitution');physical+=1
 t=(30,28,136,192,172,15,55,99);validate(t);c=cover_polynomial(t);g,G,f,w=contents(c);P=bin_polynomial(t)
 require(c[1:]==[459,2889,5310,4125,1560,300,28,1],'strict physical cover vector')
 require(P==[0,0,210,0,214,0,34,0,1]and(g,G)==(9,27),'strict physical bin/universal content')
 hD=h_subsets(wrapper(t))[-1];require(hD==7155 and gcd(hD,f[0])==27,'two actual sharp-divisor contexts')
 require([value(P,d)for d in range(1,5)]==[459,6696,50571,262944],'exact physical bin evaluations')
 print(json.dumps({'status':'PASS','central_newton_columns':columns,'signed_formal_controls':formal,'rooted_contexts':contexts,'factorial_controls':divchecks,'actual_clone_counts':clones,'matching_bank':B,'input_gcd_controls':inputs,'actual_bank_counts':physical,'strict_example':{'rows':t,'bin_profile':P,'unordered_covers':c[1:],'ordered_covers':f,'weighted_covers':w,'bin_content':g,'universal_content':G,'identity_count':f[0],'triangle_count':hD}},sort_keys=True,indent=2))
if __name__=='__main__':run()
