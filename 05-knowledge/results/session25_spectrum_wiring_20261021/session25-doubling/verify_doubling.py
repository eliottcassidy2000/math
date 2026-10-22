import math,json

def require(v):
 if not v:raise RuntimeError('verification failed')
def phi(n):return sum(math.gcd(a,n)==1 for a in range(1,n+1))
def order(d):
 a=2%d;t=1
 while a!=1:a=a*2%d;t+=1
 return t
def valuation(a):
 if a==0:return 10**9
 s=0
 while a%2==0:s+=1;a//=2
 return s
checks=0;examples=[]
for M in range(2,1001):
 r=M-1;k=2*r+1;s=valuation(M);u=M//2**s
 cycles=set()
 for a in range(M):
  if a:
   c=k-2*a;b=0 if c==r else (k-c if c>r else r-c)
   require(b==2*a%M);checks+=1
  seen={};seq=[];x=a
  while x not in seen:seen[x]=len(seq);seq.append(x);x=2*x%M
  cycle=seq[seen[x]:];cycles.add(tuple(sorted(cycle)))
  require(seen[x]==max(0,s-valuation(a)));checks+=1
 C=len(cycles)-1;formula=sum(phi(d)//order(d) for d in range(3,u+1,2) if u%d==0)
 require(C==formula);require((C==0)==(u==1));checks+=2
 if M in [2,3,4,5,6,8,21,63,341]:examples.append({'M':M,'n':4*M-3,'nonzero_cycles':C,'reversal_assignments':2**(1+C)})
print(json.dumps({'status':'PASS','moduli':[2,1000],'checks':checks,'examples':examples},indent=2))
