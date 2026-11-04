"""Independent finite controls for arbitrary rooted defect-slot substitution."""
import itertools,json,math
from functools import lru_cache

def rows(n,code):
 a=[0]*n
 for z,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if code>>z&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a

def counts(a):
 n=len(a); h=[0]*(1<<n);h[0]=1;dp=[[0]*n for _ in h]
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    dp[s][v]=1 if not t else sum(dp[t][u] for u in range(n) if t>>u&1 and a[u]>>v&1)
  h[s]=sum(dp[s])
 return h

def profiles(a,maxslots):
 h=counts(a);full=len(h)-1;p=[0]*len(h);p[0]=1;out=[p[full]]
 for k in range(maxslots):
  z=[0]*len(h)
  for s in range(len(h)):
   t=s
   while True:
    z[s]+=p[s^t]*h[t]
    if not t:break
    t=(t-1)&s
  p=z;out.append(p[full])
 return out

def inflate(q,r,b):
 outside=[i for i in range(len(q)) if i!=r];n=len(outside);a=[0]*(n+len(b))
 for i,u in enumerate(outside):
  for j,v in enumerate(outside):
   if q[u]>>v&1:a[i]|=1<<j
  for j in range(len(b)):
   if q[u]>>r&1:a[i]|=1<<(n+j)
   else:a[n+j]|=1<<i
 for j,row in enumerate(b):a[n+j]|=row<<n
 return a

def slots(q,r):
 out={}
 for p in itertools.permutations(i for i in range(len(q)) if i!=r):
  m=len(p);g=set();d=set()
  if q[r]>>p[0]&1:g.add(0)
  if q[p[-1]]>>r&1:g.add(m)
  for i in range(1,m):
   if not q[p[i-1]]>>p[i]&1:d.add(i)
   if q[p[i-1]]>>r&1 and q[r]>>p[i]&1:g.add(i)
  if d<=g:out[len(g),len(d)]=out.get((len(g),len(d)),0)+1
 return out

def evaluate(tab,p):
 return sum(num*sum((-1)**j*math.comb(d,j)*p[g-j] for j in range(d+1)) for (g,d),num in tab.items())

def strong(a):
 for r in range(len(a)):
  seen=1<<r
  while True:
   nxt=seen
   for i in range(len(a)):
    if seen>>i&1:nxt|=a[i]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<len(a))-1:return False
 return True

def nonempty_profiles(a):
 out=[0]*(len(a)+1)
 for perm in itertools.permutations(range(len(a))):
  d=sum(not (a[perm[i]]>>perm[i+1]&1) for i in range(len(a)-1))
  for l in range(d+1,len(a)+1):out[l]+=math.comb(len(a)-1-d,l-1-d)
 return out

def main():
 checks=0;clonechecks=0;bounds=0;hostile=None
 modules=[rows(n,c) for n in range(1,5) for c in range(1<<(n*(n-1)//2))]
 ps=[profiles(b,4) for b in modules]
 fs=[nonempty_profiles(b) for b in modules]
 for p,f in zip(ps,fs):
  for t in range(5):
   if p[t]!=sum(math.comb(t,l)*f[l] for l in range(1,min(t,len(f)-1)+1)):raise RuntimeError("profile dictionary")
 for n in range(2,5):
  for c in range(1<<(n*(n-1)//2)):
   q=rows(n,c)
   for r in range(n):
    tab=slots(q,r)
    for b,p,f in zip(modules,ps,fs):
     child=inflate(q,r,b);actual=counts(child)[-1];pred=evaluate(tab,p)
     positive=sum(num*sum(math.comb(g-d,l-d)*f[l] for l in range(max(1,d),min(g,len(f)-1)+1)) for (g,d),num in tab.items())
     if actual!=pred or positive!=pred:raise RuntimeError(("module",n,c,r,b,actual,pred,positive))
     if strong(q)!=strong(child):raise RuntimeError('strongness')
     checks+=1
    for k in range(1,7):
     b=[sum(1<<j for j in range(i+1,k)) for i in range(k)]
     actual=counts(inflate(q,r,b))[-1];pred=evaluate(tab,[g**k for g in range(n+1)])
     if actual!=pred:raise RuntimeError('clone')
     if strong(q):
      if actual<2**k+1:raise RuntimeError('bound')
      bounds+=1
     clonechecks+=1
    if n==3 and strong(q) and r==0 and hostile is None:
     hostile={'Q_rows':q,'root':r,'slot_table':[[g,d,num] for (g,d),num in sorted(tab.items())], 'TT_k_counts':[evaluate(tab,[g**k for g in range(n+1)]) for k in range(1,7)]}
 out={'rooted_module_checks':checks,'transitive_inflation_checks':clonechecks,'strong_growth_checks':bounds,'triangle_control':hostile,'status':'PASS'}
 print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
