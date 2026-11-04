#!/usr/bin/env python3
"""Independent representation and algorithms; no producer imports."""
import itertools,json,math

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def hp(n,arcs):
 d={(1<<u,u):1 for u in range(n)}
 for S in range(1,1<<n):
  for v in range(n):
   if not S>>v&1:continue
   if S!=1<<v:d[S,v]=sum(d.get((S^(1<<v),u),0) for u in range(n)if(u,v)in arcs)
 return sum(d.get(((1<<n)-1,v),0)for v in range(n))

def scc(n,arcs):
 reach=[]
 for v in range(n):
  seen={v};todo=[v]
  while todo:
   x=todo.pop()
   for y in range(n):
    if (x,y)in arcs and y not in seen:seen.add(y);todo.append(y)
  reach.append(seen)
 return sorted({tuple(w for w in range(n)if w in reach[v]and v in reach[w])for v in range(n)})

def inflate(n,arcs,sizes):
 vs=[(i,j)for i in range(n)for j in range(sizes[i])]
 es={(u,v)for u,(i,a)in enumerate(vs)for v,(j,b)in enumerate(vs)if(i==j and a<b)or(i,j)in arcs}
 return len(vs),es

def main():
 N=10000;prime=[True]*(N+1);prime[:2]=[False,False]
 for p in range(2,math.isqrt(N)+1):
  if prime[p]:
   for m in range(p*p,N+1,p):prime[m]=False
 base=inflated=literal=0
 for n in range(1,6):
  pairs=list(itertools.combinations(range(n),2))
  for b in range(1<<len(pairs)):
   a={(u,v)if b>>i&1 else(v,u)for i,(u,v)in enumerate(pairs)}
   cs=scc(n,a);factors=[]
   for c in cs:
    ca={(c.index(u),c.index(v))for u,v in a if u in c and v in c}
    factors.append(hp(len(c),ca))
   h=hp(n,a);need(h==math.prod(factors),'SCC product')
   non=sum(len(c)>1 for c in cs)
   need(not prime[h]or non==1,'prime core filter')
   need(all(f>=3 for c,f in zip(cs,factors)if len(c)>1),'strong lower bound');base+=1
   if n<=4:
    for x in range(n):
     rs=[1]*n;rs[x]=2;m,e=inflate(n,a,rs);h=hp(m,e)
     independent=sum(all((p[i],p[i+1])in e for i in range(m-1))for p in itertools.permutations(range(m)))
     need(h==independent,'literal inflated HP');literal+=1;inflated+=1
     if non==1:
      c=next(c for c in cs if len(c)>1);ca={(c.index(u),c.index(v))for u,v in a if u in c and v in c}
      cm,ce=inflate(len(c),ca,[rs[u]for u in c]);need(h==hp(cm,ce),'padding core equality')
      need(all((1<<rs[u])+1<=h for u in c),'triangle lower bound')
     elif non==0:need(h==1,'zero cores')
     else:need(not prime[h],'composite core product')
 # Independently recompute the integer ceil-log bound and compare with a sieve-based prime count.
 samples=list(range(3,1001))+[2000,5000,10000]
 for N in samples:
  m=N//2;b=math.comb(2*m,m);q=2*m;t=0;power=1
  while b>power:t+=1;power*=q
  actual=sum(prime[p]for p in range(3,N+1)if p!=7)
  need(max(0,t-2)<=actual,'elementary prime count lower bound')
 print(json.dumps({'status':'PASS','base_tournaments':base,'independently_counted_inflations':inflated,'literal_permutation_checks':literal,'prime_bound_cutoffs':len(samples)},indent=2,sort_keys=True))
if __name__=='__main__':main()
