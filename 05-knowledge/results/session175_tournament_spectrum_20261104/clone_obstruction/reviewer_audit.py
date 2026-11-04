#!/usr/bin/env python3
"""Independent recursive endpoint and SCC/prime-density boundary checks."""
from itertools import combinations,product,permutations
from functools import lru_cache
from math import comb,factorial,log
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent

def require(x,msg):
 if not x:raise ArithmeticError(msg)
def count(a):
 n=len(a)
 @lru_cache(None)
 def end(mask,v):
  if mask==1<<v:return 1
  rem=mask^(1<<v)
  return sum(end(rem,u) for u in range(n) if rem>>u&1 and a[u][v])
 return sum(end((1<<n)-1,v) for v in range(n))
def quotient(n,bits):
 a=[[False]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):a[i][j]=bool(bits>>k&1);a[j][i]=not a[i][j]
 return a
def strong(a):
 n=len(a);r=[x[:] for x in a]
 for k in range(n):
  for i in range(n):
   for j in range(n):r[i][j]|=r[i][k] and r[k][j]
 return n==1 or all(all(row) for row in r)
def inflate(a,rs):
 labels=[(i,k) for i,r in enumerate(rs) for k in range(r)]
 return [[(x==y and k<l) or (x!=y and a[x][y]) for y,l in labels] for x,k in labels]
def ordered(a,b):
 n=len(a);m=len(b)
 return [row+[True]*m for row in a]+[[False]*n+row for row in b]
def main():
 trials=literal_checks=0;cores=[]
 for n in (3,4):
  for code in range(1<<(n*(n-1)//2)):
   a=quotient(n,code)
   if not strong(a):continue
   cores.append(a)
   for rs in product((1,2),repeat=n):
    g=inflate(a,rs);h=count(g);require(strong(g),'strongness')
    require(all(h>=2**r+1 for r in rs),'triangle parameter bound');trials+=1
    if len(g)<=5:
     direct=sum(all(g[p[i]][p[i+1]] for i in range(len(g)-1)) for p in permutations(range(len(g))))
     require(h==direct,'literal independent');literal_checks+=1
 c=[[False,True,False],[False,False,True],[True,False,False]]
 for r in range(1,10):require(count(inflate(c,(r,1,1)))==2**r+1,'inflated triangle')
 transitive=quotient(3,7)
 pieces=[[[False]],transitive,c,cores[-1]]
 for a,b in product(pieces,repeat=2):require(count(ordered(a,b))==count(a)*count(b),'ordered sum product')
 # The singleton exception is material, even though the main theorem handles it.
 require(count(inflate([[False]],(10,)))==1,'singleton boundary')
 primes=[p for p in range(2,201) if all(p%d for d in range(2,p))]
 for m in range(1,101):
  b=comb(2*m,m);require(b*(2*m+1)>=4**m,'central binomial lower')
  pi=sum(p<=2*m for p in primes);require(b<=(2*m)**pi,'central binomial upper')
  for p in primes:
   if p>2*m:break
   q=b;e=0
   while q%p==0:q//=p;e+=1
   require(p**e<=2*m,'prime power factor bound')
 print(json.dumps({'status':'PASS','independent_binary_size_blowups':trials,'literal_path_crosschecks':literal_checks,'ordered_sum_checks':16,'prime_bound_checks':100,'triangle_r':list(range(1,10)),'singleton_boundary_H':1,'author_checker_sha256':hashlib.sha256((P/'check.py').read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()
