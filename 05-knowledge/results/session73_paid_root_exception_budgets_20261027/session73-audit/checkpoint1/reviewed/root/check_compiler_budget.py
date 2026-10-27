"""Exact independent compiler-potential controls; no removable assertions."""
from itertools import permutations
import json

def check(x,msg):
 if not x: raise RuntimeError(msg)

def fas(rows):
 n=len(rows); dp=[0]*(1<<n)
 for s in range(1,1<<n):
  dp[s]=min(dp[s^(1<<v)]+(rows[v]&(s^(1<<v))).bit_count() for v in range(n) if s>>v&1)
 return dp[-1]

def rows_of(n,mask):
 a=[0]*n; k=0
 for i in range(n):
  for j in range(i+1,n):
   u,v=(i,j) if mask>>k&1 else (j,i)
   a[u]|=1<<v;k+=1
 return a

def flip(a,p):
 n=len(a); b=[((1<<n)-1)^(1<<i)^a[i] for i in range(n)]
 for u,v in zip(p,p[1:]): b[u]|=1<<v;b[v]&=~(1<<u)
 return b

def tower(k):
 if k==0:return [0],[]
 a,ts=tower(k-1);m=len(a);n=3*m;b=[0]*n
 for i in range(3):
  for v in range(m): b[i*m+v]=(a[v]<<(i*m))|(((1<<m)-1)<<(((i+1)%3)*m))
 triangles=[tuple(i*m+x for x in t) for i in range(3) for t in ts]
 triangles += [(x,m+y,2*m+(x+y)%m) for x in range(m) for y in range(m)]
 return b,triangles

cases=paths=0; maxima=[]
for n in range(1,6):
 ps=list(permutations(range(n)));tab={};mx=0
 for mask in range(1<<(n*(n-1)//2)):
  a=rows_of(n,mask);f=fas(a);tab[tuple(a)]=f;cases+=1
  direct=min(sum(not(a[p[i]]>>p[j]&1) for i in range(n) for j in range(i+1,n)) for p in ps)
  check(f==direct,'subset DP mismatch')
 for a,f in tab.items():
  for p in ps:
   if all(a[u]>>v&1 for u,v in zip(p,p[1:])):
    g=tab[tuple(flip(a,p))];check(abs(f-g)<=n-1-min(2,n//2),'complement Lipschitz failure');mx=max(mx,abs(f-g));paths+=1
 maxima.append([n,mx])
towers=[]
for k in range(5):
 a,ts=tower(k);n=len(a);seen=set()
 for t in ts:
  check(all(a[t[i]]>>t[(i+1)%3]&1 for i in range(3)),'noncyclic packed triple')
  for i,j in ((0,1),(1,2),(2,0)):
   e=tuple(sorted((t[i],t[j])));check(e not in seen,'repeated packed edge');seen.add(e)
 check(len(seen)==n*(n-1)//2,'packing misses edge')
 back=sum(not(a[i]>>j&1) for i in range(n) for j in range(i+1,n))
 check(back==len(ts)==n*(n-1)//6,'tower optimum proof certificate failure')
 if n<=9:check(fas(a)==back,'independent tower optimum failure')
 towers.append({'k':k,'n':n,'minimum_fas':back,'cycle_debt_lower':max(0,back-n+1),'weighted_complement_bill_lower':max(0,back-n+1),'primitive_complement_clock_lower_weaker':max(0,(n-6+5)//6) if n>1 else 0,'primitive_complement_clock_lower_sharp_tariff':max(0,((n-1)*(n-6)+6*(n-3)-1)//(6*(n-3))) if n>=4 else 0})
# A legal complement need not decrease a cycle budget or FAS.
a=rows_of(3,7);p=(0,1,2);check(fas(a)==0 and fas(flip(a,p))==1,'hostile missing')
print(json.dumps({'status':'PASS','labelled_tournaments':cases,'legal_path_complements':paths,'max_absolute_fas_change_by_order':maxima,'tower_certificates':towers,'scope':'all labelled tournaments n<=5; ternary tower pair-disjoint triangle certificates k<=4; exact DP tower n<=9; no reachability claim'},sort_keys=True,indent=2))
