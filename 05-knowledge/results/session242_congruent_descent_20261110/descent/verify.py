import json
from itertools import product, permutations
from functools import lru_cache

def require(b,m):
 if not b: raise RuntimeError(m)
def word_graph(w):
 a=[2,0]; r=0
 for l in w:
  n=len(a); a.append(0)
  for i in range(n):
   if l!='O' and i==r:a[n]|=1<<i
   else:a[i]|=1<<n
  if l=='J':r=n
 return a
def ribbon(n):
 a=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   u,v=(j,i) if j-i in(1,2) else(i,j)
   a[u]|=1<<v
 return a
def count(a):
 n=len(a); dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1<<n):
  for i in range(n):
   if not dp[s][i]:continue
   avail=a[i]&~s
   while avail:
    bit=avail&-avail;avail-=bit;j=bit.bit_length()-1
    dp[s|bit][j]+=dp[s][i]
 return sum(dp[-1]) if n else 1
def join(a,b):
 n=len(a);m=len(b)
 return [v|(((1<<m)-1)<<n) for v in a]+[v<<n for v in b]
words={3:'Z',5:'OZ',31:'ZOZOZ',35:'OZJOZ',49:'OJJOOZ',63:'OZZOOZ',217:'OJOOZZOZ',287:'OOOOZJOZ',343:'OOJOZJJOZ'}
graphs={h:word_graph(w) for h,w in words.items()};graphs[41]=ribbon(6)
# Coordinates: v3,v5,v31,v41,v7.
primes=(3,5,31,41,7)
def vec(h):
 out=[]
 for p in primes:
  v=0
  while h%p==0:h//=p;v+=1
  out.append(v)
 require(h==1,'outside support');return tuple(out)
vs={h:vec(h) for h in graphs};cost={h:len(g) for h,g in graphs.items()}
for h,g in graphs.items():
 n=len(g)
 require(all(not(g[i]>>i&1) for i in range(n)),'loops')
 require(all(((g[i]>>j)&1)+((g[j]>>i)&1)==1 for i in range(n) for j in range(i)),'not tournament')
 require(count(g)==h,('seed count',h))
 if n<=7:
  v=sum(all(g[p[i]]>>p[i+1]&1 for i in range(n-1)) for p in permutations(range(n)))
  require(v==h,'permutation cross-check')
@lru_cache(None)
def best(v):
 if not any(v):return 0,()
 choices=[]
 for h,s in vs.items():
  if all(x>=y for x,y in zip(v,s)):
   prev,ws=best(tuple(x-y for x,y in zip(v,s)))
   if prev is not None:choices.append((prev+cost[h],ws+(h,)))
 return min(choices) if choices else (None,())
def carry(c):
 if c==1 or c<0:return None
 q,r=divmod(c,3)
 return 11*q+(0,5,8)[r]
def formula(a,b,d,e,c):
 vals=[]
 for y,x,z,w in product(range(a//2+1),range(b+1),range(d+1),range(e+1)):
  f=carry(c-y-x-z-w)
  if f is not None:vals.append(3*a+4*b+7*d+6*e+2*y+3*x+3*z+4*w+f)
 return min(vals) if vals else None
checked=0
for a,b,d,e in product(range(3),repeat=4):
 for c in range(13):
  v=(a,b,d,e,c);val,ws=best(v);expected=formula(*v)
  require(val==expected,('cost mismatch',v,val,expected))
  N=3**a*5**b*31**d*41**e*7**c
  require((val is None)==(N in(7,21)),('coverage mismatch',v))
  if val is not None:
   require(sum(cost[h] for h in ws)==val,'cost recipe')
   prod=1
   for h in ws:prod*=h
   require(prod==N,'recipe target')
  checked+=1
periodic=0
for a,b,d,e in product(range(4),repeat=4):
 K=a//2+b+d+e
 for c in range(K+2,K+12):
  require(formula(a,b,d,e,c+3)==formula(a,b,d,e,c)+11,'periodic')
  periodic+=1
# Independent literal construction, no multiplication formula used to count output.
join_tests=[]
for hs in ((3,5),(5,41),(31,5),(49,3),(63,5),(35,35),(217,3),(287,3)):
 g=[];N=1
 for h in hs:g=join(g,graphs[h]);N*=h
 require(count(g)==N,('join count',hs));join_tests.append([list(hs),len(g),N])
# Cheap local division fails even though a legitimate path exists.
require(best(vec(35))[0] is not None and best((0,0,0,0,1))[0] is None,'hostile35')
# Arithmetic valuations are insufficient without a physical supplier for the atom217.
require(not any(all(x<=y for x,y in zip(vs[h],vec(217))) for h in graphs if h not in(217,31)),'atom217 supply')
print(json.dumps({'status':'FINITE-EXACT controls; all-order proof in PROOF.md',
 'seeds':{str(h):{'order':len(g),'word':words.get(h),'masks':g} for h,g in sorted(graphs.items())},
 'cost_formula_cases':checked,'periodicity_cases':periodic,'literal_joins':join_tests,
 'hostiles':{'forbidden_targets':[7,21],'unsafe_division':[35,5,7],
 'greedy49_vs343_at117649':{'three49_cost':24,'two343_cost':22},
 'missing_supplier':217},
 'examples':{str(N):{'target_vector':vec(N),'minimum_bank_order':best(vec(N))[0],'recipe':best(vec(N))[1]} for N in(35,49,63,217,287,343,3969,117649,31*41*49,5*31*41*7)}},indent=2,sort_keys=True))
