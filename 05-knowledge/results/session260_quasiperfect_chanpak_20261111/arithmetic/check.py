"""Exact controls; no external dependencies. Run from any working directory."""
import json,math
from itertools import combinations,permutations
from pathlib import Path

def require(x,msg):
 if not x:raise RuntimeError(msg)
def tournament(a,b):
 n=a+b+1; adj=[0]*n
 for i,j in combinations(range(n),2):
  if j==n-1 and i<a:adj[j]|=1<<i
  else:adj[i]|=1<<j
 return adj
def flip(adj,i,j):
 out=adj.copy();out[i]^=1<<j;out[j]^=1<<i;return out
def h(adj):
 n=len(adj);dp=[{} for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for v,c in dp[mask].items():
   avail=adj[v]&~mask
   while avail:
    bit=avail&-avail;avail-=bit;w=bit.bit_length()-1;d=dp[mask|bit];d[w]=d.get(w,0)+c
 return sum(dp[-1].values())
def hp(adj):
 return sum(all(adj[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in permutations(range(len(adj))))
def fan(k):
 n=k+2;a=k//2;b=k+1-a;loss={a,b,k+1};adj=[0]*n
 for i,j in combinations(range(n),2):
  if i==0 and j in loss:adj[j]|=1
  else:adj[i]|=1<<j
 return adj
def join(x,y):
 n=len(x);tail=((1<<len(y))-1)<<n
 return [a|tail for a in x]+[a<<n for a in y]
def sigma(n):return sum(d+n//d if d*d!=n else d for d in range(1,math.isqrt(n)+1) if n%d==0)
rows=[];checks=0
for a in range(2,7):
 for b in range(2,7):
  t=tournament(a,b);base=h(t);alt=h(flip(t,a-1,a));A=2**(a-1);B=2**(b-1)
  require(base==(2**a-1)*(2**b-1)+2,'base');require(alt==base+2-A-B,'flip');require(4*alt-11==(2**(a+1)-3)*(2**(b+1)-3),'chamber');checks+=1
  if len(t)<=7:require(hp(t)==base and hp(flip(t,a-1,a))==alt,'permutation cross-check')
for a in range(2,6):
 t=tournament(a,a+1);base=h(t);values={}
 for i,j in combinations(range(len(t)),2):values.setdefault(h(flip(t,i,j)),[]).append([i,j])
 losses=[base-v for v in values if v<base]
 require(base-2 not in values,'unexpected target hit')
 rows.append(dict(a=a,base=base,target=base-2,edge_count=math.comb(len(t),2),target_edges=[],min_decrease=min(losses),attainers=values[base-min(losses)]))
fans=[]
for k in range(4,11):
 value=h(fan(k));require(value==2**k-1,'fan');fans.append(dict(k=k,order=k+2,H=value))
w105=[116,45,80,101,10,84,18]
require(h(w105)==105 and hp(w105)==105,'105 explicit witness')
joined=[]
for a in range(4,6):
 t=join(fan(a),fan(a+1));value=h(t);require(value==2**a*(2**(a+1)-3)+1,'joined neighbor');joined.append(dict(a=a,order=len(t),H=value))
source_rows=[]
for n in [20,104,464,650,1952,130304,522752]:
 require(sigma(n)==2*n+2,'historical sample')
 a=0;m=n
 while m%2==0:a+=1;m//=2
 require(n+1==(2**(a+1)-1)*(sigma(m)//2),'general factor');source_rows.append(dict(n=n,a=a,m=m,half_sigma=sigma(m)//2,Chowla=n+1))
require(sigma(8000)!=16002,'prime-guard hostile')
w217=[0]*10
for i,j in combinations(range(10),2):
 if i==0 and j in {2,5,6,9}:w217[j]|=1
 else:w217[i]|=1<<j
require(h(w217)==217,'217 apex witness')
w651=join(tournament(1,1),w217)
require(h(w651)==651,'651 product witness')
result=dict(witness217=w217,witness651=w651,status='PASS' ,crossing_formula_cases=checks,neighborhoods=rows,fans=fans,witness105=w105,joined_neighbors=joined,excess_sources=source_rows,composite_guard_hostile=dict(n=8000,sigma=sigma(8000)))
print(json.dumps(result,indent=2))
