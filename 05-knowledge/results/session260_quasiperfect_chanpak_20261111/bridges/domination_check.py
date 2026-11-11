from functools import lru_cache
from itertools import product
from pathlib import Path
import json

def h(out):
 n=len(out);full=(1<<n)-1
 @lru_cache(None)
 def f(S,i):
  if S==full:return 1
  return sum(f(S|1<<j,j) for j in range(n) if not S>>j&1 and out[i]>>j&1)
 return sum(f(1<<i,i) for i in range(n))
def dom(adj):return sum(all(S>>i&1 or S&adj[i] for i in range(len(adj))) for S in range(1<<len(adj)))
def word_graph(w):
 n=len(w);adj=[0]*n;out=[0]*(n+1)
 for i in range(n):
  for j in range(i+1,n):
   out[i]|=1<<j
   if w[i]=='A' and w[j]=='B':adj[i]|=1<<j;adj[j]|=1<<i
  if w[i]=='A':out[n]|=1<<i
  else:out[i]|=1<<n
 return adj,out

def formula(w):return 1+sum(2 if j==i+1 else 2**(j-i-1) for i in range(len(w)) for j in range(i+1,len(w)) if w[i]=='A' and w[j]=='B')
def require(q,m):
 if not q:raise RuntimeError(m)
checks=0;corevals=set()
for n in range(0,9):
 for p in product('AB',repeat=n):
  w=''.join(p);adj,out=word_graph(w);v=formula(w)
  require(dom(adj)==h(out)==v,('adapter',w));checks+=1
  H=1;x=0
  for c in w:
   if c=='A':x=2*x+1
   else:H+=x+x%2;x*=2
  require(H==v,('register',w))
  if n<=6 and w.startswith('A') and w.endswith('B'):corevals.add(v)
require(13 not in corevals and 29 not in corevals,'hostile')
# Explicit even tree: its only nontrivial automorphism exchanges leaves1,2.
adj=[0]*6
for i,j in [(0,1),(0,2),(0,3),(3,4),(4,5)]:adj[i]|=1<<j;adj[j]|=1<<i
require(dom(adj)==29,'tree29')
ineq=[]
for a,b in [(2,2),(2,3),(3,3)]:
 failures=[];equals=0
 for mask in range(1<<(a*b)):
  n=a+b;out=[0]*(n+1);adj=[0]*n
  for i in range(n):
   for j in range(i+1,n):
    forward=(j<a or i>=a or bool(mask>>(i*b+j-a)&1))
    if forward:out[i]|=1<<j
    else:out[j]|=1<<i
    if i<a<=j and mask>>(i*b+j-a)&1:adj[i]|=1<<j;adj[j]|=1<<i
   if i<a:out[n]|=1<<i
   else:out[i]|=1<<n
  D=dom(adj);H=h(out)
  if H<D:failures.append(mask)
  if H==D:equals+=1
 ineq.append({'a':a,'b':b,'graphs':1<<(a*b),'H_less_d':failures,'equality_cases':equals})
res={'status':'FINITE-EXACT','all_binary_words_length_at_most_8':checks,'core_values_length_at_most_6':sorted(corevals),'even_tree29_adjacency':adj if False else [14,1,1,17,40,16],'exploratory_bipartite_inequality':ineq}
Path(__file__).with_name('domination_receipt.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
