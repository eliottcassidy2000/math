from itertools import permutations,combinations,product
from math import factorial
from functools import lru_cache
import json

def need(x,s):
 if not x:raise RuntimeError(s)

@lru_cache(None)
def compositions(n,k):
 if k==1:return ((n,),)
 return tuple((i,)+p for i in range(n+1) for p in compositions(n-i,k-1))

def half(r,root,C,left):
 V=[v for v in range(len(r)) if v!=root and v not in C and bool(r[root]>>v&1)==left];ans=0
 for perm in permutations(V):
  for lens in compositions(len(V),len(C)+1):
   j=0;ok=True
   for i,size in enumerate(lens):
    w=perm[j:j+size];j+=size
    if any(not r[u]>>v&1 for u,v in zip(w,w[1:])):ok=False;break
    if i<len(C):
     c=C[i];cO=bool(r[root]>>c&1)
     if not w:
      if cO!=left:ok=False;break
     elif not ((r[w[-1]]>>c&1) if left else (r[c]>>w[0]&1)):ok=False;break
   ans+=ok
 return ans

def runs(r,root):
 old=[v for v in range(len(r)) if v!=root];a=[0]*(len(r)+1)
 for p in permutations(old):
  for mask in range(1<<(len(old)+1)):
   w=[]
   for i in range(len(old)+1):
    if mask>>i&1:w.append(root)
    if i<len(old):w.append(p[i])
   if all(r[u]>>v&1 for u,v in zip(w,w[1:])):a[mask.bit_count()]+=1
 return a

def check(r):
 roots=covers=0
 for root in range(len(r)):
  O=[v for v in range(len(r)) if r[root]>>v&1];I=[v for v in range(len(r)) if v!=root and r[v]>>root&1]
  edges=[(u,v) for u in O for v in I if r[u]>>v&1]
  old=[v for v in range(len(r)) if v!=root];found=[]
  for nu in range(len(old)+1):
   found=[C for C in combinations(old,nu) if all(u in C or v in C for u,v in edges)]
   if found:break
  a=runs(r,root);need(a[nu+1]>0 and not any(a[nu+2:]),'literal maximum-run degree')
  for C in found:
   L=half(r,root,C,True);R=half(r,root,C,False)
   need(L>0 and R>0 and factorial(nu)*L*R==a[nu+1],'literal half-cover factorization');covers+=1
  roots+=1
 return roots,covers

graphs=roots=covers=0
for n in range(1,6):
 pairs=list(combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  r=[0]*n
  for k,(u,v) in enumerate(pairs):
   if bits>>k&1:u,v=v,u
   r[u]|=1<<v
  nr,nc=check(r);graphs+=1;roots+=nr;covers+=nc
need((graphs,roots,covers)==(1099,5405,8747),'complete labelled controls')
print(json.dumps({'status':'PASS','labelled_graphs':graphs,'literal_root_word_checks':roots,'literal_minimum_cover_factorizations':covers},indent=2,sort_keys=True))
