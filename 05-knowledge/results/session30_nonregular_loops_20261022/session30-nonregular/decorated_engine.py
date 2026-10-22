"""Independent local-matching resummation of reduced forest clusters."""
import itertools as it
from collections import Counter,defaultdict
from fractions import Fraction
from math import factorial
import numpy as np
from core_kernel import matchings

def cores(e,l):
 for b in range(max(0,1-l),2*e+l+1):
  v=b+l;E=v+e
  if v<1 or E<1:continue
  pairs=list(it.combinations_with_replacement(range(v),2));seen=set()
  for ed in it.combinations_with_replacement(pairs,E):
   d=Counter(a for edge in ed for a in edge)
   if any(d[i]<3 for i in range(b)) or any(d[i]!=1 for i in range(b,v)):continue
   reach={0}
   while True:
    rr=reach|{a for a,c in ed if c in reach}|{c for a,c in ed if a in reach}
    if rr==reach:break
    reach=rr
   if len(reach)!=v:continue
   variants=[]
   for p in it.permutations(range(b)):
    for q in it.permutations(range(b,v)):
     perm=p+q
     variants.append(tuple(sorted(tuple(sorted((perm[a],perm[c]))) for a,c in ed)))
   key=min(variants)
   if key in seen:continue
   seen.add(key)
   aut=sum(x==key for x in variants)
   for (a,c),m in Counter(key).items():aut*=factorial(m)*(2**m if a==c else 1)
   yield key,aut

# symbol 0=R+,1=R-,2=C++,3=C--,4=C+-
def symbolic_kernel(ed,aut=1):
 E=len(ed);inc=defaultdict(list)
 for i,(a,b) in enumerate(ed):inc[a].append(2*i);inc[b].append(2*i+1)
 opts=[]
 for hs in inc.values():
  oo=[]
  for pairs in matchings(hs):
   s=len(hs)-len(pairs);oo.append((pairs,(-1)**(s-1)*factorial(s-1)))
  opts.append(oo)
 ans=defaultdict(int)
 for loc in it.product(*opts):
  mu=1;base=[]
  for pairs,w in loc:mu*=w;base+=list(pairs)
  for mask in range(1<<E):
   adj=[[] for _ in range(2*E)]
   for a,b in base+[(2*j,2*j+1) for j in range(E) if not (mask>>j)&1]:adj[a].append(b);adj[b].append(a)
   groups=[];sgns=[0]*(2*E);bad=False
   for root in range(2*E):
    if sgns[root]:continue
    stack=[root];sgns[root]=1;group=[]
    while stack:
     a=stack.pop();group.append(a)
     for b in adj[a]:
      if not sgns[b]:sgns[b]=-sgns[a];stack.append(b)
    if all(len(adj[a])==2 for a in group):bad=True;break
    groups.append(group)
   if bad:continue
   gid={a:j for j,g in enumerate(groups) for a in g}
   for signs in it.product([1,-1],repeat=len(groups)):
    keys=[]
    for j in range(E):
     s=sgns[2*j]*signs[gid[2*j]];t=-sgns[2*j+1]*signs[gid[2*j+1]]
     if not(mask>>j)&1:keys.append(0 if s==1 else 1)
     elif s==t:keys.append(2 if s==1 else 3)
     else:keys.append(4)
    ans[tuple(keys)]+=mu
 return {k:Fraction(v,aut) for k,v in ans.items() if v}

def matrices(M):
 I=np.eye(len(M));Q=M@M;P=np.linalg.inv(I-M);N=np.linalg.inv(I+M);Z=np.linalg.inv(I+Q)
 return [M@P,-M@N,-Q@(I-Q)@P@P@Z,-Q@(I-Q)@N@N@Z,Q@Z]

def evaluate(ed,syms,M):
 mats=matrices(M);letters='abcdefghijklmnopqrstuvwxyz';subs=','.join(letters[a]+letters[b] for a,b in ed)+'->'
 return sum(float(c)*np.einsum(subs,*[mats[k] for k in keys],optimize=True) for keys,c in syms.items())
if __name__=='__main__':
 import json,time
 for e,l in [(-1,3),(0,2),(1,1),(2,0),(-1,2),(0,1),(1,0)]:
  cc=list(cores(e,l));print('TYPE',e,l,'COUNT',len(cc),flush=True)
  for ed,a in cc:
   t=time.time();sy=symbolic_kernel(ed,a)
   print(ed,'aut',a,'terms',len(sy),'sec',round(time.time()-t,2),flush=True)
