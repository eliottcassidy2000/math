from itertools import permutations
from collections import Counter
from math import factorial

def run(n):
 pairs=list(__import__('itertools').combinations(range(n),2)); idx={p:i for i,p in enumerate(pairs)}; M=(1<<len(pairs))-1; perms=list(permutations(range(n)))
 def arc(T,i,j):
  return ((T>>idx[min(i,j),max(i,j)])&1) == (i>j)
 def perm(T,q):
  return sum((not arc(T,q[i],q[j]))<<idx[i,j] for i,j in pairs)
 def canon(T):return min(perm(T,q) for q in perms)
 chords=[(i,j) for i,j in pairs if j-i>1]; m=len(chords)
 til=[sum(((t>>k)&1)<<idx[e] for k,e in enumerate(chords)) for t in range(1<<m)]
 cls={canon(T) for T in til}; aut={C:sum(perm(C,q)==C for q in perms) for C in cls}
 # Independent base-tiling edge enumeration.
 direct=Counter()
 for t in range(1<<(m-1)):
  C,D=sorted((canon(til[t]),canon(til[t^((1<<m)-1)])))
  S=all(((t>>k)&1)==((t>>chords.index((n-1-j,n-1-i)))&1) for k,(i,j) in enumerate(chords))
  direct[C,D,'B' if S else 'K']+=1
 # Witness formula: agreement graph is an oriented Hamilton path.
 directed=Counter()
 for C in cls:
  for D in cls:
   for q in perms:
    same=[(i,j) for i,j in pairs if arc(C,i,j)==arc(D,q[i],q[j])]
    if len(same)!=n-1:continue
    out={i:[] for i in range(n)}; indeg=[0]*n
    for i,j in same:
     u,v=(i,j) if arc(C,i,j) else (j,i);out[u].append(v);indeg[v]+=1
    starts=[i for i in range(n) if indeg[i]==0]
    if len(starts)!=1 or any(len(out[i])>1 or indeg[i]>1 for i in range(n)):continue
    p=[starts[0]]
    while out[p[-1]] and len(p)<=n:p.append(out[p[-1]][0])
    if len(p)!=n or len(set(p))!=n:continue
    r={p[i]:p[n-1-i] for i in range(n)}
    blue=all(arc(C,i,j)!=arc(C,r[i],r[j]) for i,j in pairs)
    directed[C,D,'B' if blue else 'K']+=1
 got=Counter()
 for (C,D,c),v in directed.items():
  if v%(aut[C]*aut[D])!=0: raise AssertionError("divisibility")
  d=v//(aut[C]*aut[D]);
  if d!=directed[D,C,c]//(aut[C]*aut[D]): raise AssertionError("symmetry")
  if C<=D:
   if C==D and d%2: raise AssertionError("diagonal")
   got[C,D,c]=d//(2 if C==D else 1)
 if got!=direct: raise AssertionError((n,got-direct,direct-got))
 print(f'n={n}: {len(cls)} classes, {sum(got.values())} lines; all original colored multiplicities and diagonal factors agree')
for n in (3,4,5):run(n)
