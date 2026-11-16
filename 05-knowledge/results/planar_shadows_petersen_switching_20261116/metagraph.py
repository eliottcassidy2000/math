from itertools import combinations, permutations
from collections import Counter
n=5
free=[(i,j) for i in range(n) for j in range(i+2,n)]
pairs=list(combinations(range(n),2))
def tiling(t):
    a=[[0]*n for i in range(n)]
    for i in range(n-1):a[i][i+1]=1
    for k,(i,j) in enumerate(free):a[i][j]=t>>k&1
    for i,j in pairs:a[j][i]=1-a[i][j]
    return a
def h(a):
    return sum(all(a[i][j] for i,j in zip(p,p[1:])) for p in permutations(range(n)))
def canon(a):
    return min(sum((a[p[i]][p[j]]^b)<<k for k,(i,j) in enumerate(pairs)) for p in permutations(range(n)) for b in (0,1))
T=[tiling(t) for t in range(64)]
C=[canon(a) for a in T]
V=sorted(set(C));idx={x:i for i,x in enumerate(V)}
mult=Counter(tuple(sorted((idx[C[t]],idx[C[t^63]]))) for t in range(32))
E={e for e in mult if e[0]!=e[1]};N=[set()for _ in V]
for i,j in E:N[i].add(j);N[j].add(i)
print("nodes,edges",len(V),len(E))
print("H",[h(T[C.index(v)]) for v in V]);print("degs",sorted(map(len,N)))
print("pairs",sorted(E));print("mult",sorted(mult.items()))
P=[tuple(j for j in range(10) if len(N[i])==len(N[j])) for i in range(10)]
def enum(M={}):
 if len(M)==len(V):yield M;return
 i=min((i for i in range(len(V))if i not in M),key=lambda j:len(P[j]))
 for j in P[i]:
  if j not in M.values() and all((v in N[i])==(w in N[j]) for v,w in M.items()):
   yield from enum({**M,i:j})
print("aut",list(enum()))
print("free-bit single flip",len({tuple(sorted((idx[C[t]],idx[C[t^(1<<k)]])))for t in range(64)for k in range(6)if C[t]!=C[t^(1<<k)]}))
print("all-arc single flip",len({tuple(sorted((i,j)))for i,v in enumerate(V)for x,y in pairs for a in [list(map(list,T[C.index(v)]))] for _ in [(a[x].__setitem__(y,1-a[x][y]),a[y].__setitem__(x,1-a[y][x]))] if (j:=idx[canon(a)])!=i}))
S={tuple(sorted((idx[C[t]],idx[C[t^(1<<k)]])))for t in range(64)for k in range(6)if C[t]!=C[t^(1<<k)]}
for name,z in [('complement',E),('single_flip',S)]:
    print(name,"literal K5", [e for e in combinations(range(10),5) if all(a in z for a in combinations(e,2))])
    print(name,"literal K3,3",[
     (x,y) for x in combinations(range(10),3) for y in combinations(sorted(set(range(10))-set(x)),3)
           if x<y and all(tuple(sorted((i,j))) in z for i in x for j in y)])
print("intersect",len(E&S),"complement only",sorted(E-S),"single only",sorted(S-E))
def minor(g,k,conn):
  n=10
  def rec(i,parts):
    if len(parts)>k or len(parts)+n-i<k:return
    if i==n:
      if len(parts)!=k:return
      h=set()
      for a,b in combinations(range(k),2):
        if any(tuple(sorted(e)) in g for e in __import__('itertools').product(parts[a],parts[b])):h.add((a,b))
      if not conn(h):return
      for x in parts:
        r=set(x[:1])
        while True:
          more={i for i in x if any(tuple(sorted((i,j))) in g for j in r)}
          if not more-r:break
          r|=more
        if len(r)!=len(x):return
      yield parts
    else:
      # A spanning minor is enough in a connected graph: remaining vertices
      # can be assigned to a nearest branch set without losing edges.
      for j in range(len(parts)):
        p=parts.copy();p[j]=(*p[j],i)
        yield from rec(i+1,p)
      if len(parts)<k:yield from rec(i+1,parts+[(i,)])
  return next(rec(0,[]),None)
print("single-flip K5 minor",minor(S,5,lambda h:len(h)==10))
print("single-flip K3,3 minor",minor(S,6,lambda h:any(all(tuple(sorted((i,j)))in h for i in q for j in set(range(6))-set(q))for q in combinations(range(6),3))))
