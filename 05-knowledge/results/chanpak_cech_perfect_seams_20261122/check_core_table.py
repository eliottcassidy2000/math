from itertools import *
from collections import *
V=tuple(range(4))
names={ (0,1,2,3):'chain',(0,2,2,2):'sink+3cycle',(1,1,1,3):'source+3cycle',(1,1,2,2):'strong core' }
seen=set()
for B in product(*[{(i,j),(j,i)} for i,j in combinations(V,2)]):
 A=set(B); d=[sum((i,j) in A for j in V) for i in V]; D=tuple(sorted(d))
 if D in seen: continue
 seen.add(D)
 P=[p for p in permutations(V) if all(e in A for e in zip(p,p[1:]))]
 def cyc(S):return all(sum((a,b) in A for b in S)==1 for a in S)
 C=sum(cyc(S) for S in combinations(V,3))
 S=[sum(p[0]==i for p in P) for i in V]; E=[sum(p[-1]==i for p in P) for i in V]
 print(names[D],d,sorted(A),'C',C,'S',S,'E',E)
 out=defaultdict(list)
 for a,b in product(V,V):
  t=C+d[a]+3-d[b]+(a==b)
  f=S[a]+E[b]+sum((a,c) in A and (c,b) in A for c in V)+2*((a,b) in A)+3*(a==b)
  q=(a==b and cyc(set(V)-{a}))+sum((a,c) in A and (e,b) in A and len({a,b,c,e})==4 for c,e in product(V,V))
  out[t].append((a,b,f,q, t+f+2*q))
 for t,O in sorted(out.items()):print(t, sorted(set(o[-1] for o in O)))
