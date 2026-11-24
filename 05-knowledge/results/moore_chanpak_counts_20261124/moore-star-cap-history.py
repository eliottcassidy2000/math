from collections import Counter
from itertools import combinations

def hs(i,j):
 side,a=divmod(i,25); other,b=divmod(j,25)
 x,a=divmod(a,5);y,b=divmod(b,5)
 if side==other:
  return x==y and (a-b)%5 in ((1,4) if side==0 else (2,3))
 if side: x,a,y,b=y,b,x,a
 return (b-a-x*y)%5==0

P=list(combinations(range(5),2))
cases=[
(2,lambda i,j:(i-j)%5 in (1,4),tuple(range(5))),
(3,lambda i,j:not set(P[i])&set(P[j]),(0,1,2,6,8,9)),
(7,hs,(3,15,17,18,21,22,28,38,44,47)),
]
def ck(ok,m):
 if not ok:raise ValueError(m)
def h(A):
 n=len(A);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for s in range(1,1<<n):
  for i,v in enumerate(D[s]):
   w=A[i]&~s
   while v and w:
    b=w&-w;w-=b
    D[s|b][b.bit_length()-1]+=v
 return sum(D[-1])
for d,adj,S in cases:
 n=d*d+1
 G=[{j for j in range(n) if adj(i,j)} for i in range(n)]
 ck(all(len(g)==d for g in G),d)
 ck(all(len(G[i]&G[j])==int(j not in G[i]) for i,j in combinations(range(n),2)),d)
 cnt=Counter(len(g&set(S)) for g in G)
 ck(len(S)==d+3 and 1 not in cnt and max(cnt)==2,(d,cnt))
 cyc=[set(S)]
 comps=[]
 seen=set()
 for i in S:
  if i in seen:continue
  C={i}; todo=[i];seen.add(i)
  while todo:
   u=todo.pop()
   for v in G[u]&set(S)-seen:
    seen.add(v);C.add(v);todo.append(v)
  ck(all(len(G[u]&C)==2 for u in C),C)
  comps.append(len(C))
 if n<=10:
  print('ordinary path words and edge-reversal completion',
   h([sum(1<<j for j in g) for g in G]),
   h([sum(1<<j for j in range(n) if (i>j if j in G[i] else i<j)) for i in range(n)]))
 print(d,len(S),sorted(cnt.items()),sorted(comps))
