from itertools import permutations
from collections import Counter
def tour(n,rule):
 return [{j for j in range(n) if i!=j and rule(i,j)} for i in range(n)]

X=tour(6,lambda i,j:(i<j)!= (tuple(sorted((i,j))) in {(0,5),(1,3)}))
P=tour(7,lambda i,j:(j-i)%7 in {1,2,4})
def by_dp(adj):
 n=len(adj);dp=[[0]*n for _ in range(1<<n)];h=[1]+[0]*((1<<n)-1)
 for s in range(1,1<<n):
  for i in range(n):
   if s>>i&1:
    t=s^(1<<i)
    dp[s][i]=1 if not t else sum(dp[t][j] for j in range(n) if i in adj[j])
  h[s]=sum(dp[s])
 f=(1<<n)-1
 return h[-1],sum(h[s]*h[f^s] for s in range(1,f))
def direct(adj):
 n=len(adj)
 obj=[[],[]]
 for p in permutations(range(n)):
  failures=[i for i in range(n-1) if p[i+1] not in adj[p[i]]]
  if not failures:
   obj[0].append((p,))
   obj[1].extend((p[:i],p[i:]) for i in range(1,n))
  if len(failures)==1:
   i=failures[0]+1
   obj[1].append((p[:i],p[i:]))
 return obj

def wrap(adj,a):
 n=len(adj);S=range(n,n+a);b=n+a
 return [row|set(S) for row in adj]+[set(range(i+1,b+1)) for i in S]+[set(range(n))]
for adj,expect in [(X,(27,234)),(P,(189,1638))]:
 assert by_dp(adj)==expect==tuple(map(len,direct(adj)))
 print("core profile:",expect)
for family,exps in [(X,[315,837,1881]),(P,[2205,5859,13167])]:
 for a,expected in enumerate(exps,1):
  assert by_dp(wrap(family,a))[0]==expected
 print("wrap:",exps)
for k,O in enumerate(direct(P),1):
 reps={tuple(tuple((v-obj[0][0])%7 for v in slot) for slot in obj) for obj in O}
 assert len(reps)*7==len(O)
 print("normalized object count",k,len(reps))

base=['CAb','AbC','bCA','CAbC','AbCA','CAbCA','bC']
arcs={('C','A'),('A','b'),('b','C')}
for w in base:
 for S in range(1<<w.count('A')):
  i=0;out=''
  for c in w:
   if c!='A':out+=c
   else:
    if S>>i&1:out+=c
    i+=1
  assert out.count('C')==w.count('C')
  print(w,'→',out,all(e in arcs for e in zip(out,out[1:])))

def all_cycles(adj):
 n=len(adj)
 for l in range(3,n+1,2):
  for p in permutations(range(n),l):
   if p[0]==min(p) and all(p[(i+1)%l] in adj[p[i]] for i in range(l)):
    yield p
cy=list(all_cycles(P))
def pack(C,i=0,used=frozenset()):
 yield ()
 for j in range(i,len(C)):
  q=C[j]
  if not set(q)&used:
   for tail in pack(C,j+1,used|set(q)):
    for c in (0,1):yield ((q,c),)+tail
def rot(p,k):
 q=tuple((x+k)%7 for x in p)
 j=q.index(min(q))
 return q[j:]+q[:j]
O=set(tuple(sorted(p)) for p in pack(cy))
def T(p,k):
 return tuple(sorted((rot(c,k),color) for c,color in p))
fixed=[p for p in O if T(p,1)==p]
assert len(O)==189 and len(fixed)==7
assert len({min(T(p,i) for i in range(7)) for p in O-set(fixed)})==26
print("colored cycle packs: 7 fixed, 26 free orbits")
