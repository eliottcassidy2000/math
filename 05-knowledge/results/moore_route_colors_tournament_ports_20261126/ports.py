from itertools import combinations
from collections import Counter
from pathlib import Path

def hs():
 g=[set() for _ in range(50)]
 def edge(a,b):g[a].add(b);g[b].add(a)
 for i in range(5):
  for x in range(5):
   edge(5*i+x,5*i+(x+1)%5)
   edge(25+5*i+x,25+5*i+(x+2)%5)
   for j in range(5):edge(5*i+x,25+5*j+(x+i*j)%5)
 return g

def petersen():
 s=list(combinations(range(5),2))
 return [{j for j,t in enumerate(s) if not set(r)&set(t)} for r in s]

def cycle(n):return [{(i-1)%n,(i+1)%n} for i in range(n)]

def check_ports(G,d,W):
 n=len(G);P=[[-1]*n for _ in range(d)]
 edges=[(x,y) for x,y in combinations(range(n),2) if y in G[x]]
 for (x,y),k in zip(edges,map(int,W),strict=True):
  assert 0<=k<d and P[k][x]==P[k][y]==-1
  P[k][x]=y;P[k][y]=x
 assert all(sorted(p)==list(range(n)) for p in P)
 assert all(p[p[x]]==x and p[x]!=x for p in P for x in range(n))
 words=[list(range(n)),*P,*[[b[a[x]] for x in range(n)] for a in P for b in P if a is not b]]
 assert len(words)==n
 assert all(len(set(r))==n for r in zip(*words))
 print("Hoffman–Singleton",len(words),"route words; Latin on all",n,"roots")

def transpositions(G):
 n=len(G);d=len(G[0])
 assert n==d*d+1 and all(len(w)==d for w in G)
 for u,v in combinations(range(n),2):
  assert len(G[u]&G[v])==(0 if u in G[v] else 1)
  p=list(range(n));p[u],p[v]=v,u
  changed=sum({p[x] for x in G[p[w]]}!=G[w] for w in range(n))
  assert changed==2*d
 print("transposition controls",n,n*(n-1)//2,"all touch",2*d)


def petersen_obstruction():
 G=petersen()
 def match(todo):
  if not todo:
   yield []
   return
  v=min(todo)
  for w in sorted(G[v]&todo):
   for M in match(todo-{v,w}):yield [(v,w),*M]
 count=0
 for M in match(set(range(10))):
  H=[N.copy() for N in G]
  for u,v in M:H[u].remove(v);H[v].remove(u)
  left=set(range(10));cycles=[]
  while left:
   root=min(left);walk=[root];left.remove(root)
   while len(w:=H[walk[-1]]&left):
    t=min(w);walk.append(t);left.remove(t)
   cycles.append(len(walk))
  assert sorted(cycles)==[5,5]
  count+=1
 assert count==6
 print("Petersen: all six perfect matchings leave two odd 5-cycles")

if __name__=='__main__':
 g=hs();check_ports(g,7,Path(__file__).with_name('hscolors.txt').read_text().strip())
 petersen_obstruction()
 for G in (cycle(5),petersen(),g):transpositions(G)
