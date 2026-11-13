from itertools import combinations
class G(dict):
 def add_nodes_from(self,vs):
  for v in vs:self[v]=set()
 def add_edge(self,u,v):self[u].add(v);self[v].add(u)
 def number_of_edges(self):return sum(map(len,self.values()))//2
def clique(g,c,k):
 if k==0:return []
 for x in tuple(c):
  c.remove(x)
  if len(n:=c&g[x])>=k-1 and (r:=clique(g,n,k-1)) is not None:return [x]+r
 return None
def atlas(n):
 e=[(i,j) for i,j in combinations(range(n),2) if j>i+1]
 full=(1<<len(e))-1
 refl=[e.index((n-1-j,n-1-i)) for i,j in e]
 cls={}
 masks={}
 for m in range(1<<len(e)):
  if m in cls:continue
  A=[0]*n
  for i in range(n-1):A[i]|=1<<(i+1)
  for k,(i,j) in enumerate(e):
   if m>>k&1:A[j]|=1<<i
   else:A[i]|=1<<j
  def dfs(P,S):
   if not S:
    v=sum((not (A[P[i]]>>P[j]&1))<<k for k,(i,j) in enumerate(e))
    cls[v]=m;w=sum((v>>i&1)<<j for i,j in enumerate(refl));cls[w]=m
    return
   rem=S if not P else S&A[P[-1]]
   while rem:
    bit=rem&-rem;rem-=bit;dfs(P+[bit.bit_length()-1],S^bit)
  dfs([], (1<<n)-1)
  masks[m]=A
 B=G();M=G();K=G()
 M.add_nodes_from(masks);B.add_nodes_from(masks);K.add_nodes_from(masks)
 for m in range(1<<len(e)):
  if m>=m^full:continue
  u,v=cls[m],cls[m^full]
  if u!=v:
   M.add_edge(u,v)
   if all((m>>i&1)==(m>>j&1) for i,j in enumerate(refl)):B.add_edge(u,v)
   else:K.add_edge(u,v)
 return M,B,K,masks
def witness(G):
 k5=clique(G,set(G),5)
 k33=None
 adj={x:set(G[x]) for x in G}
 for a,b,c in combinations(G,3):
  if len(V:=adj[a]&adj[b]&adj[c])>=3:
   k33=((a,b,c),tuple(sorted(V)[:3]));break
 return k5,k33
def pet(g,ind):
 P=list(combinations(range(5),2)); A=[{j for j,z in enumerate(P) if not set(w)&set(z)} for w in P]
 U={v for v in g if len(g[v])>=3};F={}
 def extend(U):
  if len(F)==10:return [F[i] for i in range(10)]
  u=max((u for u in range(10) if u not in F),key=lambda u:len(A[u]&F.keys()))
  C=U.copy()
  for v,x in F.items():
   if v in A[u]:C.intersection_update(g[x])
   elif ind:C.difference_update(g[x])
  for x in C:
   F[u]=x
   if (r:=extend(U-{x})) is not None:return r
   del F[u]
  return None
 return extend(U)

if __name__=='__main__':
 for n in range(4,8):
  m,b,k,_=atlas(n)
  for name,g in [('all',m),('blue',b),('black',k)]:
   print(n,name,len(g),g.number_of_edges(),witness(g))
   print('petersen subgraph',pet(g,False),flush=True)
   print('petersen induced',pet(g,True),flush=True)
