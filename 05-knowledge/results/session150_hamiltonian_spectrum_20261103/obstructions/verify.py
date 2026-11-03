from itertools import combinations, permutations

def graph(n,mask):
 a=[[0]*n for _ in range(n)]
 for i,(u,v) in enumerate(combinations(range(n),2)):
  if mask>>i&1:u,v=v,u
  a[u][v]=1
 return a

def count(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if d[s][v]:
    for w in range(n):
     if not (s>>w)&1 and a[v][w]:d[s|1<<w][w]+=d[s][v]
 return sum(d[-1])

def reach(a):
 n=len(a);r=[row[:] for row in a]
 for i in range(n):r[i][i]=1
 for k in range(n):
  for i in range(n):
   for j in range(n):r[i][j] |= r[i][k] and r[k][j]
 return r

def sub(q,v,t):
 m=len(q);n=len(t); old=[u for u in range(m) if u!=v];a=[[0]*(n+m-1) for _ in range(n+m-1)]
 for i in range(n):
  for j in range(n):a[i][j]=t[i][j]
 for i,u in enumerate(old,n):
  for j,w in enumerate(old,n):a[i][j]=q[u][w]
  for x in range(n):a[i][x]=q[u][v];a[x][i]=q[v][u]
 return a

single=strong=cyclechecks=0
for m in range(1,5):
 for mask in range(1<<(m*(m-1)//2)):
  q=graph(m,mask);r=reach(q);hq=count(q)
  for v in range(m):
   component=[u for u in range(m) if r[v][u] and r[u][v]]
   if len(component)==1:
    for n in range(1,4):
     for tm in range(1<<(n*(n-1)//2)):
      t=graph(n,tm)
      if count(sub(q,v,t))!=hq*count(t): raise AssertionError((m,mask,v,n,tm))
      single+=1
   else:
    rest=[u for u in component if u!=v]
    if not any(all(q[a][b] for a,b in zip((v,)+p,p+(v,))) for p in permutations(rest)): raise AssertionError((m,mask,v,"cycle"))
    cyclechecks+=1
    for n in range(2,6):
     t=graph(n,0)
     if count(sub(q,v,t))<2**n-2: raise AssertionError((m,mask,v,n,"growth"))
     strong+=1
print('singleton-product assertions:',single)
print('strong-component cycle witnesses:',cyclechecks)
print('strong-component exponential lower-bound assertions:',strong)
