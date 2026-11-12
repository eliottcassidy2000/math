from collections import deque
class Flow:
 def __init__(self,n):self.g=[[] for _ in range(n)]
 def add(self,a,b,c):
  self.g[a].append([b,c,len(self.g[b])]);self.g[b].append([a,0,len(self.g[a])-1])
 def cut(self,s,t):
  total=0
  while True:
   level=[-1]*len(self.g);level[s]=0;q=deque([s])
   while q:
    a=q.popleft()
    for b,c,r in self.g[a]:
     if c and level[b]<0:level[b]=level[a]+1;q.append(b)
   if level[t]<0:return total,{i for i,d in enumerate(level) if d>=0}
   it=[0]*len(self.g)
   def dfs(a,f):
    if a==t:return f
    while it[a]<len(self.g[a]):
     e=self.g[a][it[a]];b,c,r=e
     if c and level[b]==level[a]+1:
      z=dfs(b,min(f,c))
      if z:e[1]-=z;self.g[b][r][1]+=z;return z
     it[a]+=1
    return 0
   while True:
    z=dfs(s,sum(e[1] for e in self.g[s]))
    if not z:break
    total+=z
