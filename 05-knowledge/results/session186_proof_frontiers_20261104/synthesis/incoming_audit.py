"""Independent literal checks for the incoming direct-query addendum."""
from itertools import combinations,permutations
import json

def graph(n,m):
 a=[[False]*n for _ in range(n)]
 for k,(u,v) in enumerate(combinations(range(n),2)):
  a[u][v]=bool(m>>k&1);a[v][u]=not a[u][v]
 return a

def count(a):
 return sum(all(a[u][v] for u,v in zip(p,p[1:])) for p in permutations(range(len(a))))
def ext(a,v):
 n=len(a);b=[r+[False,False] for r in a]+[[False]*(n+2) for _ in range(2)]
 b[n][n+1]=True
 for i in range(n):b[i][n]=(i==v);b[n][i]=(i!=v);b[n+1][i]=True
 return b
p1=p2=0
for n in range(1,5):
 for m in range(1<<(n*(n-1)//2)):
  a=graph(n,m);h=count(a);e=[0]*n;c=[[0]*n for _ in range(n)]
  for p in permutations(range(n)):
   if all(a[u][v] for u,v in zip(p,p[1:])):e[p[-1]]+=1
   for k in range(1,n):
    if all(a[p[t]][p[t+1]] for t in range(n-1) if t!=k-1):c[p[k-1]][p[k]]+=1
  f=[count(ext(a,v)) for v in range(n)]
  for v in range(n):
   if f[v]!=h+2*(e[v]+sum(c[v])):raise RuntimeError('P1')
   p1+=1
  for u,v in combinations(range(n),2):
   if not all(a[u][w]==a[v][w] for w in range(n) if w not in (u,v)):continue
   if not a[u][v]:u,v=v,u
   ids=[w for w in range(n) if w!=v];q=[[a[x][y] for y in ids] for x in ids];z=ids.index(u)
   eq=0;row=0
   for p in permutations(range(n-1)):
    if p[-1]==z and all(q[x][y] for x,y in zip(p,p[1:])):eq+=1
    for k in range(1,n-1):
     if p[k-1]==z and all(q[p[t]][p[t+1]] for t in range(n-2) if t!=k-1):row+=1
   if f[v]-f[u]!=2*(eq+row) or f[v]-f[u]<2*count(q):raise RuntimeError('P2')
   p2+=1
seeds=((16,5,33,39,46,3),(16,5,49,39,42,3),(18,4,49,39,42,3));word=[5,10,6,11,7,0,4,3,8,1,2,9]
for s in seeds:
 a=[[False]*12 for _ in range(12)]
 for u,v in combinations(range(10),2):
  old_u=u if u<5 else 5;old_v=v if v<5 else 5
  a[u][v]=(u<v) if old_u==old_v else bool(s[old_u]>>old_v&1);a[v][u]=not a[u][v]
 for v in range(10):
  for port,mask in [(10,32),(11,64)]:a[v][port]=bool(mask>>v&1);a[port][v]=not a[v][port]
 a[10][11]=True
 if sorted(word)!=list(range(12)) or not all(a[u][v] for u,v in zip(word,word[1:])):raise RuntimeError('word')
 runs=sum(5<=v<10 and (i==0 or not(5<=word[i-1]<10)) for i,v in enumerate(word))
 if runs!=5:raise RuntimeError('runs')
print(json.dumps({'status':'PASS','literal_P1_vertices':p1,'literal_twin_pairs':p2,'five_run_words':3,'compiled_count_arithmetic':6*(841+6241)+72+72+32},sort_keys=True,indent=2))
