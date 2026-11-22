from itertools import combinations,permutations
from sys import argv

def reads(a):
    n=len(a); full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    for i in range(n):dp[1<<i][i]=1
    for S in range(1,full):
        for i,v in enumerate(dp[S]):
            if not v:continue
            r=a[i]&~S
            while r:
                j=(r&-r).bit_length()-1;r&=r-1
                dp[S|1<<j][j]+=v
    p=[sum(r) for r in dp];p[0]=1
    return p[-1],sum(p[S]*p[full^S] for S in range(full+1))

def strong(a):
    n=len(a);full=(1<<n)-1
    for i in range(n):
        S=1<<i
        while True:
            T=S
            for j in range(n):
                if S>>j&1:T|=a[j]
            if S==T:break
            S=T
        if S!=full:return False
    return True
n=int(argv[1]) if len(argv)>1 else 6
E=list(combinations(range(n),2)); N=1<<len(E)
G=[-1]*N; R=[]; reps=[]
def rows(m):
 a=[0]*n
 for t,(i,j) in enumerate(E):
  if m>>t&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return tuple(a)
def betteredge(m,q,i,j):
 return (m,q,(i,j))

for m in range(N):
 if G[m]!=-1:continue
 a=rows(m); c=len(R)
 R.append(a);reps.append(m)
 for pi in permutations(range(n)):
  G[sum((1<<t) for t,(i,j) in enumerate(E) if a[pi[i]]>>pi[j]&1)]=c
Rst={i for i,a in enumerate(R) if strong(a)}
adj=[set() for _ in R]; cert={}
for m in range(N):
 if G[m] not in Rst:continue
 for t,(i,j) in enumerate(E):
  q=m^(1<<t)
  if G[q] in Rst and G[m]!=G[q]:
   adj[G[m]].add(G[q]);cert.setdefault((G[m],G[q]),(m,q,(i,j)))
print('strong-classes',len(Rst))
D={(0,1),(0,2),(1,2),(1,3),(2,3),(0,4),(0,5),(4,5),(4,6),(5,6),(3,6)}
h=[reads(a)[0] for a in R]
best=[-1,None,0]
def match(cur,unused,units,tolls):
 if units+(11-units-tolls)<=best[0]:return
 if not unused:
  best[0]=units;best[1]=cur
  print("newbest", best[0], [reps[cur[v]] for v in range(7)], [h[cur[v]] for v in range(7)], flush=True)
  return
 v=max(unused,key=lambda v:(sum(tuple(sorted((u,v))) in D for u in cur),sum(v in x for x in D)))
 for q in sorted(Rst-set(cur.values())):
  if not all((a in adj[q])==(tuple(sorted((u,v))) in D) for u,a in cur.items()):continue
  new=[abs(h[a]-h[q])==2 for u,a in cur.items() if a in adj[q]]
  match({**cur,v:q},unused-{v},units+sum(new),tolls+len(new)-sum(new))
match({},set(range(7)),0,0)
r=best[1]
print('spindle',r)
if r:
 print('masks',[reps[r[v]] for v in range(7)])
 print('heights',[reads(R[r[v]])[0] for v in range(7)])
 print('witnesses',[(i,j,cert[r[i],r[j]]) for i,j in sorted(D)])
