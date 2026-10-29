"""Direct permutation enumeration for signed switching symmetries; no producer imports."""
from itertools import combinations, permutations
from collections import Counter
import json

def need(x,msg):
 if not x:raise RuntimeError(msg)
def mat(n,b):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):a[i][j]=1 if b>>k&1 else -1;a[j][i]=-a[i][j]
 return a
def mapsign(a,p,e):
 n=len(a);d=[1]+[e*a[p[0]][p[i]]*a[0][i] for i in range(1,n)]
 return all(a[p[i]][p[j]]==e*d[i]*d[j]*a[i][j] for i,j in combinations(range(n),2))
def v2(x):
 r=1
 while x%2==0:r*=2;x//=2
 return r
summary=[]
for n in range(3,6):
 hist=Counter()
 for b in range(1<<(n*(n-1)//2)):
  a=mat(n,b);gp=[];gm=[]
  for p in permutations(range(n)):
   if mapsign(a,p,1):gp.append(p)
   if mapsign(a,p,-1):gm.append(p)
  need(not set(gp)&set(gm),'sign is not unique n>=3')
  need(n%v2(len(gp))==0,'ordinary Sylow bound')
  need((2*n)%v2(len(gp)+len(gm))==0,'extended Sylow bound')
  norm=[[a[i][j]*(a[0][i] if i else 1)*(a[0][j] if j else 1) for j in range(n)] for i in range(n)]
  for e,G in [(1,gp),(-1,gm)]:
   fixed={p for p in G if p[0]==0}
   direct={(0,)+p for p in permutations(range(1,n)) if all(norm[p[i-1]][p[j-1]]==e*norm[i][j] for i,j in combinations(range(1,n),2))}
   need(fixed==direct,'root stabilizer signed deletion identity')
  hist[(len(gp),len(gp)+len(gm))]+=1
 summary.append({'n':n,'labelled_tournaments':1<<(n*(n-1)//2),'order_histogram':{str(k):v for k,v in sorted(hist.items())}})
# Explicit sharp case: root 0 dominates cyclic triangle 1->2->3->1.
a=[[0,1,1,1],[-1,0,1,-1],[-1,-1,0,1],[-1,1,-1,0]]
g=[p for p in permutations(range(4)) if mapsign(a,p,1)];ge=[p for p in permutations(range(4)) if mapsign(a,p,1) or mapsign(a,p,-1)]
even=[p for p in permutations(range(4)) if sum(p[i]>p[j] for i,j in combinations(range(4),2))%2==0]
need(g==even and len(ge)==24,'A4/S4 exact sharp witness')
print(json.dumps({'status':'FINITE-EXACT','universe':summary,'sharp_four_vertex_groups':['A4','S4']},indent=2,sort_keys=True))
