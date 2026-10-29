"""Independent generator closure certificate for order-eight switching groups."""
from itertools import combinations,permutations
from collections import Counter
import json

def need(x,msg):
 if not x:raise RuntimeError(msg)
def comp(a,b):return tuple(a[b[i]] for i in range(len(a)))
def closure(gens,n):
 G={tuple(range(n))};front=list(G)
 while front:
  p=front.pop()
  for g in gens:
   h=comp(g,p)
   if h not in G:G.add(h);front.append(h)
 return G
def mapsign(a,p,e):
 d=[1]+[e*a[p[0]][p[i]]*a[0][i] for i in range(1,len(a))]
 return all(a[p[i]][p[j]]==e*d[i]*d[j]*a[i][j] for i,j in combinations(range(len(a)),2))
q=7;n=8;sq={1,2,4};a=[[0]*n for _ in range(n)]
for i in range(1,n):a[0][i]=1;a[i][0]=-1
for i,j in combinations(range(1,n),2):a[i][j]=1 if (j-i)%q in sq else -1;a[j][i]=-a[i][j]
# Index 0=infinity, index x+1=x. Elementary PSL generators x+1 and -1/x.
t=(0,)+tuple(1+(x+1)%q for x in range(q))
s=(1,0)+tuple(1+(-pow(x,-1,q))%q for x in range(1,q))
scale=(0,)+tuple(1+(3*x)%q for x in range(q))
G=closure([t,s],n);E=closure([t,s,scale],n)
need(len(G)==168 and len(E)==336,'generator closure sizes')
pos={p for p in permutations(range(n)) if mapsign(a,p,1)}
neg={p for p in permutations(range(n)) if mapsign(a,p,-1)}
need(pos==G and pos|neg==E and not pos&neg,'full signed switching identification')
# Independently count conjugacy classes via centralizers, then a union-size test.
counts=Counter(sum(comp(g,h)==comp(h,g) for h in G) for g in G)
classes=[]
for centralizer,total in counts.items():
 size=len(G)//centralizer;need(total%size==0,'class multiplicity');classes.extend([size]*(total//size))
classes.sort();need(classes==[1,21,24,24,42,56],'class size certificate')
possible={1}
for size in classes[1:]:possible|={s+size for s in list(possible)}
normal=[s for s in sorted(possible) if len(G)%s==0];need(normal==[1,168],'simplicity by normal class-union orders')
print(json.dumps({'status':'FINITE-EXACT','positive_group':168,'extended_group':336,'generator_model':'translation, negative inversion, nonsquare dilation','conjugacy_class_sizes':classes,'possible_normal_orders':normal},sort_keys=True,indent=2))
