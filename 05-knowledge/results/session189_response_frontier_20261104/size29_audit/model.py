from pathlib import Path
import itertools,json
D=Path(__file__).resolve().parent;P=D/'input'
FIVE_A=[4,1,2,7,15];FIVE_B=[16,1,3,7,14];FOURS=[[2,4,1,7],[10,12,9,0],[6,4,8,3]];CYCLE=[2,4,1]
def need(x,m):
 if not x:raise RuntimeError(m)
def profile(r):
 n=len(r);q=[0]*n
 from math import comb
 for p in itertools.permutations(range(n)):
  bad=sum(not(r[p[i]]>>p[i+1]&1) for i in range(n-1))
  for h in range(bad+1,n+1):q[h-1]+=comb(n-1-bad,h-1-bad)
 return q

def inflate(q,locals,a,c):
 k=len(q);cells=[];extra=3*k
 for i in range(k):
  cell=list(range(3*i,3*i+3))
  for _ in range(len(locals.get(i,CYCLE))-3):cell.append(extra);extra+=1
  cells.append(cell)
 rows=[0]*(extra+2);x=extra;y=extra+1
 for i,cell in enumerate(cells):
  local=locals.get(i,CYCLE)
  for u,v in enumerate(cell):
   for h,w in enumerate(cell):
    if local[u]>>h&1:rows[v]|=1<<w
   for j,other in enumerate(cells):
    if q[i]>>j&1:rows[v]|=sum(1<<w for w in other)
   for port,mask in [(x,a),(y,c)]:
    if mask>>i&1:rows[v]|=1<<port
    else:rows[port]|=1<<v
 rows[x]|=1<<y
 return rows,cells

def ambient(q,a,c):
 rows=[0]*31
 for i in range(9):
  for v in range(3*i,3*i+3):
   rows[v]|=1<<(3*i+(v-3*i+1)%3)
   for j in range(9):
    if q[i]>>j&1:rows[v]|=7<<(3*j)
   rows[v]|=1<<27;rows[28]|=1<<v
   for port,mask in [(29,a),(30,c)]:
    if mask>>i&1:rows[v]|=1<<port
    else:rows[port]|=1<<v
 rows[27]|=1<<28;rows[29]|=(1<<30)|(1<<27);rows[30]|=1<<27;rows[28]|=(1<<29)|(1<<30)
 return rows,[list(range(3*i,3*i+3)) for i in range(9)]+[[27],[28]]

def compress(rows,cells):
 # Only disjoint physical directed triangles which are genuine modules are
 # compressed. Every other physical vertex remains individually named.
 tri=[];single=[]
 for cell in cells:
  remaining=set(cell)
  for subset in itertools.combinations(cell,3):
   if not set(subset)<=remaining:continue
   if not all(sum(rows[v]>>w&1 for w in subset)==1 for v in subset):continue
   if not all(len({rows[w]>>v&1 for v in subset})==1 for w in cell if w not in subset):continue
   tri.append(list(subset));remaining-=set(subset)
  single += [[v] for v in sorted(remaining)]
 blocks=tri+single+[[len(rows)-2],[len(rows)-1]];out=[0]*len(blocks)
 need(sorted(v for b in blocks for v in b)==list(range(len(rows))),'physical partition')
 for i,part in enumerate(blocks):
  for j,other in enumerate(blocks):
   if i==j:continue
   edge={rows[v]>>w&1 for v in part for w in other};need(len(edge)==1,'rotation block has uniform exterior')
   if edge=={1}:out[i]|=1<<j
 return len(tri),out
