from collections import Counter
from itertools import permutations
import json, random
BULK=-1

def canon(parts):return tuple(sorted(parts))
def introduce(states,u,A):
 out=Counter()
 for parts,num in states.items():
  left=[-1]+[i for i,(s,e)in enumerate(parts)if e==BULK or e not in A]
  right=[-1]+[i for i,(s,e)in enumerate(parts)if s!=BULK and s in A]
  for i in left:
   for j in right:
    if i==j and i!=-1:continue
    new=[part for z,part in enumerate(parts)if z!=i and z!=j]
    new.append((parts[i][0]if i!=-1 else u,parts[j][1]if j!=-1 else u))
    out[canon(new)]+=num
 return out

def forget(states,u,p):
 out=Counter()
 for parts,num in states.items():
  new=tuple((BULK if s==u else s,BULK if e==u else e)for s,e in parts)
  if sum(s==BULK for s,e in new)<=p:out[canon(new)]+=num
 return out

def read(states,r):
 h=s=L=T=C=0
 for parts,num in states.items():
  q=len(parts)
  if q==1:
   h+=num
   if parts[0][0]==r:s+=num
  if q in(1,2)and any(e==r for a,e in parts):L+=num
  if q<=2:C+=num*(2 if q in(1,2)else 1)
  if q in(1,2,3)and any(a==r for a,e in parts):T+=num*([1,2,2][q-1])
 return h,s,L,h+L-s,C,T

def literal(rows,r):
 n=len(rows);M=1<<n;F=M-1
 dp=[[0]*n for _ in range(M)];rp=[[0]*n for _ in range(M)]
 for i in range(n):dp[1<<i][i]=rp[1<<i][i]=1
 H=[0]*M;H[0]=1
 for mask in range(1,M):
  for i in range(n):
   if mask>>i&1:
    for j in range(n):
     if not(mask>>j&1):
      if rows[i]>>j&1:dp[mask|1<<j][j]+=dp[mask][i]
      if rows[j]>>i&1:rp[mask|1<<j][j]+=rp[mask][i]
  H[mask]=sum(dp[mask])
 L=sum(dp[mask][r]*H[F^mask]for mask in range(M)if mask>>r&1)
 C=sum(H[mask]*H[F^mask]for mask in range(M))
 BB=[0]*M
 for mask in range(M):
  sub=mask
  while True:
   BB[mask]+=H[sub]*H[mask^sub]
   if not sub:break
   sub=(sub-1)&mask
 T=sum(rp[mask][r]*BB[F^mask]for mask in range(M)if mask>>r&1)
 return H[F],rp[F][r],L,H[F]+L-rp[F][r],C,T

def full(rows,p=3):
 states=Counter({():1})
 for u in range(len(rows)):
  states=introduce(states,u,{v for v in range(u)if rows[u]>>v&1})
 for u in range(1,len(rows)):states=forget(states,u,p)
 return states

def need(c,m):
 if not c:raise RuntimeError(m)
checks=0
for n in range(1,6):
 for bits in range(1<<(n*(n-1)//2)):
  rows=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if bits>>k&1:rows[i]|=1<<j
    else:rows[j]|=1<<i
    k+=1
  states=full(rows)
  need(read(states,0)==literal(rows,0),'complete graph transfer')
  checks+=1
random.seed(203)
random_checks=0;max_patterns=0
for W in (2,3,4):
 for trial in range(30):
  p=3;states=Counter({():1});active={};rows=[];root=None
  for vertex in range(9):
   if len(active)==W:
    u=random.choice([u for u in active if u!=root]);states=forget(states,u,p);del active[u]
   u=min(set(range(W))-set(active))
   A={z for z in active if random.randrange(2)}
   row=0
   for z in range(vertex):
    if z in {active[j]for j in A}:row|=1<<z
    else:rows[z]|=1<<vertex
   rows.append(row);states=introduce(states,u,A);active[u]=vertex
   if root is None:root=u
   max_patterns=max(max_patterns,len(states))
  for u in list(active):
   if u!=root:states=forget(states,u,p);del active[u]
  need(read(states,root)==literal(rows,active[root]),'reused register transfer')
  random_checks+=1
print(json.dumps({'status':'PASS','all_labeled_graphs_through5':checks,'random_width_2_3_4_order9_words':random_checks,'largest_reached_pattern_dictionary':max_patterns,'readouts':['h','s','L','B','C2','T_root_start_2ordinary'],'scope':'Finite controls; general bounded-width proof separate.'},indent=2))
