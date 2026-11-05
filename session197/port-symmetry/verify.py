#!/usr/bin/env python3
from itertools import permutations,product
import json

def require(v,msg):
 if not v:raise RuntimeError(msg)

def tournament():
 rows=[0]*9
 for a in range(9):
  for b in range(9):
   A,i=divmod(a,3);B,j=divmod(b,3)
   if (A==B and (j-i)%3==1) or (A!=B and (B-A)%3==1):rows[a]|=1<<b
 return rows

def auto(rows,p):
 return all(bool(rows[a]>>b&1)==bool(rows[p[a]]>>p[b]&1) for a in range(len(rows)) for b in range(len(rows)))

def rotation(k,n):
 p=list(range(n))
 for j in range(3):p[3*k+j]=3*k+(j+1)%3
 return p

def birth(rows,ports,side,S,drop):
 n=len(rows);chosen={ports[i] for i in range(2) if S>>i&1};out=rows[:]+[0]
 for v in range(n):
  win=(v in chosen) if side=='R' else (v not in chosen)
  if win:out[n]|=1<<v
  else:out[v]|=1<<n
 newports=ports+[n];del newports[drop];return out,newports

def count_h(rows):
 n=len(rows);D=[[0]*n for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1:
    rem=mask^(1<<v)
    if rem:D[mask][v]=sum(D[rem][u] for u in range(n) if rem>>u&1 and rows[u]>>v&1)
 return sum(D[-1])

def strong(rows):
 n=len(rows)
 for r in range(n):
  seen=1<<r;last=-1
  while seen!=last:
   last=seen
   for v in range(n):
    if seen>>v&1:seen|=rows[v]
  if seen!=(1<<n)-1:return False
 return True

T=tournament();require(T==[58,60,57,464,480,456,135,263,71],'explicit lexicographic graph')
require(strong(T),'strong input');h=count_h(T)
literal=sum(all(T[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(9)))
require(h==literal==3159,'independent H3159 count')
for k in range(3):
 p=rotation(k,9);require(auto(T,p),'module rotation');require([p[p[p[v]]] for v in range(9)]==list(range(9)),'rotation order3')
first=second=0;counted=0;ops=list(product(('R','L'),range(4),range(3)))
for ports in permutations(range(9),2):
 ports=list(ports);empty=[k for k in range(3) if all(v//3!=k for v in ports)]
 require(empty,'unexposed module')
 for side,S,d in ops:
  out,ps=birth(T,ports,side,S,d)
  for k in empty:
   p=rotation(k,10);require(auto(out,p) and all(p[v]==v for v in ps),'one-step persistent stabilizer')
  first+=1
  if tuple(ports) in ((0,1),(0,3),(0,6)):
   require(count_h(out)%3**len(empty)==0,'subgroup divisibility control');counted+=1
  # Every next side, mask and forgetting choice is checked at adjacency level.
  k=empty[0]
  for side2,S2,d2 in ops:
   final,qs=birth(out,ps,side2,S2,d2);p=rotation(k,11)
   require(auto(final,p) and all(p[v]==v for v in qs),'two-step persistent automorphism');second+=1
# A concrete contrast: the prior H49 core has no nonidentity automorphism.
U=[122,112,67,102,108,4,32]
auts=sum(auto(U,p) for p in permutations(range(7)));require(auts==1,'rigid H49 contrast')
print(json.dumps({'rows':T,'H':h,'strong':True,'ordered_initial_port_pairs':72,'first_birth_checks':first,'second_birth_checks':second,'Hamilton_divisibility_controls':counted,'initial_same_module_divisor':9,'initial_distinct_module_divisor':3,'H49_automorphism_count':auts},sort_keys=True,indent=2))
