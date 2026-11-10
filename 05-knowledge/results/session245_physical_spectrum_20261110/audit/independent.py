from itertools import permutations,product
from math import gcd
import json

def check(v,m):
 if not v:raise ValueError(m)
def paths(a):
 n=len(a);out=[[] for _ in range(1<<n)];out[0]=[()]
 for r in range(1,n+1):
  for p in permutations(range(n),r):
   if all(a[x]>>y&1 for x,y in zip(p,p[1:])):out[sum(1<<x for x in p)].append(p)
 return out
def inter(a,c):
 ps=paths(a);f=len(ps)-1
 h=len(ps[f]);q=A=G=T=0
 for m in range(f+1):
  for p in ps[m]:
   for r in ps[f^m]:
    q+=1;u=bool(r and c>>r[0]&1);v=bool(p and c>>p[-1]&1)
    A+=u;G+=v;T+=u*v
 return h,q,A,G,T,sum(bool(p and c>>p[0]&1) for p in ps[f]),sum(bool(p and c>>p[-1]&1) for p in ps[f])
def closure(a,c):
 n=len(a);return [x|(1<<n)|((1<<(n+1)) if not(c>>i&1) else 0) for i,x in enumerate(a)]+[1<<(n+1),c]
def strong(a):
 for i in range(len(a)):
  seen={i}
  while True:
   new=seen|{j for x in seen for j in range(len(a)) if a[x]>>j&1}
   if new==seen:break
   seen=new
  if len(seen)!=len(a):return False
 return True
rows=[]
for a,cs,expected in [([18,8,3,5,14],(3,6),[(13,104,32,41,6,5,6),(13,104,32,44,8,4,6)]),([32,17,3,7,13,30],(17,33),[(37,352,96,102,12,12,10),(37,352,96,96,32,10,10)])]:
 for c,e in zip(cs,expected):
  z=inter(a,c);check(z==e,'interface');g=closure(a,c);p=paths(g);f=len(p)-1;r=len(g)-1
  U=sum(len(p[m])*sum(bool(t and t[0]==r) for t in p[f^m]) for m in range(f+1))
  D=len(p[f^(1<<r)]);E=sum(t[-1]==r for t in p[f])
  check((len(p[f]),E,U,D)==(z[0]+2*z[2]-z[4],z[0],z[0]+2*z[2],z[0]),'closure sidecars')
  rows.append(dict(core=a,cut=c,interface=z,closed=g,state=[len(p[f]),E,U,D],strong=strong(g)))
# Independent arbitrary-signature gap compiler, including k=0 boundary.
def graph(k,arcs,w):
 a=arcs[:]+[0]*len(w)
 for i,sig in enumerate(w,k):
  for j in range(k):
   if sig>>j&1:a[i]|=1<<j
   else:a[j]|=1<<i
  for j in range(i+1,len(a)):a[i]|=1<<j
 return a
def count(k,arcs,w):
 total=0
 a=graph(k,arcs,w)
 for p in permutations(range(k)):
  for labels in product(range(k+1),repeat=len(w)):
   path=[]
   for gap in range(k+1):
    path.extend(k+i for i,x in enumerate(labels) if x==gap)
    if gap<k:path.append(p[gap])
   total+=all(a[x]>>y&1 for x,y in zip(path,path[1:]))
 return total
cases=0
for k in range(4):
 pairs=list((i,j) for i in range(k) for j in range(i+1,k))
 for code in range(1<<len(pairs)):
  arcs=[0]*k
  for b,(i,j) in enumerate(pairs):
   x,y=(i,j) if code>>b&1 else (j,i);arcs[x]|=1<<y
  for w in product(range(1<<k),repeat=2):
   a=graph(k,arcs,w);check(count(k,arcs,w)==len(paths(a)[-1]),'arbitrary apex exactness');cases+=1
  # Full annihilator recurrence via counts for each possible letter and a fixed context.
  for letter in range(1<<k):
   vals=[count(k,arcs,(0,)+(letter,)*m+((1<<k)-1,)) for m in range(k+3)]
   coeff=[1]
   for root in range(k+2):
    nxt=[0]*(len(coeff)+1)
    for i,x in enumerate(coeff):nxt[i]-=root*x;nxt[i+1]+=x
    coeff=nxt
   check(sum(x*y for x,y in zip(coeff,vals))==0,'annihilator arbitrary signature')
print(json.dumps(dict(status='PASS',literal_cells=rows,arbitrary_apex_cases=cases),indent=2))
