import itertools as it
from collections import Counter
import sympy as s
def isomaps(A,B):
 da=[sum(r) for r in A]; db=[sum(r) for r in B]
 order=sorted(range(8),key=lambda i:sum(d==da[i] for d in db))
 def rec(p,used):
  if len(p)==8:yield p.copy();return
  i=order[len(p)]
  for j in range(8):
   if j in used or da[i]!=db[j]:continue
   if all(A[i][u]==B[j][v] and A[u][i]==B[v][j] for u,v in p.items()):
    p[i]=j;yield from rec(p,used|{j});del p[i]
 yield from rec({},set())
# Alternate rank certificate: restriction to sum-zero, all 1033 odd tournaments n<=5.
count=0
for n in [1,3,5]:
 B=s.zeros(n,n-1)
 for i in range(n-1):B[i,i]=1;B[n-1,i]=-1
 edges=list(it.combinations(range(n),2))
 for bits in it.product([-1,1],repeat=len(edges)):
  M=s.zeros(n)
  for (i,j),v in zip(edges,bits):M[i,j]=v;M[j,i]=-v
  R=B.T*M*B
  assert R.det()%2==1
  count+=1
print('EXACT all odd tournaments n=1,3,5 sum-zero restriction determinants odd:',count)
# Invariant Cayley masks from direct Boolean constancy constraints on the28 ± pairs.
sets=[]
# Independent count using orbits on signed assignments via direct closure of each residue.
orbits=[];left=set(range(1,57))
while left:
 a=left.pop();o={a};b=7*a%57
 while b!=a:o.add(b);b=7*b%57
 left-=o;orbits.append(o)
assert Counter(map(len,orbits))=={1:2,3:18}
opts=[];used=set()
for o in orbits:
 if min(o) in used:continue
 neg={57-a for a in o};used|=o|neg;opts.append((o,neg))
for choices in it.product([0,1],repeat=10):
 S=set().union(*(opts[i][c] for i,c in enumerate(choices)))
 assert len(S)==28 and S|{57-a for a in S}==set(range(1,57))
 assert {7*a%57 for a in S}==S
 sets.append(tuple(sorted(S)))
assert len(set(sets))==1024
print('EXACT distinct labelled Frobenius-invariant57 Cayley rules:1024')
# Conference: independent constraint-backtracking isomorphisms, brute Hamilton paths, anchored DP cycles.
n=8;C=s.zeros(n)
for i in range(n):
 for j in range(n):
  if i!=j:C[i,j]=(1 if i==0 else -1 if j==0 else 1 if (j-i)%7 in (1,2,4) else -1)
assert C*C==-7*s.eye(8)
x=s.symbols('x');target=x*(x-3)*(x*x+x+2)**3
reps={};sizes=Counter();sources=[];sinks=[]
for mask in range(128):
 d=[1]+[1-2*((mask>>i)&1) for i in range(7)];D=s.diag(*d);M=D*C*D
 assert s.eye(n)-M*M==8*s.eye(n)
 A=(s.ones(n)-s.eye(n)+M)/2
 assert s.expand(A.charpoly(x).as_expr()-target)==0
 G=[[bool(A[i,j]) for j in range(n)] for i in range(n)]
 deg=tuple(sorted(map(sum,G)))
 if deg not in reps:reps[deg]=(mask,G)
 else:assert next(isomaps(G,reps[deg][1]),None) is not None
 sizes[deg]+=1
 if 7 in deg:sources.append(mask)
 if 0 in deg:sinks.append(mask)
assert sorted(sizes.values())==[8,8,56,56]
for deg,(mask,G) in reps.items():
 A=G
 H=sum(all(A[p[k]][p[k+1]] for k in range(7)) for p in it.permutations(range(n)))
 cyc={}
 for size in [3,5,7]:
  for vs in it.combinations(range(n),size):
   start=min(vs);remaining=set(vs)-{start}
   dp={(frozenset([start]),start):1}
   for r in range(1,size):
    nd={}
    for (seen,last),value in dp.items():
     for j in set(vs)-seen:
      if A[last][j]:key=(seen|{j},j);nd[key]=nd.get(key,0)+value
    dp=nd
   cyc[frozenset(vs)]=sum(value for (seen,last),value in dp.items() if A[last][start])
 c=[sum(v for k,v in cyc.items() if len(k)==r) for r in [3,5,7]]
 tri=[(k,v) for k,v in cyc.items() if len(k)==3]
 c33=sum(v*w for (k,v),(l,w) in it.combinations(tri,2) if not k&l)
 c35=sum(v*w for k,v in tri for l,w in cyc.items() if len(l)==5 and not k&l)
 assert H==1+2*sum(c)+4*(c33+c35)
 expected=(14,42,24,7,0) if H==189 else (14,42,48,19,12)
 assert tuple(c+[c33,c35])==expected
 aut=sum(1 for _ in isomaps(G,G))
 assert aut==(21 if H==189 else 3)
 print('EXACT representative',mask,'size',sizes[deg],'H',H,'cycles',c+[c33,c35],'Aut',aut)
for code in [set(sources),{127^m for m in sinks}]:
 assert len(code)==8 and all(a^b in code for a in code for b in code)
 lines=[{i for i in range(7) if not(m>>i&1)} for m in code if m]
 assert len(lines)==7 and all(len(L)==3 for L in lines)
 assert Counter(tuple(sorted(pair)) for L in lines for pair in it.combinations(L,2))==Counter(it.combinations(range(7),2))
print('EXACT source simplex:',sources,'complemented sink simplex:',sorted(127^m for m in sinks))
assert set(sources)!={127^m for m in sinks}
print('ALL NEW-ADDITION INDEPENDENT CHECKS PASS')
