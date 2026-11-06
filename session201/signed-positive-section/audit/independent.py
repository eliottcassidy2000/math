import itertools,json,hashlib
from pathlib import Path

def require(x,msg):
 if not x: raise RuntimeError(msg)
# Independent positive-coordinate transitions, not candidate raw functions.
def nxt(y):
 a,f,e,p=y
 return [(a+e,a+2*f+e+2*p,0,e+p),(2*a+f+e+p,f,e+p,0),(f+p,0,a+e+p,e)]
def raw(y):
 a,f,e,p=y
 return a+2*e+p,e,a+f+2*e+2*p,e+p
cap=1000
seen={(0,0,0,1):''}
levels=[(0,0,0,1)]
for y in levels:
 for c,z in zip('OZJ',nxt(y)):
  if raw(z)[2]<=cap and z not in seen:seen[z]=seen[y]+c;levels.append(z)
S={raw(y):w for y,w in seen.items()};pair={(1,1)}|{(s[0],s[2]) for s in S}
sliced=[]
for s in S:
 D,E,U,d=s;m=U-D;n=6*U-7*D
 if (E,d)!=(2*m,3*m):continue
 if n==m:require((n,n) in pair,'exceptional actual parent')
 else: require(n>0 and m>=3*n and (n,(m-n)//2) in pair,'canonical actual parent')
 sliced.append(s)
cores=sorted(s for s in S if s[2]<=11 and s[1]>0)
require(not any(s[0]==7 for s in cores),'D7 core')
# Matrix multiplication independent of state iteration.
O=[[1,0,0,0],[0,0,0,0],[0,0,2,0],[0,0,0,1]]
Z=[[1,-1,1,0],[0,0,0,1],[0,0,2,-1],[0,0,0,1]]
J=[[1,-1,1,0],[1,-1,0,0],[1,0,1,0],[1,0,0,0]]
Q=[[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]]
def mul(A,B):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*B)] for r in A]
def tr(A):return list(map(list,zip(*A)))
W=O
for i in range(4):W=mul(J,W)
require(W==[[5,0,12,0],[2,0,4,0],[6,0,14,0],[3,0,6,0]],'W matrix')
for A,B in [(O,Z),(Z,O),(J,J)]:require(mul(tr(A),Q)==mul(Q,B),'adjoint')
require(mul(Z,J)==mul(J,O),'chronological JZ=OJ')
# Literal paths by permutations, covers by ordered partition plus permutations.
def literal(word,singleton=False):
 N=1 if singleton else 2;edges=set() if singleton else {(0,1)};root=0
 for c in word:
  for v in range(N):edges.add((N,v) if c!='O' and v==root else (v,N))
  if c=='J':root=N
  N+=1
 def paths(vs):return [p for p in itertools.permutations(vs) if all(e in edges for e in zip(p,p[1:]))]
 vs=set(range(N));full=paths(vs);deletion=paths(vs-{root});U=0
 for mask in range(1<<N):
  part={i for i in range(N) if mask>>i&1}
  if root not in part:continue
  U+=sum(p[0]==root for p in paths(part))*len(paths(vs-part))
 response=(len(full),sum(p[-1]==root for p in full),U,len(deletion))
 start=[p for p in full if p[0]==root]
 cyc=sum((p[-1],root) in edges for p in start)
 return response,len(start),cyc
checks=[]
for word,singleton,parent in [('OJJJJ',True,(1,1)),('OJJJJ',False,(1,2)),('ZOJJJJ',False,(3,3)),('OZOJJJJ',False,(5,7)),('ZJOJJJJ',False,(5,6))]:
 r,st,cy=literal(word,singleton);D,U=parent
 require(r==tuple(sum(a*b for a,b in zip(row,(D,0,U,0))) for row in W),'literal profile')
 require(st==cy==1,'literal starts/cycles')
 checks.append({'word':word,'singleton':singleton,'profile':r,'starts':st,'cycles':cy})
print(json.dumps({'status':'PASS','independent_method':'positive-coordinate orbit and permutation paths/covers','cap':cap,'source_profiles':len(S),'slice_profiles':len(sliced),'slice':sorted(sliced),'core_U11_positive_endpoint':cores,'literal_checks':checks},indent=2))
