"""Independent exact square-zero ring matrices, no producer imports."""
from itertools import permutations
import json

def require(c,label):
 if not c:raise RuntimeError(label)
def add(a,b):
 c=a.copy()
 for m,v in b.items():
  c[m]=c.get(m,0)+v
  if not c[m]:del c[m]
 return c
def neg(a):return {m:-v for m,v in a.items()}
def mul(a,b):
 c={}
 for m,v in a.items():
  for n,w in b.items():
   if not m&n:c[m|n]=c.get(m|n,0)+v*w
 return {m:v for m,v in c.items() if v}
def zero(n,m=None):return [[{} for j in range(n if m is None else m)] for i in range(n)]
def eye(n):return [[{0:1} if i==j else {} for j in range(n)] for i in range(n)]
def plus(A,B):return [[add(a,b) for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def minus(A,B):return plus(A,[[neg(a) for a in r] for r in B])
def product(A,B):
 C=zero(len(A),len(B[0]))
 for i in range(len(A)):
  for j in range(len(B[0])):
   for k in range(len(B)):C[i][j]=add(C[i][j],mul(A[i][k],B[k][j]))
 return C
def inv_unit(A,variables):
 I=eye(len(A)); N=minus(I,A);P=I;S=I
 for k in range(1,variables+1):P=product(P,N);S=plus(S,P)
 require(product(P,N)==zero(len(A)),'nilpotence')
 require(product(A,S)==I,'inverse certificate')
 return S
def ax(A,offset=0):return [[{1<<(offset+j):a} if a else {} for j,a in enumerate(row)] for row in A]
def paths(A):
 n=len(A);R=inv_unit(minus(eye(n),ax(A)),n)
 return [[mul({1<<i:1},v) for v in row] for i,row in enumerate(R)]
checks=[]
for A in [[[0,1],[0,0]],[[0,1,0],[0,0,1],[1,0,0]],[[0,1,1],[0,0,1],[0,0,0]]]:
 n=len(A);N=2*n
 AD=[[0]*N for _ in range(N)]
 for i in range(n):
  for j in range(n):
   AD[i][j]=A[i][j];AD[i][j+n]=A[i][j]+int(i==j)
   AD[i+n][j]=A[i][j];AD[i+n][j+n]=A[j][i]
 B=minus(eye(N),ax(AD));full=inv_unit(B,N)
 U=[r[:n] for r in B[:n]];V=[r[n:] for r in B[:n]]
 W=[r[:n] for r in B[n:]];Z=[r[n:] for r in B[n:]]
 Zi=inv_unit(Z,N);Si=inv_unit(minus(U,product(product(V,Zi),W)),N)
 reconstructed=[Si,[[neg(v) for v in row] for row in product(product(Si,V),Zi)],[[neg(v) for v in row] for row in product(product(Zi,W),Si)],plus(Zi,product(product(product(product(Zi,W),Si),V),Zi))]
 blocks=[[r[:n] for r in full[:n]],[r[n:] for r in full[:n]],[r[:n] for r in full[n:]],[r[n:] for r in full[n:]]]
 require(reconstructed==blocks,'four Schur inverse blocks')
 P=paths(AD);expected=zero(N)
 for k in range(1,N+1):
  for p in permutations(range(N),k):
   if all(AD[a][b] for a,b in zip(p,p[1:])):
    m=sum(1<<v for v in p); expected[p[0]][p[-1]]=add(expected[p[0]][p[-1]],{m:1})
 require(P==expected,'every labelled port coefficient equals direct simple paths')
 phi={0:1}
 for row in P:
  for entry in row:phi=add(phi,entry)
 ocf={0:1}
 for i in range(N):ocf=mul(ocf,{0:1,1<<i:1})
 for k in range(3,N+1,2):
  for cycle in permutations(range(N),k):
   if cycle[0]!=min(cycle):continue
   if all(AD[cycle[i]][cycle[(i+1)%k]] for i in range(k)):
    ocf=mul(ocf,{0:1,sum(1<<i for i in cycle):2})
 require(phi==ocf,'all squarefree coefficients equal odd-cycle collection expansion')
 checks.append({'odd_cycle_polynomial':True,'old_n':n,'double_n':N,'all_four_schur_blocks':True,'every_port_coefficient':True})
print(json.dumps(checks,indent=2))
