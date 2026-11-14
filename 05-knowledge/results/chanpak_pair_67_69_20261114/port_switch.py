from itertools import combinations, product
from collections import deque

BANK={
  1:(2,'0010',1,0),
  115:(5,'0001010000110000110011110',4,4),
  163:(5,'0010110000010001110001110',0,4),
  251:(5,'0000110010110001010001110',4,4),
  256:(5,'0000110010110001010001110',1,4),
  258:(5,'0000110010110001010001110',0,1),
  259:(5,'0001110001110000110000110',0,1),
}
ROWS={258:(31,46,36),251:(33,42,36),259:(33,46,34),163:(29,48,40),256:(33,48,34),115:(37,42,32)}
PAIRS={(258,251):72,(259,163):73,(256,115):74}

def make(n,code):
 return [[int(code[n*i+j]) for j in range(n)]for i in range(n)]
def chain(A,u,i):
 m,s,w,v=BANK[i]
 B=make(m,s);n=len(A)
 C=[r+[1]*m for r in A]+[[0]*n+r for r in B]
 C[u][n+w]=0;C[n+w][u]=1
 return C,n+v
def profile(A):
 n=len(A); D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1 and (m:=mask^(1<<v)):
    D[mask][v]=sum(D[m][u]*A[u][v] for u in range(n))
 return [1]+[sum(row) for row in D[1:]]
def pointed(A,r):
 n=len(A); F=(1<<n)-1; H=profile(A);D=[[0]*n for _ in range(1<<n)];D[1<<r][r]=1
 for m in range(1<<n):
  for v in range(n):
   if D[m][v]:
    for w in range(n):
     if not m>>w&1 and A[v][w]:
      D[m|(1<<w)][w]+=D[m][v]
 x=[0,0]
 for m,d in enumerate(D): x[m.bit_count()%2]+=sum(d)*H[F^m]
 return H[-1],*x
def mod_source(B):
 n=len(B)
 return [[0,1]+[0]*n,[0,0]+[1]*n]+[[1,0]+row for row in B]
def cert(A,r,i):
 W,s=chain(A,r,i);W,s=chain(W,s,1)
 row=ROWS[i]; p=pointed(A,r)
 assert profile(W)[-1]==sum(x*y for x,y in zip(row,p))
 return profile(W)[-1]
def tournaments(n):
 ij=list(combinations(range(n),2))
 for bits in product((0,1),repeat=len(ij)):
  A=[[0]*n for _ in range(n)]
  for (i,j),c in zip(ij,bits):A[i][j]=c;A[j][i]=1-c
  yield A
for n in range(1,4):
 for B in tournaments(n):
  b=profile(B); h=b[-1]; F=(1<<n)-1
  Q=sum(b[s]*b[F^s] for s in range(1<<n))
  T=mod_source(B)
  assert pointed(T,0)==(h+Q,Q//2,h+Q//2)
  for (i,j),v in PAIRS.items():
   assert cert(T,0,i)==v*Q+67*h
   assert cert(T,0,j)==v*Q+69*h
for i,j in PAIRS:
 assert cert([[0]],0,i)==67;assert cert([[0]],0,j)==69

for m in range(1,8):
 B=mod_source([[int(i<j) for j in range(m)] for i in range(m)])
 b=profile(B); F=(1<<len(B))-1
 assert b[-1]==2**m+1
 assert sum(b[i]*b[F^i] for i in range(F+1))==2**(m+2)+2*(3**m-1)

N=69; initial=(1,1,1,1); reached={initial}; terminal=set(); q=deque([initial])
while q:
 H,E,L,h=q.popleft()
 for nxt in [(H,0,2*L,h),(H+L-E,h,2*L-h,h),(H+L-E,H-E,H+L,H)]:
  a,b,c,d=nxt
  if a>N: continue
  if c-b>N-a: terminal.add(nxt);continue
  if nxt not in reached: reached.add(nxt);q.append(nxt)
assert len(reached)==116 and len(terminal)==176
assert all(v[0]!=67 for v in reached|terminal)
assert (69,1,73,1) in reached|terminal
print("PASS: independent subset Hamiltonian recounts, three paired reader slopes, and the complete root-word prune")
