"""Squarefree path algebra; fixed physical names."""
def matrix(R):
 n=len(R);N=1<<n;G=[[[0]*N for _ in R] for _ in R]
 for s in range(n):
  d=[[0]*n for _ in range(N)];d[1<<s][s]=1
  for M in range(1,N):
   for v,z in enumerate(d[M]):
    if not z:continue
    G[s][v][M]=z;A=R[v]&~M
    while A:
     b=A&-A;A-=b;d[M|b][b.bit_length()-1]+=z
 return G

def plus(*seq):return [sum(v) for v in zip(*seq)]
def minus(A,B):return [x-y for x,y in zip(A,B)]
def product(A,B):
 N=len(A);D=[0]*N
 for M in range(N):
  S=M
  while True:
   D[M]+=A[S]*B[M^S]
   if S==0:break
   S=(S-1)&M
 return D

def full(A,B):return sum(A[s]*B[len(A)-1-s] for s in range(len(A)))
def kernel(R,out,inn,i,j):
 n=len(R);N=1<<n;G=matrix(R);E=[plus(*(G[s][v] for s in range(n))) for v in range(n)];S=[plus(*G[v]) for v in range(n)];zero=[0]*N
 alpha=plus(zero,*(E[v] for v in inn));alpha[0]=1
 beta=plus(zero,*(S[v] for v in out));beta[0]=1
 ell=[plus(zero,*(G[v][w] for w in inn)) for v in range(n)];rr=[plus(zero,*(G[w][v] for w in out)) for v in range(n)];tau=plus(zero,*(G[u][v] for u in out for v in inn));tau=[2*x for x in tau];tau[0]=1
 ab=product(alpha,beta)
 term1=2*full(ab,minus(G[i][j],G[j][i]))
 es=minus(product(E[i],S[j]),product(E[j],S[i]));term2=-full(tau,es)
 term3=full(alpha,minus(product(rr[i],S[j]),product(rr[j],S[i])))
 term4=full(beta,minus(product(E[i],ell[j]),product(E[j],ell[i])))
 return [term1,term2,term3,term4]
