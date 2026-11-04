from itertools import permutations,combinations


def need(b,s):
 if not b:raise RuntimeError(s)

def mul(a,b):
 N=len(a);out=[0]*N
 for s in range(N):
  t=s
  while True:
   out[s]+=a[t]*b[s^t]
   if not t:break
   t=(t-1)&s
 return out

def add(a,b,c=1):return [x+c*y for x,y in zip(a,b)]

def paths(rows):
 n=len(rows);N=1<<n;G=[[[0]*N for _ in rows]for _ in rows]
 # Endpoint subset recurrence, independent from producer's permutation enumerator.
 for start in range(n):
  dp=[[0]*n for _ in range(N)];dp[1<<start][start]=1
  for mask in range(N):
   for end,v in enumerate(dp[mask]):
    if not v:continue
    G[start][end][mask]+=v;avail=rows[end]&~mask
    while avail:
     b=avail&-avail;avail-=b;dp[mask|b][b.bit_length()-1]+=v
 return G

def outer(rows,into):
 n=len(rows);N=1<<n;G=paths(rows);S=[[sum(G[i][j][s]for j in range(n))for s in range(N)]for i in range(n)];E=[[sum(G[j][i][s]for j in range(n))for s in range(N)]for i in range(n)]
 I=[i for i in range(n)if into>>i&1];O=[i for i in range(n)if not(into>>i&1)];one=[1]+[0]*(N-1);alpha=[one[s]+sum(E[i][s]for i in I)for s in range(N)];beta=[one[s]+sum(S[i][s]for i in O)for s in range(N)];tau=[sum(G[u][v][s]for u in O for v in I)for s in range(N)];ell=[[sum(G[i][v][s]for v in I)for s in range(N)]for i in range(n)];rho=[[sum(G[u][i][s]for u in O)for s in range(N)]for i in range(n)];return alpha,beta,tau,S,E,ell,rho

def peval(p,tau):
 out=[0]*len(tau);power=[1]+[0]*(len(tau)-1)
 for a in p:out=[x+a*y for x,y in zip(out,power)];power=mul(power,tau)
 return out

def derivative(p):return[(j+1)*p[j+1]for j in range(len(p)-1)]
def kernel(rows,p,q,r):
 n=len(rows);L=[0]*(n-1);R=[0]*(n-1)
 for v in permutations(range(n)):
  for cut in range(1,n):
   pol=[1]
   for j in range(n-1):
    if j+1==cut:continue
    a=(rows[v[j]]>>v[j+1])&1;pol=[a*pol[0]]+[a*pol[k]+pol[k-1]for k in range(1,len(pol))]+[pol[-1]]
   if v[-1]==r and((v[cut-1],v[cut])in[(p,q),(q,p)]):L=[a+b for a,b in zip(L,pol)]
   if v[0]==r and((v[cut-1],v[cut])in[(p,q),(q,p)]):R=[a+b for a,b in zip(R,pol)]
 return L,R

def tt(order):return[sum(1<<v for v in order[order.index(u)+1:])for u in range(len(order))]
def capacity(rows):
 n=len(rows);C=[[0]*n for _ in range(n)]
 for v in permutations(range(n)):
  bad=[j for j in range(n-1)if not(rows[v[j]]>>v[j+1]&1)]
  if len(bad)>1:continue
  for cut in (bad if bad else range(n-1)):C[v[cut]][v[cut+1]]+=1
 return C

def sympoly(rows,z,G):
 n=len(rows);out=[0]*(n-1)
 for v in permutations(range(n)):
  for cut in range(1,n):
   if not((v[cut-1]==z and v[cut]in G)or(v[cut]==z and v[cut-1]in G)):continue
   pol=[1]
   for j in range(n-1):
    if j+1==cut:continue
    a=(rows[v[j]]>>v[j+1])&1;pol=[a*pol[0]]+[a*pol[k]+pol[k-1]for k in range(1,len(pol))]+[pol[-1]]
   out=[a+b for a,b in zip(out,pol)]
 return out

def replacement_formula(A,B,rows,into,p,q,r,i):
 L,R=kernel(A,p,q,r);l,rr=kernel(B,p,q,r);P=[x-y for x,y in zip(L,l)];Q=[x-y for x,y in zip(R,rr)];alpha,beta,tau,S,E,ell,rho=outer(rows,into);x=add(mul(mul(alpha,S[i]),peval(P,tau)),mul(mul(beta,E[i]),peval(Q,tau)));y=mul(mul(alpha,beta),add(mul(ell[i],peval(derivative(P),tau)),mul(rho[i],peval(derivative(Q),tau))));return add(x,y)[-1],P,Q
