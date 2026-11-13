from itertools import combinations, product
def mat(n,c):
 A=[0]*n
 for i,(u,v) in enumerate(combinations(range(n),2)):
  if c>>i&1:A[u]|=1<<v
  else:A[v]|=1<<u
 return A
def stats(A,cuts=None):
 n=len(A);N=(1<<n)-1
 s=[[0]*n for _ in range(N+1)];e=[[0]*n for _ in range(N+1)]
 for i in range(n):s[1<<i][i]=e[1<<i][i]=1
 for S in range(1,N+1):
  if S&(S-1)==0:continue
  for i in range(n):
   if S>>i&1:
    T=S^(1<<i)
    s[S][i]=sum(s[T][j] for j in range(n) if A[i]>>j&1)
    e[S][i]=sum(e[T][j] for j in range(n) if A[j]>>i&1)
 h=[1]+list(map(sum,s[1:]))
 ans={}
 for B in (range(N+1) if cuts is None else cuts):
  start=sum(s[N][i] for i in range(n) if B>>i&1)
  end=sum(e[N][i] for i in range(n) if not B>>i&1)
  sstar=sum(h[N^U]*sum(s[U][i] for i in range(n) if B>>i&1) for U in range(1,N+1))
  estar=sum(h[N^U]*sum(e[U][i] for i in range(n) if not B>>i&1) for U in range(1,N+1))
  split=sum(sum(s[U][i] for i in range(n) if B>>i&1)*sum(e[N^U][i] for i in range(n) if not B>>i&1) for U in range(1,N))
  ans[B]=((start,sstar),(end,estar,split))
 return h[N],sum(h[U]*h[N^U] for U in range(N+1)),ans
def formula(L):
 H=0
 for i,(h,q,s,ss,e,ee,b) in enumerate(L):
  others=1
  for j in range(len(L)):
   if j!=i:others*=L[j][0]
  H+=others*(b+(s if i==0 else 0)+(e if i==len(L)-1 else 0))
  for j in range(i+1):
   if j==i and i+1<len(L):
    out=1
    for k in range(len(L)):
     if k not in (i,i+1):out*=L[k][0]
    H+=out*e*L[i+1][2]
   elif j<i:
    out=1
    for k in range(len(L)):
     if k<j or k>i:out*=L[k][0]
     if j<k<i:out*=L[k][1]
    H+=out*ee*L[j][3]
 return H
from itertools import product,combinations
P=[((56,17,3,6,12,30,63),14),((52,33,2,7,14,28,63),40),((52,33,2,7,14,28,63),17)]
def get(A,c):
 h,q,d=stats(A,[c]);(s,S),(e,E,B)=d[c];return h,q,s,S,e,E,B
def pad(A):return tuple(a|(1<<len(A)) for a in A)+(0,)
def comp(v):
 h,q,s,S,e,E,B=v
 return h,2*q,s,2*S,h,E+q,B+S
L=[]
for A,c in P:
 u=get(A,c); w=get(pad(pad(A)),c)
 assert comp(comp(u))==w,(u,w)
 L.append(w)
print(L)
def val(x):
 k=0
 while x%11==0:x//=11;k+=1
 return k,x%11
for d in range(1,5):
 D=list(product(range(3),repeat=d))
 C=[formula([L[1]]+[L[i] for i in w[::-1]]) for w in D]
 assert C==sorted(set(C))
 for i,j in combinations(range(len(C)),2):
  v,s=val(C[j]-C[i]);z=next(t for t in range(d) if D[i][t]!=D[j][t])
  assert v==z
  signed = pow(s,5,11)==1
  assert signed== (pow((L[D[j][z]][5]-L[D[i][z]][5])%11,5,11)!=1)
 print(d,len(C),min(C),max(C))
for k in range(1,4):
 mod=11**k; h,q,_,S,_,_,_=L[0]
 C=set()
 for r in [0,1]:
  E=L[r][5]; c=S*E*pow((1-h*pow(q,-1,mod))%mod,-2,mod)%mod
  for n in range(k+2,k+2+5*11**(k-1)):
   assert ((L[r][4]+n*L[r][6])*pow(h,n-1,mod) + S*E*sum((n-1-g)*pow(q,g,mod)*pow(h,n-2-g,mod) for g in range(n-1)))%mod==c*pow(q,n-2,mod)%mod
   C.add(c*pow(q,n-2,mod)%mod)
 assert len(C)==10*11**(k-1)
