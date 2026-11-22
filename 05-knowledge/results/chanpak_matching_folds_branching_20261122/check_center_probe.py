from itertools import product
from collections import Counter
HOSTS = [
 ((58,4,1,6,14,30),(2,5,0,4,3,1)),
 ((38,8,2,5,15,30),(3,5,4,0,2,1))
]
def H(t):
 n=len(t); d=[[0]*n for _ in range(2**n)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,2**n):
  for i in range(n):
   E=t[i]&~S
   while E:
    z=E&-E; E-=z
    d[S|z][z.bit_length()-1]+=d[S][i]
 return sum(d[-1])
def extend(t,s,c):
 z=len(t); M=sum(1<<i for i in c)
 N=tuple(x|((1<<z) if not(M>>i&1) else 0) for i,x in enumerate(t))+(M,)
 S=s+(z,)
 if not all((N[i]>>j&1)==(N[S[j]]>>S[i]&1) for i in range(z+1) for j in range(z+1)):
  raise ValueError("not anti-invariant")
 return N
if __name__=='__main__':
 for t,s in HOSTS:
  pairs=[(i,s[i]) for i in range(6) if i<s[i]]
  results=[]
  for c in product(*pairs):
   N=extend(t,s,c); results.append((c,H(N)))
  print(pairs, results,sorted(Counter(h for c,h in results).items()))

# Explicit combinatorial endpoint controls for the graded law.
def endpoint_tables(t):
 from itertools import permutations
 n=len(t); end=[[0]*n for _ in range(1<<n)]; start=[[0]*n for _ in range(1<<n)]
 for k in range(1,n+1):
  for p in permutations(range(n),k):
   if all(t[a]>>b&1 for a,b in zip(p,p[1:])):
    S=sum(1<<i for i in p); end[S][p[-1]]+=1; start[S][p[0]]+=1
 return end,start
def subs(S):
 U=S
 while True:
  yield U
  if not U:
   break
  U=(U-1)&S
def profiles(t,r):
 e,s=endpoint_tables(t); h=[1]+[sum(q) for q in e[1:]]
 d=[1]+[0]*(len(h)-1); tables=[d]
 for _ in range(r):
  d=[sum(h[U]*d[S^U] for U in subs(S)) for S in range(len(h))]
  tables.append(d)
 return tables,e,s
def R(t,sigma,r,D,E,S):
 M=(1<<len(t))-1
 return sum(E[A][i]*S[B][sigma[i]]*D[r-1][M^A^B]
            for A in range(1,M) for B in subs(M^A) if B
            for i in range(len(t)))
if __name__=='__main__':
 for index,(t,sig) in enumerate(HOSTS):
  d,e,s=profiles(t,6)
  if any(d[r][-1] != r**6+8*r**4+8*r*r for r in range(7)):
   raise RuntimeError("host profile mismatch")
  Ra,Rb=((38,14),(30,10))[index]
  if any(R(t,sig,r,d,e,s) != 6*r**4+Ra*r*r+Rb for r in range(1,6)):
   raise RuntimeError("reflected seam polynomial mismatch")
  normals=[profiles(extend(t,sig,c),3)[0]
           for c in product(*( (i,sig[i]) for i in range(6) if i<sig[i]))]
  check=[]
  for r in (1,2,3):
   Rr=R(t,sig,r,d,e,s)
   average_numerator=sum(p[r][-1] for p in normals)
   if 4*average_numerator != len(normals)*r*(2*d[r][-1]+d[r+1][-1]+d[r-1][-1]+Rr):
    raise RuntimeError("graded average mismatch")
   check.append(Rr)
  print("reflected seam grading:",check)
