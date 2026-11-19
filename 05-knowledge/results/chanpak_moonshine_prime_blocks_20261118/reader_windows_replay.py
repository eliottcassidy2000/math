from itertools import combinations, permutations
from math import comb

def rows(n,code):
 A=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if not(code>>k&1):i,j=j,i
  A[i]|=1<<j
 return A

def coeff(A):
 n=len(A);out=[0]*((n+1)//2+1)
 for p in permutations(range(1,n)):
  for m in range(1,1<<n):
   r=m.bit_count()
   if r>=len(out) or (m&(m<<1)):continue
   w=[]
   for i in range(n):
    if m>>i&1:w.append(0)
    if i<n-1:w.append(p[i])
   if all(A[i]>>j&1 for i,j in zip(w,w[1:])):out[r]+=1
 return out[1:]

def inflate(A,m):
 q=[0]*m+list(range(1,len(A)))
 return [sum(1<<j for j,b in enumerate(q) if (i<j if a==b else A[a]>>b&1)) for i,a in enumerate(q)]

def h(A):
 n=len(A);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S,row in enumerate(D):
  for i,d in enumerate(row):
   w=A[i]&~S
   if d:
    while w:
     b=w&-w;w^=b;D[S|b][b.bit_length()-1]+=d
 return sum(D[-1])

for n,codes,start,k in ((6,[318,315,335,343],11,6),(7,[1403,1399,11454,5374,1679],29,18)):
 for i,z in enumerate(codes):
  A=rows(n,z); c=coeff(A)
  expected=[start+2*i,k]+[0]*(len(c)-2)
  if c!=expected:raise ValueError((z,c,expected))
  for m in range(1,7):
   if h(inflate(A,m))!=start+2*i+k*((1<<m)-2):raise ValueError((z,m))
  print(n,z,A,c)
print("PASS: 54 literal independently counted supports, plus exact slot coefficients")
