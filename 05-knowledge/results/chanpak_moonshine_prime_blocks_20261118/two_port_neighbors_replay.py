from itertools import combinations, product
from functools import lru_cache

def req(b,info):
 if not b:raise ValueError(info)
def rows(n,code):
 A=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if not code>>k&1:i,j=j,i
  A[i]|=1<<j
 return A
def visits(A):
 n=len(A);F=(1<<(n-2))-1
 @lru_cache(None)
 def f(s,a,b,tail):
  W=A[tail] if tail<n else (1<<n)-1
  v={}
  def add(w):
   for p,c in w.items():v[p]=v.get(p,0)+c
  if s==F and a and b:v[a,b]=1
  if W&1 and a<n-1:add(f(s,a+1,b,0))
  if W&2 and b<n-1:add(f(s,a,b+1,1))
  for j in range(2,n):
   if W>>j&1 and not s>>(j-2)&1:add(f(s|1<<(j-2),a,b,j))
  return v
 return f(0,0,0,n)
def inflate(A,m,k):
 labels=[0]*m+[1]*k+list(range(2,len(A)))
 return [sum(1<<j for j,b in enumerate(labels) if (i<j if a==b else A[a]>>b&1)) for i,a in enumerate(labels)]
def h(A):
 n=len(A);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S,row in enumerate(D):
  for i,d in enumerate(row):
   if not d:continue
   w=A[i]&~S
   while w:
    z=w&-w;w-=z;D[S|z][z.bit_length()-1]+=d
 return sum(D[-1])
x=[rows(7,z) for z in (22775,21446)]
y=[[80,121,227,5,12,25,56,123],[96,121,203,33,13,20,56,123],[240,249,203,161,140,20,184,32]]
z=[rows(7,p) for p in (22921,20463)]
def relabel(T,order):
 return [sum(((T[i]>>j)&1)<<k for k,j in enumerate(order)) for i in order]
a=[54,100,64,23,70,92,9]; b=a.copy();b[3]^=1<<4;b[4]^=1<<3
w=[relabel(T,[1,2,0,3,4,5,6]) for T in (a,b)]
tests=0
for bank,primary,rest,form in [
 (x,[31,33],{(1,2):14,(2,1):18,(2,2):8},lambda m,k: 8*2**(m+k)+2**(m+1)-2**(k+1)-1),
 (y,[177,179,181],{(1,2):88,(2,1):179,(2,2):88,(3,1):48,(3,2):24},lambda m,k: 16*2**(m+k)+3*2**m+24*2**k*3**m-16*2**k-5),
 (w,[89,91],{(1,2):44,(2,1):89,(2,2):44,(3,1):24,(3,2):12},lambda m,k: 2**m+2**k*(8*2**m+12*3**m-8)-1),
 (z,[15,17],{(1,2):6,(2,1):9,(2,2):6},lambda m,k: 3*(2**m-1)*(2**(k+1)-1)+6),
]:
 for i,A in enumerate(bank):
  req(visits(A)=={(1,1):primary[i],**rest},(A,visits(A)))
  for m,k in product(range(1,4),repeat=2):
   got=h(inflate(A,m,k)); req(got==form(m,k)+2*i,(m,k,got));tests+=1
  print(A,visits(A))
print(tests,'literal inflation checks')
