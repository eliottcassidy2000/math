from itertools import permutations
from math import factorial

def require(v,msg):
 if not v:raise RuntimeError(msg)

def wrap(rows):
 n=len(rows)
 return [r|(1<<(n+1)) for r in rows]+[(1<<n)-1,1<<n]

def operator(p):
 q=[0]*(len(p)+2)
 for k,a in enumerate(p):
  q[k+2]+=a;q[k+1]+=(2*k+1)*a;q[k]+=(k*k+2*k)*a
  if k:q[k-1]+=k*(k-1)*a
 return q

def count_h(rows):
 n=len(rows);f=[[0]*n for _ in range(1<<n)]
 for v in range(n):f[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1:
    rem=mask^(1<<v)
    if rem:f[mask][v]=sum(f[rem][u] for u in range(n) if rem>>u&1 and rows[u]>>v&1)
 return sum(f[-1]) if n else 1

def cover_polynomial(rows):
 n=len(rows)
 if not n:return [1]
 out=[0]*(n+1)
 for word in permutations(range(n)):
  bad=sum((not(rows[word[i]]>>word[i+1]&1))<<i for i in range(n-1))
  for cut in range(1<<(n-1)):
   if not(bad&~cut):out[1+cut.bit_count()]+=1
 for k in range(n+1):
  require(out[k]%factorial(k)==0,'unordered cover factorial')
  out[k]//=factorial(k)
 return out

def det(a):
 a=[r[:] for r in a];n=len(a);last=1;sgn=1
 for k in range(n-1):
  if a[k][k]==0:
   j=next((j for j in range(k+1,n) if a[j][k]),None)
   if j is None:return 0
   a[k],a[j]=a[j],a[k];sgn=-sgn
  v=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    num=a[i][j]*v-a[i][k]*a[k][j];require(num%last==0,'exact determinant division');a[i][j]=num//last
  for i in range(k+1,n):a[i][k]=0
  last=v
 return sgn*a[-1][-1]

def transitive(n):return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]

def strong(rows):
 n=len(rows)
 for v in range(n):
  seen=1<<v;old=-1
  while old!=seen:
   old=seen
   for u in range(n):
    if seen>>u&1:seen|=rows[u]
  if seen!=(1<<n)-1:return False
 return True

def feedback_number(rows):
 n=len(rows);largest=0
 for mask in range(1<<n):
  k=mask.bit_count()
  if k<=largest:continue
  degrees=sorted((rows[v]&mask).bit_count() for v in range(n) if mask>>v&1)
  if degrees==list(range(k)):largest=k
 return n-largest
