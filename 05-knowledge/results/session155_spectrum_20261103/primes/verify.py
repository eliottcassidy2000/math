from itertools import permutations, product

def fixture(m,mask):
 n=m+2; A=[[0]*n for _ in range(n)]
 for i in range(m):
  for j in range(i+1,m): A[i][j]=1
 for i in range(m):
  for a in range(2):
   x=m+a
   if mask>>(2*i+a)&1: A[x][i]=1
   else: A[i][x]=1
 x,y=m,m+1
 if mask>>(2*m)&1:A[x][y]=1
 else:A[y][x]=1
 return A

def pdp(A):
 n=len(A);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   for j in range(n):
    if not (mask>>j)&1 and A[i][j]: d[mask|1<<j][j]+=d[mask][i]
 return sum(d[-1])

def brute(A):
 return sum(all(A[p[i]][p[i+1]] for i in range(len(p)-1)) for p in permutations(range(len(A))))

def transitive(A,excepts=()):
 v=[i for i in range(len(A)) if i not in excepts]; return all(not(A[i][j] and A[j][k] and A[k][i]) for i in v for j in v for k in v)

def strong(A):
 n=len(A); R=[sum(1<<j for j in range(n) if A[i][j])|1<<i for i in range(n)]
 for k in range(n):
  for i in range(n):
   if R[i]>>k&1:R[i]|=R[k]
 return all(r==(1<<n)-1 for r in R)

def wordcount(w):
 h,a,p=1,0,0
 for bit in w:
  if bit:a,p=2*a+1,1
  else:h,a,p=h+a+p,2*a,0
 return h

for n in range(2,8):
 vals=set()
 for b in product((0,1), repeat=n-2):
  word=(1,)+b+(0,); h=wordcount(word); vals.add(h)
  # Independent vertex-level subset DP for every strong word.
  A=[[0]*(n+1) for _ in range(n+1)]
  for i in range(n):
   for j in range(i+1,n):A[i][j]=1
   if word[i]:A[n][i]=1
   else:A[i][n]=1
  if pdp(A)!=h or not strong(A):raise SystemExit('word verification failed')
 print('one-apex length',n,'values',sorted(vals))
for m,mask,H in [(3,6,13),(4,14,29),(5,202,35),(5,106,49)]:
 A=fixture(m,mask); P=pdp(A);Q=brute(A); ok=strong(A); f2=transitive(A,(m,m+1)) and all(not transitive(A,(i,)) for i in range(m+2))
 print('H',H,'order',m+2,'mask',mask,'rows',[sum(1<<j for j,a in enumerate(row) if a) for row in A],'DP',P,'brute',Q,'strong',ok,'FVS exactly 2',f2)
 if (P,Q,ok,f2)!=(H,H,True,True):raise SystemExit('FAIL')
