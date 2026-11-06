import itertools,json

def require(value):
 if not value: raise AssertionError("verification failed")

def count(a):
 n=len(a); dp=[[0]*n for _ in range(1<<n)]
 for i in range(n): dp[1<<i][i]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    dp[s][v]+=sum(dp[s^(1<<v)][u] for u in range(n) if s>>u&1 and a[u][v])
 return sum(dp[-1])
def mat(n,mask):
 a=[[0]*n for _ in range(n)]; bit=0
 for i in range(n):
  for j in range(i+1,n):
   x= 1 if j==i+1 else 1-((mask>>bit)&1)
   if j!=i+1:bit+=1
   a[i][j]=x;a[j][i]=1-x
 return a
def complement(a):
 n=len(a);return [[0 if i==j else a[i][j] if abs(i-j)==1 else 1-a[i][j] for j in range(n)]for i in range(n)]
def pad(a,k):
 n=len(a);b=[[int(i<j) for j in range(n+2*k)]for i in range(n+2*k)]
 for i in range(n):
  for j in range(n):b[i+k][j+k]=a[i][j]
 return b
rows=[]
for n in range(1,5):
 for mask in range(1<<((n-1)*(n-2)//2)):
  t=mat(n,mask); u=complement(t);h=count(t);g=count(u)
  for k in (1,2,3):
   t1=pad(t,k);u1=complement(t1);hh=count(t1);gg=count(u1)
   require(hh==h)
   require(gg>=g+2)
   if k==3:require(gg>=4*g+1)
   rows.append([n,mask,k,h,g,gg])

def auts(a):
 groups={}
 for v,row in enumerate(a):groups.setdefault(sum(row),[]).append(v)
 cells=list(groups.values());out=[]
 for images in itertools.product(*(itertools.permutations(c) for c in cells)):
  p=list(range(len(a)))
  for cell,img in zip(cells,images):
   for v,w in zip(cell,img):p[v]=w
  if all(a[i][j]==a[p[i]][p[j]] for i in range(len(a)) for j in range(len(a))):out.append(p)
 return out
for n in range(1,6):
 for mask in range(1<<((n-1)*(n-2)//2)):
  t=mat(n,mask);u=complement(t);t1=pad(t,3);u1=complement(t1)
  require(len(auts(t1))==len(auts(t)))
  require(len(auts(u1))==sum(p[0]==0 and p[n-1]==n-1 for p in auts(u)))

print(json.dumps({'status':'PASS','universe':'all fixed-path tilings at orders 1..4, padding k=1,2,3','count_cases':len(rows),'automorphism_cases':76,'automorphism_universe':'all fixed-path tilings at orders 1..5, padding k=3','rows':rows},indent=2))
