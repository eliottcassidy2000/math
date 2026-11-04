"""Standalone exact checks for a three-state, three-run unit gate."""
import itertools,json,math

def need(x,msg):
 if not x:raise RuntimeError(msg)
def seed(mask):
 a=[0]*6
 for bit,(i,j) in enumerate(itertools.combinations(range(5),2)):
  if mask>>bit&1:a[i]|=1<<j
  else:a[j]|=1<<i
 for i in range(5):
  if i<2:a[5]|=1<<i
  else:a[i]|=1<<5
 return a
def profile(a,v=5):
 out=[0]*6
 for p in itertools.permutations([i for i in range(len(a)) if i!=v]):
  for mask in range(1,1<<len(a)):
   good=True
   for j in range(len(a)):
    if mask>>j&1:
     if j and not(a[p[j-1]]>>v&1):good=False;break
     if j<len(p) and not(a[v]>>p[j]&1):good=False;break
    elif 0<j<len(p) and not(a[p[j-1]]>>p[j]&1):good=False;break
   if good:out[mask.bit_count()-1]+=1
 return out
def inflate(a,r):
 b=[0]*(5+r)
 for i in range(5):
  b[i]=a[i]&31
  for j in range(5,5+r):
   if a[i]>>5&1:b[i]|=1<<j
   else:b[j]|=1<<i
 for i in range(5,5+r):
  for j in range(i+1,5+r):b[i]|=1<<j
 return b
def hp(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   c=dp[mask][i]
   if not c:continue
   opts=a[i]&~mask
   while opts:
    b=opts&-opts;opts-=b;j=b.bit_length()-1;dp[mask|b][j]+=c
 return sum(dp[-1])
def strong(a):
 for x in range(len(a)):
  seen=1<<x;old=0
  while old!=seen:
   old=seen
   for i in range(len(a)):
    if seen>>i&1:seen|=a[i]
  if seen!=(1<<len(a))-1:return False
 return True
def qcover(a,k):
 # independent color assignment, empty colors excluded, internal HP subset count
 n=len(a);cache={}
 for s in range(1,1<<n):
  ids=[i for i in range(n) if s>>i&1]
  sub=[sum(1<<j for j,v in enumerate(ids) if a[u]>>v&1) for u in ids];cache[s]=hp(sub)
 ans=0
 for word in itertools.product(range(k),repeat=n):
  sets=[sum(1<<i for i,x in enumerate(word) if x==c) for c in range(k)]
  if all(sets):ans+=math.prod(cache[s] for s in sets)
 return ans
def substitute(a,b):
 n=len(b);o=list(b)+[0]*5
 for i in range(5):
  for j in range(5):
   if a[i]>>j&1:o[n+i]|=1<<(n+j)
  for j in range(n):
   if a[i]>>5&1:o[n+i]|=1<<j
   else:o[j]|=1<<(n+i)
 return o
def run():
 seeds=[seed(x) for x in (24,280,281)];profiles=[]
 for j,a in enumerate(seeds):
  p=profile(a);need(p==[27+2*j,21,4,0,0,0],'skeleton table');profiles.append(p);need(strong(a),'seed strong')
 need([(seeds[j][i]^seeds[j+1][i]).bit_count() for j in range(2) for i in range(6)].count(1)==4,'edge differences')
 records=[]
 for r in range(1,8):
  center=4*3**r+9*2**r
  counts=[]
  for j,a in enumerate(seeds):
   b=inflate(a,r);v=hp(b);need(v==center-3+2*j,('inflation',r,j,v));need(strong(b),'inflation strong');counts.append(v)
  records.append({'r':r,'counts':counts})
 controls=0
 for n in range(1,5):
  for mask in range(1<<(n*(n-1)//2)):
   b=[0]*n
   for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
    if mask>>bit&1:b[i]|=1<<j
    else:b[j]|=1<<i
   h=hp(b);q2=qcover(b,2);q3=qcover(b,3)
   for j,a in enumerate(seeds):
    got=hp(substitute(a,b));need(got==(27+2*j)*h+21*q2+4*q3,'module identity');need(strong(substitute(a,b)),'module strong');controls+=1
 return {'seeds':seeds,'profiles':profiles,'clock_controls':records,'arbitrary_module_controls':controls,'status':'exact checks pass; no spectrum-completeness claim'}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
