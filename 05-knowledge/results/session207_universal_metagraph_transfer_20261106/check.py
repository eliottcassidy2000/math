import itertools,json

def ck(v,msg):
 if not v:raise AssertionError(msg)
def tournament(n,mask,path=False):
 t=[[0]*n for _ in range(n)];z=0
 for i in range(n):
  for j in range(i+1,n):
   bit=0 if path and j==i+1 else (mask>>z)&1
   if not(path and j==i+1):z+=1
   t[i][j]=1-bit;t[j][i]=bit
 return t
def hp(t):
 n=len(t)
 if not n:return 1
 out=[sum(t[i][j]<<j for j in range(n))for i in range(n)]
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for v,c in enumerate(dp[s]):
   if c:
    z=out[v]&~s
    while z:
     bit=z&-z;w=bit.bit_length()-1;z-=bit;dp[s|bit][w]+=c
 return sum(dp[-1])
def induced(t,vs):return [[t[i][j]for j in vs]for i in vs]
def left(t,q,a):
 n=len(t);r=[[int(j==i+1 or i>=j+2)for j in range(a)]for i in range(a)]
 s=[[0]*(n+a)for _ in range(n+a)]
 for i in range(a):
  for j in range(a):s[i][j]=r[i][j]
 for i in range(n):
  for j in range(n):s[a+i][a+j]=t[i][j]
  for j in range(a):s[a+i][j]=1
 if a:s[a-1][a+q]=1;s[a+q][a-1]=0
 return s
def weights(t,q):
 n=len(t);w=[0,0]
 for k in range(n):
  for tail in itertools.permutations([v for v in range(n)if v!=q],k):
   p=(q,)+tail
   if all(t[p[i]][p[i+1]]for i in range(len(p)-1)):
    w[len(p)%2]+=hp(induced(t,[v for v in range(n)if v not in p]))
 return w
def complement(t):return [[0 if i==j else t[i][j]if abs(i-j)==1 else 1-t[i][j]for j in range(len(t))]for i in range(len(t))]
def padded_target(t,a,b):
 n=len(t);N=n+a+b;s=[[int(i<j)for j in range(N)]for i in range(N)]
 for i in range(n):
  for j in range(n):s[a+i][a+j]=t[i][j]
 return complement(s)
def auts(t):
 cells={}
 for i,r in enumerate(t):cells.setdefault(sum(r),[]).append(i)
 cells=list(cells.values());out=[]
 for imgs in itertools.product(*(itertools.permutations(c)for c in cells)):
  p=list(range(len(t)))
  for c,img in zip(cells,imgs):
   for v,w in zip(c,img):p[v]=w
  if all(t[i][j]==t[p[i]][p[j]]for i in range(len(t))for j in range(len(t))):out.append(p)
 return out
pointed=0;tail_cases=0;one_sided_group_cases=0
for n in range(1,5):
 for mask in range(1<<(n*(n-1)//2)):
  t=tournament(n,mask)
  for q in range(n):
   we,wo=weights(t,q);f=[hp(t),hp(t)+2*we,hp(t)+2*we+2*wo]
   for a in range(3,6):f.append(sum(f[-3:]))
   for a in range(6):ck(hp(left(t,q,a))==f[a],('pointed tail',n,mask,q,a));tail_cases+=1
   ck(len(auts(left(t,q,3)))==sum(p[q]==q for p in auts(t)),('one-sided stabilizer',n,mask,q));one_sided_group_cases+=1
   pointed+=1
rows=[];aut_cases=0
for n in range(1,5):
 for mask in range(1<<((n-1)*(n-2)//2)):
  t=tournament(n,mask,True);u=complement(t);g={}
  for a in range(5):
   for b in range(5):
    v=padded_target(t,a,b);g[a,b]=hp(v)
    if a>=3:ck(g[a,b]==sum(g[a-i,b]for i in(1,2,3)),('left recurrence',n,mask,a,b))
    if b>=3:ck(g[a,b]==sum(g[a,b-i]for i in(1,2,3)),('right recurrence',n,mask,a,b))
    if a>=3 and b>=3:
     ck(len(auts(v))==sum(p[0]==0 and p[n-1]==n-1 for p in auts(u)),('stabilizer',n,mask,a,b));aut_cases+=1
    rows.append([n,mask,a,b,g[a,b]])
# Positive nontrivial stabilizer control: source -> C3 -> sink.
u=tournament(5,0);u[1][3]=0;u[3][1]=1;t=complement(u)
ck(len(auts(u))==3,'C3 endpoint stabilizer input')
nontrivial=[]
for a,b in ((3,3),(3,4),(4,3),(4,4)):
 v=padded_target(t,a,b);order=len(auts(v));ck(order==3,'nontrivial stabilized group');nontrivial.append([a,b,hp(v),order])
print(json.dumps({'status':'PASS','arbitrary_pointed_cores':pointed,'tail_cases':tail_cases,'one_sided_stabilizer_cases':one_sided_group_cases,'two_sided_cases':len(rows),'stabilizer_cases':aut_cases,'nontrivial_stabilizer_controls':nontrivial,'universes':{'pointed':'all labelled tournaments n1..4, every vertex q, tail a0..5','two_sided':'all fixed-path tilings n1..4, tails a,b0..4','stabilizer':'all those cores, a,b in {3,4}'},'rows':rows},indent=2))
