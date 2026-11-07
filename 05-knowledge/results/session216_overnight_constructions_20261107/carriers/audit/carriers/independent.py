from pathlib import Path
import sys,json,itertools,random,importlib.util
P=Path(__file__).parent;sys.path.insert(0,str(P/'frozen'))
import probe as prod

def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return c

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def coeff(rows,v):
 # A transition between external vertices is either direct, or a single
 # module visit. This contracts placeholders before the subset recurrence.
 ex=[u for u in range(len(rows)) if u!=v];n=len(ex);N=1<<n;root=ex.index(0)
 def paths(start_root=False):
  dp=[[[] for _ in ex] for S in range(N)]
  for j,u in enumerate(ex):
   if not start_root or j==root:dp[1<<j][j]=[1,int(not start_root and rows[v]>>u&1)]
  for S in range(1,N):
   for j,u in enumerate(ex):
    if not dp[S][j]:continue
    for l,w in enumerate(ex):
     if S>>l&1:continue
     weight=[int(rows[u]>>w&1),int(bool(rows[u]>>v&1) and bool(rows[v]>>w&1))]
     dp[S|1<<l][l]=add(dp[S|1<<l][l],mul(dp[S][j],weight))
  allp=[[] for _ in range(N)];endp=[[] for _ in range(N)]
  allp[0]=[] if start_root else [1,1]
  for S in range(1,N):
   for j,u in enumerate(ex):
    allp[S]=add(allp[S],mul(dp[S][j],[1,int(bool(rows[u]>>v&1))])) if dp[S][j] else allp[S]
   endp[S]=dp[S][root]
  return allp,endp
 hp,ep=paths();sp,_=paths(True)
 def conv(a,b,S):
  z=[];A=S
  while True:
   if a[A] and b[S^A]:z=add(z,mul(a[A],b[S^A]))
   if not A:break
   A=(A-1)&S
  return z
 F=N-1;D=F^(1<<root)
 values=[hp[F],sp[F],ep[F],conv(ep,hp,F),conv(sp,hp,F),hp[D],conv(hp,hp,F),conv(hp,hp,D)]
 return [(x+[0]*(len(rows)+2))[1:len(rows)+2] for x in values]

def seeds():
 out=[]
 for X,js,r,a,b in prod.PAIRS:out.append((X,prod.remark(a,r),prod.remark(b,r)))
 for X,a,r,b,t,w in [(4225,(2,20,41,67,77,27,167,63),6,(6,12,24,49,195,23,47,111),4,'JO'),(4357,(10,20,65,22,37,15,59),3,(2,48,11,67,45,77,23),2,'ZJO')]:
  a,b=prod.remark(a,r),prod.remark(b,t)
  for c in w:a,b=prod.graph_move(a,c),prod.graph_move(b,c)
  out.append((X,a,b))
 return out

def profile(a):
 # Independent forward extension, with root-start counts tracked separately.
 n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)];sd=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for j in range(n):dp[1<<j][j]=1
 sd[1][0]=1
 for S in range(1,N):
  H[S]=sum(dp[S])
  for j in range(n):
   for u in range(n):
    if not S>>u&1 and a[j]>>u&1:
     dp[S|1<<u][u]+=dp[S][j];sd[S|1<<u][u]+=sd[S][j]
 F=N-1
 start=[sum(x) for x in sd]
 def cov(S):
  total=0;A=S
  while True:
   total+=H[A]*H[S^A]
   if not A:return total
   A=(A-1)&S
 L=sum(dp[A][0]*H[F^A] for A in range(N));R=sum(start[A]*H[F^A] for A in range(N))
 return (H[F],start[F],dp[F][0],L,R,H[F^1],cov(F),cov(F^1))

def substitute(rows,v,M):
 ex=[i for i in range(len(rows)) if i!=v];labels=ex+list(range(len(rows),len(rows)+len(M)));q=len(ex);out=[]
 for i,u in enumerate(labels):
  row=0
  for j,w in enumerate(labels):
   if i==j:continue
   b=(M[i-q]>>(j-q)&1) if i>=q and j>=q else rows[v if i>=q else u]>>(v if j>=q else w)&1
   row|=b<<j
  out.append(row)
 return out

def covers(M,K):
 n=len(M);H=[0]*(1<<n)
 for S in range(1,1<<n):
  vertices=[j for j in range(n) if S>>j&1]
  H[S]=sum(all(M[a]>>b&1 for a,b in zip(p,p[1:])) for p in itertools.permutations(vertices))
 dp=[1]+[0]*((1<<n)-1);res=[]
 for k in range(K):
  nxt=[0]*(1<<n)
  for S in range(1<<n):
   A=S
   while A:nxt[S]+=H[A]*dp[S^A];A=(A-1)&S
  dp=nxt;res.append(dp[-1])
 return res

symbolic=0;matches=[];literal=0
for X,a,b in seeds():
 banks=[]
 for rows in (a,b):
  bank=[]
  for v in range(1,len(rows)):
   C=coeff(rows,v);expected=prod.polynomials(rows,v)
   if C!=expected:raise RuntimeError(('coefficients',X,v,C,expected))
   symbolic+=1;s=C[1]
   if any(s[1:]) or not s[0]%3 or s[0]>26:continue
   j=s[0]//3;B=[h+l-z for h,l,z in zip(C[0],C[3],C[1])];B[0]+=2*j
   key=[B,C[2],add(C[4],C[7]),C[5]];bank.append((key,v,C))
  banks.append(bank)
 for ka,va,ca in banks[0]:
  for kb,vb,cb in banks[1]:
   if ka!=kb:continue
   matches.append((X,va,vb))
   for rows,v,C in ((a,va,ca),(b,vb,cb)):
    for M in ((0,),(2,0),(6,4,0),(2,4,1)):
     c=covers(M,len(C[0]));pred=tuple(sum(x*y for x,y in zip(col,c)) for col in C)
     actual=profile(substitute(rows,v,M))
     if actual!=pred:raise RuntimeError(('literal',X,v,M,actual,pred))
     literal+=1
if matches!=[(759,1,1),(821,1,1),(4225,9,9),(4357,9,9)]:raise RuntimeError(matches)

# Physically test right-amalgamation with ordinary random rooted tournaments.
def random_graph(n,rng):
 a=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   x,y=(i,j) if rng.randrange(2) else (j,i);a[x]|=1<<y
 return a

def join(U,W):
 n=len(U);m=len(W);a=list(U)+[0]*(m-1)
 for i in range(m):
  for j in range(m):
   if W[i]>>j&1:a[0 if i==0 else n+i-1]|=1<<(0 if j==0 else n+j-1)
 for i in range(1,n):
  for j in range(n,n+m-1):a[i]|=1<<j
 return a
rng=random.Random(216);rightchecks=0;singleton=0
for z in range(160):
 U=random_graph(rng.randrange(2,6),rng);W=random_graph(rng.randrange(1,5),rng)
 h,s,e,L,R,d,C,b=profile(U);hw,sw,ew,Lw,Rw,q,Cw,bw=profile(W)
 pp=profile(join(U,W));B=h+L-s;V=R+b
 actual=(pp[0]+pp[3]-pp[1],pp[2],pp[4]+pp[7],pp[5],pp[1])
 if len(W)==1:expected=(B,e,V,d,s);singleton+=1
 else:expected=(d*hw+(V-d)*Lw+(e-d)*sw+(B-e-V+d)*q,d*ew,bw*V+d*(Rw-bw),d*q,s*q)
 if actual!=expected:raise RuntimeError(('right law',U,W,actual,expected))
 rightchecks+=1
print(json.dumps({'symbolic_positions':symbolic,'matches':matches,'independent_literal_profiles':literal,'right_product_literal_checks':rightchecks,'singleton_identity_checks':singleton},indent=2))
