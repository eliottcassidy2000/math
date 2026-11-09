"""Standalone exact audit; Python standard library, no external inputs or asserts."""
import itertools as it, json, math, hashlib
from pathlib import Path

def require(p,s):
 if not p: raise RuntimeError(s)
def mm(a,b):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]
def mv(a,v):return [sum(x*y for x,y in zip(r,v)) for r in a]
def sub(a,b):return [[x-y for x,y in zip(r,s)] for r,s in zip(a,b)]
def zero(a):return all(x==0 for r in a for x in r)
def cm(w):
 m=[[0]*6 for _ in range(6)]
 for j in range(6):
  p,r=divmod(j,2)
  for (pp,rr),v in w.items():
   if not r and rr==2:continue
   m[2*(pp or p)+int(bool(r or rr))][j]+=v
 return m
def categories(adj,mask):
 n=len(adj);w={c:0 for c in it.product(range(3),repeat=2) if c!=(0,0)}
 for p in it.permutations(range(n)):
  for k in range(n+1):
   a,b=p[:k],p[k:]
   if any(not(adj[u]>>v&1) for seq in (a,b) for u,v in zip(seq,seq[1:])):continue
   pp=(2 if mask>>a[-1]&1 else 1) if a else 0
   rr=(1 if mask>>b[0]&1 else 2) if b else 0
   w[pp,rr]+=1
 return w
def local(adj,mask):return cm(categories(adj,mask))
def data(adj):
 n=len(adj);N=1<<n;h=[1]+[0]*(N-1);end=[[0]*n for _ in range(N)];first=[[0]*n for _ in range(N)]
 for s in range(1,N):
  for v in range(n):
   if s>>v&1:
    r=s^(1<<v)
    end[s][v]=sum(end[r][u] for u in range(n) if adj[u]>>v&1) if r else 1
    first[s][v]=sum(first[r][u] for u in range(n) if adj[v]>>u&1) if r else 1
  h[s]=sum(end[s])
 alpha=[sum(first[s][v]*h[(N-1)^s] for s in range(1,N)) for v in range(n)]
 return h[-1],sum(h[s]*h[(N-1)^s] for s in range(N)),alpha,end[-1],first[-1]
def graph(blocks,masks):
 n=sum(map(len,blocks));a=[0]*(n+1);off=0
 for B,m in zip(blocks,masks):
  for u,w in enumerate(B):
   a[off+u]=(w<<off)|sum(1<<j for j in range(off+len(B),n))
   if m>>u&1:a[n]|=1<<(off+u)
   else:a[off+u]|=1<<n
  off+=len(B)
 return a
def hamilton(a):
 n=len(a);dp=[{} for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for u,c in dp[s].items():
   todo=a[u]&~s
   while todo:
    b=todo&-todo;v=b.bit_length()-1;todo-=b;d=dp[s|b];d[v]=d.get(v,0)+c
 return sum(dp[-1].values())
def strong(a):
 for i in range(len(a)):
  s=1<<i
  while True:
   z=s
   for j in range(len(a)):
    if s>>j&1:z|=a[j]
   if z==s:break
   s=z
  if s!=(1<<len(a))-1:return False
 return True
W=[[0,0],[1,0],[0,0],[0,1],[0,0],[-1,-1]]
V=[r+[int(i==2)-int(i==4)] for i,r in enumerate(W)]
def inW(a):return all(a[i][j]==0 for i in (0,2,4) for j in range(len(a[0]))) and all(sum(a[i][j] for i in (1,3,5))==0 for j in range(len(a[0])))
def inV(a):return all(a[0][j]==0 and a[2][j]+a[4][j]==0 and a[1][j]+a[3][j]+a[5][j]==0 for j in range(len(a[0])))
e=[0,0,0,1,0,-1];rho=[1,0,1,0,1,0];rankone=[[a*b for b in rho] for a in e]
# Eight category generators prove the strengthened annihilation universally.
for c in it.product(range(3),repeat=2):
 if c!=(0,0):require(inW(mm(cm({c:1}),V)),('M V in W',c))
Fs=[]
for k in range(4):
 x,y,a,b=[int(j==k) for j in range(4)]
 F=cm({(0,1):x,(0,2):-x,(1,0):y,(2,0):-y,(1,1):a,(1,2):b,(2,1):-x-a,(2,2):x-b})
 require(inV(F) and zero(mm(F,W)),'defect image and kernel');Fs.append(F)
for i,F in enumerate(Fs):
 for j,G in enumerate(Fs):
  require(mm(F,G)==[[int(i==0 and j==1)*z for z in r] for r in rankone],'adjacent product')
  for K in Fs:require(zero(mm(F,mm(G,K))),'third order annihilation')
B=[40,1,3,6,15,30];H,Q,alph,ends,starts=data(B)
require((H,Q,alph)==(19,256,[21,20,29,35,73,59]),'H19 data')
fib={}
for m in range(64):fib.setdefault(sum(a for i,a in enumerate(alph) if m>>i&1),[]).append(m)
doubles=sorted((a,m) for a,m in fib.items() if len(m)>1)
require(doubles==[(94,[17,40]),(108,[24,38]),(114,[19,42]),(123,[21,44]),(129,[25,39]),(143,[23,46])],'complete fixed alpha fibers')
S0,S1=local([0],0),local([0],1);v=[1,0,0,0,0,0];ell=[1,1,1,1,0,0]
def scalar(x):return sum(a*b for a,b in zip(ell,x))
def evalword(ms):
 x=mv(S1,v)
 for m in ms:x=mv(m,x)
 return scalar(mv(S0,x))
mat={m:local(B,m) for _,pair in doubles for m in pair}
pairs=[]
for a,(m,n) in doubles:
 x=sum(starts[i]*(((n>>i)&1)-((m>>i)&1)) for i in range(6))
 y=-sum(ends[i]*(((n>>i)&1)-((m>>i)&1)) for i in range(6))
 F=sub(mat[n],mat[m]);require(mm(F,F)==[[x*y*z for z in r] for r in rankone],'physical self interaction')
 pairs.append({'alpha':a,'masks':[m,n],'x':x,'y':y})
# Mixed menu positions, no separators. Check ALL choices and nearest-neighbor law.
reg=[]
for t in range(1,10):
 ps=[pairs[i%len(pairs)] for i in range(t)]
 base=[mat[p['masks'][0]] for p in ps];d0=evalword(base);linear=[]
 for i,p in enumerate(ps):
  z=list(base);z[i]=mat[p['masks'][1]];linear.append(evalword(z)-d0)
 edge=[ps[i]['y']*ps[i+1]['x']*19**(t-2) for i in range(t-1)]
 vals=[]
 for bs in it.product((0,1),repeat=t):
  d=evalword([mat[p['masks'][b]] for p,b in zip(ps,bs)])
  predicted=d0+sum(b*l for b,l in zip(bs,linear))+sum(bs[i]*bs[i+1]*edge[i] for i in range(t-1))
  require(d==predicted,'nearest neighbor formula');require(d>=256**t,'all-mask entropy lower bound');vals.append(d)
 # Independent O(t k^2) extremal dynamic programming.
 lo={0:d0,1:d0+linear[0]};hi=dict(lo)
 for i in range(1,t):
  lo={b:min(z+b*linear[i]+a*b*edge[i-1] for a,z in lo.items()) for b in (0,1)}
  hi={b:max(z+b*linear[i]+a*b*edge[i-1] for a,z in hi.items()) for b in (0,1)}
 require(min(lo.values())==min(vals) and max(hi.values())==max(vals),'chain extrema')
 reg.append({'t':t,'choices':2**t,'distinct':len(set(vals)),'minimum':min(vals),'maximum':max(vals)})
# Literal whole-graph independent path DP for genuinely interacting adjacent defects.
literal=[]
for bs in it.product((0,1),repeat=2):
 masks=[(24,38)[b] for b in bs];G=graph([[0],B,B,[0]],[1]+masks+[0]);d=hamilton(G)
 require(d==evalword([mat[m] for m in masks]) and strong(G),'literal interacting graph')
 literal.append({'bits':bs,'D':d,'adj':G,'sidecar':[19**2,19**2+2*256**2+2*108*(256+19),19**2]})
require(literal[3]['D']-literal[2]['D']-literal[1]['D']+literal[0]['D']==20,'physical quadratic coefficient20')
# Explicit entropy injection: all 256 ordered covers, all 64 root masks.
covers=[]
for perm in it.permutations(range(6)):
 for cut in range(7):
  L,R=perm[:cut],perm[cut:]
  if all(B[u]>>v&1 for seq in (L,R) for u,v in zip(seq,seq[1:])):covers.append((L,R))
require(len(covers)==256,'ordered cover count for injection')
for mask in range(64):
 G=graph([[0],B,[0]],[1,mask,0]);paths=set()
 for L,R in covers:
  path=tuple(u+1 for u in L)+(7,8,0)+tuple(u+1 for u in R)
  require(len(set(path))==9 and all(G[u]>>v&1 for u,v in zip(path,path[1:])),'literal injected path')
  paths.add(path)
 require(len(paths)==256,'literal injection uniqueness')
# Noncoprime physical carry hostile, exact masks reconstructed from subset weights.
C=[16,1,3,7,14];h,q,a,_,_=data(C);require((h,q,a)==(9,84,[9,10,16,25,15]),'carry module')
lookup={sum(a[i] for i in range(5) if m>>i&1):m for m in range(32)}
carry=[]
for alphas in ((10,56),(16,0)):
 masks=[lookup[z] for z in alphas];G=graph([[0],C,C,[0]],[1]+masks+[0]);d=hamilton(G)
 require(strong(G),'carry strong graph')
 U=h*h+2*q*q+2*(q*alphas[0]+h*alphas[1]);carry.append({'alphas':alphas,'masks':masks,'D':d,'sidecar':[h*h,U,h*h],'adj':G})
require(carry[0]['sidecar']==carry[1]['sidecar'] and carry[0]['D']!=carry[1]['D'],'physical carry equal sidecar distinct scalar')
# Exhaust exact carry equivalence on full bounded digit alphabets, not only realizable subsets.
carry_tests=0
for H,Q,t in ((1,2,4),(2,6,3),(3,12,3),(4,10,3),(6,15,3),(9,84,2),(19,256,2)):
 g=math.gcd(H,Q);h,q=H//g,Q//g
 # Enumerate differences, sufficient since every delta in range is difference of digits.
 for delta in it.product(range(-(Q-H),Q-H+1),repeat=t):
  eq=sum(d*h**i*q**(t-1-i) for i,d in enumerate(delta))==0
  cs=[];c=0;valid=True
  for d in delta:
   numerator=d+q*c
   if numerator%h:valid=False;break
   c=numerator//h;cs.append(c)
  valid=valid and c==0
  require(eq==valid,'carry iff')
  if eq:require(all(abs(c)<g for c in cs),'carry alphabet bound')
  carry_tests+=1
print(json.dumps({'passed':True,'H19_module':{'adj':B,'H':19,'Q':256,'alpha':alph,'end':ends,'start':starts},'all_nonsingleton_alpha_fibers':pairs,'quadratic_regression':reg,'literal_adjacent_witnesses':literal,'noncoprime_carry_witnesses':carry,'carry_difference_words_tested':carry_tests,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'proof_sha256':hashlib.sha256(Path(__file__).with_name('PROOF.md').read_bytes()).hexdigest()},indent=2))
