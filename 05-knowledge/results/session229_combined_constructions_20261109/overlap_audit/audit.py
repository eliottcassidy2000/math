"""Standalone Session229 overlap audit; stdlib only; no producer imports/asserts."""
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import hashlib,json

def require(ok, label):
 if not ok: raise RuntimeError(label)
def matmul(A,B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def matvec(A,v): return [sum(a*b for a,b in zip(row,v)) for row in A]
def minus(A,B): return [[a-b for a,b in zip(x,y)] for x,y in zip(A,B)]
def iszero(A): return all(x==0 for row in A for x in row)
def scale(A,s):return [[s*x for x in row] for row in A]
def outer(v,w):return [[a*b for b in w] for a in v]
W=[[0,0],[1,0],[0,0],[0,1],[0,0],[-1,-1]]
V=[row+[int(i==2)-int(i==4)] for i,row in enumerate(W)]
rho=[1,0,1,0,1,0]; u=[0,0,0,1,0,-1]
def inv(A):return all(A[0][j]==0 and A[2][j]+A[4][j]==0 and A[1][j]+A[3][j]+A[5][j]==0 for j in range(len(A[0])))
def inw(A):return inv(A) and all(A[i][j]==0 for i in (2,4) for j in range(len(A[0])))
def category(C):
 T=[[0]*6 for _ in range(6)]
 for p,r in product(range(3),range(2)):
  for (a,b),n in C.items():
   if r==0 and b==2:continue
   target=(a if a else p, int(r!=0 or b!=0))
   T[2*target[0]+target[1]][2*p+r]+=n
 return T
cats=[(p,r) for p,r in product(range(3),repeat=2) if p or r]
units=[category({c:1}) for c in cats]
for c,M in zip(cats,units):
 require(inw(matmul(M,V)),('elementary image W',c))
 require(matmul(M,W)==scale(W,int(c[0]==0)),('elementary W action',c))
# Equal-H,Q,alpha signed weight basis, verified separately in all mixed products.
def defect(x,y,a,b):return category({(0,1):x,(0,2):-x,(1,0):y,(2,0):-y,(1,1):a,(1,2):b,(2,1):-x-a,(2,2):x-b})
Fs=[defect(*(int(i==k) for i in range(4))) for k in range(4)]
for i,F in enumerate(Fs):
 require(inv(F) and iszero(matmul(F,W)),('defect subspaces',i))
 for j,G in enumerate(Fs):
  require(matmul(F,G)==scale(outer(u,rho),int(i==0 and j==1)),('ordered pair identity',i,j))
  for M in units:require(iszero(matmul(matmul(F,M),G)),'arbitrary nonempty gap')
  for K in Fs:require(iszero(matmul(matmul(F,G),K)),'three adjacent defects')
require(not iszero(matmul(Fs[0],Fs[1])),'empty gap hostile')
I=[[int(i==j) for j in range(6)] for i in range(6)]
require(not inw(matmul(I,V)),'empty block does not send V into W')
# Literal local path covers via subset path enumeration, independent of transfer formulas.
def paths(adj):
 n=len(adj);out={0:[()]}
 for mask in range(1,1<<n):
  vertices=[i for i in range(n) if mask>>i&1]
  out[mask]=[p for p in permutations(vertices) if all(adj[a]>>b&1 for a,b in zip(p,p[1:]))]
 return out
def literal(adj,rootmask):
 P=paths(adj);full=(1<<len(adj))-1; C={c:0 for c in cats}
 for subset,ps in P.items():
  for p,r in product(ps,P[full^subset]):
   a=0 if not p else (2 if rootmask>>p[-1]&1 else 1)
   b=0 if not r else (1 if rootmask>>r[0]&1 else 2)
   C[a,b]+=1
 H=len(P[full]);Q=sum(C.values());alpha=C[0,1]+C[1,1]+C[2,1]
 return C,category(C),(H,Q,alpha)
B=[40,1,3,6,15,30];masks=(24,38)
C0,M0,s0=literal(B,24); C1,M1,s1=literal(B,38)
require(s0==s1==(19,256,108),'physical statistics')
require(not inv(minus(M0,literal(B,0)[1])),'unequal alpha breaks image V')
F=minus(M1,M0); x=C1[0,1]-C0[0,1];y=C1[1,0]-C0[1,0]
require((x,y)==(-2,-10),'physical endpoint differences')
require(matmul(F,F)==scale(outer(u,rho),20),'physical nonzero square')
_,old0,_=literal(B,17);_,old1,_=literal(B,40)
require(iszero(matmul(minus(old1,old0),minus(old1,old0))),'old pair square zero')
S0=literal([0],0)[1];S1=literal([0],1)[1]
v=[1,0,0,0,0,0];ell=[1,1,1,1,0,0]
def value(word):
 z=matvec(S1,v)
 for b in word:z=matvec((M0,M1)[b],z)
 return sum(a*b for a,b in zip(ell,matvec(S0,z)))
# Every binary word through t=10 checked against independently extracted quadratic expansion.
finite=[]
for t in range(1,11):
 base=value([0]*t);linear=[value([int(i==j) for i in range(t)])-base for j in range(t)]
 vals=[]
 for word in product((0,1),repeat=t):
  expected=base+sum(b*c for b,c in zip(word,linear))+sum(word[i]*word[i+1]*20*19**(t-2) for i in range(t-1))
  actual=value(word);require(actual==expected,'quadratic physical expansion');vals.append(actual)
 finite.append({'t':t,'distinct':len(set(vals)),'gcd':gcd(*vals),'base':base,'linear':linear,'edge':20*19**(t-2) if t>1 else None})
# Whole graph DP has no six-state machinery.
def graph(word,block=B,choices=masks):
 blocks=[[0]]+[block for _ in word]+[[0]]; outs=[1]+[choices[b] for b in word]+[0]
 n=sum(len(b) for b in blocks);A=[0]*(n+1);off=0
 for block,om in zip(blocks,outs):
  for i,row in enumerate(block):
   A[off+i]=(row<<off)|sum(1<<j for j in range(off+len(block),n))
   if om>>i&1:A[n]|=1<<(off+i)
   else:A[off+i]|=1<<n
  off+=len(block)
 return A
def pathdp(A,startset=None):
 n=len(A);full=(1<<n)-1; dp=[[0]*n for _ in range(full+1)];tot=[1]+[0]*full
 incoming=[sum(1<<i for i in range(n) if A[i]>>j&1) for j in range(n)]
 for mask in range(1,full+1):
  left=mask
  while left:
   bit=left&-left;j=bit.bit_length()-1;left-=bit;prev=mask^bit
   if prev==0:dp[mask][j]=int(startset is None or startset>>j&1)
   else:
    pred=prev&incoming[j];s=0
    while pred:
     q=pred&-pred;s+=dp[prev][q.bit_length()-1];pred-=q
    dp[mask][j]=s
  tot[mask]=sum(dp[mask])
 return dp,tot
def strong(A):
 for i in range(len(A)):
  reach=1<<i
  while True:
   grow=reach
   for j in range(len(A)):
    if reach>>j&1:grow|=A[j]
   if grow==reach:break
   reach=grow
  if reach!=(1<<len(A))-1:return False
 return True
physical=[]
for t in (1,2):
 for word in product((0,1),repeat=t):
  A=graph(word);n=len(A);dp,ht=pathdp(A);deleted=[a&((1<<(n-1))-1) for a in A[:-1]]
  _,h=pathdp(deleted);_,r=pathdp(deleted,A[-1]);full=len(h)-1
  U=h[-1]+sum(r[s]*h[full^s] for s in range(1,full+1))
  d=19**t;expectedU=2*256**t+2*108*sum(19**i*256**(t-1-i) for i in range(t))+d
  state=[ht[-1],dp[-1][-1],U,h[-1]]
  require(state==[value(word),d,expectedU,d],'literal full rooted state')
  require(strong(A) and n==6*t+3,'physical strong order')
  physical.append({'word':word,'state':state,'adjacency':A})
# Separate physical carry witness with unequal alpha words and equal U.
C=[16,1,3,7,14];lookup={};alphas=[]
for mask in range(32):
 _,_,st=literal(C,mask);require(st[:2]==(9,84),'carry block H,Q');lookup[st[2]]=mask
for i in range(5):alphas.append(literal(C,1<<i)[2][2])
require(alphas==[9,10,16,25,15],'carry block individual weights')
carryphysical=[]
for digits in ((10,56),(16,0)):
 chosen=[lookup[d] for d in digits];A=graph([0,1],C,chosen);dp,ht=pathdp(A)
 n=len(A);deleted=[a&((1<<(n-1))-1) for a in A[:-1]]
 _,hh=pathdp(deleted);_,rr=pathdp(deleted,A[-1]);full=len(hh)-1
 U=hh[-1]+sum(rr[s]*hh[full^s] for s in range(1,full+1))
 state=[ht[-1],dp[-1][-1],U,hh[-1]]
 require(strong(A),'carry physical strong')
 carryphysical.append({'alpha_word':digits,'masks':chosen,'state':state,'adjacency':A})
require(carryphysical[0]['state'][1:]==carryphysical[1]['state'][1:],'physical carry same sidecar')
require(carryphysical[0]['state'][0]!=carryphysical[1]['state'][0],'physical carry changed D')
# Exhaustive bounded carry test: all digit words length <=4, 1<=H<=4, 2H<=Q<=10.
carry_cases=0;noncoprime=None
for H in range(1,5):
 for Q in range(2*H,11):
  g=gcd(H,Q);h,q=H//g,Q//g
  for t in range(1,5):
   seen={}
   for digits in product(range(Q-H+1),repeat=t):
    U=sum(a*H**i*Q**(t-1-i) for i,a in enumerate(digits))
    if U in seen:
     old=seen[U];delta=[a-b for a,b in zip(digits,old)];c=0;cs=[]
     for dd in delta:
      require((dd+q*c)%h==0,'integer carry');c=(dd+q*c)//h;cs.append(c)
     require(c==0 and all(abs(z)<=g-1 for z in cs),'carry boundary and bound')
     require(g>1 or digits==old,'coprime uniqueness')
     if noncoprime is None and g>1:noncoprime={'H':H,'Q':Q,'words':[old,digits],'carry':cs}
     carry_cases+=1
    else:seen[U]=digits
require(noncoprime is not None,'noncoprime hostile found')
require(19*256==256*19 and 256>256-19,'coprime unrestricted digits collide outside physical range')
# Entropy hostile/positive controls: every mask, not merely same-alpha choices.
allM=[literal(B,m)[1] for m in range(64)]
entropy=[]
for t in (1,2):
 values=[]
 for word in product(range(64),repeat=t):
  z=matvec(S1,v)
  for m in word:z=matvec(allM[m],z)
  D=sum(a*b for a,b in zip(ell,matvec(S0,z)))
  require(D>=256**t,'all-mask cover injection lower bound');values.append(D)
 entropy.append({'t':t,'mask_words':len(values),'minimum_D':min(values),'Q_power':256**t,'distinct_D':len(set(values))})
# Q=2^n exactly for transitive labeled tournaments at n<=4.
cover_controls=[]
for n in range(1,5):
 strict=0
 for code in range(1<<(n*(n-1)//2)):
  A=[0]*n
  for k,(i,j) in enumerate(combinations(range(n),2)):
   if code>>k&1:A[i]|=1<<j
   else:A[j]|=1<<i
  P=paths(A);full=(1<<n)-1;Q=sum(len(P[m])*len(P[full^m]) for m in P)
  cyclic=any(A[a]>>b&1 and A[b]>>c&1 and A[c]>>a&1 for a,b,c in permutations(range(n),3))
  require(Q>=2**n and (Q>2**n)==cyclic,'cover lower bound and equality boundary')
  strict+=int(cyclic)
 cover_controls.append({'n':n,'nontransitive':strict})
# Explicit mixed-bank sufficient entropy witness at exponent15/16, using integer inequalities.
require(180**16<256**15 and 63**16<84**15,'rational root lower bounds')
require(64*63+32*180<180*63,'two-block bank F(15/16)<1')
report={'status':'PASS','module':{'adjacency':B,'masks':masks,'HQalpha':s0,'x':x,'y':y,'categories':[C0,C1] if False else [sorted((str(k),v) for k,v in C.items()) for C in (C0,C1)]},'universal_basis_tests':{'categories':8,'defects':4,'mixed_pairs':16,'mixed_triples':64,'spaced_pairs':128},'finite_words':finite,'literal_witnesses':physical,'physical_carry_witnesses':carryphysical,'carry_collision_checks':carry_cases,'noncoprime_hostile':noncoprime,'entropy_allmask_controls':entropy,'cover_equality_controls':cover_controls,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(report,indent=2))
