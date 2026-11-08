"""Independent audit. Python 3 standard library only. Assertions are not used."""
import itertools as it, json, math, hashlib
from pathlib import Path

def check(x,msg):
 if not x: raise RuntimeError(msg)
def mm(a,b):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]
def mv(a,v):return [sum(x*y for x,y in zip(r,v)) for r in a]
def sub(a,b,h):return [[x-h*y for x,y in zip(r,s)] for r,s in zip(a,b)]
def zero(a):return all(not x for r in a for x in r)
def local(adj,out):
 n=len(adj); m=[[0]*6 for _ in range(6)]
 for perm in it.permutations(range(n)):
  for cut in range(n+1):
   p,r=perm[:cut],perm[cut:]
   if any(not(adj[a]>>b&1) for seq in (p,r) for a,b in zip(seq,seq[1:])):continue
   for j in range(6):
    pp,rr=divmod(j,2)
    if not rr and r and r[0] not in out:continue
    np=(2 if p[-1] in out else 1) if p else pp
    nr=int(bool(rr or r));m[2*np+nr][j]+=1
 return m

def stats(adj):
 n=len(adj);full=(1<<n)-1;end=[[0]*n for _ in range(full+1)];h=[1]+[0]*full
 for mask in range(1,full+1):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v)
    end[mask][v]=sum(end[rest][u] for u in range(n) if adj[u]>>v&1) if rest else 1
  h[mask]=sum(end[mask])
 # Compute first endpoint counts independently by reversing arcs.
 rev=[sum(1<<j for j in range(n) if adj[j]>>i&1) for i in range(n)]
 first=[[0]*n for _ in range(full+1)]
 for mask in range(1,full+1):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v);first[mask][v]=sum(first[rest][u] for u in range(n) if adj[v]>>u&1) if rest else 1
 alpha=[sum(first[s][v]*h[full^s] for s in range(1,full+1)) for v in range(n)]
 return h, end, first, alpha

def module_stats(adj):
 h,e,f,a=stats(adj);full=len(h)-1
 return h[-1],sum(h[s]*h[full^s] for s in range(full+1)),a

def catmatrix(weights):
 m=[[0]*6 for _ in range(6)]
 for j in range(6):
  p,r=divmod(j,2)
  for (pp,rr),value in weights.items():
   if not r and rr==2:continue
   m[2*(pp or p)+int(bool(r or rr))][j]+=value
 return m
W=[[0,0],[1,0],[0,0],[0,1],[0,0],[-1,-1]]
V=[r+[int(i==2)-int(i==4)] for i,r in enumerate(W)]
def in_v(C):return all(C[0][j]==0 and C[2][j]+C[4][j]==0 and C[1][j]+C[3][j]+C[5][j]==0 for j in range(len(C[0])))
def in_w(C):return in_v(C) and not any(C[i][j] for i in (2,4) for j in range(len(C[0])))
for pp,rr in it.product(range(3),repeat=2):
 if pp or rr:check(in_v(mm(catmatrix({(pp,rr):1}),V)),f'category V invariant {(pp,rr)}')
U=[]
for k in range(4):
 x,y,a,b=[int(j==k) for j in range(4)]
 F=catmatrix({(0,1):x,(0,2):-x,(1,0):y,(2,0):-y,(1,1):a,(1,2):b,(2,1):-x-a,(2,2):x-b})
 check(in_v(F),'universal image V');check(zero(mm(F,W)),'universal kills W');U.append(F)
S=[local([0],set() if b==0 else {0}) for b in range(2)]
for s in S:check(in_w(mm(s,V)),'singleton V to W')
check(any(not zero(mm(f,g)) for f in U for g in U),'hostile empty spacer fails universally')
B=[40,1,3,6,15,30];h,q,alpha=module_stats(B)
check((h,q,alpha)==(19,256,[21,20,29,35,73,59]),'literal block statistics')
masks=[17,40];check([sum(alpha[i] for i in range(6) if m>>i&1) for m in masks]==[94,94],'same alpha')
M=[local(B,{i for i in range(6) if m>>i&1}) for m in masks]
F=sub(M[1],M[0],1);check(not zero(F),'different transfers');check(in_v(F) and zero(mm(F,W)),'physical defect')
for m in M:check(mm(m,W)==[[19*z for z in r] for r in W],'physical W eigenvalue')
I=[[int(i==j) for j in range(6)] for i in range(6)]
A=mm(S[1],M[0]);pol=mm(mm(mm(A,sub(A,I,19)),sub(A,I,19)),sub(A,I,512));check(zero(pol),'exact matrix recurrence')
v=[1,0,0,0,0,0];ell=[1,1,1,1,0,0]
def scalar(v):return sum(x*y for x,y in zip(ell,v))
def div(a,b):check(a%b==0,'integer expression');return a//b
def baseline(t):return div(19950*512**t+(13862*t-2289)*19**t,15979)
def weights(t):return [div(-(8400*512**i*19**(t-1-i)+(474 if i==t-1 else 1460)*19**(t-1)),493) for i in range(t)]
x=mv(S[1],v);initial=[]
for i in range(4):
 b=scalar(mv(S[0],mv(M[0],x)));g=scalar(mv(S[1],mv(F,x)));k=scalar(mv(S[0],mv(F,x)))
 check(b==baseline(i+1),'baseline recurrence initials');check(g==div(-(8400*512**i+1460*19**i),493),'interior initials');check(k==div(-(8400*512**i+474*19**i),493),'last initials');initial.append([b,g,k]);x=mv(A,x)
# Symbolic characteristic roots and repeated root suffice for all formula sequences.
p=lambda z:z*(z-19)**2*(z-512)
check(p(19)==p(512)==0,'recurrence roots')
# Derivative at repeated root, by exact product rule.
check((19-19)**2*(19-512)+19*2*(19-19)*(19-512)+19*(19-19)**2==0,'repeated root')
check(8400*(512-20*19)-19*(20*1460-474)>0,'all-order superincreasing inequality')
check(math.gcd(493,19*5)==1 and 8400%19!=0 and 512%19!=0 and 474%5!=0 and 19%5!=0,'gcd prime exclusions')
check(8400%4==0 and 474%4==2 and 19%2==1 and 493%2==1,'gcd two-adic certificate')
check(19950%2==0 and 13862%2==0 and 2289%2==1 and 15979%2==1,'baseline odd certificate')
check(math.gcd(653,635)==1,'t1 gcd')
order4=[]
for code in range(64):
 a=[0]*4
 for k,(u,w) in enumerate(it.combinations(range(4),2)):
  if code>>k&1:a[u]|=1<<w
  else:a[w]|=1<<u
 h,q,al=module_stats(a)
 if (h,q)==(5,32):
  sums=[sum(al[i] for i in range(4) if mask>>i&1) for mask in range(16)]
  check(len(set(sums))==16,'order4 no same-alpha masks');order4.append(a)
check(len(order4)==24,'order4 labeled count')
def graph(bits):
 blocks=[[0]];outs=[1]
 for i,b in enumerate(bits):
  blocks.append(B);outs.append(masks[b]);blocks.append([0]);outs.append(int(i<len(bits)-1))
 n=sum(map(len,blocks));a=[0]*(n+1);off=0
 for block,out in zip(blocks,outs):
  for u,mask in enumerate(block):
   a[off+u]=(mask<<off)|sum(1<<j for j in range(off+len(block),n))
   if out>>u&1:a[n]|=1<<(off+u)
   else:a[off+u]|=1<<n
  off+=len(block)
 return a

def strong(a):
 for v in range(len(a)):
  reach=1<<v
  while True:
   new=reach
   for i in range(len(a)):
    if reach>>i&1:new|=a[i]
   if new==reach:break
   reach=new
  if reach!=(1<<len(a))-1:return False
 return True

def literalstate(a):
 h,e,f,al=stats(a);h2,_,_,al2=stats([x&((1<<(len(a)-1))-1) for x in a[:-1]])
 r=len(a)-1
 return [h[-1],e[-1][r],h2[-1]+sum(al2[i] for i in range(r) if a[r]>>i&1),h2[-1]]
witnesses=[]
for t in (1,2):
 vals=[]
 for bits in it.product((0,1),repeat=t):
  a=graph(bits);state=literalstate(a);expected=[baseline(t)+sum(b*w for b,w in zip(bits,weights(t))),19**t,div(700*512**t-207*19**t,493),19**t]
  check(state==expected,'literal whole graph state');check(strong(a),'literal strong');check(len(a)==7*t+2,'literal order')
  witnesses.append({'bits':bits,'adj':a,'state':state});vals.append(state[0])
 check(len(set(vals))==2**t and math.gcd(*vals)==1,'literal fiber')
for t in range(1,11):
 vals=[]
 for bits in it.product((0,1),repeat=t):
  x=mv(S[1],v)
  for i,b in enumerate(bits):x=mv(S[1 if i<t-1 else 0],mv(M[b],x))
  d=scalar(x);check(d==baseline(t)+sum(b*w for b,w in zip(bits,weights(t))),'all choice transfer regression');vals.append(d)
 check(len(set(vals))==2**t and math.gcd(*vals)==1,'finite fiber regression')
# A deliberately changed root mask must change U; no unnoticed fixed-alpha assumption.
a=graph([0]);a[-1]^=1<<1;a[1]^=1<<(len(a)-1)
check(literalstate(a)[2]!=witnesses[0]['state'][2],'hostile altered root mask changes sidecar')
print(json.dumps({'passed':True,'module':{'adj':B,'H':19,'Q':256,'alpha':alpha,'root_masks':masks},'matrices':M,'defect':F,'recurrence_initials':initial,'order4_injective_blocks':len(order4),'literal_witnesses':witnesses,'finite_transfer_t_max':10,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'proof_sha256':hashlib.sha256(Path(__file__).with_name('PROOF.md').read_bytes()).hexdigest()},indent=2))
