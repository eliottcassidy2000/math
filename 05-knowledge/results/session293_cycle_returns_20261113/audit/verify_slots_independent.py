#!/usr/bin/env python3
"""Independent exact finite-slot and arbitrary-middle substitution controls."""
import json
from itertools import permutations,combinations
from functools import lru_cache
checks=0

def require(q,label):
 global checks
 checks+=1
 if not q:raise RuntimeError(label)
def tt(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def tournaments(n):
 pairs=list(combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  r=[0]*n
  for j,(a,b) in enumerate(pairs):
   if mask>>j&1:r[a]|=1<<b
   else:r[b]|=1<<a
  yield tuple(r)
@lru_cache(None)
def h(r,mask):
 if mask==0:return 1
 d={(1<<i,i):1 for i in range(len(r)) if mask>>i&1}
 for _ in range(1,mask.bit_count()):
  nd={}
  for (used,last),num in d.items():
   avail=r[last]&mask&~used
   while avail:
    bit=avail&-avail;avail-=bit;j=bit.bit_length()-1
    key=(used|bit,j);nd[key]=nd.get(key,0)+num
  d=nd
 return sum(d.values())
def H(r):return h(r,(1<<len(r))-1)
def covers(r,k,mask):
 if k==0:return int(mask==0)
 if mask.bit_count()<k:return 0
 ans=0;sub=mask
 while sub:
  ans+=h(r,sub)*covers(r,k-1,mask^sub);sub=(sub-1)&mask
 return ans

def slotpoly(ext,A):
 # Each chosen gap represents precisely one nonempty path run of M.
 n=len(ext);B=set(range(n))-set(A);result=[0]*(n+2)
 for p in permutations(range(n)):
  legal=[g for g in range(n+1) if (g==0 or p[g-1] in A) and (g==n or p[g] in B)]
  for mask in range(1<<len(legal)):
   occupied={g for j,g in enumerate(legal) if mask>>j&1}
   # Literal adjacency check after insertion of run symbols.
   seq=[]
   for g in range(n+1):
    if g in occupied:seq.append(None)
    if g<n:seq.append(p[g])
   ok=True
   for x,y in zip(seq,seq[1:]):
    if x is None:ok=y in B
    elif y is None:ok=x in A
    else:ok=bool(ext[x]>>y&1)
    if not ok:break
   if ok:result[len(occupied)]+=1
 return result

def substitute(ext,A,M):
 # Retain old exterior labels, append middle block labels.
 n=len(ext);r=list(ext);shifted=[]
 for i in range(n):
  if i in A:r[i]|=((1<<len(M))-1)<<n
 for row in M:shifted.append((row<<n)|sum(1<<i for i in range(n) if i not in A))
 return tuple(r+shifted)
def exterior(bits):
 r=list(tt(6))+[0]
 for i,b in enumerate(bits):
  if b=='1':r[6]|=1<<i
  else:r[i]|=1<<6
 return tuple(r)
X=exterior('101101');Y=exterior('011110')
px=slotpoly(X,set(range(3)));py=slotpoly(Y,set(range(3)))
require(px==[15,31,15]+[0]*6,'X slot coefficients')
require(py==[17,33,15]+[0]*6,'Y slot coefficients')
cores=[r for n in range(1,4) for r in tournaments(n)]+list(tournaments(4))[::7]
rows=[]
for C in cores:
 full=(1<<len(C))-1;hc=H(C);Z=sum(h(C,s)*h(C,full^s) for s in range(full+1))
 cx=H(substitute(X,{0,1,2},C));cy=H(substitute(Y,{0,1,2},C))
 require(cx==hc+15*Z and cy==3*hc+15*Z,'direct arbitrary core scaffold')
 require(Z==2*hc+covers(C,2,full),'two-cover Z identity')
 rows.append([list(C),hc,Z,cx,cy])
# All exterior tournaments through4, both nonempty uniform sides with A initial half;
# representative arbitrary middle TT2/C3, plus empty boundary.
general=0
for n in range(2,5):
 for ext in tournaments(n):
  A=set(range(n//2));p=slotpoly(ext,A)
  amask=sum(1<<i for i in A);bmask=((1<<n)-1)^amask
  HA=h(ext,amask);HB=h(ext,bmask)
  E=[(a,b) for a in A for b in range(n) if b not in A and ext[b]>>a&1]
  WE=sum(h(ext,amask^(1<<a))*h(ext,bmask^(1<<b)) for a,b in E)
  require(p[1]>=HA*HB,'one-slot baseline lower bound')
  require(p[2]>=len(E)*HA*HB,'two-slot feedback lower bound')
  require(not any(p[min(len(A),n-len(A))+2:]),'slot upper bound')
  for M in [(),tt(1),tt(2),(2,4,1)]:
   value=sum(c*covers(M,k,(1<<len(M))-1) for k,c in enumerate(p))
   require(H(substitute(ext,A,M))==value,'general uniform seam exact compiler')
   if M:
    lower=H(substitute(ext,A,tt(len(M))))
    require(value>=lower,'transitive seam lower bound')
    full=(1<<len(M))-1
    Wodd=sum(h(M,S)*h(M,full^S) for S in range(full+1) if S.bit_count()%2)
    require(value>=HA*HB*H(M)+2*WE*Wodd,'general feedback odd-cycle lower bound')
    require(value>=HA*HB+(2**len(M))*WE,'universal activated length budget')
   general+=1
print(json.dumps({'status':'PASS','checks':checks,'X_slots':px,'Y_slots':py,'arbitrary_core_cases':len(cores),'general_cases':general,'scaffold_rows':rows},indent=2,sort_keys=True))
