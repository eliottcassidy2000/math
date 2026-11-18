from unit_gain_replay import state,forward,H,K
from itertools import combinations

def candidate(M,g):
 out=[]
 for h in range(1,M//3+1,2):
  if (M-2*h-g)%3:continue
  eta=(M-3*h)//2;d=(M-2*h-g)//3
  if 0<=eta<=d<=h and d>=1:out.append((h,d-eta,h+eta,d))
 return out
for A,r,M,g,p in [((8,1,3,6),0,15,2,(5,1,5,1)),((20,1,2,7,14),2,31,-2,(9,3,11,5))]:
 assert state(A,r)==p and p in candidate(M,g)
 assert H(forward(A,r,[0,1]))==M
 assert H(forward(A,r,[1,1]))==M+g
print('target-selected physical controls: PASS')
checked=0
for n in range(2,6):
 pairs=list(combinations(range(n),2))
 for code in range(1<<len(pairs)):
  A=[0]*n
  for k,(i,j) in enumerate(pairs):
   if code>>k&1:A[i]|=1<<j
   else:A[j]|=1<<i
  for r in range(n):
   if A[r].bit_count()!=1:continue
   h,e,t,d=state(A,r);eta=t-h
   assert e==d-eta and 0<=eta<=d<=h
   assert H(forward(A,r,[0,1]))==3*h+2*eta
   assert H(forward(A,r,[1,1]))==4*h+4*eta-3*d
   checked+=1
for M in range(3,301,2):
 for g in(-2,2):
  expected=[]
  for h in range(1,M//3+1,2):
   if 5*h<max(M+2*g,M-g) or (h-g+M)%3:continue
   eta=(M-3*h)//2;d=(M-2*h-g)//3
   if d>=1:expected.append((h,d-eta,h+eta,d))
  assert candidate(M,g)==expected
for w in range(1,65536):assert (K(w)+(w&1)==1)==(w==1)
print('sole-exit profiles and exact inverse constraints: PASS',checked)
