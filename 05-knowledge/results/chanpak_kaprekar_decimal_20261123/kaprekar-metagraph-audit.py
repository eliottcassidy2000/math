from itertools import combinations, permutations
from collections import Counter
from functools import lru_cache

def H(rows):
 n=len(rows);L=1<<n;d=[[0]*n for _ in range(L)]
 for i in range(n):d[1<<i][i]=1
 for m,row in enumerate(d):
  for i,x in enumerate(row):
   if not x:continue
   v=rows[i]&~m
   while v:
    b=v&-v;v-=b;d[m|b][b.bit_length()-1]+=x
 return sum(d[-1])

def strong(R):
 n=len(R)
 for v in range(n):
  u=1<<v
  for _ in R:
   for i,r in enumerate(R):
    if u>>i&1:u|=r
  if u!=(1<<n)-1:return False
 return True

def anchor(n):
 edges=[(i,j) for i,j in combinations(range(n),2) if j>i+1]
 for mask in range(1<<len(edges)):
  R=[1<<(i+1) if i+1<n else 0 for i in range(n)]
  for k,(i,j) in enumerate(edges):
   if mask>>k&1:R[i]|=1<<j
   else:R[j]|=1<<i
  yield tuple(R)

def flip(R,e):
 i,j=e; R=list(R);R[i]^=1<<j;R[j]^=1<<i;return tuple(R)

def iso(R,S):
 n=len(R)
 return any(all((R[i]>>j&1)==(S[P[i]]>>P[j]&1) for i,j in combinations(range(n),2))
            for P in permutations(range(n)))

def ordinal(blocks):
 R=[]
 for A in blocks:
  n=len(R);R=[r|(((1<<len(A))-1)<<n) for r in R]+[r<<n for r in A]
 return tuple(R)

def low(n):
 nine=[A for A in anchor(5) if H(A)==9]
 assert len(nine)==18
 cores=[(6,4,8,19,7),(2,4,9,19,7)]
 assert not iso(*cores)
 assert all(any(iso(A,B) for B in cores) for A in nine)
 for A in cores:
  assert H(A)==9 and strong(A)
  for u in range(n-4):yield ordinal([(0,)]*u+[A]+[(0,)]*(n-5-u))
 T=(2,4,1)
 for u in range(n-5):
  for v in range(n-5-u):
   yield ordinal([(0,)]*u+[T]+[(0,)]*v+[T]+[(0,)]*(n-6-u-v))

if __name__=='__main__':
 assert max(map(H,anchor(4)))==5
 assert max(map(H,anchor(5)))==15
 assert not any(H(A)==63 for A in anchor(7))
 follows={}
 for n in (8,9):
  E=list(combinations(range(n),2)); lows=list(low(n)); flows=[]
  for A in lows:
   assert H(A)==9
   for e in E:
    B=flip(A,e)
    if H(B)==81 and strong(B):flows.append((A,B,e))
  assert (len(lows),len(flows))=={8:(14,12),9:(20,13)}[n]
  for A,B,e in flows: assert not any(H(flip(B,f))==63 for f in E)
  follows[n]=flows
  print(n,len(lows),len(flows))
 R=(486,484,488,499,487,448,384,256,0)
 reads=[H(R)]
 for e in ((2,8),(2,4),(7,8)):
  R=flip(R,e);reads.append(H(R))
 assert reads==[9,81,115,63]
 print(reads)
 # All possible height-27 targets on eight vertices. Explicitly preserve
 # duplicate Hamilton-anchored presentations to avoid a quotient dependency.
 A=(6,4,8,19,7); B=(2,4,9,19,7); C=(2,4,1); Z=(0,)
 P=[ordinal(x) for K in (A,B) for x in ((K,C),(C,K))]
 for n,N in ((6,18),(7,117)):
  anchored=[R for R in anchor(n) if H(R)==27]
  assert len(anchored)==N
  for R in anchored:
   if not strong(R):continue
   for u in range(9-n):P.append(ordinal([Z]*u+[R]+[Z]*(8-n-u)))
 assert len(P)==220
 for R in P:
  assert len(R)==8 and H(R)==27
  assert not any(H(flip(R,e))==63 for e in combinations(range(8),2))
 print('no 63↔27, 220 complete order-eight presentations')
 # Optimum loop and labelled closure.
 A=(486,484,488,499,487,448,384,256,0)
 R=A; reads=[H(R)]
 for e in ((2,8),(2,4),(7,8),(6,8),(2,4),(6,8),(2,8)):
  R=flip(R,e);reads.append(H(R))
 p=(0,1,2,3,4,5,6,8,7)
 assert all((A[i]>>j&1)==(R[p[i]]>>p[j]&1) for i,j in combinations(range(9),2))
 assert reads==[9,81,115,63,37,27,45,9]
 assert R!=A and flip(R,(7,8))==A
 print('optimal unlabelled cycle',reads,'closing relabelling',p)
