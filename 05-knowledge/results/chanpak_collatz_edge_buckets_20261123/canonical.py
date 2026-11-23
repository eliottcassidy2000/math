from itertools import combinations,permutations,product
from pathlib import Path
import sys
def canonical(A,e=None):
 n=len(A)
 def C(P):
  return sum(((A[P[i]]>>P[j])&1)<<k for k,(i,j) in enumerate(combinations(range(n),2)))
 def by_scores(R):
  d={}
  for v in R:d.setdefault(A[v].bit_count(),[]).append(v)
  return (sum(W,()) for W in product(*(permutations(d[k]) for k in sorted(d))))
 if e is None:return min(C(P) for P in by_scores(range(n)))
 return min(C(E+P) for E in (e,e[::-1]) for P in by_scores(set(range(n))-set(e)))
def flip(A,e):
 a,b=e;B=list(A);B[a]^=1<<b;B[b]^=1<<a;return tuple(B)
def strong(A):
 n=len(A)
 for v in range(n):
  X=1<<v
  while True:
   Y=X
   for j in range(n):
    if X>>j&1:Y|=A[j]
   if Y==X:break
   X=Y
  if X!=((1<<n)-1):return False
 return True
B63=set(); C1=set();C2=set();S1=set();S2=set()
for line in Path(sys.argv[1]).read_text().splitlines():
 s=line.split()
 if s[0]=='HIT':
  A=tuple(map(int,s[6:14]))
  assert strong(A)
  B63.add(canonical(A))
  for j in range(14,len(s),3):
   assert s[j]=='L'
   e=tuple(map(int,s[j+1:j+3]))
   C1.add(canonical(A,e))
   S1.add(strong(flip(A,e)))
 if s[0]=='DIRECT61':
  A=tuple(map(int,s[3:11]))
  for j in range(11,len(s),2):
   e=tuple(map(int,s[j:j+2]))
   C2.add(canonical(A,e))
   S2.add(strong(flip(A,e)))
assert (B63,C1,C2,S1,S2)==({3145732,2359556},{969,1004},{449,2044},{True},{True})
print("canonical order8 H63 cores",sorted(B63))
print("canonical marked 63 toward 61",sorted(C1),"all strong")
print("canonical marked 61 toward 65",sorted(C2),"all strong")
def H(A):
 d=[[0]*len(A) for _ in range(1<<len(A))]
 for i in range(len(A)):d[1<<i][i]=1
 for mask in range(1,len(d)):
  for i,h in enumerate(d[mask]):
   if h:
    for j in range(len(A)):
     if A[i]>>j&1 and not mask>>j&1:d[mask|1<<j][j]+=h
 return sum(d[-1])
A=(158,28,24,16,32,79,159,62)
E=(1,5);F=(5,7)
cube=[A,flip(A,E),flip(A,F),flip(flip(A,E),F)]
assert list(map(H,cube))==[61,63,65,69] and all(map(strong,cube))
print("literal strong four corners: 61,63,65,69")

