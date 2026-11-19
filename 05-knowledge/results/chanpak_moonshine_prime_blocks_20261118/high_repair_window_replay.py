from itertools import combinations
from collections import defaultdict

def h(out):
 n=len(out);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S,row in enumerate(D):
  for i,val in enumerate(row):
   todo=out[i]&~S
   if not val:continue
   while todo:
    b=todo&-todo;D[S|b][b.bit_length()-1]+=val;todo^=b
 return sum(D[-1])
def out(n,z):
 A=[0]*n
 for j,(a,b) in enumerate(combinations(range(n),2)):
  if z>>j&1:A[a]|=1<<b
  else:A[b]|=1<<a
 return A

W=[([0,57,35,85,5,89,23],2,[33,48],[129,131,133,135]),([0,153,51,85,225,203,135,13],16,[15,46],[593,595,597,599]),([0,113,67,375,165,197,273,335,55],8,[201,216],[1483,1485,1487,1489]),([0,29,785,565,545,327,543,895,603,35],8,[137,152],[1721,1723,1725,1727])]
for base,d,states,expected in W:
 n=len(base)
 assert all(((base[i]>>j&1)+(base[j]>>i&1))==1 for i,j in combinations(range(n),2))
 def ev(s):
  A=base.copy()
  for i in range(n-1):
   if s>>i&1:A[0]^=1<<(i+1);A[i+1]^=1
  return h(A)
 pairs=[(ev(s),ev(s^d)) for s in states]
 assert all(b-a==2 for a,b in pairs)
 assert (states[0]^states[1]).bit_count()==2
 assert not((states[0]^states[1])&d)
 assert sorted(x for pair in pairs for x in pair)==expected
 print(n,pairs,'repairbit',d,'commonmove',states[0]^states[1])
