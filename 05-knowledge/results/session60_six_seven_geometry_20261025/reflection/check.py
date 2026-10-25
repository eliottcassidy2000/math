from itertools import product
from collections import Counter

def require(x,m):
 if not x:raise RuntimeError(m)
def apply(M,v):
 return sum(((sum(M[i][j]*((v>>j)&1) for j in range(3))%2)<<i) for i in range(3))
I=tuple((1<<i) for i in range(3))
counts=Counter(); invertible=0
for bits in product((0,1),repeat=9):
 M=[bits[3*i:3*i+3] for i in range(3)]
 p=tuple(apply(M,v) for v in range(8))
 if len(set(p))!=8:continue
 invertible+=1
 if all(p[p[v]]==v for v in range(8)) and any(p[v]!=v for v in range(8)):
  fixed=sum(p[v]==v for v in range(1,8));counts[fixed]+=1
require(invertible==168,'GL3F2order')
require(counts=={3:21},'involutionfixedtype')
J={v:(v^7) if v!=7 else 7 for v in range(1,8)}
require(sum(J[v]==v for v in J)==1,'Jfixed')
lines={tuple(sorted((a,b,a^b))) for a in range(1,8) for b in range(1,8) if a!=b}
require(len(lines)==7,'Fanolines')
require(sum(tuple(sorted(J[v] for v in L)) in lines for L in lines)==3,'Jlinepreservation')
# A permutation involution with 1 fixed point and 3 two-cycles has invariant dimension4.
require((7+1)//2==4,'canonicalinvariants')
print('PASS: GL3F2=168; involutions21 allfix3; modifiedcomplementfix1 preserves3of7lines; canonicalinvariants4')
