def require(x):
 if not x: raise AssertionError("transfer check failed")
from collections import Counter
from functools import lru_cache
C=((0,0,0,0,4),(3,0,0,0,0),(2,5,0,0,0),(2,5,5,0,0),(0,2,2,3,0))
@lru_cache(None)
def transfer(n, fields):
 if n==0: return {(0,0):1}
 plus=transfer(n-1, tuple(fields[i]+C[i][n-1] for i in range(n-1)))
 minus=transfer(n-1, tuple(fields[i]-C[n-1][i] for i in range(n-1)))
 out=Counter(plus)
 for (k,e),v in minus.items():out[k+1,e+fields[-1]+sum(C[n-1][:n-1])]+=v
 return dict(out)
Z=transfer(5,(0,)*5)
direct=Counter()
for s in range(32):
 k=s.bit_count(); e=sum(C[i][j] for i in range(5) for j in range(5) if s>>i&1 and not s>>j&1)
 direct[k,e]+=1
require(Z==direct)
require(sum(Z.values())==32)
P=Counter()
for (k,e),v in Z.items():
 if k not in (0,5):P[9+2*e]+=v
require(sorted(P)==[15,17,19,23,25,27,29,33,37,41])
require([P[x] for x in sorted(P)]==[2,4,2,6,2,2,2,4,4,2])
B={2,3,5,7,11,13,17,19,23,27}
print('layers (u exponent, q exponent, multiplicity):',sorted((k,e,v) for (k,e),v in Z.items()))
print('proper response coefficients:',sorted(P.items()))
print('common exponents and their ear coefficients:',sorted((b,P[b]) for b in B if b in P))
print('PASS transfer equals all 32 direct cut weights')
