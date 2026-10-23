"""Small independent synthesis normalization tests, not packet proof replays."""
from fractions import Fraction
from math import factorial, isqrt
import json

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def v2(n):
    return (n & -n).bit_length()-1

checks = 0
for n in range(3, 2002, 2):
    a=(n-1).bit_length()-1
    B=min(n-v2(n+1),n+1-a)
    repunit=(n & (n+1)) == 0
    piece=n-v2(n+1) if repunit else n+1-a
    require(B==piece, f'binary piece at {n}')
    transitive=(n-v2(n+1)==B)
    m=n+1
    family=repunit or (m % 3 == 0 and m//3 >= 2 and ((m//3)&(m//3-1))==0)
    require(transitive==family, f'transitive family at {n}')
    checks+=2

require(Fraction(226,315)<Fraction(8,9)<Fraction(1657897,1801800), 'gap order')
require(Fraction(8,9)*Fraction(1,2)==Fraction(4,9),'soft gap')
checks+=2

# Construct Eulerian polynomials by their elementary coefficient recurrence.
A=[1]
polys={1:A}
for n in range(2,6):
    A=[(k+1)*(A[k] if k<len(A) else 0)+(n-k)*(A[k-1] if k else 0) for k in range(n)]
    polys[n]=A
C=[0]*5
for i,x in enumerate(polys[3]):
    for j,y in enumerate([1,-2,1]): C[i+j]+=2*x*y
hostile=[x+3*y for x,y in zip(polys[5],C)]
require(hostile==[7,38,30,38,7], 'hostile coefficients')
require(sum(hostile)==factorial(5), 'hostile total')
checks+=2

primes=[2,3,5,7,11,13]
exponents=[5,7,3,4,1,2]
import itertools
hist={}
basins={True:0,False:0}
N=1
for p,e in zip(primes,exponents): N*=p**e

def height(b):
    q=N;k=0
    while q % b == 0:
        q//=b;k+=1
    return k
for es in itertools.product(*(range(e+1) for e in exponents)):
    if not any(es): continue
    h=min(E//e for E,e in zip(exponents,es) if e)
    hist[h]=hist.get(h,0)+1
    b=h; seen=set()
    while b>1 and b not in seen:
        seen.add(b);b=height(b)
    basins[b>1]+=1
require(hist=={1:5616,2:120,3:16,4:4,5:2,7:1}, 'height polynomial')
require(basins=={True:143,False:5616}, 'height basins')
for b,h in [(4,2),(2,5),(5,3),(3,7),(7,4),(13,2)]:
    require(height(b)==h, 'height cycle')
checks+=8

# Enumerate all plane monomial charges through rigidity threshold.
for degree in range(71):
    for a in range(degree+1):
        b=degree-a;q=a-b
        for target in [-1,1]:
            if (q-target)%8==0 and (q-target)%9==0:
                require(q==target,'degree70 alias')
                checks+=1
require((-71-1)%8==0 and (-71-1)%9==0 and -71!=1,'sharp shear')
checks+=1
print(json.dumps({'status':'PASS','independent_checks':checks,'binary_odd_orders':[3,2001], 'scope':'Small arithmetic and normalization tests only; packet proofs and audits reviewed separately.'},indent=2))
