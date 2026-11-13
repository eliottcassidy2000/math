from itertools import permutations
from math import isqrt

def req(ok):
    if not ok:raise RuntimeError('replay mismatch')

def ext(A,U):
    n=len(A);return tuple(a+((not (U>>i&1))<<n) for i,a in enumerate(A))+(U,)

def valid(A,p): return all(A[i]>>j&1 for i,j in zip(p,p[1:]))

def coeff(A):
    n=len(A); ans=[]
    for j in range(1,min(sum(A[i]&1 for i in range(n)),A[0].bit_count())+2):
        ans.append(sum(valid(A,p) for p in set(permutations(tuple(range(1,n))+(0,)*j))))
    return ans

A=(24,9,3,4,14);B=(44,33,2,6,15,28)
req(coeff(A)==[13,12,2] and coeff(B)==[25,24,6])
for Q,u,want in ((A,27,[29,24,4]),(B,4,[43,40,10]),(B,10,[75,72,18]),(B,12,[53,48,12])):
    req(coeff(ext(Q,u))==want)
Z=ext(B,10)
for i in (2,3,6):
    req(tuple(j for j in (0,1,4,5) if Z[i]>>j&1)==(1,))
req(valid(Z,(3,2,6,3)))
for j,want in ((1,[13,10,2]),(2,[6,7,2])):
    vals=[sum(valid(B,p) for p in set(permutations((0,)*i+(2,)*j+(1,4,5)))) for i in (1,2,3)]
    req(vals==want)
print('four refined quotient equalizers, exact triple transfer, and two-column clock verified')
