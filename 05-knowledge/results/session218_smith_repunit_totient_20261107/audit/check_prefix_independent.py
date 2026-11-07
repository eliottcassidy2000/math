"""Independent complete totient prefix via direct factorization, exact strong counts."""
from math import comb,isqrt

def phi(n):
    v=n; t=n; p=2
    while p*p<=t:
        if t%p==0:
            v=v//p*(p-1)
            while t%p==0: t//=p
        p+=1
    if t>1: v=v//t*(t-1)
    return v
bound=2*76**2
image={phi(n) for n in range(1,bound+1)}
missing=[m for m in range(2,77,2) if m not in image]
assert missing==[14,26,34,38,50,62,68,74,76]
# Exact integer recurrence, with no modular reductions inside the recurrence.
strong=[0]*124
for n in range(1,124):
    strong[n]=2**comb(n,2)-sum(comb(n,j)*strong[j]*2**comb(n-j,2) for j in range(1,n))
assert strong[1:5]==[1,0,2,24]
p=251
assert all(p%d for d in range(2,isqrt(p)+1))
N=2**50-1
assert N%p==0 and N>p
assert ((N-1)%(p-1))+1==123
D=strong[123]-2**61
assert D%p==229
print('PASS complete even-nontotient prefix through76 via direct-factorization phi, n<=11552:',missing)
print('PASS exact-integer strong recurrence: D(123) mod251 =229; 251 divides2^50-1; reduced index123.')
print('Given the cited quasi-periodicity theorem, nonzero229 excludes exponent50 from the Camion exponent set.')
