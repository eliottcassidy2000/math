from functools import cache
from math import gcd
from fractions import Fraction

@cache
def P(n,a,b,x):
    if not n:
        return Fraction(1)
    return a*sum(x[i]*b**i*P(n-d,a,b,tuple(b**gcd(d,e)*x[j]
                 for j,e in enumerate(range(1,n-d+1,2))))
                 for i,d in enumerate(range(1,n+1,2)))/n

def at1(n,a,b):
    return P(n,a,b,(1,)*((n+1)//2))

A=[at1(n,1,2) for n in range(9)]
B=[at1(n//2,2**(n%2),4) for n in range(9)]
if A!=[1,1,1,2,4,12,56,456,6880] or B!=[1,1,1,2,2,8,12,88,176]:
    raise ValueError(A,B)
print('ordinary:',A,'self-converse:',B)
