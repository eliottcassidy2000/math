from itertools import product
from math import gcd
from fractions import Fraction
S=(0,8,1,0)
def mul(x,y):
    a,b,c,d=x;e,f,g,h=y
    return ((a*e+b*g)%9,(a*f+b*h)%9,(c*e+d*g)%9,(c*f+d*h)%9)
def neg(m):return tuple(-x%9 for x in m)
def can(m):return min(m,neg(m))
I=can((1,0,0,1))
G={can((a,b,c,d)) for a,b,c,d in product(range(9),repeat=4) if (a*d-b*c)%9==1}
assert len(G)==324
assert len([a for a in G if can(mul(a,S))==can(mul(S,a))])==12
assert len([a for a in G if can(mul(a,a))==I and a!=I])==27

a,b=1465,2322
D=2**b-3**a
x=Fraction(3**a-2**a,D)
assert gcd(D,63*72)==1 and D%63==19 and D%72==37
assert x//1==47 and x.denominator>1
y=x
for i in range(b):
    assert y.denominator%2 and (y.numerator%2)==(i<a)
    assert y>=x and y.denominator>1
    y=(3*y+1)/2 if i<a else y/2
assert y==x
print('G9 projective centralizer 12; 27 involutions; exact primitive wheel-passing rational cycle checked')
