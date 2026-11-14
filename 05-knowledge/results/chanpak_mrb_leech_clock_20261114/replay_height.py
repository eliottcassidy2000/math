from fractions import Fraction as F
from functools import cache

num,den=203,500
@cache
def up(x):
    n,d=x.numerator,x.denominator
    scale=10**7
    lo,hi=0,scale
    goal=scale**den*d**num
    while hi-lo>1:
        m=(lo+hi)//2
        if m**den*n**num>=goal: hi=m
        else: lo=m
    assert hi**den*n**num>=goal
    return F(hi,scale)

z=up(F(2))
w=[up(F(3)),sum(up(F(1+2**r)) for r in range(2,33))
                  +up(F(2**33))/(1-z)]
P=(F(2,3),F(4,5))
def M(e,f,g):
    return [[w[j]*up(e*(1-p)*(1-q)+f*(p*(1-q)+(1-p)*q)+g*p*q)
             for j,q in enumerate(P)] for p in P]
v=[F(1000),F(991)]
u=v.copy()
for mat in reversed([M(9,11,13),M(3,3,4),M(17,20,24),M(5,6,7),M(57,68,81)]):
    u=[sum(a*b for a,b in zip(row,u)) for row in mat]
assert all(a<=F(99,100)*b for a,b in zip(u,v))
print('two width types, rational upper cycle ratios:',
      tuple(round(float(a/b),9) for a,b in zip(u,v)))
