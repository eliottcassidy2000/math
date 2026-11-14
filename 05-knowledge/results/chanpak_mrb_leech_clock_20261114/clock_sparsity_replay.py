from fractions import Fraction as F
from math import prod
from clock_replay import K, QS

def power_bounds(x,s,den=10**8):
    x=F(x); s=F(s)
    a,b=x.numerator,x.denominator
    n,d=s.numerator,s.denominator
    bound=den**d*b**n
    xscale=a**n
    lo,hi=0,den
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**d*xscale>=bound: hi=mid
        else: lo=mid
    return F(lo,den),F(hi,den)

def upper(q):
    e,f,g=K[q]
    s=F(203,500)
    C=F(f+2*g,3)
    def ub(x):return power_bounds(x,s)[1]
    return sum(ub(F(e+2*f,3)+C*2**r) for r in range(1,33)) + ub(C*2**33)/(1-ub(2))

def lower(q):
    _,f,g=K[q]
    s=F(401,1000)
    return sum(power_bounds(f+g*2**r,s)[0] for r in range(1,33))

u,l=prod(map(upper,QS)),prod(map(lower,QS))
assert u<F(99,100) and l>F(101,100)
assert sum(QS)==21
print(f'ordered clock weights: at .406 < {float(u):.9f}; at .401 > {float(l):.9f}')

def inject_weight(q):
    _,f,g=K[q]
    s=F(49,125)
    start=2 if q==3 else 1
    return sum(power_bounds(f+g*2**r,s)[0] for r in range(start,start+31))
z=prod(map(inject_weight,QS))
assert z>F(102,100)
print(f'vector-decodable cycle alphabet weight > {float(z):.9f}')
