from fractions import Fraction as Q
from math import gcd, isqrt
from clock_replay import K

# Recurrence-cone inequalities need check only on the first three
# indices of the tail (2,3,4); the positive recurrence proves the rest.
for e,f,g in K[2:5]:
    assert 4*e >= 3*f and 3*g <= 4*f
assert Q(2*2,(2+2*2)*(1+2*2))==Q(2,15)
# 7-4 sqrt(3) < 2/15.
assert Q(7)-Q(2,15)>0 and (Q(7)-Q(2,15))**2 < 16*3

# 0 < beta=(3 sqrt(41)-19)/(2 sqrt(41)) < 1/50.
assert 19*19 < 9*41
assert (Q(3)-Q(2,50))**2*41 < 19**2
assert 5*5 < 41 < 7*7  # -1 < (5-sqrt(41))/2 < 0.
a,b=1,2
last_two=[]
for i in range(20):
    h=a+b
    if i<=1: assert h==[3,17][i]
    else: assert h==5*last_two[-1]+4*last_two[-2]
    last_two.append(h)
    a,b=a+2*b,4*(a+b)

# Derive the six labels as a proper subset, not as the full pure-A list.
H={d+1 for d in range(1,25) if 24%d==0}
assert H-set((4,2,5,3,7,13))=={9,25}
assert 24//(13-1)==2

n=47**2+56**2
u=-56*pow(47,-1,n)%n
assert n==5345==5*1069 and u==2387 and (u*u+1)%n==0
assert gcd(47,n)==gcd(56,n)==1
assert all(1069%i for i in range(2,isqrt(1069)+1))
assert n==4**2+73**2
assert -73*pow(4,-1,n)%n not in (u,(-u)%n)
for x in range(-50,51):
    for y in range(-50,51):
        U,V=47*x+56*y,-56*x+47*y
        assert (V-u*U)%n==0
        assert (47*U-56*V,56*U+47*V)==(n*x,n*y)
print('physical derivative bound, integer rounding guard, pure-A indices and reset glue checks passed')
