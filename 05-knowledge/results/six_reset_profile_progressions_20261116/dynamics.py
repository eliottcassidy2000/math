from fractions import Fraction as F
from math import gcd

def order(a,n):
    assert gcd(a,n)==1
    x=a%n;t=1
    while x!=1:
        x=x*a%n;t+=1
    return t
assert (order(2,3),order(2,7),order(2,9),order(2,63))==(2,3,6,6)
assert 2**6-1==63==7*9 and 756==12*63
f=lambda x:x*x-F(29,16)
T=[F(-7,4),F(5,4),F(-1,4)]
assert [f(x) for x in T]==T[1:]+T[:1]
assert f(F(1,4))==T[0] and f(F(-5,4))==T[2]
m=8*T[0]*T[1]*T[2]
assert m==F(35,8) and m+1==F(43,8)
Q=[]
a,b=F(-1),F(1,12)
for j in range(6):
    Q.append((a,b))
    a,b=a*a+33*b*b-F(71,48),2*a*b
assert len(set(Q))==6 and Q[0]==(a,b)
assert Q[3:]==[(a,-b) for a,b in Q[:3]]
assert len({Q[i+1][0]-Q[i][0]**2 for i in range(3)})==2
print('Mersenne/residue clocks and exact quadratic three- and six-step checks pass')
