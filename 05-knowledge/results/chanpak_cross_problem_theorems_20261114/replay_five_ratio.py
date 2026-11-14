from fractions import Fraction as F
from math import gcd

def ratio(j):
    return (4*F(2)**j+2*(F(3)**j-1))/(1+F(2)**j)

assert ratio(-3)*ratio(5)*ratio(4)*(ratio(2)/ratio(-2))**2 == -2**15

for mod in range(3,500,2):
    if gcd(mod,3*5*7*11*17)>1:
        continue
    period=1
    while pow(2,period,mod)!=1 or pow(3,period,mod)!=1:
        period+=1
    def unit_ratio(j):
        j=j if j>0 else j%period or period
        h=(1+pow(2,j,mod))%mod
        q=(4*pow(2,j,mod)+2*(pow(3,j,mod)-1))%mod
        assert gcd(h*q,mod)==1
        return q*pow(h,-1,mod)%mod
    assert (pow(2,-15,mod)*unit_ratio(-3)*unit_ratio(5)*unit_ratio(4)
            *pow(unit_ratio(2),2,mod)*pow(unit_ratio(-2),-2,mod))%mod==mod-1

for m,a,s,mod in [(4,1,0,15),(16,3,0,51),(3,2,1,85)]:
    h=(1+2**m)**a
    q=(2**(m+2)+2*(3**m-1))**a*2**s
    assert (h+2*q)%mod==0 and gcd(h*q,mod)==1

assert pow(2,16,255)==1 and pow(3,17,255)==3
for m in range(1,17):
    h=1+2**m
    q=2**(m+2)+2*(3**m-1)
    assert (gcd(h*q,255)==1)==(m==16)
    if m==16:
        assert q*pow(h,-1,85)%85==2
assert all(pow(2,k,85)!=84 for k in range(8))
print("five-ratio, positive-lift, and three-prime compatibility controls passed")
