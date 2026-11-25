from math import factorial, isqrt
from fractions import Fraction as F
W,D0,mu,sig,T = 140,F(4149,100),F(115,4),F(365625,8),34
levels = [(F(1,d-1)-F(1,d),(d-mu)**2,57+F(d*(57-d),57)) for d in range(T+1,58)]
points={F(0)}
for _,eta,_ in levels:
    points.update(eta*h for h in range(1,int(sig//eta)+1))
assert sig < 4*T**2*(W-D0)
minimum=None
for V in sorted(points):
    x=(sig-V)*(W-D0)*10**10
    S=W*mu-F(isqrt(x.numerator//x.denominator)+1,10**5)
    rhs=S/57+sum(coef*max(0,S-loss*(V//eta)) for coef,eta,loss in levels)
    if minimum is None or minimum>rhs: minimum=rhs
assert minimum>D0
q,m,b,E=70,1555,1317,52929
alpha=W/(W-D0)
x=D0/(2*(q-1))
beta=1-x+x*x/2-x*x*x/6
num=(21**(21*m)*69**(69*2*q)*factorial(m)*b**E*alpha.numerator**W*beta.numerator)**22
den=(22**(20*m)*70**(68*2*q)*factorial(q)**2*(b+1)**E*alpha.denominator**W*beta.denominator)**22*(2*m)
assert (den<<(18997*22)) < num < (den<<(18998*22))
assert len(points)==6269
print(len(points),'variance breakpoints; uniform strict exponent 18997')
