from math import factorial, isqrt
from fractions import Fraction as F
def minsq(odd,even,budget):
    u,v=divmod(budget,odd+even)
    return odd+even*4*u*u+odd*(4*u*u+4*u)+(4+8*u)*min(v,even)+(8+8*u)*max(0,v-even)
def cert_cuts():
    n,k,r,q=3250,1625,23,78
    near=0
    for a in range(r,k+1):
        for b in range(max(a,k+q-a+1),k+1):
            d=r*(a+b-k-q)
            z=115*a*b-2*n*d
            if z>=0 and z*z>=225*a*b*(n-a)*(n-b): continue
            near+=1
            c=d-1
            for x,y in ((a,b),(b,a)):
                rem=n-x
                if 57*y<c or (115*y-2*c)*x-2*c*rem<=2*x*rem: continue
                if minsq(0,x,c)+minsq(y,rem-y,57*y-c)>225*y+4*y*y: break
            else: raise AssertionError((a,b,d))
    return near
assert 35*1547*2>65*1625
assert cert_cuts()==99

W,D0,mu,sig,T = 156,F(606,10),F(115,4),F(365625,8),35
levels = [(F(1,d-1)-F(1,d),(d-mu)**2,F(d*(57-d),57)) for d in range(T+1,58)]
points={F(0)}
for _,eta,_ in levels:
    points.update(eta*h for h in range(1,int(sig//eta)+1))
assert F(214)/(2*T*9)+(F(1,T)-F(1,57))*W*mu/(2*(W-D0))<1
assert sig<214**2 and W-D0>81
def sqrt_up(x, scale=10**5):
    x=F(x)*scale**2
    return F(isqrt(x.numerator//x.denominator)+1,scale)
minimum=None
for V in sorted(points):
    S=W*mu-sqrt_up((sig-V)*(W-D0))
    rhs=S/57
    for coef,eta,loss in levels:
        h=V//eta
        incoming=min(57*h,W+F(h*(h-1),2),sqrt_up((W-D0)*h*(56+h)))
        rhs+=coef*max(0,S-loss*h-incoming)
    if minimum is None or minimum>rhs: minimum=rhs
assert minimum>D0
q,m,b,E=78,1547,1309,52929
alpha=W/(W-D0)
x=D0/(2*(q-1))
beta=1-x+x*x/2-x*x*x/6
num=(22**(22*m)*77**(77*2*q)*factorial(m)*b**E*alpha.numerator**W*beta.numerator)
den=(23**(21*m)*78**(76*2*q)*factorial(q)**2*(b+1)**E*alpha.denominator**W*beta.denominator)
assert den<<19025 < num < den<<19026
assert len(points)==5127
print('99 near-cuts, 5127 variance breakpoints, uniform exponent 19025')
