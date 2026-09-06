"""Exact failure of coefficient positivity for the minimal h4 amplitude interpolant.

This tests the degree-three interpolant on the positive square-root branch,
not every possible higher-degree or multi-window amplitude certificate.
"""
from fractions import Fraction as F
from math import comb,factorial,isqrt
import sys
sys.stdout.reconfigure(newline='\n')
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise ArithmeticError(label)


def add(*rows):return sum(a for a,b in rows),sum(b for a,b in rows)
def neg(I):return -I[1],-I[0]


def mul(A,B):
    vals=[a*b for a in A for b in B]
    return min(vals),max(vals)


def div(A,B):
    need(not B[0]<=0<=B[1],'interval denominator excludes zero')
    return mul(A,(1/B[1],1/B[0]))


def horner(row,I):
    v=(F(0),F(0))
    for c in reversed(row):v=add(mul(v,I),(c,c))
    return v


def sqrtI(I):
    M=10**20
    a=isqrt(I[0].numerator*M*M//I[0].denominator)
    b=isqrt(I[1].numerator*M*M//I[1].denominator)+1
    lo,hi=F(a,M),F(b,M)
    need(0<lo and lo*lo<=I[0]<=I[1]<=hi*hi,'rational square-root enclosure')
    return lo,hi


def main():
    h,g=4,14
    P=[F(1,11),-10,84,-60,1]
    beta={j:comb(13-2*j,1+j) for j in range(-1,5)}
    W=[];T=[]
    for j in range(-1,9):
        W.append(F((-1)**(j%2)*comb(28,2+2*j)*sum(beta[a]*beta.get(j-a,0) for a in beta)))
        T.append(F((-1)**(j%2)*factorial(28),factorial(2+j)*factorial(24-3*j)*factorial(2+2*j)))
    first=[F((-1)**j*factorial(14),factorial(1+j)*factorial(12-3*j)*factorial(1+2*j)) for j in range(5)]
    need([v/first[-1] for v in first]==P,'actual monic first polynomial')
    for mass,z in ((14,1),(28,2)):
        actual={}
        for a in range(mass+1):
            for b in range(mass-a+1):
                c=mass-a-b
                if -27*a+b+15*c:continue
                need((c-z)%2==0,'literal channel congruence')
                actual[(c-z)//2]=factorial(mass)//(factorial(a)*factorial(b)*factorial(c))
        expected={j:int((-1)**j*first[j]) for j in range(5)} if mass==14 else {j:int((-1)**(j%2)*T[j+1]) for j in range(-1,9)}
        need(actual==expected,'literal first/doubled charge fibers with full carry')
    Is=[(F(8419,849544),F(11993,1210189)),(F(259526,2155711),F(291249,2419213)),
        (F(11376199,8744207),F(1124341,864214)),(F(120947883,2065060),F(91906768,1569213))]
    for I in Is:
        need(horner(P,(I[0],I[0]))[0]*horner(P,(I[1],I[1]))[0]<0,'opposite-sign exact root bracket')
        need(horner(T,I)[1]<0 and horner(W,I)[1]<0,'actual and hit responses remain negative')
    need(all(Is[i][1]<Is[i+1][0] for i in range(3)),'four disjoint roots exhaust the quartic')
    Z=[sqrtI(I) for I in Is]
    total=(F(0),F(0))
    for i in range(4):
        denominator=(F(1),F(1))
        for j in range(4):
            if j!=i:denominator=mul(denominator,add(Z[i],neg(Z[j])))
        # The z^2 coefficient of the ith monic cubic Lagrange numerator.
        weight=div(neg(add(*(Z[j] for j in range(4) if j!=i))),denominator)
        total=add(total,mul(div(horner(T,Is[i]),horner(W,Is[i])),weight))
    need(total[1]<0,'strict negative quadratic coefficient of the unique cubic amplitude interpolant')
    scaled=[v*10**6 for v in total]
    lower=scaled[0].numerator//scaled[0].denominator
    upper=-((-scaled[1].numerator)//scaled[1].denominator)
    need((lower,upper)==(-344095,-344094),'frozen rational millionth enclosure')
    print('ACTUAL support=(-27,1,15) firstmass14 doubledmass28 P=s^4-60s^3+84s^2-10s+1/11')
    print('MINIMAL_AMPLITUDE_INTERPOLANT quadratic_coefficient_enclosure=',lower,'/1000000,',upper,'/1000000')
    for i,I in enumerate(Is):print('ROOT',i+1,','.join(map(str,I)))
    print('SCOPE coefficient positivity of the unique degree<4 positive-branch interpolant fails; higher-degree or multi-window certificates remain open')
    print('PASS',GATES,'always-active gates; exact path convolution, literal charge fibers, rational interval interpolation')


if __name__=='__main__':main()
