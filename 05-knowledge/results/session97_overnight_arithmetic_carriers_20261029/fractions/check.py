"""Exact fruit cubic, symmetry-cover and rational-reconstruction certificates."""
from fractions import Fraction as F
from math import gcd, lcm
from itertools import permutations

def require(ok, msg):
    if not ok: raise ArithmeticError(msg)

def add(P,Q,A=109,B=224):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and y==-v:return None
    m=(3*x*x+2*A*x+B)/(2*y) if P==Q else (v-y)/(u-x)
    z=m*m-A-x-u
    return z,m*(x-z)-y

def mul(n,P):
    if n<0:return mul(-n,(P[0],-P[1]))
    Q=None
    while n:
        if n&1:Q=add(Q,P)
        P=add(P,P);n//=2
    return Q

def fruit(P):
    if P is None:return (F(1),F(-1),F(0))
    x,y=P
    return 56-x+y,56-x-y,-56-12*x

def primitive(v):
    m=lcm(*(x.denominator for x in v));a=[int(x*m) for x in v]
    g=gcd(*a);a=[x//g for x in a]
    if sum(a)<0:a=[-x for x in a]
    return tuple(a)

def audit(P):
    x,y=P;require(y*y==x**3+109*x*x+224*x,'curve')
    a,b,c=fruit(P);s=a+b+c;q=a*b+a*c+b*c;r=a*b*c
    require(s**3-6*s*q+7*r==0,'cubic')
    if not (a+b)*(a+c)*(b+c):return
    require(a/(b+c)+b/(a+c)+c/(a+b)==4,'original denominators')
    d=6*(a+b)-c
    require(d!=0 and (-28*(a+b+2*c)/d,364*(a-b)/d)==P,'inverse')
    q/=s*s;w=7*(a-b)*(b-c)*(c-a)/s**3
    require(w*w==(q+1)*(1+29*q-196*q*q),'quotient')
    Z=-196*(q+1);Y=-196*w
    require(Y*Y==Z**3+421*Z*Z+43904*Z,'isogenous curve')
    require(Z==x*(x-56)**2/(x-4)**2 and Y==(x-56)*(x*x+44*x+224)*y/(x-4)**3,'isogeny map')
    roots=[a/s,b/s,c/s]
    require(all(t**3-t*t+q*t-(6*q-1)/7==0 for t in roots),'root guard')
    pos=all(t>0 for t in roots)
    require(pos==(q>F(1,6) and 1+29*q-196*q*q>0),'positive chamber')
    for z in permutations(roots):
        require(sum(z)==1 and z[0]*z[1]+z[0]*z[2]+z[1]*z[2]==q,'S3 quotient')

def main():
    G=(F(-4),F(28));T=(F(56),F(728))
    require(mul(6,T) is None and mul(3,T)==(0,0),'torsion')
    first=None
    for n in range(1,21):
        for k in range(6):
            P=add(mul(n,G),mul(k,T));audit(P)
        p=primitive(fruit(mul(n,G)))
        if min(p)>0 and first is None:first=n
    require(first==9,'first tested positive multiple')
    p=primitive(fruit(mul(9,G)))
    require(tuple(sorted(map(lambda x:len(str(x)),p)))==(79,80,81),'digits')
    # Rational-root theorem: primitive cubic 7t^3-7t^2+1 has candidates +-1,+-1/7.
    require(all(7*t**3-7*t*t+1!=0 for t in [F(1),F(-1),F(1,7),F(-1,7)]),'hostile irreducibility')
    require((0+1)*(1+29*0-196*0)==1,'hostile discriminant square')
    # No repeated rational coordinates: t=1 or 14t²-11t+1=0; t=1 is boundary.
    require(8**2<65<9**2,'irrational branch values')
    H=(F(-196),F(196));H3=add(add(H,H,421,43904),H,421,43904)
    require(H3==(F(-225),F(210)),'3H=phi(G) return identity')
    R=None
    for _ in range(35):R=add(R,H,421,43904)
    Z,Y=R;q=-Z/196-1;w=-Y/196
    require(Y*Y==Z**3+421*Z*Z+43904*Z,'35H membership')
    require(q>F(1,6) and 1+29*q-196*q*q>0,'positive-real hostile')
    require(w*w==(q+1)*(1+29*q-196*q*q),'positive hostile square discriminant')
    require(q.numerator%2==0 and q.denominator%2==1,'local obstruction parity')
    require(all((t**3+t*t+1)%2 for t in range(2)),'no local root')
    print('35H: positive real chamber and rational orientation, but no Q_2 lift')
    print('PASS: 120 curve/map/quotient/chamber checks; 9G digits 79,80,81; nonsplit rational cover hostile')
    print('9G primitive triple:',p)
if __name__=='__main__':main()
