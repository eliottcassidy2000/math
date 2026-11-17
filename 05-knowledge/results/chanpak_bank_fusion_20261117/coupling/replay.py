from fractions import Fraction as Q
from random import Random
from itertools import product

def ck(v):
    if not v:raise ArithmeticError('control failed')
def TT(n):
    return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def join(A,B):
    x=((1<<len(B))-1)<<len(A)
    return tuple(y|x for y in A)+tuple(y<<len(A) for y in B)
def bits(u):
    return bin(u)[2:] if u else ''
def make(A,r,u):
    w=bits(u);k=len(w);B=TT(k);S=join(B,A);i=len(S)
    tag=sum(int(b)<<j for j,b in enumerate(w))|(((1<<len(A))-1)^(1<<r))<<k
    return tuple(z|((1-(tag>>j&1))<<i) for j,z in enumerate(S))+(tag,),i
def contrast(a):
    return tuple(sum(((a[j]>>i)&1)<<j for j in range(len(a))) for i in range(len(a)))
def paths(a):
    n=len(a);U=1<<n;d=[[0]*n for _ in range(U)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,U):
        for i,z in enumerate(d[s]):
            if not z:continue
            b=a[i]&~s
            while b:
                x=b&-b;b-=x;d[s|x][x.bit_length()-1]+=z
    return [1]+[sum(b) for b in d[1:]],d
def profile(a,r):
    s,d=paths(a);_,e=paths(contrast(a));U=len(s)-1
    return (s[U],sum(d[j][r]*s[U^j] for j in range(U+1)),e[U][r],d[U][r])
def F(u):
    v=1
    while u:
        x,b=divmod(u,2)
        if not b:v+=x+x%2
        u=x
    return v
def compile(sig,u):
    h,R,S,e=sig;k=u.bit_length();w=bits(u);rc=0
    for b in reversed(w):rc=(rc<<1)+(1-int(b))
    return (F(u)*h+(u+1)*R-(1-u%2)*S, (1<<k)*R+h*(1+rc), h if u else h-S, e)
C=(4,1,2);g0=profile(C,0);z=0
for seq in product(range(5),repeat=3):
    a,r,sig=C,0,g0
    for u in seq:
        a,r=make(a,r,u);sig=compile(sig,u)
    ck(profile(a,r)==sig and sig[-1]==1)
    z+=1
print('PASS',z,'independent literal subset-DP migration chains')
s=Q(1,12)+Q(2,32)
for k in range(2,9):
    b=2**(k+3)+2*(3**k+1);c=4*(3**(k-1)+1)
    s+=2**(k-2)*(1/Q(b)+2/Q(b+c)+1/Q(b+2*c))
    ck(Q(28,25)**100>b+2*c)
ck(Q(51,100)**100 > Q(1,2)**99)
ck((Q(17,25)/2)**100 > Q(1,3)**99)
upper=s*Q(28,25)+Q(51,100)+Q(51,100)*Q(17,25)**9/(1-Q(17,25))
ck(upper<Q(99,100))
print('PASS exact rational weight bound for exponent 99/100:',float(upper))
strong_upper=upper-Q(51,100)
ck(strong_upper<Q(477084,1000000))
ck((strong_upper+Q(51,100))*Q(51,49)<Q(103,100))
ck(Q(1007,1000)**1000>1000)
ck(Q(1007,1000)*(Q(987084,1000000)+Q(103,100000))<1)
print('PASS exact weak-stage frequency certificate')
