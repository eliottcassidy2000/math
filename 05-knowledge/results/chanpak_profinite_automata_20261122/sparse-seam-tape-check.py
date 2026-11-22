from math import gcd
from random import Random
rng=Random(26)
def f(y):
    v=0
    while y:
        if y&1:
            u=y+1;y>>=(u&-u).bit_length()-1
        else:
            u=y>>((y&-y).bit_length()-1)
            v+=y-u+1;y=u
    return v
def order2(m):
    if m==1:return 1
    L=1;q=2%m
    while q!=1:L+=1;q=2*q%m
    return L
def rev(x,r):
    z=0
    for j in range(r):z=2*z+(x&1);x>>=1
    return z
def reset(word,r,m,ya,f2,fm):
    U=1<<r;M=m*U;d=order2(m);L=((r+d-1)//d)*d
    ya%=m;f2%=U;fm%=m
    y0=int(word,2);y=y0%M;g=f(y0)%M
    pieces=[]
    def app(w):
        nonlocal y,g
        for ch in w:
            if ch=='0':g=(g+y+(y&1))%M;y=2*y%M
            else:y=(2*y+1)%M
        pieces.append(w)
    app('0'*L)
    delta=(g-f2)%U
    for s in range(1,r):
        if delta>>s&1:app('1'+'0'*s+'1'+'0'*L)
    assert g%U==f2
    P='0'*(L+1);Q='1'+'0'*L;T=P*(L-1)+Q
    app(T*((ya-y)%m))
    t=(-1-ya)%m
    b=(fm-g)%m
    app(T*t)
    assert y%m==(-1)%m
    app(Q*b)
    app(T*((-t)%m))
    assert y%U==0 and y%m==ya and g%U==f2 and g%m==fm
    return word+''.join(pieces)
def general(word,r,m,a,t):
    U=1<<r;M=U*m;b=a%U
    ym=((a-b)*pow(U,-1,m))%m if m>1 else 0
    zmod=U*(ym*pow(U,-1,m)%m) if m>1 else 0
    k=(U-1-rev(b,r))*zmod+f(b)
    s=reset(word,r,m,ym,(t-k)%U,(t-k)%m)
    return s+format(b,f"0{r}b")
for r,m in ((1,1),(3,1),(4,3),(4,5),(3,7),(5,9),(4,25),(5,35),(5,71)):
    U=1<<r;M=U*m
    for i in range(3):
        a=rng.randrange(M);t=2*rng.randrange(M//2)
        w=general('1010',r,m,a,t)
        x=int(w,2)
        assert x%M==a and f(x)%M==t
        assert '11' not in w[:-r]
# Fixed strong chart tracked along a single nested tape; its literal counts are inherited.
w='10';n=10**18+3
for r,m in ((1,1),(3,5),(5,25),(6,75),(7,225)):
    U=1<<r;M=U*m
    target=(n-149)%M
    # Bezout choice on odd component, which may mix 5 and other primes.
    a=(-2*target)%m;b=(13*target)%m
    dz=((n-149)*pow(25,-1,U))%U
    nxt=reset(w,r,m,a,dz,b)
    assert nxt.startswith(w) and '11' not in nxt
    x=int(nxt,2)
    assert (149+162*x+25*f(x))%M==n%M
    w=nxt
for n in range(4,30):
    assert f(1<<n)-f((1<<n)+2)==(1<<(n-2))-2
    assert (f(4<<n)-f(5<<n))%4==2
for a in range(2,20):
    for b in range(2,20):
        assert f(((1<<a)-1)<<b)==(1<<(a+b))-(1<<a)-(1<<b)+2
for a in range(1,20):
    for b in range(1,20):
        y=int('1'+'0'*a+'1'+'0'*b,2)
        assert f(y)==(1<<(a+b+1))-(1<<a)+(1<<b)
        assert (f(y)==(1<<(y.bit_length()-1)))==(a==b)
print('sparse programming, unrestricted suffixes, nested literal seam phases, and binary controls pass')
