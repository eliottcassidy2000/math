from math import gcd
from replay_prefix_inverse import count, F

def need(x):
    if not x: raise AssertionError("prefix phase replay failed")

def order(m):
    if m==1:return 1
    z=2%m;o=1
    while z!=1:z=2*z%m;o+=1
    return o

def move(word,mod=0,state=(0,0)):
    u,g=state
    for ch in word:
        if ch=="A":u*=2
        elif ch=="T":g+=u;u=4*u+1
        else:raise AssertionError(ch)
        if mod:u%=mod;g%=mod
    return u,g

def iinv(a,m):
    return pow(a,-1,m) if m>1 else 0

def rooted_word(M,c,k):
    if M<1 or k<1 or (not M%2 and not c%2):
        raise ValueError("bad phase")
    m=M;b=0
    while not m%2:b+=1;m//=2
    a=max(b-1,0);D2=1<<a
    t=1<<k;Ak=568*t-370;C0=185*t-87
    o=order(m);j=a+(-a-4)%o
    reset="A"*j+"T"
    pre=""
    if m>1:
        d=gcd(50,m)
        u0=((c-C0)*iinv(2*Ak,d))%d
        g0=((C0+2*Ak*u0-c)//d*iinv(50//d,m//d))%(m//d)
        z=pow(2,j,m)
        u0=((u0-1)*iinv(4*z,m))%m
        g0=(g0-z*u0)%m
        D="A"*((-2)%o)+"T"
        S="A"*((-1)%o)+D+"A"
        V=D+D+S*(m-1)
        ell=(4*(g0-u0*(u0-1)*iinv(8,m)))%m
        pre=V*ell+D*u0
        need(move(pre,m)==(u0,g0))
    word=pre+reset
    u,g=move(word,D2)
    need(u==1%D2)
    if a:
        want=(Ak-(c-C0)//2)*iinv(25,D2)%D2
        shifts=(want-g)*iinv(m,D2)%D2
        for i in range(a):
            if shifts>>i&1:
                d=a+(-a-i-4)%o
                loop="A"*i+"T"+"A"*d+"T"
                word+=loop*m
    u,g=move(word,M)
    need((C0+2*Ak*u-50*g-c)%M==0)
    return word

def main():
    for M in (1,2,3,4,5,7,9,12,14,16,25,32,40,49,56,75,80,112):
        for c in range(M):
            if M%2 or c%2:rooted_word(M,c,3)
    for M,c,k in [(12,5,1),(112,95,2),(25,13,4),(16,3,1)]:
        word=rooted_word(M,c,k)
        u,g=move(word)
        r=2*u
        if r.bit_length()<500:
            need(2*g==r+1-F(r) and (count(4*u,k)-c)%M==0)
    print("Physical word phase compiler for odd and mixed moduli: PASS")
if __name__=='__main__':
    main()
