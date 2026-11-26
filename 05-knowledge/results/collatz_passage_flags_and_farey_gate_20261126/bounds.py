"""Rational certificate for bounded expanding-packet compression in Collatz."""
from fractions import Fraction as F

def log_interval(x,n=100):
    t=(F(x)-1)/(F(x)+1)
    assert 0<=t<1
    s=2*sum(t**(2*k+1)/(2*k+1) for k in range(n))
    return s, s+2*t**(2*n+1)/((2*n+1)*(1-t*t))

l2,h2=log_interval(2)
l3,h3=log_interval(3)
lo,hi=l3/h2,h3/l2
gamma_base=F(317,200)
mu=F(83,200)
N=1<<71
assert hi<gamma_base and F(133,10)*F(46057,100000)<9*l2
assert (4*N)**200>2**83*(3*N+1)**200

def first_farey(lo,hi):
    a,b,c,d=0,1,1,0
    steps=0
    while True:
        assert c*b-a*d==1
        p,q=a+c,b+d
        if F(p,q)<lo:a,b=p,q
        elif F(p,q)>hi:c,d=p,q
        else:return p,q,steps,(a,b),(c,d)
        steps+=1

P,K,steps,left,right=first_farey(lo,hi+1/(3*N*l2))
t=K.bit_length()-1
_,ln_remainder=log_interval(F(K,2**t))
log_K=t+ln_remainder/l2

def ceil(a):return -((-a.numerator)//a.denominator)
def base_data(B):
    assert 3**B<N
    return B+ceil(gamma_base*B),3**B
def bound(H,r,beta):
    gamma=gamma_base*beta
    assert gamma>1
    return (beta+(gamma-1)/mu)*H*(gamma**r-1)/(gamma-1)+r/mu

def maximum(R,H):
    r=0
    while bound(H,r+1,R)<K:r+=1
    return r, ceil(bound(H,r,R))

def coefficients(w):
    a,b,B=1,1,0
    for e in w:a,b,B=3*a,b*2**e,3*B+b
    return a,b,B

def period(n):
    w=[];x=n
    while True:
        z=3*x+1;e=(abs(z)&-abs(z)).bit_length()-1
        w.append(e);x=z//2**e
        if x==n:return tuple(w)

def valuations(parities):
    ix=[i for i,b in enumerate(parities) if b=='1']+[len(parities)]
    assert ix[0]==0
    return tuple(y-x for x,y in zip(ix,ix[1:]))

def height_lib(W):
    from math import gcd
    # Account for any internal cut by at most B and the true anchor denominators.
    B=max(sum(w) for w in W)
    L=C=0;M=F(0)
    for w in W:
        p,q,b=coefficients(w);assert p>q
        d=p-q;g=gcd(d,b);d//=g;b//=g
        L=max(L,(d-1).bit_length())
        C=max(C,(b+d-1)//d)
        M=max(M,F(len(w)*p,p-q))
    assert C<N
    return B+L,M
def main():
    print("First relaxed admissible denominator",P,K,steps,left,right)
    assert P-1 < K*lo < K*hi < P
    print("Least-counterexample first coefficient drop excluded through shortcut step",P-1)
    for B in [11,22,32,44]:
        e,C=base_data(B)
        H=e+9+F(143,10)*log_K
        for beta in [F(1),F(2,3),F(7,11)]:
            print("B, beta, excluded r, upper bound at barrier",B,str(beta),*maximum(beta,H))
    # Additional reciprocal-sum accounting for three inherited motifs and all rotations.
    wb=valuations("1111111110011011010000")
    assert (len(wb),sum(wb))==(14,22)
    for named in [("negative -5", (period(-5),)),("negative -17",(period(-17),)),
                  ("both signed anchors",(period(-5),period(-17))),
                  ("-17 and 22-step expanding reset",(period(-17),wb))]:
        name,W=named
        W=[w[i:]+w[:i] for w in W for i in range(len(w))]
        e,M=height_lib(W);beta=max(F(len(w),sum(w)) for w in W)
        r,upper=maximum(beta,e+9+F(143,10)*log_K)
        while True:
            r1=r+1
            # ceil log2 of the rational constant, bounded by integer powers.
            z=r1*(M+5);j=0
            while 2**j<z:j+=1
            H=e+9+min(F(143,10)*log_K,j+F(133,10)*log_K)
            upper1=bound(H,r1,beta)
            if upper1>=K:break
            r,upper=r1,ceil(upper1)
        print(name,"beta",str(beta),"excluded r",r,"max bound",upper,"reciprocal price",M)
    # Genuine positive cylinders: expanding packet (1,2) with arbitrarily many rises/falls.
    for r in range(1,33):
        x=2*8**r-5
        for t in range(r):
            assert x==-5+2*9**t*8**(r-t)
            h=3*x+1;assert h%4==2
            y=h//2;h=3*y+1;assert h%8==4
            x=h//4;assert x<y
    print("32 expanding-cylinder controls OK")

if __name__=='__main__':main()
