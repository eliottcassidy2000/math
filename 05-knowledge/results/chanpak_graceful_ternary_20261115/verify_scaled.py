from math import isqrt
from random import Random

def p(x):
    w=0
    while x:
        x,b=divmod(x,2)
        if not b:w+=x+(x&1)
    return w

def h(x,y,j):
    return 432*x*y+129*x+258*y+15*(x&1)+30*(y&1)+45*(p(x)+p(y))+135+2*j

def inverse(N,G,upper):
    if G<30 or G%30 or not N&1:return None
    t,r=divmod(G,90)
    d=2*t+bool(r)
    j=2*N%3
    b=864*d+516
    c=258*d+135+2*j
    if N<c:return None
    v=(isqrt(b*b+3456*(N-c))-b)//1728
    if v<1 or (r and (v&1)!=(r==30)):return None
    x,y=(2*(v+d),v) if upper else (2*v,v+d)
    return (x,y,j) if h(x,y,j)==N else None

def need(b):
    if not b:raise ValueError('scaled switch mismatch')

if __name__ == "__main__":
    rng=Random(20261115)
    for _ in range(1600):
        d=rng.randint(1,1<<160)
        v=rng.randint(1,1<<160)
        j=rng.randrange(3)
        u=v+d
        a,b=h(2*u,v,j),h(2*v,u,j)
        G=a-b
        need(G==45*d+15*((u&1)-(v&1)) and G%30==0)
        need(inverse(a,G,True)==(2*u,v,j) and inverse(b,G,False)==(2*v,u,j))
        need(inverse(a+6,G,True) is None and inverse(b+6,G,False) is None)

    for D in range(30,1000,30):
        expected={}
        for v in range(1,100):
            d=2*(D//90)+bool(D%90)
            for j in range(3):
                a,b=h(2*(v+d),v,j),h(2*v,v+d,j)
                if a-b==D:
                    expected[a,True]=(2*(v+d),v,j)
                    expected[b,False]=(2*v,v+d,j)
        for (N,t),w in expected.items():
            need(inverse(N,D,t)==w)
    print("1600 large scaled swaps, near-miss rejection, and 33 short gap lanes")
