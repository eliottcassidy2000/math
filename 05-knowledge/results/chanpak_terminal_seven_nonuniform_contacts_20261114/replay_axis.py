from functools import cache
from itertools import product
from math import isqrt

@cache
def C(q):
    if q<3: return ((1,1,1),(1,2,2),(3,3,4))[q]
    return tuple(map(sum,zip(C(q-3),C(q-2),C(q-1))))
def act(q,x,y):
    e,f,g=C(q)
    return e*x+f*y,f*x+g*y

@cache
def word(a,b):
    """Exact inverse on reset-free physical states; return one gap word."""
    if not (a&1 and not b&1 and 1<=a and 2*a<=b):return None
    if a==1:return (b.bit_length()-1,) if not b&(b-1) else None
    q=1
    while True:
        e,f,g=C(q)
        if e+2*f>a:return None
        d=e*g-f*f
        if d:
            B=b
            r=0
            while not B&1:
                B//=2
                r+=1
                aa,bb=g*a-f*B,e*B-f*a
                if aa%d==bb%d==0:
                    P=word(aa//d,bb//d)
                    if P is not None:return P+(0,)*q+(r,)
        q+=1

def ceil_div(a,b):return -((-a)//b)

def seven(p):
    """Decide 7p for terminal q=7, for p any odd prime distinct from seven."""
    r=1
    while True:
        y=1<<r;r+=1
        if 193+230*y>7*p:return None
        if y%7==4:
            A,B=(57+68*y)//7,(68+81*y)//7
            for a in range(p*pow(A,-1,2*B)%(2*B),p//(A+2*B)+1,2*B):
                W=word(a,(p-A*a)//B)
                if W is not None:return W+(0,)*7+(r-1,)
        else:
            l=81*p//(81*y+68)+1;l+=l&1
            for v in range(l,230*p//(230*y+193)+1,2):
                a,b=(81*y+68)*v-81*p,68*p-(68*y+57)*v
                W=word(a,b)
                if W is not None:return W+(0,)*7+(r-1,)

def evaluate(w):
    a,b=1,1<<w[0];q=0
    for r in w[1:]:
        if not r:q+=1;continue
        a,b=act(q,a,b);b<<=r;q=0
    return a+b

def dot(X,Y):return sum(x*y for x,y in zip(X,Y))
def probe(steps):
    a=b=1
    for q,y in steps:a,b=act(q,a,y*b)
    return a,b
def prime(p):return p>1 and all(p%i for i in range(2,isqrt(p)+1))

if __name__=='__main__':
    assert act(7,1,4)==(329,392)
    assert act(11,1,2)==(2209,2632)
    for p in (271,3631):
        assert prime(p) and evaluate(seven(p))==7*p
    assert evaluate((2,)+(0,)*7+(6,))==7*3631
    assert all(seven(p) is None for p in (3,5,11,13,17))
    for q in range(90):
        if q>=2:assert act(q,1,0)==act(q-2,1,2)
        if q>=1:assert act(q,0,1)==act(q-1,1,1)
    checks=[
        (((8,2),), (423939,423941)),
        (((12,2),), (4851811,4851809)),
        (((2,8),(4,2)), (1209689,1209691)),
        (((1,4),(2,8),(4,2)),(11985689,11985691)),
    ]
    for steps,p in checks:
        R=probe(steps)
        assert p==(dot((291,758),R),dot((179,852),R))
    for q in (10,14):
        for r in range(1,20):assert abs(dot((56,-47),probe(((q,1<<r),))))==1
    print('literal decoder, both seven-axis directions, boundary lift and all special probes replayed')
