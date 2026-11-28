"""Exact CRT steering for moving-root Chan–Pak registers (standard library only)."""
from functools import lru_cache
from math import gcd
from random import Random
from itertools import product

I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
RAW = {
 'P':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
 'K':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
 'M':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))
}

def mul(a,b,m):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(4))%m
                       for j in range(4)) for i in range(4))

def mat(v,a,m):
    return tuple(sum(x*y for x,y in zip(row,v))%m for row in a)

def chain(m,*steps):
    a=I; count=0
    for b,n in steps:
        a=mul(b,a,m);count+=n
    return a,count

def power(m,w,j):
    a,_=w; s=I; t=j
    assert j>=0
    while t:
        if t&1:s=mul(a,s,m)
        a=mul(a,a,m)
        t>>=1
    return s,w[1]*j

def word(m,s):return chain(m,*((RAW[x],1) for x in s))
def apply(m,v,w):return mat(v,w[0],m)

def factor(m):
    a=3; out=[]
    while a*a<=m:
        if m%a==0:
            q=1
            while m%a==0:m//=a;q*=a
            out.append((a,q))
        a+=2
    if m>1:out.append((m,m))
    return out

def crt(m,items):
    return sum(v*(m//q)*pow(m//q,-1,q) for q,v in items)%m

def ord2(m):
    x=2%m;n=1
    while x!=1:x=x*2%m;n+=1
    return n

def vec(x,y,g,m):
    return ((x+y+g)%m,g%m,(x+g)%m,g%m)

def xyz(v,m):
    h,e,u,g=v
    assert e==g
    return (u-g)%m,(h-u)%m,g

@lru_cache(None)
def steps(m):
    K=word(m,'K'); order=ord2(m)
    E=chain(m,word(m,'PK'),power(m,K,(-2)%order))
    F=chain(m,word(m,'PPK'),power(m,K,(-3)%order))
    G=chain(m,power(m,K,(-1)%order),E,K)
    Gi=chain(m,power(m,K,(-1)%order),power(m,E,m-1),K)
    C=chain(m,E,G,power(m,E,m-1),Gi)
    J=chain(m,F,power(m,E,(-3*pow(2,-1,m))%m))
    return E,C,J

def evaluate(m,v,actions):
    return apply(m,v,chain(m,*actions))

def finish_y(m,v,y,parts):
    E,C,J=steps(m)
    x,Y,g=xyz(v,m)
    out={C:[],J:[]}
    for p,q in parts:
        op=J if p==3 else C
        x1,y1,g1=xyz(apply(m,v,op),m)
        assert (x1,g1)==(x,g)
        tt=(y-Y)*pow((y1-Y)%q,-1,q)%q
        out[op].append((q,tt))
    return chain(m,*(power(m,op,crt(m,rs)) for op,rs in out.items()))

def fill_unit(m,v,x,y):
    x0,_,g=xyz(v,m)
    assert gcd(g,m)==1
    E,_,_=steps(m)
    S=power(m,E,(2*(x-x0)*pow(g,-1,m))%m)
    T=finish_y(m,apply(m,v,S),y,factor(m))
    w=chain(m,S,T)
    assert apply(m,v,w)==vec(x,y,g,m)
    return w

def hensel_xunit(p,q,g,x,t,y):
    j=0;a=1
    inv4=pow(4,-1,q);inv8=pow(8,-1,q)
    # End y when starting at x0=x-j*g/2 in the target g slice.
    def f(j):return (-x*inv4+g*pow(2,-1,q)-t-j*x+g*j*inv8+g*j*(j+1)*inv4-y)%q
    while a<q:
        assert f(j)%a==0
        j+=a*(-(f(j)//a)*pow(-x,-1,p)%p)
        a*=p
    assert f(j)%q==0
    return j

def primitive(m,x,y,g):
    seed=(1,1,1,1)
    assert gcd(gcd(gcd(x,y),g),m)==1
    if gcd(g,m)==1:
        a=fill_unit(m,seed,0,g-1) # h=g on the g=1 slice
        b=chain(m,a,word(m,'MK'))
        return chain(m,b,fill_unit(m,apply(m,seed,b),x,y))
    # Choose a single unit t for the supplying slice before MPK.
    choices=[]; parts=factor(m)
    i2=pow(2,-1,m);i4=pow(4,-1,m)
    for p,q in parts:
        tv=(-y-x*i4+g*i2)%q if g%p==0 and x%p==0 else 1
        assert gcd(tv,q)==1
        choices.append((q,tv))
    t=crt(m,choices)
    js=[]
    live=[]
    for p,q in parts:
        if g%p!=0:live.append((p,q));j=0
        elif x%p==0:j=0
        else:j=hensel_xunit(p,q,g,x,t,y)
        js.append((q,j))
    j=crt(m,js)
    x0=(x-j*g*i2)%m
    u0=((x0-2*g)*i4)%m
    a=fill_unit(m,seed,0,t-1)
    b=chain(m,a,word(m,'MK'))
    c=fill_unit(m,apply(m,seed,b),u0-t,g-u0)
    d=chain(m,b,c,word(m,'MPK'),power(m,steps(m)[0],j))
    v=apply(m,seed,d)
    assert xyz(v,m)[0]==x
    w=chain(m,d,finish_y(m,v,y,live))
    assert apply(m,seed,w)==vec(x,y,g,m)
    return w

def steer(m,x,y,g):
    seed=(1,1,1,1)
    delta=gcd(gcd(gcd(x,y),g),m)
    if delta==1:return primitive(m,x,y,g)
    a=fill_unit(m,seed,-1,delta)
    b=chain(m,a,word(m,'MK'))
    assert apply(m,seed,b)==(delta%m,)*4
    q=m//delta
    if q==1:return b
    # The matrix of a word reduced modulo q has the same action on delta multiples modulo m.
    c=primitive(q,x//delta,y//delta,g//delta)
    w=chain(m,b,c)
    assert apply(m,seed,w)==vec(x,y,g,m)
    return w

def dyadic(h,u,g):
    v=u+g; s=(v&-v).bit_length()-1
    assert s>=1 and g%2==1 and u%2==1
    small=v>>s
    return s-1,h-(small-g if s==1 else v//2),small

def forward_block(h,u,g,r):
    return h+(u-g if r==0 else (1<<r)*u),(1<<(r+1))*u-g,g

if __name__=='__main__':
    rng=Random(20261128)
    for m in [3,5,9,15,21,27,45,63,75,81,105,135,225,315,675,945,3465,10395]:
        v=product(range(m),repeat=3) if m<=9 else (
            tuple(rng.randrange(m) for _ in range(3)) for _ in range(100))
        count=0; depth=0
        for x,y,g in v:
            w=steer(m,x,y,g)
            V=vec(x,y,g,m)
            assert apply(m,(1,1,1,1),w)==V
            h,e,u,g0=V
            M=(h-e+u,h-e,h+u,h)
            P=(h,0,2*u,g0)
            assert apply(m,(1,1,1,1),chain(m,w,word(m,'M')))==tuple(x%m for x in M)
            assert apply(m,(1,1,1,1),chain(m,w,word(m,'P')))==tuple(x%m for x in P)
            count+=1;depth=max(depth,w[1])
        print('CRT',m,count,'largest expanded witness',depth,flush=True)
    for h in range(1,40,2):
        for g in range(1,h+1,2):
            for u in range(h,60,2):
                for r in range(7):
                    hh,uu,_=forward_block(h,u,g,r)
                    assert dyadic(hh,uu,g)==(r,h,u)
    for h in range(1,1000,2):
        u=2*h+3
        assert dyadic(h,u,1)==(0,-1,h+2)
    print('Dyadic fold and negative-predecessor controls pass')
