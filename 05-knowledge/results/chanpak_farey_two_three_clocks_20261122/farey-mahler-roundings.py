"""Exact inverse on the graph of *all* 3x/2 rounding scripts at each depth."""
from argparse import ArgumentParser
from functools import cache
from math import gcd, isqrt

CHARTS = {4:(2465,2293760,1344270,712944,1740,986,423385),
          8:(37265,144834560,84839910,57602544,24360,14906,33838015)}
N61 = (1<<60)*((1<<61)-1)-1

def f(n):
    z = 0
    while n:
        n,b=divmod(n,2)
        if not b: z+=n+(n&1)
    return z

def height(C,x,y):
    h,q,c,t,e,r,k=C
    return k+q*x*y+t*x+c*y+r*(x&1)+e*(y&1)+h*(f(x)+f(y))

def ceil(a,b): return -(-a//b)

def bounds(x,d):
    lo=hi=x
    for _ in range(d):
        lo=3*lo//2
        hi=(3*hi+1)//2
    return lo,hi

@cache
def words(b):
    # All padded reversals and literal f reads, indexed by low suffix.
    if b==0: return ((0,0),)
    prev=words(b-1)
    return tuple((rev+((j&1)<<(b-1)),F+(0 if j&1 else ((j>>1)+((j>>1)&1))))
                 for j in range(1<<b) for rev,F in (prev[j>>1],))

def row(h,s,M,l,u,p):
    if M<0 or u<l: return
    l=max(l,ceil(M,s+h));u=min(u,M//s)
    if u<l:return
    m=2*h
    g=gcd(s,m)
    if M%g:return
    m//=g
    a=M//g*pow(s//g,-1,m)%m
    if m%2==0:
        if a%2!=p:return
    else:
        if a%2!=p:a+=m
        m*=2
    y0=l+(a-l)%m
    if y0>u:return
    scale=max(2,isqrt(M*h)//s)
    B=1<<(scale.bit_length()-1)
    if (u-y0)//m < B//2:
        for y in range(y0,u+1,m):
            if s*y+h*f(y)==M:yield y
        return
    bb=B.bit_length()-1
    d=gcd(m,B)
    pm=m//d
    inv=pow(B//d,-1,pm)
    for j in range(a%d,B,d):
        if (j&1)!=p:continue
        rev,F=words(bb)[j]
        R=M-s*j-h*F
        if R<0:break
        A=s*B+h*(B-1-rev)
        lo=max(0,ceil(l-j,B),ceil(R,A+h))
        hi=min((u-j)//B,R//A)
        b0=(a-j)//d*inv%pm
        for z in range(lo+(b0-lo)%pm,hi+1,pm):
            if A*z+h*(f(z)+(j<B//2)*(z&1))==R:yield B*z+j

def scan(C,N,maxdepth):
    h,q,c,t,e,r,k=C
    for sw in (False,True):
        T,Cc,R,E=(c,t,e,r) if sw else (t,c,r,e)
        A=B=1
        for depth in range(maxdepth+1):
            AA=q*A;BB=T*B+Cc*A;CC=(N-k)*B
            if CC>=0:
                z=(isqrt(BB*BB+4*AA*CC)-BB)//(2*AA)
                for x in range(max(1,z-1),z+2):
                    lo,hi=bounds(x,depth)
                    s=Cc+q*x
                    base=N-k-T*x-R*(x&1)-h*f(x)
                    for p in (0,1):
                        for y in row(h,s,base-E*p,lo,hi,p):
                            pair=(y,x) if sw else (x,y)
                            assert height(C,*pair)==N
                            yield depth,sw,*pair
            A*=3;B*=2

def deterministic(C,N):
    """All depths in both deterministic floor/ceil sheets, both orientations."""
    _,q,c,t,_,_,k=C
    xbound=isqrt(max(0,N-k)//q)
    for sw in (False,True):
        T,Cc=(c,t) if sw else (t,c)
        def affine(x,y): return k+T*x+Cc*y+q*x*y
        for rnd in (0,1):
            xmin=1 if rnd else 2
            test=xmin
            d=0
            hits=[]
            while affine(xmin,test)<=N:
                def Y(x):
                    for _ in range(d):x=(3*x+rnd)//2
                    return x
                lo,hi=xmin,xbound+1
                while hi-lo>1:
                    m=(hi+lo)//2
                    if affine(m,Y(m))<=N:lo=m
                    else:hi=m
                y=Y(lo)
                pair=(y,lo) if sw else (lo,y)
                if height(C,*pair)==N:hits.append((d,*pair))
                test=(3*test+rnd)//2
                d+=1
            if not rnd and height(C,1,1)==N:hits.append(("fixed",1,1))
            yield sw,rnd,d,hits

def check():
    assert all(rev==int(f"{a:010b}"[::-1],2) and ff==f(a)
               for a,(rev,ff) in enumerate(words(10)))
    for C in CHARTS.values():
        for x in range(1,30):
            for d in range(10):
                lo,hi=bounds(x,d)
                y=(lo+hi)//2
                for sw in (False,True):
                    a,b=(y,x) if sw else (x,y)
                    assert (d,sw,a,b) in set(scan(C,height(C,a,b),d))
    # Plant an elongated positive row that uses the low-suffix branch.
    h,s,y=2465,17,(1<<24)+12345
    M=s*y+h*f(y)
    assert y in set(row(h,s,M,1,1<<28,y&1))
    print("rounding shell, both orientations and split-row controls PASS")

if __name__=="__main__":
    a=ArgumentParser(description=__doc__)
    a.add_argument("--depth",type=int,default=72)
    a.add_argument("--check",action="store_true")
    a.add_argument("--deterministic",action="store_true")
    b=a.parse_args()
    if b.check:check()
    elif b.deterministic:
        for j,C in CHARTS.items():
            for sw,rnd,depth,hits in deterministic(C,N61):
                print("spine",j,"swapped",int(sw),"round up",rnd,"depths tested",depth,"hits",hits)
    else:
        for j,C in CHARTS.items():
            print("spine",j,"depths 0 through",b.depth,"hits",list(scan(C,N61,b.depth)),flush=True)
