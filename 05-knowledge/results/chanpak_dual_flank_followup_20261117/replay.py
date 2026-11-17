from functools import cache
from math import gcd
D=(54,32,2,119,6,20,55); mark=26

def check(a,b):
    if a!=b:raise AssertionError((a,b))

def TT(n):
    return tuple((1<<n)-(1<<(j+1)) for j in range(n))

def extreme(m):
    r=list(TT(m+2));r[0]^=1<<(m+1);r[-1]=1
    return tuple(r)

def pack(m,q):
    return q**(m+2)+q*((q+1)**m-(q-1)**m)

def ending(a):
    n=len(a);N=1<<n;H=[1]+[0]*(N-1);D=[[0]*n for _ in range(N)]
    for j in range(n):D[1<<j][j]=1
    for m in range(1,N):
        H[m]=sum(D[m])
        for j,v in enumerate(D[m]):
            if not v:continue
            z=a[j]&~m
            while z:
                b=z&-z;z-=b
                D[m|b][b.bit_length()-1]+=v
    return H,D

def reverse(a):
    return tuple(sum(((x>>i)&1)<<j for j,x in enumerate(a)) for i in range(len(a)))

def split(h):
    P=[]
    for m in range(len(h)):
        s=m;c=0
        while True:
            c+=h[s]*h[m^s]
            if not s:break
            s=(s-1)&m
        P.append(c)
    return P

@cache
def small(a,q,k,t):
    h,e=ending(a);_,s=ending(reverse(a))
    p=h if k==1 else split(h)
    f=len(h)-1
    out=[p[f],0,0,0,0]
    for m in range(f+1):
        A=t**m.bit_count()*p[f^m]
        out[1]+=h[m]*A
        out[2]+=A*sum(c for j,c in enumerate(s[m]) if q>>j&1)
        out[3]+=A*sum(c for j,c in enumerate(e[m]) if not(q>>j&1))
    v=[i for i in range(len(a)) if q>>i&1]
    for i in v:
        n=len(a);d=[[0]*n for _ in range(1<<n)];d[1<<i][i]=1
        for m in range(1,1<<n):
            out[4]+=sum(c for j,c in enumerate(d[m]) if not(q>>j&1))*t**m.bit_count()*p[f^m]
            for j,x in enumerate(d[m]):
                if not x:continue
                z=a[j]&~m
                while z:
                    b=z&-z;z-=b;d[m|b][b.bit_length()-1]+=x
    return tuple(out)

def word(a,k,t):
    n=len(a);e=k+t
    L=sum(k**i*t*e**(n-1-i) for i,x in enumerate(a) if x=='1')
    R=sum(k**(n-1-i)*t*e**i for i,x in enumerate(a) if x=='0')
    Z=sum(k**(i+n-1-j)*t*t*e**(j-i-1) for i,x in enumerate(a) if x=='1'
                                     for j in range(i+1,n) if a[j]=='0')
    return k**n,e**n,L,R,Z

def compose(a,b):
    p,c,L,R,Z=a;s,d,l,r,z=b
    return p*s,c*d,L*d+p*l,R*s+c*r,Z*s+p*z+L*r

def reflect(a):
    p,c,L,R,Z=a
    return p,c,R,L,Z

def join(a,b):
    n=len(a)
    return tuple(r|(((1<<len(b))-1)<<n) for r in a)+tuple(r<<n for r in b)

def blocks(u,w,e,j):
    def fw(x):
        bits=bin(x)[2:]
        return TT(len(bits)), sum(int(b)<<i for i,b in enumerate(bits))
    def glue(x,y):
        a,m=x;b,n=y
        return join(a,b),m|(n<<len(a))
    return glue(fw(u),e),(D,mark),glue(fw(w),j)

def literal(u,w,e,j,flip):
    X,C,Y=blocks(u,w,e,j)
    if flip: X,Y=Y,X
    a,m=X;b,n=C;c,t=Y
    T=join(join(a,b),reverse(c))
    M=m|(n<<len(a))|((((1<<len(c))-1)^t)<<(len(a)+len(b)))
    N=len(T)
    return tuple(r|((1-(M>>i&1))<<N) for i,r in enumerate(T))+(M,)

def readings(u,w,e,j,flip):
    if flip: u,w,e,j=w,u,j,e
    out=[]
    for k in (1,2):
        ps=[]
        for t in (1,-1):
            xx=compose(word(bin(u)[2:],k,t),small(*e,k,t))
            yy=compose(word(bin(w)[2:],k,t),small(*j,k,t))
            cc=small(D,mark,k,t)
            ps.append(compose(compose(xx,cc),reflect(yy)))
        out.append(k*(ps[0][0]+ps[0][4]+ps[1][4]))
    return tuple(out)

def profile(e):
    a,m=e
    return (ending(a)[0][-1],*small(a,m,2,1)[:3],
            sum(x for i,x in enumerate(ending(a)[1][-1]) if m>>i&1))

def prime(n):
    from math import isqrt
    if n<2 or n%2==0:return n==2
    return all(n%d for d in range(3,isqrt(n)+1,2))

def qread(a,slot):
    A=tuple(r&((1<<slot)-1) for r in a[:slot])
    H,E=ending(A);_,S=ending(reverse(A))
    F=(1<<slot)-1
    W=sum(H[F^m]*S[m][8] for m in range(1,F+1) if m.bit_count()>=2 and (m>>8&1))
    assert a[slot]==256 and all(r>>slot&1 for i,r in enumerate(a[:slot]) if i!=8)
    return W,H[F]-E[F][8]+H[F^256]-W

def lucas(n, g, factors):
    from math import gcd,prod
    check(prod(p**k for p,k in factors),n-1)
    check(pow(g,n-1,n),1)
    for p,_ in factors:
        check(prime(p),True)
        check(gcd(pow(g,(n-1)//p,n)-1,n),1)

def main():
    e=((22,12,56,49,34,3),57); j=((10,20,25,18,1),21)
    check(profile(e),(43,424,2259,1230,30)); check(profile(j),(13,104,471,222,9))
    check([readings(6,7,e,j,b) for b in (0,1)],[(570968953,3203502193724),(570968959,3203502193724)])
    q=(270,140,264,272,167,463,31,333,82)
    p=(4,5,7,6,0,1,2,3,8)
    a=[0]*10
    for i,v in enumerate(p):
        a[v]=sum(1<<w for w in p[i+1:] if (v,w)!=(4,8))+(512 if v!=8 else 16)
    a[9]=256
    check(qread(tuple(a),9),(128,1))
    for x,y in ((1,7),(1,8),(3,4),(6,8),(4,6)):
        a[x]^=1<<y;a[y]^=1<<x
    check(tuple(r&511 for r in a[:9]),q)
    check(qread(tuple(a),9),(896,7))
    for f in (0,1):
        h,P=readings(6,7,e,j,f);N=512*P+h
        check(N,1640193694155641+6*f)
        check(7*N,11481355859089487+42*f)
        check(prime(N),True)
    e=(TT(3),6); j=((10,12,9,0),1)
    check([profile(x) for x in (e,j)],[(1,8,27,10,1),(3,24,99,21,0)])
    for k in range(8):
        u,w=(5*4**k-2)//3,(4*4**k-1)//3
        a,b=readings(u,w,e,j,0),readings(u,w,e,j,1)
        check(b,(a[0]+6,a[1]))
        if k==0:
            check(a,(93859,58530716))
            for flip,z in enumerate((a,b)):
                H,_=ending(literal(u,w,e,j,flip))
                check((H[-1],split(H)[-1]),z)
    e=((22,60,56,49,32,1),24);j=(e[0],36)
    check([profile(x) for x in (e,j)],[(19,256,1611,404,10),(19,256,1611,404,7)])
    for u in (1,2,3,4,10,48,100):
        a,b=readings(u,u,e,j,0),readings(u,u,e,j,1)
        check(b,(a[0]+114,a[1]))
        if u==48:
            check(a,(39155020907,2982526817837240))
            check(prime(a[0]) and prime(b[0]),True)
    a,b=(readings(39,39,e,j,k) for k in (0,1))
    check(a,(25924835857,2843189815014472))
    check(b,(a[0]+114,a[1]))
    N=256*a[1]+a[0]
    check(N,727856618568540689)
    lucas(N,3,((2,4),(214639,1),(398207,1),(532241,1)))
    lucas(N+114,5,((2,1),(11,1),(107,1),(457,1),(491,1),(1377976099,1)))
    check((7*N,7*(N+114)),(5094996329979784823,5094996329979785621))
    print("Two-prime sevenfold square, equal-length pump, and common-address prime pair verified.")

if __name__=='__main__':main()
