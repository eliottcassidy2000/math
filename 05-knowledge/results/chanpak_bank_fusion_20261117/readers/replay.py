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

def read(m,u,w,k,turn=False,contact=1):
    T=[]
    for t in (1,-1):
        A=compose(word(bin(u)[2:],k,t),small(extreme(m),contact,k,t))
        B=word(bin(w)[2:],k,t)
        C=small(D,mark,k,t)
        src=compose(compose(B,C),reflect(A)) if turn else compose(compose(A,C),reflect(B))
        T.append(src)
    return k*(T[0][0]+T[0][4]+T[1][4])

def joint(a,b):
    n=len(a)
    return tuple(r|(((1<<len(b))-1)<<n) for r in a)+tuple(r<<n for r in b)

def literal(m,u,w,turn,contact=1):
    def fw(x):
        a=bin(x)[2:];n=len(a);return TT(n),sum(int(b)<<i for i,b in enumerate(a))
    def glue(a,b):
        x,s=a;y,t=b;return joint(x,y),s+(t<<len(x))
    def rf(a):
        x,s=a;return reverse(x),((1<<len(x))-1)^s
    U=glue(fw(u),(extreme(m),contact));W=fw(w);C=(D,mark)
    X=glue(glue(W if turn else U,C),rf(U if turn else W))
    a,q=X;n=len(a)
    return tuple(r|((1-(q>>i&1))<<n) for i,r in enumerate(a))+(q,)

def carry_follow(c,a,b,p,q):
    for A,B in zip(a,b):
        z=c-p*int(B)+q*int(A);check(z%3,0);c=2*z//3
    return c

def addresses(k):
    a='1'+'1001'*k+'1'
    b='0'+'0000'*k+'000001'
    return int(a[::-1],2),int(b[::-1],2)

def check_ids(m,u,w,contact=1):
    U=bin(u)[2:];W=bin(w)[2:]
    x=word(U,2,1);y=word(W,2,1)
    a=compose(x,small(extreme(m),contact,2,1))
    return a[2]*y[0]-a[0]*y[2]

def lift(a,k):
    n=len(a)-1
    B=TT(k);new=((1<<k)-1)<<n
    U=tuple((x&((1<<n)-1))|(new if x>>n&1 else 0) for x in a[:-1])
    return U+tuple((x<<n)|a[-1] for x in B)

def strong(a):
    for i in range(len(a)):
        r=1<<i
        for _ in a:
            for j,x in enumerate(a):
                if r>>j&1:r|=x
        if r!=(1<<len(a))-1:return False
    return True

def reader():
    a=[750,16,18,6,41,78,158,318,127,510]
    flips=[(1,8),(2,6),(3,6),(5,7),(3,7),(0,7),(4,7)]
    prof=[(292,1),(791,117),(1150,147),(1307,195),(1315,273),(1315,335),(2017,19),(2044,7)]
    for i,(alpha,beta) in enumerate(prof):
        check([j for j in range(9) if a[j]>>9&1],[0])
        check(strong(a),True)
        check(ending(lift(a,1))[0][-1],2*alpha+beta)
        check(ending(lift(a,2))[0][-1],4*alpha+beta)
        if i<len(flips):
            j,k=flips[i];a[j]^=1<<k;a[k]^=1<<j
            A,B=prof[i+1]
            check(A>=alpha and 2*(A-alpha)+B-beta>0,True)
    check(prof[-1],tuple(7*x for x in prof[0]))
    return prof

def prime_certificate(n):
    from math import isqrt
    factors=(2,3,8779,829097,14119373)
    check(n-1,2**3*3*8779*829097*14119373)
    for p in factors:
        check(all(p%d for d in range(2,isqrt(p)+1)),True)
    check(pow(5,n-1,n),1)
    for p in factors:check(gcd(pow(5,(n-1)//p,n)-1,n),1)

def main():
    check((pack(6,2),pack(6,3)),(1712,18657))
    for m in range(1,8):
        X=small(extreme(m),1,2,1)
        check((X[0],X[1],X[2]),(pack(m,2),pack(m,3),3**(m+1)-2**m))
    check((small(D,mark,2,1)[2:4],small(D,mark,2,-1)[2:4]),((1397,1641),(-109,-109)))
    for m,u,w in ((1,1,2),(2,1,2),(1,3,2),(1,1,4)):
        for t in (False,True):
            a=literal(m,u,w,t);j=ending(a)[0]
            check((read(m,u,w,1,t),read(m,u,w,2,t)),(j[-1],split(j)[-1]))
        check(read(m,u,w,2,True)-read(m,u,w,2),-488*check_ids(m,u,w))
    p,q,L=small(extreme(2),12,2,1)[:3]
    check((p,q,L),(32,117,36))
    check(sum(ending(extreme(2))[1][-1][j] for j in (2,3)),3)
    check(carry_follow(2*L,'1','0',p,q),126)
    check(carry_follow(126,'1001','0000',p,q),126)
    check(carry_follow(126,'100000','000001',p,q),0)
    for k in range(12):
        u,w=addresses(k)
        check((u.bit_length(),w.bit_length()),(2+4*k,7+4*k))
        check(check_ids(2,u,w,12),0)
        h,p=read(2,u,w,1,contact=12),read(2,u,w,2,contact=12)
        check((read(2,u,w,1,True,12),read(2,u,w,2,True,12)),(h+6,p))
        if k==0:
            print('k=0: u,w,h,P2,lower factor-seven rails:',u,w,h,p,68*p+5*h,34*p+5*h,22*p+5*h)
    prof=reader()
    u,w=addresses(0);h,P=read(2,u,w,1,contact=12),read(2,u,w,2,contact=12)
    n=292*2**21*P+h
    check(n,2466476862736151977);prime_certificate(n)
    print('padding 21 prime and physical sevenfold:',n,7*n)
    print('12 pumped levels; independent raw checks on eight small hosts; exact carry-loop certificate')
if __name__=='__main__':main()
