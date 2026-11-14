from math import lcm,gcd

def req(ok):
    if not ok:raise RuntimeError('replay mismatch')

def factors(n):
    out=[]
    p=2
    while p*p<=n:
        if n%p==0:
            e=0
            while n%p==0:n//=p;e+=1
            out.append((p,e))
        p+=1
    if n>1:out.append((n,1))
    return out

def matmul(A,B,m):
    return [[sum(A[i][k]*B[k][j] for k in range(3))%m for j in range(3)] for i in range(3)]

def W(n,a,b,m):
    # sum_{k<n}(n-k) a^(n-k-1)b^k, valid without division
    A=[[a,a,1],[0,a,1],[0,0,b]]
    Z=[[1,0,0],[0,1,0],[0,0,1]]
    while n:
        if n&1:Z=matmul(Z,A,m)
        n//=2;A=matmul(A,A,m)
    return Z[0][2]

def Q(d,n,z,m):
    ans=pow(z,d*n+1,m)
    for e in (-1,1):
        w=W(n,pow(z,d,m),pow(z+e,d,m),m)
        for j in range(1,d):ans+=w*pow(z,d-j,m)*pow(z+e,j-1,m)
    return ans%m


def template(d,n,z,m):return Q(d,n,z,m)

def ordered_product(L,m):
    x=1
    for v in L:x=x*v%m
    return x

def inv(t,m):return pow(t,-1,m)

def basef(d,p,x):
    if p==5:
        return (((3-x)*inv(2,p)) if d==3 else ((x*x+x+1)*inv(8,p)))%p
    return (x+4)*inv(5,p)%p

def pairbanks(p):
    d=3 if p==5 else 2
    A={}
    for x in range(1,p):
        y=inv(x,p);q=basef(d,p,x)*basef(d,p,y)%p
        if q:A[q]=[(d,x),(d,y)]
    if p==5:A[3]=[(3,2),(2,3)]
    P={1:[]}
    for _ in range(3):
        P={q*a%p: L+B for q,L in P.items() for a,B in A.items()}
    req(len(P)==p-1)
    return P

def hsolve(d,p,e,P,x):
    n=0
    for t in range(e):
        step=P*p**t; m=p**(t+1)
        L=[n+step*j for j in range(p) if Q(d,n+step*j,1,m)==x%m]
        req(len(L)==1); n=L[0]
    return n

def packet(p,e,earlier,future,h,q):
    d=3 if p==5 else 2
    m=p**e;h%=m;q%=m
    step=lcm(p-1,*[r**(f+2)*(r-1) for r,f in earlier])
    req(gcd(step,p)==1)
    f=lambda x:basef(d,p,x)
    if q%p==0:
        b=1;x=(-4 if p!=5 else 3)%p
        req(f(x)==0 and ((x*x-b)%p)!=0)
    else:
        b=next(i for i in range(1,p) if pow(i,(p-1)//2,p)==p-1)
        x=next(i for i in range(1,p) if f(i)*f(b*inv(i,p)%p)%p)
    ht=h*inv(b,m)%m
    if not ht%p:
        remaining=[(2,ht)]
    else:
        u=next(i for i in range(1,p) if f(i)*f(ht*inv(i,p)%p)%p)
        remaining=[(d,u),(d,ht*inv(u,m)%m)]
    if q%p:
        q0=f(x)*f(b*inv(x,p)%p)%p
        for dd,v in remaining:q0=q0*basef(dd,p,v%p)%p
        corrections=pairbanks(p)[q*inv(q0,p)%p]
        for j in range(0,len(corrections),2):
            (dd,u),(ff,_)=corrections[j:j+2]
            remaining.extend([(dd,u),(ff,inv(u,m))])
    fixed=[]
    qq=1
    for dd,v in remaining:
        n=hsolve(dd,p,e,step,v)
        fixed.append((dd,n))
        qq=qq*Q(dd,n,2,m)*inv(2,m)%m
    req(qq%p!=0)
    for t in range(1,e):
        mod=p**(t+1)
        opts=[]
        for j in range(p):
            u=x+p**t*j
            v=b*inv(u,m)%m
            a=hsolve(d,p,e,step,u); c=hsolve(d,p,e,step,v)
            if (qq*Q(d,a,2,mod)*Q(d,c,2,mod)*inv(4,mod)-q)%mod==0:opts.append(u)
        req(len(opts)==1);x=opts[0]
    fixed+=[(d,hsolve(d,p,e,step,x)),(d,hsolve(d,p,e,step,b*inv(x,m)%m))]
    F=[]
    for dd,n in fixed:
        L=step*m
        for y in future:
            L=lcm(L,y-1)
            req(gcd(L,y)==1)
            n=next(n+L*j for j in range(y) if Q(dd,n+L*j,1,y)*Q(dd,n+L*j,2,y)%y)
            L*=y
        if not n:n=L
        F.append((dd,n))
    a=ordered_product([Q(v,n,1,m) for v,n in F],m)
    b=ordered_product([Q(v,n,2,m)*inv(2,m) for v,n in F],m)
    req((a,b)==(h,q))
    return F

def compile_(spec,H,Q2):
    M=ordered_product([p**e for p,e in spec],10**100)
    h=q=1; output=[];padding=0;prev=[]
    for i,(p,e) in enumerate(spec):
        m=p**e
        F=packet(p,e,prev,[q for q,_ in spec[i+1:]],H*inv(h,m),Q2*inv(q,m))
        prev.append((p,e))
        period=lcm(*[(r-1)*r**(f-1) for r,f in prev])
        pad=(-len(F))%period
        h=h*ordered_product([Q(d,n,1,M) for d,n in F],M)%M
        q=q*ordered_product([Q(d,n,2,M) for d,n in F],M)*pow(2,pad,M)%M
        output+=F;padding+=pad
        req((h% m,q% m)==(H%m,Q2%m))
    req((h,q)==(H%M,Q2%M))
    return output,padding

for specs in ([(5,2)],[(7,2)],[(17,2)],[(5,1),(7,1),(11,1),(19,1)],[(5,2),(11,2)]):
    p=specs[0][0]
    for H,Q2 in ((0,0),(0,1),(1,0),(2,3),(17,16),(p,1),(5*p,7*p+1),(41,0)):
        compile_(specs,H,Q2)
print('joint prescribed H and Q(2) over prime powers and mixed odd moduli checked')

def helper(p,e,allp):
    m=p**e
    step=lcm(p-1,*[q*(q-1) for q in allp if q<p])
    if p==7:
        n=10;req(step==6 and Q(2,n,3,7)==0)
        start=1
    else:n=0;start=0
    for t in range(start,e):
        mod=p**(t+1)
        choices=[n+j*step*p**t for j in range(p) if Q(2,n+j*step*p**t,3,mod)==0]
        req(len(choices)==1);n=choices[0]
    step*=m
    for q in allp:
        if q<=p:continue
        step=lcm(step,q-1)
        n=next(n+j*step for j in range(q) if Q(2,n+j*step,1,q)*Q(2,n+j*step,2,q)%q)
        step*=q
    req(n>0 and Q(2,n,3,m)==0)
    req(all(Q(2,n,1,q)*Q(2,n,2,q)%q for q in allp))
    return 2,n

def strong_target(spec,A,B):
    from math import prod
    M=prod(p**e for p,e in spec)
    V=[helper(p,e,[q for q,_ in spec]) for p,e in spec]
    u=ordered_product([Q(d,n,1,M) for d,n in V],M)
    v=ordered_product([Q(d,n,2,M) for d,n in V],M)
    req(ordered_product([Q(d,n,3,M) for d,n in V],M)==0)
    q=(B+2*A)*inv(6,M)%M
    h=(A-q)%M
    T,pad=compile_(spec,h*inv(u,M),q*inv(v,M))
    U=V+T
    x=ordered_product([Q(d,n,1,M) for d,n in U],M)
    y=ordered_product([Q(d,n,2,M) for d,n in U],M)*pow(2,pad,M)%M
    z=ordered_product([Q(d,n,3,M) for d,n in U],M)*pow(3,pad,M)%M
    req(((x+y)%M,(4*y+2*z-2*x)%M)==(A%M,B%M))

for spec in ([(7,3)],[(11,2)],[(7,2),(11,2),(19,1)]):
    for A,B in ((0,0),(0,1),(1,0),(137,303),(1001,2025)):
        strong_target(spec,A,B)
print('strong joint wrap controls including exceptional 7-prime phase pass')

def strong5_packet(e,A,B,future):
    p=5;m=p**e
    q=(B+2*A)*inv(8,p)%p;h=(A-q)%p
    d=3;f=lambda x:basef(3,5,x)
    if q==0:b,x=1,3
    else:
        b=2
        x=next(i for i in range(1,p) if f(i)*f(b*inv(i,p)%p)%p)
    ht=h*inv(b,p)%p
    # The first letter below is the free H coordinate.
    if not ht:dd,u,static=2,0,[]
    else:
        dd=3
        c=next(i for i in range(1,p) if f(i)*f(ht*inv(i,p)%p)%p)
        u=ht*inv(c,p)%p;static=[(3,c)]
    if q:
        q0=f(x)*f(b*inv(x,p)%p)*basef(dd,p,u)%p
        for ff,v in static:q0=q0*basef(ff,p,v)%p
        corrections=pairbanks(p)[q*inv(q0,p)%p]
        for j in range(0,len(corrections),2):
            (ff,v),(gg,_)=corrections[j:j+2]
            static.extend([(ff,v),(gg,inv(v,m))])
    fix=[(ff,hsolve(ff,p,e,4,v)) for ff,v in static]
    def model(u,x):
        return fix+[(dd,hsolve(dd,p,e,4,u)),(3,hsolve(3,p,e,4,x)),(3,hsolve(3,p,e,4,b*inv(x,m)%m))]
    def wrapped(F,mod):
        a,c,k=[ordered_product([Q(gg,n,z,mod)*inv(z,mod)%mod for gg,n in F],mod) for z in (1,2,3)]
        return ((a+c)%mod,(4*c+2*k-2*a)%mod)
    req(wrapped(model(u,x),5)==(A%5,B%5))
    for t in range(1,e):
        mod=5**(t+1)
        opts=[(v,w) for v in (u+j*5**t for j in range(5)) for w in (x+j*5**t for j in range(5)) if wrapped(model(v,w),mod)==(A%mod,B%mod)]
        req(len(opts)==1);u,x=opts[0]
    F=[]
    for gg,n in model(u,x):
        D=4*m
        for r in future:
            D=lcm(D,r-1)
            n=next(n+D*j for j in range(r) if Q(gg,n+D*j,1,r)*Q(gg,n+D*j,2,r)%r);D*=r
        if not n:n=D
        F.append((gg,n))
    req(wrapped(F,m)==(A%m,B%m))
    return F,(-len(F))%((p-1)*p**(e-1))

def strong_full(spec,A,B):
    from math import prod
    M=prod(p**e for p,e in spec);e0=spec[0][1];req(spec[0][0]==5)
    prev=[spec[0]]
    U,pad=strong5_packet(e0,A,B,[p for p,_ in spec[1:]])
    def actual(z):
        return ordered_product([Q(d,n,z,M) for d,n in U],M)*pow(z,pad,M)%M
    for i,(p,e) in enumerate(spec[1:],1):
        m=p**e
        prior=lcm(*[r**(f+2)*(r-1) for r,f in prev]);D=lcm(p-1,prior)
        if p==7:
            n=next(prior*j for j in range(42) if (prior*j)%6==4 and (prior*j)%7==3)
            start=1;req(Q(2,n,3,p)==0)
        else:n=0;start=0
        for t in range(start,e):
            mod=p**(t+1)
            opts=[n+j*D*p**t for j in range(p) if Q(2,n+j*D*p**t,3,mod)==0]
            req(len(opts)==1);n=opts[0]
        D*=m
        for q,_ in spec[i+1:]:
            D=lcm(D,q-1)
            n=next(n+j*D for j in range(q) if Q(2,n+j*D,1,q)*Q(2,n+j*D,2,q)%q)
            D*=q
        q=(B+2*A)*inv(6,m)%m;h=(A-q)%m
        # The killer's normalized Q2 has factor 1/2, since all letters will be padded jointly.
        F=packet(p,e,prev,[q for q,_ in spec[i+1:]],h*inv(actual(1)*Q(2,n,1,m),m),q*inv(actual(2)*Q(2,n,2,m)*inv(2,m),m))
        F.append((2,n))
        prev.append((p,e))
        phi=lcm(*[(r-1)*r**(f-1) for r,f in prev])
        pad+=(-len(F))%phi;U+=F
    h,q,z=[actual(i) for i in (1,2,3)]
    req(((h+q)%M,(4*q+2*z-2*h)%M)==(A%M,B%M))
for spec in ([(5,2)],[(5,3)],[(5,2),(7,2)],[(5,1),(7,1),(11,1),(19,1)]):
    for A,B in ((0,0),(0,1),(1,0),(137,303),(1001,2025)):
        strong_full(spec,A,B)
print('strong wrapped 5-adic tangent and all-prime joint controls pass')
