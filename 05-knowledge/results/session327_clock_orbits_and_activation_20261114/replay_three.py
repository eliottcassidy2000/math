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

def profile(d,n,M):
    h=Q(d,n,1,3*M)
    q=Q(d,n,2,3*M)*pow(2,-(d*n+1),3*M)%(3*M)
    req((q-h)%3==0)
    return h%M, ((q-h)//3)%M

def product(a,b,M):
    h,d=a; u,v=b
    return h*u%M,(h*v+u*d+3*d*v)%M

def length(x,mod,parity,floor):
    while x%2!=parity or x<floor:
        x+=mod
    return x

def compile3(h,d,e,floor=None):
    if e==0:
        return []
    M=3**e; floor=2*e+20 if floor is None else floor
    h%=M;d%=M
    if h%3==0 and d%3==0:
        n0=length(8,18,0,floor)
        g,j=profile(3,n0,3*M)
        req(g%3==j%3==0 and g%9!=0)
        u,v=g//3,j//3
        H=(h//3)*pow(u,-1,M//3)%(M//3)
        D=(d//3-v*H)*pow(u+3*v,-1,M//3)%(M//3)
        return [(3,n0)]+compile3(H,D,e-1,floor)
    unit=h%3!=0
    da=2 if unit else 3
    ca=0 if unit else 2
    def eval_(n,m,M):
        return product(profile(da,n,M),profile(3,m,M),M)
    n=next(n for i in range(3) for n in [length(ca+3*i,9,not unit,floor)]
           if any(eval_(n,length(j,3,0,floor),9)[0]% (3 if unit else 9)==h%(3 if unit else 9)
                  and eval_(n,length(j,3,0,floor),9)[1]%3==d%3 for j in range(3)))
    m=next(m for j in range(3) for m in [length(j,3,0,floor)]
           if eval_(n,m,9)[0]%(3 if unit else 9)==h%(3 if unit else 9) and eval_(n,m,9)[1]%3==d%3)
    for t in range(1,e):
        step=3**t
        for i in range(3):
            for j in range(3):
                N,L=n+6*step*i,m+2*step*j
                H,D=eval_(N,L,9*step)
                if H%(step*(3 if unit else 9))==h%(step*(3 if unit else 9)) and D%(3*step)==d%(3*step):
                    break
            else:
                continue
            break
        else:
            raise RuntimeError((h,d,e,t,n,m))
        n,m=N,L
    req(eval_(n,m,M)==(h,d))
    return [(da,n),(3,m)]

def wrap_compile(A,C,e,order_parity):
    # Input C is normalized output q, and must match A mod 3.
    M=3**e; D=(C-A)//3
    epsilon=1 if order_parity else -1
    # A=(3epsilon-1)h+9epsilon*d, D'=(1/2-3epsilon)d + slopes*h.
    a=3*epsilon-1; b=9*epsilon
    c=(-pow(4,-1,M)*3 if epsilon==1 else 7*pow(4,-1,M))%M
    f=(pow(2,-1,M)-3*epsilon)%M
    invdet=pow((a*f-b*c)%M,-1,M)
    h=(f*A-b*D)*invdet%M
    d=(a*D-c*A)*invdet%M
    L=compile3(h,d,e)
    n=sum(k*j+1 for k,j in L)
    # In Z/3^e, 2^(3^(e-1)) = -1 for e>=1.
    period=2*3**(e-1)
    n_final=(n+e+period-1)//period*period
    if epsilon==-1:
        n_final+=period//2
    pad=n_final-n
    req(pad>=e and n_final%2==1-order_parity)
    internal_h,internal_q,internal_r=1,pow(2,pad,M),pow(3,pad,M)
    for k,j in L:
        internal_h=internal_h*Q(k,j,1,M)%M
        internal_q=internal_q*Q(k,j,2,M)%M
        internal_r=internal_r*Q(k,j,3,M)%M
    actual_h=(-internal_h+3*internal_q)%M
    actual_q=(-6*internal_h+4*internal_q+10*internal_r)%M
    req(actual_h==A%M and actual_q*pow(2,-n_final-3,M)%M==C%M)
    return L,pad

def actual_slot(rows,z):
    n=len(rows); N=1<<n
    paths=[[0]*n for _ in range(N)]
    H=[0]*N
    for S in range(1,N):
        for j in range(n):
            if S>>j&1:
                rest=S^(1<<j)
                paths[S][j]=1 if not rest else sum(paths[rest][i] for i in range(n) if rows[i]>>j&1)
        H[S]=sum(paths[S])
    H[0]=1
    cnt=[1]+[0]*(N-1)
    for k in range(z):
        nc=[0]*N
        for S in range(N):
            U=S
            while True:
                nc[S]+=H[U]*cnt[S^U]
                if not U: break
                U=(U-1)&S
        cnt=nc
    return cnt[-1]

def literal_wrap(rows,m):
    n=len(rows); end=n+m
    return [x+((1<<m)-1)*(1<<n) for x in rows]+[
        ((1<<(end+1))-1)^((1<<(i+1))-1) for i in range(n,end)]+[(1<<n)-1]

if __name__=='__main__':
    import random
    rng=random.Random(327)
    for rows in ((0,), (2,0), (2,4,1), (8,1,3,6), (10,28,25,16,1)):
        for m in range(1,4):
            for z in range(1,5):
                h,q,r=(actual_slot(rows,t) for t in (z-1,z,z+1))
                expect=z**(m+1)*q + z*((z+1)**m-z**m)*(r-q)+z*((z-1)**m-z**m)*(h-q)
                req(expect==actual_slot(literal_wrap(rows,m),z))
    for e in range(1,5):
        M=3**e
        cases=[(h,d) for h in range(M) for d in range(M)] if e<=2 else [(rng.randrange(M),rng.randrange(M)) for _ in range(120)]
        cases+=[(0,0),(0,1),(0,2),(3,0),(3,1),(3,2),(9,3)]
        for h,d in cases:
            a=(1,0)
            for k,j in compile3(h,d,e):
                a=product(a,profile(k,j,M),M)
            req(a==(h%M,d%M))
            for parity in (0,1):
                wrap_compile(h,h+3*d,e,parity)
    print("p=3 complete compatible joint digits, both parities, and literal wrap agree")
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


def adjust_future(d,n,L,future):
    extra5=(11 in future and d==3)
    if extra5 and L%5==0:
        req(n%5!=2)
    for p in sorted(set(future)|({5} if extra5 and L%5 else set())):
        L=lcm(L,p-1)
        req(gcd(L,p)==1)
        n=next(n+L*j for j in range(p) if (not(extra5 and p==5) or (n+L*j)%5!=2)
              and Q(d,n+L*j,1,p)*Q(d,n+L*j,2,p)%p)
        L*=p
    return n or L

def mixed_compile(A,B,spec,parity):
    M=1
    for p,e in spec:
        M*=p**e
    E=next((e for p,e in spec if p==3),0)
    S=3**E
    P=[(p,e) for p,e in spec if p!=3]
    eps=1 if parity==1 else -1
    pad0=(1-parity)
    L=[]
    if E:
        h=A*inv(2,S)%S;q=B*inv(4*eps,S)%S
        req((q-h)%3==0)
        L0=compile3(h,(q-h)//3,E)
        L=[(d,adjust_future(d,n,2*3**(E+2),[p for p,e in P]),1) for d,n in L0]
        N0=sum(d*n+1 for d,n,k in L)
        period=2*3**(E-1)
        pad0=((-N0 + (0 if eps==1 else period//2))-(E%period))%period+E
    pads=pad0
    old=[(3,E)] if E else []
    def evaluate(L,pad,z,N):
        return ordered_product([pow(Q(d,n,z,N),k,N) for d,n,k in L]+[pow(z,pad,N)],N)
    for i,(p,e) in enumerate(P):
        N=p**e
        h=A*inv(2*evaluate(L,pads,1,N),N)%N
        q=B*inv(2*evaluate(L,pads,2,N),N)%N
        future=[r for r,f in P[i+1:]]
        F=[(d,n,1) for d,n in packet(p,e,old,future,h,q)]
        old.append((p,e))
        if E:
            for d in {d for d,n,k in F}:
                g,j=profile(d,2*3**(E+2),S)
                req(g==1)
                kappa=(g+3*j)%S
                order=3**(E-1)
                while order>1 and pow(kappa,order//3,S)==1:
                    order//=3
                count=sum(k for a,n,k in F if d==a)
                times=(-count)%order
                if times:
                    P0=lcm(*[r**(f+2)*(r-1) for r,f in old])
                    n=adjust_future(d,P0,P0,future)
                    F.append((d,n,times))
        L+=F
        per=lcm(*[r**(f-1)*(r-1) for r,f in old])
        pads+=(-sum(k for d,n,k in F))%per
    n=sum(k*(d*a+1) for d,a,k in L)+pads
    period=lcm(*[(p-1)*p**(e-1) for p,e in spec],2)
    m=max(1,(E+period-1)//period)*period
    h,q,r=(evaluate(L,pads,z,M) for z in (1,2,3))
    t=pow(2,m,M);s=pow(3,m,M)
    A1=((3-t)*h+(t-1)*q)%M
    B1=(2*t*q + 2*(s-t)*(r-q) + 2*(1-t)*(h-q))%M
    req((A1,B1)==(A%M,B%M) and ((n+m+1)%2)==parity)
    return L,pads,m,n+m+1

if __name__=='__main__':
    for spec in ([(3,1)],[(3,4)],[(3,3),(5,2)],[(3,5),(5,1),(11,1)],[(3,3),(7,1),(11,1)],[(3,2),(5,1),(11,1),(19,1)],[(5,1),(7,1)]):
        M=1
        for p,e in spec:
            M*=p**e
        pairs=[(0,0),(0,3),(3,0),(9,9)]
        pairs.extend((rng.randrange(M),rng.randrange(M)) for _ in range(8))
        for A,B in pairs:
            for parity in (0,1):
                b=B if spec[0][0]!=3 else B+(pow(-1,parity)*A-B)%3
                mixed_compile(A,b,spec,parity)
    print("mixed odd moduli: modular wrapper and explicit CRT literal lengths passed")
