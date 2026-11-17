from random import Random

def F(u):
    a,z=0,1
    for w in bin(u)[2:]:
        b=int(w)
        if not b:z+=a+a%2
        a=2*a+b
    return z

def rho_num(k,u):
    return int(format(u,'0'+str(k)+'b')[::-1],2)

def update(h,R,S,k,u):
    if not k:return h+R-S,h+R,h-S
    c=(1<<k)-rho_num(k,u)
    return F(u)*h+(u+1)*R-(1-u%2)*S,(1<<k)*R+c*h,(h if u.bit_length()==k else 0)

def grow(A,root,k,u):
    n=len(A);s=n+k
    out=list(A)+[0]*(k+1)
    out[root]|=1<<s
    out[s]=((1<<n)-1)^(1<<root)
    for j in range(k):
        i=n+j
        out[i]=((1<<n)-1)+sum(1<<v for v in range(i+1,s))
        if u>>(k-1-j)&1:out[s]|=1<<i
        else:out[i]|=1<<s
    return tuple(out),s

def table(A):
    n=len(A)
    D=[[0]*n for _ in range(1<<n)]
    for j in range(n):D[1<<j][j]=1
    for S in range(1,1<<n):
        for i,z in enumerate(D[S]):
            if not z:continue
            p=A[i]&~S
            while p:
                bit=p&-p;p-=bit
                D[S|bit][bit.bit_length()-1]+=z
    return D

def state(A,r):
    D=table(A);M=len(D)-1
    h=[1]+[sum(row) for row in D[1:]]
    AT=tuple(sum(((row>>i)&1)<<j for j,row in enumerate(A)) for i in range(len(A)))
    S=table(AT)[M][r]
    return (h[M],sum(D[V][r]*h[M^V] for V in range(1,M+1)),S),D[M][r],sum(h[V]*h[M^V] for V in range(M+1))

def run(word):
    A,r=(2,4,1),0
    v=(3,3,1)
    for k,u in word:A,r=grow(A,r,k,u);v=update(*v,k,u)
    d,e,p=state(A,r)
    if d!=v or e!=1:raise ArithmeticError(('literal',word,d,e,v))
    return v,p,A,r

def ck(x,y):
    if x!=y: raise ArithmeticError((x,y))

def line(z):
    q=2**z
    t=1+q//4
    k=2*z*t
    u=(2**k-1)//(q+1)-1
    K=2*q-2
    v=(4*2**K-6)//5
    return q,k,u,K,v

def power_inverse(N,z):
    q,k,u,K,v=line(z)
    m=F(u)+2**k-1
    b=u+1
    hf,rf,_=update(0,-1,0,K,v)
    H0=(187*hf+3736*rf)
    hf,rf,_=update(1,q+1,0,K,v)
    H1=(187*hf+3736*rf)
    h,r=divmod(N-H0,H1)
    if r or h<1:return None
    X=(m-1)*h-b
    A=m-1-b
    v,r=divmod(X,A)
    if r or v<1:return None
    j=0
    while v>1:
        v,r=divmod(v,m)
        if r:return None
        j+=1
    return j

def edge(Q,q):
    i=q.bit_length()-1
    t=1+Q//2-q//4
    if t<1:return None
    k=2*i*t;b=(2**k-1)//(q+1);u=b-1
    return k,u,F(u)+(Q+1)*b,b

def finish(q,h):
    K=2*q-2;v=(4*2**K-6)//5
    hf,rf,_=update(h,(q+1)*h-1,0,K,v)
    return 187*hf+3736*rf

def inverse(N):
    from functools import cache
    @cache
    def back(q,h):
        if h==1:return (q,)
        if h<1:return None
        i=q.bit_length()-1;Q=q//2 if q>4 else 4
        while Q<=q//2+(h.bit_length()-1)//i:
            E=edge(Q,q)
            if E is not None:
                k,u,m,b=E;old,r=divmod(h+b,m)
                if old>0 and r==0:
                    w=back(Q,old)
                    if w is not None:return w+(q,)
            Q*=2
        return None
    q=4
    while 2*q-3<=N.bit_length()-1:
        H0=finish(q,0);H1=finish(q,1)-H0
        h,r=divmod(N-H0,H1)
        if r==0:
            w=back(q,h)
            if w is not None:return w
        q*=2
    return None

def read_word(w):
    h=1
    for Q,q in zip(w,w[1:]):
        k,u,m,b=edge(Q,q);h=m*h-b
    return finish(w[-1],h)

def main():
    for z in range(2,11):
        q,k,u,K,v=line(z)
        a=q+1
        for t in range(1,1+q//4+3):
            L=2*z*t
            U=(2**L-1)//a-1
            ck(a*(F(U)-1)-2**L+rho_num(L,U),2*(t-1)-q//2)
        ck((10*F(v)-8*(2**K-rho_num(K,v)),10*(v+1)-8*2**K),(2*a,-2))
        lo=[]
        h,R,S=1,q,0
        for j in range(7):
            ck(a*h-R,1);ck(S,0)
            hf,rf,s=update(h,R,S,K,v)
            ck(5*hf-4*rf,1)
            n=187*hf+3736*rf
            ck(197*hf+3728*rf,n+2)
            ck(power_inverse(n,z),j);ck(power_inverse(n+1,z),None)
            lo.append(n)
            h,R,S=update(h,R,S,k,u)
        if z<=4:print('z, line, widths, first lower endpoint:',z,a,k,K,lo[0])
    q,k,u,K,v=line(2)
    T=(6,4,0);r=2
    T,r=grow(T,r,k,u);ck(state(T,r)[0],update(1,q,0,k,u))
    T=(6,4,0);r=2
    T,r=grow(T,r,K,v);ck(state(T,r)[0][:2],(241,301))
    from random import Random
    R=Random(98)
    for _ in range(115):
        w=[2**R.randrange(2,6)]
        h=1
        for _ in range(R.randrange(6)):
            Q=w[-1];q=2**R.randrange(2,Q.bit_length()+1)
            k,u,m,b=edge(Q,q);old=h
            h=m*h-b
            ck(update(old,(Q+1)*old-1,0,k,u),(h,(q+1)*h-1,0))
            ck(m&-m,q)
            w.append(q)
        N=read_word(w)
        ww=inverse(N)
        if ww is None:raise ArithmeticError(('missed', w))
        ck(read_word(ww),N)
    print('115 arbitrary changing-line target decodes.')
    print('All signed-line widths and 63 exact target inverses; two independent literal DP checks.')

if __name__=='__main__':main()
