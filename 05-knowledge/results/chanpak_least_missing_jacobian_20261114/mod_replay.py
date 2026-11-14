from itertools import product

def ck(x):
    if not x: raise AssertionError("modular compiler control")

I=(1,0,0,0) # (P,H)->(uP+b,H+cP+d)
def comp(F,G,n):
    u,b,c,d=F; v,w,e,f=G
    return (u*v%n,(u*w+b)%n,(e+c*v)%n,(f+c*w+d)%n)
def power(F,k,n):
    G=I
    while k:
        if k&1:G=comp(F,G,n)
        F=comp(F,F,n);k//=2
    return G
def apply(F,ph,n):
    u,b,c,d=F;p,h=ph
    return (u*p+b)%n,(h+c*p+d)%n
def b(j,n):
    q=pow(2,j,n)
    return (2*q%n,(2*q-2)%n,q,q)
def ord2(m):
    t=1;q=2%m
    while q!=1:
        q=q*2%m;t+=1
    return t
def larg(j,t,a):
    return a+(j-a)%t
def compile_residue(a,m,N,P_odd=None):
    n=(1<<a)*m
    p,h=2%n,3%n
    for j in range(1,a):
        if ((h^N)>>j)&1:
            p,h=apply(b(j,n),(p,h),n)
    p,h=apply(b(a,n),(p,h),n)
    ck((h-N)%(1<<a)==0 and (p+2)%(1<<a)==0)
    if m==1:return p,h
    t=ord2(m)
    U=b(larg(-1,t,a),n);V=b(larg(0,t,a),n)
    Ui=power(U,m-1,n)
    R=comp(V,Ui,n);Ri=power(R,t-1,n)
    S=comp(comp(comp(R,U,n),Ri,n),power(U,(m-2)%m,n),n)
    C=comp(comp(comp(S,U,n),power(S,m-1,n),n),Ui,n)
    ck(apply(C,(p,h),n)[0]%m==p%m)
    step=3*pow(4,-1,m)%m
    ck((apply(C,(p,h),n)[1]-h)%m==step)
    if P_odd is not None:
        p,h=apply(power(U,(p-P_odd)%m,n),(p,h),n)
    originalp=p%m
    if m%3==0:
        L=larg(1,t,a)
        D=comp(U,power((2,0,1,0),L,n),n)
        W=power(D,1<<max(1,a-1),n)
        for _ in range(2):
            if h%3==N%3:break
            p,h=apply(W,(p,h),n)
            k=(p-originalp)%m
            ck(k%3==0)
            p,h=apply(power(U,k,n),(p,h),n)
        ck(h%3==N%3 and p%m==originalp)
    d=3 if m%3==0 else 1
    r=(N-h)%m
    ck(r%d==0)
    k=0 if m==d else (r//d)*pow(step//d,-1,m//d)%(m//d)
    p,h=apply(power(C,k,n),(p,h),n)
    ck(p%m==originalp and (p+2)%(1<<a)==0 and (h-N)%n==0)
    return p,h

def compile_joint(a,m,P,N):
    n=(1<<a)*m
    w=[(P>>i)&1 for i in reversed(range(a))]
    F=I; bias=0
    for bit in w:
        F=comp((2,1,0,0) if bit else (2,0,1,bias),F,n)
        bias=bit
    u,b,c,d=F
    p_odd=(P-b)*pow(u,-1,m)%m if m>1 else 0
    p0=((p_odd+2)*pow(1<<a,-1,m)%m)*(1<<a)-2 if m>1 else -2
    h0=N-c*p0-d
    ck(h0%2==1)
    ans=apply(F,compile_residue(a,m,h0,p_odd),n)
    ck(ans==(P%n,N%n))

for a in range(1,8):
    for m in range(1,60,2):
        n=(1<<a)*m
        targets=range(1,n,2) if n<=320 else [1,7,13,21,67,n-1]
        for N in targets:
            compile_residue(a,m,N)

for a in range(1,7):
    for m in range(1,32,2):
        for P in (0,2,6,(1<<a)*m-2):
            for N in (1,7,13,21,67):
                compile_joint(a,m,P,N)

def full(w):
    p,s,h=0,0,1
    for x in w:
        if not x:h+=p+s
        p,s=2*p+x,x
    return h

for x in (7,13,21):
    for n in range(2,1+(x-1).bit_length()+1):
        for w in product((0,1),repeat=n-2):
            ck(full((1,)+w+(0,))!=x)
print('Mixed modulus compiler verified for odd parts through 59, binary depths through 7; exact small hole controls pass.')

from functools import cache
def hp(T):
    n=len(T)
    @cache
    def start(S,i):
        R=S^(1<<i)
        return sum(start(R,j) for j in range(n) if T[i]&R&(1<<j)) if R else 1
    ans=sum(start((1<<n)-1,i) for i in range(n))
    start.cache_clear()
    return ans
def painted(w):
    n=len(w)
    return tuple(((1<<n)-(2<<i))|((1<<n) if not b else 0)
        for i,b in enumerate(w))+(sum(b<<i for i,b in enumerate(w)),)
for n in range(2,8):
    for w in product((0,1),repeat=n):
        ck(full(w)==hp(painted(w)))
ck(hp((18,8,3,5,14))==13)
print('All literal one-apex words through seven spine vertices; explicit full-spectrum 13 witness.')
