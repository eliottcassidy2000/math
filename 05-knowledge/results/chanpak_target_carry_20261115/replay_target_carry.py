from math import gcd, isqrt
from random import Random

B,Q,T,S,h = 2855,5760,4034,3914,45
LEVELS={2:5,3:8,4:6,5:10,6:9,8:8}

def need(t, *details):
    if not t: raise AssertionError(details)

def p(x):
    z=0
    while x:
        x,b=divmod(x,2)
        if b==0:z+=x+(x&1)
    return z

def prev(x,r):
    v=0
    for _ in range(r):
        v=(v<<1)|(x&1)
        x>>=1
    return v

def table(M,L):
    D=1<<L
    out=[0]
    for x in range(1,2*M*D):
        u,b=divmod(x,2)
        out.append((out[u]+((u+(u&1))//2 if b==0 else 0))%M)
    return out

def oddpart(M):
    return M//(M&-M)

def branch(s,D):
    # Largest segment of the block; ties go to the no-carry state.
    return (0,s+1,s,0) if s+1>=D-s-1 else (s+1,D,s+D,1)

def certify(M,L):
    D=1<<L
    m=oddpart(M)
    out=table(M,L)
    goal=(1<<(M*m))-1
    for s in range(D):
        lo,hi,c,_=branch(s,D)
        for u in range(0,len(out),D):
            for v in range(0,len(out),D):
                hit=hit1=goal
                for b in range(lo,hi):
                    idx=((out[u+b]+out[v+c-b])%M)*m+(u+b)%m
                    hit &= ~(1<<idx)
                    idx1=((out[u+b]+out[v+c-b]+(b&1))%M)*m+(u+b)%m
                    hit1 &= ~(1<<idx1)
                    if hit==hit1==0:break
                need(hit==hit1==0,M,L,s,u//D,v//D,hit,hit1)
    print('all fixed-sum blocks',M,L,D*(2*M)**2)

def place(N,M):
    L=LEVELS[M]
    D=1<<L
    m=oddpart(M)
    d=gcd(45,M)
    out=table(M,L)
    a=isqrt((N-B)//Q)+2
    r=max(a.bit_length()//4, (M&-M).bit_length())
    q=1<<r
    # The target chooses the full high sum, not a prescribed all-ones suffix.
    while True:
        u,v=divmod(a,q)
        if (B+2*S*q*u+(S+T)*v-N)%d==0: break
        a+=1
    K=2*u
    s=K%D
    lo,hi,c,carry=branch(s,D)
    common=2*p(v)+(q-1-prev(v,r))*K
    need(K%2==0)
    jfirst=(u+D-1)//D
    jlast=K//D-1
    need(jlast>jfirst)
    def witness(j):
        vprefix=K//D-j-carry
        need(vprefix>=0)
        x0=D*j; y0=D*vprefix
        # High address U=0 mod m. All remaining ordinary terms
        # and the shared low-word term are then fixed mod 2M.
        F0=B+Q*v*(q*K+v)+T*v+S*(q*K+v)+45*common
        debt=(N-F0)//2
        need((N-F0)%(2*d)==0)
        mod=M//d
        z=(debt//d)*pow(45//d,-1,mod)%mod if mod>1 else 0
        for b in range(lo,hi):
            U=x0+b; V=y0+c-b
            if U%m==0 and (out[(j%(2*M))*D+b]+out[(vprefix%(2*M))*D+c-b]+(b&1)*(v<q//2)+
                         (p(j)+p(vprefix))//2-
                         (p(j%(2*M))+p(vprefix%(2*M)))//2)%mod==z:
                x,y=q*U+v,q*V+v
                F=B+Q*x*y+T*x+S*y+45*(p(x)+p(y))
                need((F-N)%(2*M)==0 and x+y==2*a,N,M)
                return F
        raise AssertionError((N,M,j))
    left,right=jfirst,jlast
    A,Z=witness(left),witness(right)
    need(A>=N>Z,N,M,A,Z)
    while right-left>1:
        mid=(left+right)//2
        value=witness(mid)
        if value>=N:left,A=mid,value
        else:right,Z=mid,value
    return A,Z

if __name__=='__main__':
    rng=Random(1225)
    for M,L in LEVELS.items():
        certify(M,L)
        for n in (101,127,251,383,719):
            N=(1<<n)+rng.getrandbits(n)|1
            A,Z=place(N,M)
            # Generous fixed-modulus upper bound, through all tested scales.
            need(max(A-N,N-Z)<=10_000_000*(1<<L)*(1<<((3*n+7)//8)),M,n)
    print('30 target-directed brackets, no high-sum residue alignment')
