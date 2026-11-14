"""Construct physical terminal parameters in a prescribed 2,3,5,7 congruence class."""
from math import gcd
from itertools import product
from replay_terminals import C0,SOURCES,jac
M=((0,1,0),(0,0,1),(1,1,1))

def multiply(A,B,m):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3))%m
                       for j in range(3)) for i in range(3))

def Cmod(q,m):
    out=((1,0,0),(0,1,0),(0,0,1)); base=M
    while q:
        if q&1:out=multiply(out,base,m)
        base=multiply(base,base,m);q>>=1
    return tuple(sum(out[0][i]*C0[i][j] for i in range(3))%m for j in range(3))

def gate(q,r,m):
    e,f,g=Cmod(q,m);y=pow(2,r,m)
    return (e+y*f)%m,(f+y*g)%m

def crt(constraints):
    x,m=0,1
    for y,n in constraints:
        g=gcd(m,n)
        assert (y-x)%g==0
        x+=m*((y-x)//g*pow(m//g,-1,n//g)%(n//g))
        m*=n//g;x%=m
    return x,m

def vector_lift(p,L,R,q,r,s,u,v):
    # Initial (q,r) fixed mod (L,R), correct mod p.
    for k in range(1,s):
        modulus=p**(k+1)
        hit=next((Q,S) for i,j in product(range(p),repeat=2)
                 for Q,S in [(q+i*L*p**(k-1),r+j*R*p**(k-1))]
                 if gate(Q,S,modulus)==(u%modulus,v%modulus))
        q,r=hit
    return q,L*p**(s-1),r,R*p**(s-1)

def compile_target(n,A=5,B=3,C=3,D=2):
    # H/7 will be congruent to odd n modulo 2^A*3^B*5^C*7^D.
    assert n%2 and min(A,B,C,D)>=0
    N=7*n
    a,T,u,P=next(source for source in SOURCES if source[0]%16==N%16)
    s=max(B,2)
    def eval3(q,r,t,m):
        U,V=gate(q,r,m)
        return (a*U+(P<<t)*V)%m
    if N%3:
        t=4;r3=3
        q3=next(q for q in range(13) if eval3(q,r3,t,3)==N%3
                and eval3(q+13,r3,t,9)!=eval3(q,r3,t,9))
        for j in range(1,s):
            q3=next(q for i in range(3) for q in [q3+i*13*3**(j-1)]
                    if eval3(q,r3,t,3**(j+1))==N%3**(j+1))
    else:
        r3=3
        q3,t=next((q,t) for q,t in product(range(13),range(4,10))
                 if eval3(q,r3,t,9)==N%9 and eval3(q,r3+2,t,9)!=N%9)
        for j in range(2,s):
            r3=next(r for i in range(3) for r in [r3+i*2*3**(j-1)]
                    if eval3(q3,r,t,3**(j+1))==N%3**(j+1))
    b=P<<t
    q0=next(q for q in (12,28,44) if q%3==q3%3)
    # 7: fix second divided coefficient zero, since a is a unit.
    m7=7**D
    U7=7*(n*pow(a,-1,m7)%m7)
    q7,mq7,r7,mr7=vector_lift(7,48,3,q0,0,D+1,U7,0)
    # 5: lift from any odd-r nonsingular representative of the coefficient vector.
    s5=max(C,1); m5=5**s5
    if a%5:U5,V5=N*pow(a,-1,m5)%m5,0
    else:U5,V5=0,N*pow(b,-1,m5)%m5
    q0_5,r0_5=next((q,r) for q,r in product(range(31),(1,3))
        if gate(q,r,5)==(U5%5,V5%5)
        for du,dv,eu,ev in [jac(5,31,4,q,r)] if (du*ev-dv*eu)%5)
    q5,mq5,r5,mr5=vector_lift(5,31,4,q0_5,r0_5,s5,U5,V5)
    # 2: take r >= A so width contributes 0 mod 2^A.
    q2=q0
    for j in range(4,A):
        m=1<<(j+1)
        q2=next(q for q in (q2,q2+48*2**(j-4))
                for e,f,_ in [Cmod(q,m)] if (a*e+b*f)%m==N%m)
    q,_=crt([(q2,16*2**max(A-4,0)),(q3,13*3**(s-1)),(q5,mq5),(q7,mq7)])
    r,mr=crt([(r3,2*3**(s-1)),(r5,mr5),(r7,mr7)])
    while r<max(A,5):r+=mr
    m=7*2**A*3**B*5**C*7**D
    U,V=gate(q,r,m)
    assert (a*U+b*V-N)%m==0
    return T,u,t,q,r

if __name__=='__main__':
    for n in range(-37,40,2):
        for exponents in ((0,0,0,0),(3,6,4,2),(11,2,3,5),(8,8,2,1)):
            compile_target(n,*exponents)
    for n in (1,3,7,9,25733,180131):
        print(n,compile_target(n,8,4,3,3))
