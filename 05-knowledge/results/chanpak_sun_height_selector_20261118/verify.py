from itertools import permutations
minus=((2,4,1,7,15),4)
plus=((2,20,9,3,13),13)
small={(1,-1):(8,1,67,38,15,23,59),(1,1):(4,65,18,7,11,31,61),
       (2,-1):(4,17,162,7,13,27,63,123),(2,1):(66,48,3,7,13,29,62)}

def cols(a):
    n=len(a);f=[[0]*n for _ in range(1<<n)]
    for j in range(n):f[1<<j][j]=1
    for S in range(1,1<<n):
        for j,v in enumerate(f[S]):
            U=a[j]&~S
            while U and v:
                b=U&-U;U-=b;f[S|b][b.bit_length()-1]+=v
    return [1]+[sum(q) for q in f[1:]],f
def apex(a,q):
    n=len(a)
    return tuple(v|((1-(q>>i&1))<<n) for i,v in enumerate(a))+(q,)
def TT(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def ord(a,b):return tuple(v|(((1<<len(b))-1)<<len(a)) for v in a)+tuple(v<<len(a) for v in b)
def sig(a,q):
    h,f=cols(a);M=len(h)-1
    R=sum(h[M^S]*sum(v for j,v in enumerate(row) if not(q>>j&1)) for S,row in enumerate(f))
    return h[-1],R,cols(apex(a,q))[0][-1]
def build(m,sign):
    if m<1 or sign not in (-1,1):raise ValueError('m≥1 and sign ±1')
    if m<3:return small[m,sign]
    a,q=minus if sign<0 else plus
    return apex(ord(TT(m),a),3|(q<<m))
def strong(T):
    n=len(T)
    for v in range(n):
        R=1<<v
        for _ in T:
            for j in range(n):
                if R>>j&1:R|=T[j]
        if R!=(1<<n)-1:return False
    return True


from math import isqrt
from random import Random

def ck(ok,*ctx):
    if not ok:raise RuntimeError(ctx)
def F(u):
    f=1
    while u:
        b=u//2
        if not u&1:f+=b+(b&1)
        u=b
    return f
def H(k,u):
    h=12*(1<<k)-1
    return (36*3**k-h)*u+h*F(u)+2*h-5
def reverse(v,r):
    return int(format(v,'0%db'%r)[::-1],2) if r else 0
def lift(N,L):
    if not N&1 or L<1:raise ValueError('odd target and L>=1')
    u=1
    for b in range(1,L):
        if (N-H(0,u))>>b&1:
            r=u.bit_length()
            v=u+(1<<(r+b))
            ck(H(0,v)-H(0,u)==(1<<b)*(36*(1<<r)-11*reverse(u,r)))
            u=v
    return u
def divceil(x,m):return -((-x)//m)
def select(N):
    if N<1 or not N&1:return None
    W,z,k=36,1,0
    while 24*z<=N:
        h=12*(1<<k)-1
        R=W-h
        A0,r=R,0
        while A0*A0<N:
            A0*=2
            r+=1
        gen=range(1,1<<r,2) if r else (0,)
        inv=pow(A0,-1,h)
        for s in gen:
            C=(1<<r)-1-reverse(s,r)
            g=int(bool(r) and s<(1<<(max(0,r-1))))
            A=A0+h*C
            D=N-(2*h-5)-R*s-h*(F(s)-1)
            lo=max(0,divceil(D-2*h,A+h))
            hi=(D-h)//A
            residue=D*inv%h
            for q in range(lo+(residue-lo)%h,hi+1,h):
                u=(q<<r)+s
                if u>0 and u&1 and H(k,u)==N:return k,u
        k+=1
        W*=3
        z*=3
    return None
def host(k,u):
    D,q=plus
    for _ in range(k):D=apex(D,q)
    w=bin(u)[2:]
    return apex(ord(TT(len(w)),D),sum(int(x)<<i for i,x in enumerate(w))|(q<<len(w)))

def main():
    witnesses={49:(110,100,32,6,79,88,12),51:(110,100,96,6,79,88,8),
        41:(24,37,17,38,10,21),61:(38,116,120,115,65,16,33)}
    for target,T in witnesses.items():
        ck(cols(T)[0][-1]==target,target,cols(T)[0][-1])
    ck(49*51==2499 and 41*61==2501 and 24*26==624)
    for L in range(1,11):
        for N in range(1,1<<L,2):
            u=lift(N,L)
            ck(u.bit_length()<=L*(L+1)//2 and H(0,u)%(1<<L)==N)
    for u in range(1,64,2):
        T=host(0,u)
        ck(strong(T) and cols(T)[0][-1]==H(0,u))
    known=set()
    for k in range(10):
        for u in range(1,2001,2):
            X=H(k,u)
            if 24*u*3**k>10000:break
            known.add(X)
    for N in range(1,10001,2):
        a=select(N)
        ck(bool(a)==(N in known),N,a)
        if a:
            ck(H(*a)==N)
    rng=Random(26)
    for k in range(9):
        for _ in range(8):
            u=2*rng.randrange(40000)+1
            N=H(k,u)
            a=select(N)
            ck(bool(a) and H(*a)==N)
    for N in (25,49,51,239,241,343,239**2,110443,1000000007):
        print(N,select(N))
    print('PASS: exact compact dyadic lifts and complete absolute inverse controls')

if __name__=='__main__':main()
