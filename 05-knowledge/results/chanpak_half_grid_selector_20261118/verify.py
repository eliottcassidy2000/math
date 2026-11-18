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

from collections import deque
from math import gcd
from random import Random

MENU = [(23,3),(37,1),(1,3),(43,1),(23,1),(19,3),(11,1),(5,1),
        (89,1),(7,3),(13,1),(5,3),(47,3),(7,1),(1,1),(35,1)]

def ck(ok,*msg):
    if not ok:
        raise RuntimeError(msg)
def F(u):
    f=1
    while u:
        q=u>>1
        if not u&1:
            f+=q+(q&1)
        u=q
    return f
def S(u):
    return u-F(u)
def Hmod(u,k,M):
    return (36*(u%M)*pow(3,k,M)+(12*pow(2,k,M)-1)*(2-S(u))-5)%M
def menu(N,L,min_k):
    L=max(5,L)
    v,j=MENU[(N%32)//2]
    y=1<<(L-2)
    r=((N-S(v)+7)//4)*pow(9*v,-1,y)%y
    period=1<<(L-4)
    k=next(k for k in range(period) if pow(3,k,y)==r)
    c=max(min_k,L)
    k+=period*max(0,(c-k+period-1)//period)
    return v,k,L
def short_select(N,L,M):
    v,k,L=menu(N,L,M)
    if not M:
        return v,k
    r=max(3,v.bit_length())
    B=int(format(v,'0%db'%r)[::-1],2)
    y=3**M
    target=(2-(N+5)*pow(12*pow(2,k,y)-1,-1,y))%y
    d=L+1+(L+1)%2
    c=1 if (target-S(v))%3 else (1<<d)+1
    mult=(target-S(v)-S(c))*pow(B*c,-1,y)%y
    period=2*3**(M-1)
    b=next(b for b in range(period) if pow(2,b,y)==mult)
    b+=period*max(0,(L-b+period-1)//period)
    u=(c<<(b+r))+v
    ck(S(u)==S(v)+S(c)+B*(c<<b))
    return u,k

def word_for(M,goal):
    start=(2%M,1%M);Q=deque([start]);prev={start:(None,'')}
    while Q:
        a,b=cur=Q.popleft()
        if cur==goal:
            path=[]
            while cur!=start:
                cur,w=prev[cur]
                path.append(w)
            return ''.join(reversed(path))
        for w,u,z in [('1',2*a,b),('01',4*a-2,b+a),('001',8*a-6,b+3*a-2)]:
            y=(u%M,z%M)
            if y not in prev:
                prev[y]=(cur,w);Q.append(y)
    raise RuntimeError((M,goal))

def controller(N,L,M):
    v,k,L=menu(N,L,0)
    if M==1:return v,k
    H=12*2**k-1
    R=36*3**k-H
    target=(N-2*H+5)%M
    U=next(u for u in range(M) if (target-R*u)%gcd(H,M)==0)
    m=M//gcd(H,M)
    V=((target-R*U)//gcd(H,M))*pow(H//gcd(H,M),-1,m)%m if m!=1 else 0
    suffix='1'*L+format(v,'0%db'%L)
    r=L
    while pow(2,r+3,M)!=1:r+=1
    # Invert the fixed suffix via its triangular affine matrix.
    a,b=1,0
    for ch in suffix:
        b+=a+(a&1) if ch=='0' else 0
        a=2*a+int(ch)
    ev,gv=3,0
    for ch in suffix:
        gv+=ev+(ev&1) if ch=='0' else 0
        ev=2*ev+int(ch)
    s=pow(2,len(suffix),M)
    z=(U-a)*pow(s,-1,M)%M   # starting odd value = 1+z
    if z&1:z+=M              # use an even difference to preserve parity
    g=(gv-b)//2
    pre=((2+z)%M,(V-b-g*z)%M)
    w='1'+word_for(M,pre)
    before=F(int(w,2))
    # Dyadic fixed suffix drift, evaluated after guard without W.
    drift=F(int('1'+suffix,2))-F(1)
    # guard erases prior low digits, so this drift is fixed mod 2^L.
    t=((before+drift-F(v))//2)*pow(M,-1,1<<(L-1))%(1<<(L-1))
    u=int(w+('1'*r+'001')*(M*t)+suffix,2)
    ck(u%(1<<L)==v%(1<<L) and F(u)%(1<<L)==F(v)%(1<<L),N,L,M)
    ck(u%M==U and F(u)%M==V,N,L,M)
    return u,k

def literal(u,D,q):
    w=bin(u)[2:];n=len(w)
    return apex(ord(TT(n),D),sum(int(c)<<i for i,c in enumerate(w))|(q<<n))

def main():
    D,q=plus
    for k in range(6):
        h,R,K=sig(D,q);A,E=cols(D)
        z=sum((-1)**i.bit_count()*A[(len(A)-1)^i]*
              sum(n for j,n in enumerate(row) if not q>>j&1)
              for i,row in enumerate(E))
        ck((h,h+R,K,z)==(12*2**k-1,36*3**k,24*2**k-1,7-12*2**k))
        for u in (1,3,5,7,9,19,89):
            g=literal(u,D,q)
            ck(strong(g) and cols(g)[0][-1]==36*u*3**k+h*(2-S(u))-5)
        D=apex(D,q)
    for i,(v,j) in enumerate(MENU):
        rr=max(3,v.bit_length())
        ck((36*v*j+S(v)-7)%32==2*i+1 and int(format(v,'0%db'%rr)[::-1],2)%3)
    rng=Random(23)
    for _ in range(192):
        N=2*rng.randrange(1,100000)+1
        L,M=rng.randint(5,12),rng.randint(0,7)
        u,k=short_select(N,L,M)
        md=2**L*3**M
        ck(Hmod(u,k,md)==N%md)
    for M in (1,3,5,7,9,15,21,25,35,49,63):
        for _ in range(7):
            N=2*rng.randrange(100000)+1
            L=rng.randint(5,7)
            u,k=controller(N,L,M)
            md=2**L*M
            ck(Hmod(u,k,md)==N%md,N,L,M)
    print("PASS: 42 direct strong hosts, signed clone profiles, 192 short and 77 mixed-modulus targets.")

if __name__=='__main__':
    main()
