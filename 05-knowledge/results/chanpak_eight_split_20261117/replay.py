from itertools import product
from math import isqrt

def arc(a,i,j):return a[i]>>j&1

def strong(a):
    n=len(a)
    for i in range(n):
        if arc(a,i,i):raise ValueError('loop')
        for j in range(i):
            if arc(a,i,j)+arc(a,j,i)!=1:raise ValueError('not tournament')
        reached=1<<i
        for _ in a:
            for j,v in enumerate(a):
                if reached>>j&1:reached|=v
        if reached!=(1<<n)-1:return False
    return True

def substitute(R,B):
    b=len(R)-1;n=len(B)
    return tuple(
        ((R[i]&((1<<b)-1))|(((1<<n)-1)<<b if arc(R,i,b) else 0))
        for i in range(b)
    )+tuple((v<<b)|R[b] for v in B)

def add(R,Q):
    n=len(R);N=1<<n
    old=tuple(a|2*N|((0 if Q>>i&1 else 1)*N) for i,a in enumerate(R))
    return old+(Q,N)

def end_tables(a):
    n=len(a);N=1<<n
    d=[[0]*n for _ in range(N)]
    h=[0]*N;h[0]=1
    for j in range(n):d[1<<j][j]=1
    for m in range(1,N):
        h[m]=sum(d[m])
        for j,v in enumerate(d[m]):
            if not v:continue
            s=a[j]&~m
            while s:
                b=s&-s;s-=b;d[m|b][b.bit_length()-1]+=v
    return h,d

def profile(a,Q=0):
    inv=[sum(arc(a,j,i)<<j for j in range(len(a))) for i in range(len(a))]
    h,ends=end_tables(a)
    _,starts=end_tables(inv);f=len(h)-1
    E=[sum(v for i,v in enumerate(row) if Q>>i&1) for row in ends]
    S=[sum(v for i,v in enumerate(row) if Q>>i&1) for row in starts]
    W=sum(S[m]*h[f^m] for m in range(f+1))
    U=sum(S[m]*E[f^m] for m in range(f+1))
    return h[f],sum(h[m]*h[f^m] for m in range(f+1)),W,U

def is_prime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))

TT=lambda n:tuple((1<<n)-(2<<i) for i in range(n))
C3=(2,4,1)
H5=(2,4,9,3)
H223=(2,40,11,49,71,213,143,31)
H599=(96,21,25,51,225,198,142,15)

def ordinal(a,b):
    return tuple(x|(((1<<len(b))-1)<<len(a)) for x in a)+tuple(x<<len(a) for x in b)

DATA=(
  ((38,20,32,7,77,26,47),3,34,5),
  ((116,57,34,181,36,0,62,119),34,22,5),
)
for i,(R,Q,a,b) in enumerate(DATA):
    marked=len(R)-1
    T=add(R,Q)
    for rows in (R,T):
        if sum(arc(rows,j,marked) for j in range(len(rows)))!=1:raise ValueError('nonunique entry')
    if strong(R)!=(i==0) or not strong(T):raise ValueError('strength')
    for B in (TT(1),TT(2),TT(3),C3,H5):
        h,P,*_=profile(B)
        R1=substitute(R,B);T1=add(R1,Q)
        J,_,W,U=profile(R1,Q)
        H=profile(T1)[0]
        if J!=a*P+b*h or not strong(T1) or H!=7*J or 2*W-U!=6*J:raise ValueError((J,H))
        print(i,len(B),(h,P),J,W,U,W-U,H)
a,b=DATA[0][2:]
for B,j,p in ((TT(1),0,73),(TT(1),2,277),(TT(1),4,1093),(H223,3,1115227),(H599,3,2231219)):
    h,P,*_=profile(B); P*=2**j
    if a*P+b*h!=p or not is_prime(p):raise ValueError(p)
    print('prime',h,P,p,7*p)

h,p,*_=profile(H5)
if 34*p*32+5*h!=34841 or not is_prime(34841) or not is_prime(4357) or not is_prime(181):
    raise ValueError("extra prime examples")
if 34*16*8+5!=4357 or 22*8+5!=181:
    raise ValueError("transitive arithmetic")
R,Q,_,_=DATA[1]
def induced(a,s):
    return tuple(sum(arc(a,i,j)<<k for k,j in enumerate(s)) for i in s)
if not strong(induced(R,[i for i in range(len(R)) if i!=5])) or R[5]!=0:
    raise ValueError("second rail terminal")
