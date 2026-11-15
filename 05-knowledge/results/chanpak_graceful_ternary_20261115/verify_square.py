from itertools import product

def p(x):
    w=0
    while x:
        x,r=divmod(x,2)
        if not r:w+=x+(x&1)
    return w

def h(x,y,j):
    return 432*x*y+129*x+258*y+15*(x&1)+30*(y&1)+45*(p(x)+p(y))+135+2*j

ALPHABET=((6,17),(10,15),(15,10),(17,6))
BY_A={a*a:(i,b) for i,(a,b) in enumerate(ALPHABET)}
ROOT=(1,3)
Q=324

def extend(seed,a,b,c=18):
    x,y=seed; S=x+y
    u,v=a*a*S-x,b*b*S-y
    if (u+x,x+y,y+v,u+v)!=(a*a*S,S,b*b*S,c*c*S):
        raise ValueError("absorption identities")
    return u,v

def decode_pair(x,y):
    rev=[]
    while (x,y)!=ROOT:
        S=x+y
        if x<1 or y<1 or S%Q:return None
        T=S//Q
        A=(x+T-1)//T
        if A not in BY_A:return None
        i,b=BY_A[A]
        rev.append(i)
        x,y=A*T-x,b*b*T-y
    return tuple(reversed(rev))

def select(N):
    j=2*N%3
    if h(*ROOT,j)==N:return ()
    S=sum(ROOT)*Q
    while N>700*S*S:S*=Q
    if N<35*S*S:return None
    L,R=1,2*(S//3)
    def read(i):
        t,flip=divmod(i,2)
        if flip: return t,S-t
        return S-t,t
    # Indices 2,3 are the first positive-minimum pair; 2(S//3)+1 is last.
    L,R=2,R+1
    while L<R:
        m=(L+R)//2
        if h(*read(m),j)<N:L=m+1
        else:R=m
    xy=read(L)
    return decode_pair(*xy) if h(*xy,j)==N else None

checks=0
for d in range(7):
    seen=set()
    for w in product(range(4),repeat=d):
        xy=ROOT
        for i in w:xy=extend(xy,*ALPHABET[i])
        if xy in seen or decode_pair(*xy)!=w:
            raise ValueError("word collision")
        seen.add(xy)
        for j in range(3):
            if select(h(*xy,j))!=w:
                raise ValueError("scalar inverse")
            checks+=1
A=(5,5,7);B=(29,29,41)
for x in range(1,22):
    z=(x,25-x)
    if extend(extend(z,*A),*B)!=extend(extend(z,*B),*A):
        raise ValueError("commutation")
print(checks,"exact scalar inverses through depth 6; Pell commutation checked")
