from itertools import permutations, product, combinations
from math import isqrt,gcd
from functools import cache
from collections import Counter

def need(x):
    if not x:raise ArithmeticError('literal replay failed')

def forward(n):return tuple(((1<<n)-1)^((2<<i)-1) for i in range(n))
def zig(n):return tuple(sum(1<<j for j in range(n)if j<=i-2 or j==i+1) for i in range(n))
def converse(T):
    n=len(T);f=(1<<n)-1
    return tuple(f^(1<<i)^r for i,r in enumerate(T))
def rel(T,p):return tuple(sum((T[i]>>j&1)<<k for k,j in enumerate(p))for i in p)
def mark(T):
    n=len(T); C=[tuple(i for i in range(n-1) if T[i].bit_count()==k) for k in range(n)]
    return min(rel(T,sum(x,())+(n-1,)) for x in product(*(permutations(c) for c in C if c)))
@cache
def canon(T):
    n=len(T); C=[tuple(i for i,v in enumerate(T) if v.bit_count()==k)for k in range(n)]
    return min(rel(T,sum(x,())) for x in product(*(permutations(c) for c in C if c)))
def merged(T):return min(canon(T),canon(converse(T)))
def frame_line(T,p):
    X=rel(T,p); n=len(T); f=(1<<n)-1
    need(all(X[i]>>(i+1)&1 for i in range(n-1)))
    return tuple(f^(1<<i)^x^((1<<(i-1) if i else 0)|(1<<(i+1) if i<n-1 else 0))for i,x in enumerate(X))
def flip(T,edges):
    a=list(T)
    for i,j in edges:a[i]^=1<<j;a[j]^=1<<i
    return tuple(a)
def strong(T):
    n=len(T); f=(1<<n)-1
    for v in range(n):
        z=1<<v
        for _ in T:
            for i,r in enumerate(T):
                if z>>i&1:z|=r
        if z!=f:return False
    return True

def counts(T):
    n=len(T);f=(1<<n)-1
    need(all(0<=r<=f and not(r>>i&1) for i,r in enumerate(T)))
    need(all((T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(n) for j in range(i)))
    C=[[0]*n for _ in range(f+1)]
    for i in range(n): C[1<<i][i]=1
    for S in range(1,f+1):
        for i,v in enumerate(C[S]):
            if not v:continue
            b=T[i]&~S
            while b:
                c=b&-b;b^=c;C[S|c][c.bit_length()-1]+=v
    H=[1]+[sum(L)for L in C[1:]]
    return H[-1],sum(v*H[f^S] for S,v in enumerate(H))

def comp(n,prev=0):
    if n==0:
        yield ()
        return
    for m in range(1,n+1):
        if m==prev==1:continue
        for rest in comp(n-m,m):
            yield (m,)+rest

def blocks(K):
    z=0; G=[]
    for k in K:
        G.append(tuple(range(z,z+k)));z+=k
    return G

def second(K):
    C=blocks(K)
    E=[e for A,B in zip(C,C[1:])for e in ((A[-1],B[0]),(A[0],B[-1]))]
    return flip(forward(sum(K)),E)

def three(l,m,r):
    A=1<<(l+r)
    return (1+ (1<<(m+1))*((1<<l)+(1<<r))+ 2*A*3**m -
             (l==0)*(1<<r)-(r==0)*(1<<l))

def prime(p):return p>1 and (p==2 or p%2 and all(p%d for d in range(3,isqrt(p)+1,2)))

def gate():
    H=[1,1,1]
    for n in range(3,15):H.append(sum(H[-3:]))
    b=[0]*15;b[0:2]=[1,1]
    for n in range(2,15):b[n]=sum(b[k] for k in (n-2,n-4,n-6) if k>=0)+(n in (4,5))
    for n in range(3,15):
        S=list(comp(n))
        need(len(S)==H[n] and sum(k==k[::-1] for k in S)==b[n])
        if n<=9:need(counts(zig(n))[0]==H[n])
        for K in S:
            if n>=10 and len(K)!=2:continue
            W=second(K);P=sum(blocks(K)[::-1],())
            need(rel(W,P)==frame_line(zig(n),P))
            if len(K)==3 and min(K)>=2:
                a,m,c=K;need(strong(W) and counts(W)[0]==three(a-2,m,c-2))
    for K,h,p in [((3,2,6),721,103),((6,8,6),3375617,482231),((5,14,8),4900119553,700017079)]:
        a,m,c=K
        need(h==three(a-2,m,c-2)==7*p and prime(p))
    for k in range(1,22):
        need(three(3*k-2,4*k,5*k-2)%7==(5+5*pow(2,k,7))%7!=0)

def tail(T,m,right=False):
    if not m:return T
    k=len(T)
    if not right:
        mask=((1<<m)-1)<<k
        return tuple(r|mask for r in T)+tuple(x<<k for x in forward(m))
    return tuple(r|((i==k-1)<<k) for i,r in enumerate(T))+tuple(
      (r<<k)|(((1<<k)-1)^((i==0)<<(k-1))) for i,r in enumerate(zig(m)))

def obstruction():
    E=[(i,j) for i in range(5) for j in range(i+2,5)]
    T=[merged(flip(forward(5),[e for i,e in enumerate(E) if a>>i&1])) for a in range(64)]
    lines={frozenset((u,v)) for a,u in enumerate(T) for v in (T[a^63],) if u!=v}
    need(any(all(frozenset(pair) in lines for pair in combinations(k,2))
             for k in combinations(set(T),5)))

    A=((10,20,33,38,13,19),(6,12,24,33,11,23),(6,12,40,17,7,27))
    B=((6,12,56,49,35,3),(6,40,42,49,7,17),(6,44,56,49,35,1))
    P=[[(3,5,1,4,2,0),(3,1,2,5,4,0),(3,2,5,0,1,4)],
       [(5,1,2,4,3,0),(5,0,2,4,1,3),(5,1,3,0,2,4)],
       [(5,4,2,3,0,1),(5,4,0,1,2,3),(5,0,1,2,3,4)]]
    need(all(strong(T) for T in A+B) and len(set(merged(T) for T in A+B))==6)
    need(all(canon(T)!=canon(zig(6)) and canon(converse(T))!=canon(zig(6)) for T in B))
    for i,j in product(range(3),repeat=2):
        need(mark(frame_line(A[i],P[i][j]))==mark(B[j]))
        for m in (1,2,6,7,10):
            F=P[i][j]+tuple(range(6,6+m))
            C=frame_line(tail(A[i],m),F)
            need(merged(C)==merged(tail(B[j],m,True)))
    for m in range(1,7):
        need(len({merged(tail(T,m))for T in A}|{merged(tail(T,m,True))for T in B})==6)

def gluing_examples():
    TT=((0,), (2,0), (6,4,0))
    data=(
     ((10,20,33,38,77,211,143,31),(258,4,305,327,395,347,407,47,130),(319,5080),(679,14840),9693389),
     ((32,17,11,19,69,94,143,63),(0,33,19,23,355,205,143,31,239),(217,3760),(245,8512),4199473),
     (TT[0],(512,525,537,593,67,287,935,831,543,48),(1,2),(1547,47936),7069),
     (TT[1],(2,260,25,67,331,287,167,63,205),(1,4),(721,16184),6833),
     ((2,4,257,263,271,31,319,127,163),(6,96,554,387,15,345,157,311,87,507),(119,4996),(4879,107296),71521273),
     (TT[0],(16,5,521,771,270,863,927,575,647,19),(1,2),(665,26096),3823),
     (TT[1],(4,289,66,7,271,349,411,319,13),(1,4),(371,10024),4243),
     (TT[2],(0,29,57,113,225,195,135,15),(1,8),(175,3976),3851))
    for A,B,a,b,p in data:
        need(counts(A)==a and counts(B)==b)
        u,R=a;v,S=b
        need((R-u)*(S-v)+2*u*v==7*p and prime(p))
        print('prime quotient',p,'shifted product',(R-u)*(S-v),'gaps',(R-u,S-v))
    need(69**2==4761 and 119**2==14161 and 4700**2+8211**2==9461**2 and gcd(4700,8211)==1)
    need(3543==3*1181 and 8267==7*1181 and prime(1181) and (217+8267)%3==0 and 217%3!=0)

if __name__=='__main__':
    gate();obstruction();gluing_examples()
    print('all shell laws, pointed K3,3 lifts and exact prime reads: PASS')
