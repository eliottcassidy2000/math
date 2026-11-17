from itertools import product
from math import isqrt

def dp(rows,root=None):
    n=len(rows);N=1<<n
    d=[[0]*n for _ in range(N)]
    for i in range(n) if root is None else (root,):d[1<<i][i]=1
    for s in range(1,N):
        for i,v in enumerate(d[s]):
            if not v:continue
            r=rows[i]&~s
            while r:
                b=r&-r;r-=b;d[s|b][b.bit_length()-1]+=v
    return [1]+[sum(x) for x in d[1:]]

def profile(rows,root):
    h=dp(rows);p=dp(rows,root);N=len(h)-1
    return (h[N],h[N^(1<<root)],*(sum(t*h[N^s] for s,t in enumerate(p) if s and s.bit_count()%2==k) for k in (0,1)))

def core(w,pos):
    k=2*w[pos]+w[pos+1]
    b=w[:pos]+(1,)*k+(0,)*(3-k)+w[pos+2:]
    n=len(b)
    t=[sum(1<<j for j in range(i+1,n))|((not bit)<<n) for i,bit in enumerate(b)]
    t[pos]^=1<<(pos+2);t[pos+2]^=1<<pos
    return tuple(t+[sum(bit<<j for j,bit in enumerate(b))]),n

def value(w,pos):
    n=len(w);r=n-pos-2
    k=2*w[pos]+w[pos+1]
    v=sum(d<<(n-1-i) for i,d in enumerate(w))
    p=v>>(r+2)
    L=3*(p//2)+(p%2)
    z=int(bool(r) and not w[pos+2])
    eta=(L-z+1 if k==1 else -L+z-2 if k==2 else 0)
    def weight(t):return 1<<max(0,t-1)
    B=1+2*sum(weight(j-i-1) for j in range(n) if not w[j] for i in range(j) if w[i])
    h=3*B+2*eta
    R=3*(1+v)
    D=3*(1-w[-1]) if r else 3-k
    return h,3,(R-D)//2,(R+D)//2

def append(rows,root,z):
    a=list(rows)
    for i in range(z.bit_length()):
        n=len(a);b=z>>i&1
        a=[r+((i!=root or not b)<<n) for i,r in enumerate(a)]+[(1<<root) if b else 0]
    return a

def step(p,z):
    h,g,e,o=p
    for i in range(z.bit_length()):
        if z>>i&1:h,e,o=h+2*e,e+o-g,e+o
        else:e,o=e+o,e+o
    return h,g,e,o

def f(z):
    a=1
    while z:
        if not z&1:
            u=z//2
            a+=u+(u&1)
        z//=2
    return a

def target(n,A,J,base_delta):
    # n = 3(A(z+1)+f(z)+J)+base_delta
    if (n-base_delta)%3:return []
    y=(n-base_delta)//3-J
    if y<A+1:return []
    k=max(1,isqrt(y).bit_length()-1)
    q=1<<k
    ans=[]
    def reversal(r):
        b=0
        for i in range(k): b=(b<<1)|((r>>i)&1)
        return b
    for r in range(q):
        K=(A+1)*q-1-reversal(r)
        C=A*(r+1)+f(r)-1
        low=max(0,(y-C-2+K)//(K+1))
        high=(y-C-1)//K
        for b in range(low,high+1):
            z=q*b+r
            if A*(z+1)+f(z)==y:ans.append(z)
    return ans

stable_a=((26,12,9,16,6),(14,20,24,18,1),4,4,0,-2,0)
stable_b=((62,36,40,34,14,16),(42,60,57,32,9,16),4,6,-2,0,2)

if __name__=='__main__':
    tests=0
    for n in range(2,8):
        for w in product((0,1),repeat=n):
            for pos in range(n-1):
                r,j=core(w,pos)
                assert profile(r,j)==value(w,pos),(w,pos,profile(r,j),value(w,pos))
                tests+=1
    print('literal checks of the general one-triangle law',tests)
    for a,b,x,A,J,lo,hi in (stable_a,stable_b):
        pa,pb=profile(a,x),profile(b,x)
        assert (pa[0]+2,*pa[1:])==pb
        print(pa,pb)
        for z in range(64):
            p,q=[profile(append(t,x,z),x) for t in (a,b)]
            assert p==step(pa,z) and q==step(pb,z)
            assert (p[0],q[0])==tuple(3*(A*(z+1)+f(z)+J)+j for j in (lo,hi))
        for n in (23,67,69,127,343,511,1091,16643,6017,6019,10007,10009,30011,30013,1048509,1000000007,1000000009):
            for j,s in ((lo,a),(hi,b)):
                for z in target(n,A,J,j):
                    assert step(profile(s,x),z)[0]==n
                    print('prescribed',n,'profile A',A,'address',z)
