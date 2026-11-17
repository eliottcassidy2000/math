from random import Random
from functools import cache

def ck(ok,*args):
    if not ok:raise ArithmeticError(args)

def profile(a):
    n=len(a);N=1<<n
    D=[[0]*n for _ in range(N)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,N):
        for j,w in enumerate(D[s]):
            if not w:continue
            b=a[j]&~s
            while b:
                v=b&-b;b-=v;D[s|v][v.bit_length()-1]+=w
    H=[1]+[sum(r) for r in D[1:]]
    return H[-1],sum(H[s]*H[(N-1)^s] for s in range(N))

def add(A,x,y,rev=False):
    n=len(A)
    return tuple(r|(1-(x>>i&1))<<n|(1-(y>>i&1))<<(n+1) for i,r in enumerate(A))+(x|((not rev)<<(n+1)),y|(rev<<n))

def T(m):
    return tuple((1<<m)-(1<<(i+1)) for i in range(m))

@cache
def K(q,m,i,j):
    if i==j:return q**(m-1)
    if i>j:return 0
    g=j-i-1
    return q**(i+m-1-j)*((q+1)**g-(q-1)**g)//2

def D(q,m,x,y):
    return sum(K(q,m,i,j)*((x>>i&1)*(1-(y>>j&1))-(y>>i&1)*(1-(x>>j&1)))
               for i in range(m) for j in range(i,m))

def Dmod(q,m,x,y,a):
    z=0
    for i in range(m):
        if q==2 and i>=a:break
        js=range(i,min(m,i+a+2)) if q==1 else range(max(i,m-a),m)
        for j in js:
            c=(x>>i&1)*(1-(y>>j&1))-(y>>i&1)*(1-(x>>j&1))
            if c:z+=K(q,m,i,j)*c
    return z%(1<<a)

def spine(a,H,P,m=None):
    mod=1<<a
    ck(a>=2 and P%2==0)
    if m is None:m=4*((a+2)*(a+4)+1)
    ck(m%2==0 and m>=4*((a+2)*(a+4)))
    X=Y=(((1<<(a-2))-1)<<(m-a))
    residue=P%mod
    for i in range(0,a-1,2):
        v=[sum(K(2,m,j,b) for b in (m-2,m-1)) for j in (i,i+1)]
        g=1<<(i+1)
        ck(all(w%g==0 and (w//g)%2 for w in v),a,i,v)
        mask=min(mod,g*4)-1
        d=next((d,e) for d in (-1,0,1) for e in (-1,0,1)
            if not (residue-d*v[0]-e*v[1])&mask)
        residue-=d[0]*v[0]+d[1]*v[1]
        for j,e in zip((i,i+1),d):
            if e==1:X|=1<<j
            if e==-1:Y|=1<<j
    ck(Dmod(2,m,X,Y,a)==P%mod,a,X,Y)
    addresses=[]
    for r in range(a):
        i=2*a+2+(2*a+4)*r; j=i+r+2
        X|=1<<j;Y|=1<<j
        addresses.append(i)
    mask=(Dmod(1,m,X,Y,a)-H)%mod
    for r,i in enumerate(addresses):
        if mask>>r&1:X|=1<<i
    ck(Dmod(1,m,X,Y,a)==H%mod,a,X,Y)
    ck(Dmod(2,m,X,Y,a)==P%mod,a,X,Y)
    return m,X,Y

def main():
    rng=Random(922)
    for m in range(1,10):
        for _ in range(15):
            x,y=(rng.randrange(1<<m) for _ in range(2))
            h,p=profile(add(T(m),x,y));a,b=profile(add(T(m),x,y,True))
            ck((a-h,b-p)==(2*D(1,m,x,y),4*D(2,m,x,y)),m,x,y)
            if m>=3:
                for t in range(2,m):
                    for q in (1,2):
                        ck(Dmod(q,m,x,y,t)==D(q,m,x,y)%(1<<t),q,m,x,y,t)
    for a in range(2,12):
        for H,P in [(5,-2)]+[(rng.randrange(1<<a),2*rng.randrange(1<<a)) for _ in range(3)]:
            spine(a,H,P)
    print("135 literal transitive-remainder switches; 40 simultaneous two-adic compilers." )

if __name__=="__main__":main()
