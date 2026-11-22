"""Exact joint modular counts for affine nonnegative address pairs."""
from functools import lru_cache
from itertools import product

def read(x):
    z=0
    while x:
        x,b=divmod(x,2)
        if not b:z+=x+(x&1)
    return z

def height(c,x,y):
    h,q,C,T,e,r,k=c
    return k+q*x*y+C*y+T*x+e*(y&1)+r*(x&1)+h*(read(x)+read(y))

def count(c,target,m,A,B,D,E,lo,hi):
    """Count j in [lo,hi] with H(A+B*j,D+E*j)=target mod m."""
    assert 0<=lo<=hi and min(A+B*lo,A+B*hi,D+E*lo,D+E*hi)>=0
    h,q,C,T,e,r,k=c
    # Choose enough bits for j AND both resulting addresses.
    s=max(1,hi.bit_length(),(A+B*lo).bit_length(),(A+B*hi).bit_length(),
          (D+E*lo).bit_length(),(D+E*hi).bit_length())
    p2=[pow(2,i,m) for i in range(s)]
    @lru_cache(None)
    def dp(i,a,d,L,U,x,y,W,V,l,v,F):
        if U<0:return 0
        if i==s:
            if a or d or L or U:return 0
            return (k+q*x*y+C*y+T*x+F-target)%m==0
        ans=0
        for b in (0,1):
            na,u=divmod(a+B*b,2)
            nd,w=divmod(d+E*b,2)
            nf=(F+h*(u*(W+l)+w*(V+v)) + (r*u+e*w if i==0 else 0))%m
            ans+=dp(i+1,na,nd,(L+b)//2,(U-b)//2,
                    (x+p2[i]*u)%m,(y+p2[i]*w)%m,
                    (2*W+1-u)%m,(2*V+1-w)%m,1-u,1-w,nf)
        return ans
    return dp(0,A,D,-lo,hi,0,0,0,0,0,0,0)

def main():
    import random
    rng=random.Random(23951)
    c=(2465,2293760,1344270,712944,1740,986,423385)
    for _ in range(150):
        lo=rng.randrange(10);hi=lo+rng.randrange(35)
        A=rng.randrange(100); B=rng.randrange(-12,13)
        D=rng.randrange(100); E=rng.randrange(-12,13)
        A+=max(0,-A-B*hi);D+=max(0,-D-E*hi)
        m=rng.choice((2,3,4,6,12,24,48))
        t=rng.randrange(m)
        assert count(c,t,m,A,B,D,E,lo,hi)==sum(height(c,A+B*j,D+E*j)%m==t for j in range(lo,hi+1))
    print("150 exact affine/brute-force controls passed")
    for m in (3,4,6,12,24,48):
        n=count(c,(1<<60)*((1<<61)-1)-1,m,10**15,10,3*10**15,-7,0,4095)
        print(m,n,"of",4096)
if __name__=="__main__":main()
