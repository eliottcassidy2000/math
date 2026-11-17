from functools import cache

def ck(a,b):
    if a!=b:raise ArithmeticError((a,b))
def full(n):
    return [((1<<(i-1))-1 if i>=2 else 0)|((1<<(i+1)) if i<n-1 else 0) for i in range(n)]
def edge(n):
    a=[(1<<n)-(1<<(i+1)) for i in range(n)]
    a[0]^=1<<(n-1);a[-1]=1
    return a
def flipped(a,i,j):
    a=list(a);a[i]^=1<<j;a[j]^=1<<i
    return a
def H(a):
    n=len(a);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for m,r in enumerate(D):
        for i,x in enumerate(r):
            if x:
                w=a[i]&~m
                while w:
                    z=w&-w;w-=z;D[m|z][z.bit_length()-1]+=x
    return sum(D[-1])
def paths(a):
    n=len(a);N=(1<<n)-1
    def walk(p,m):
        if m==N:
            yield p
        else:
            r=a[p[-1]]&~m
            while r:
                j=r&-r;r-=j
                yield from walk(p+(j.bit_length()-1,),m|j)
    for i in range(n):yield from walk((i,),1<<i)
@cache
def trib(n):
    return (1,1,1)[n] if n<3 else sum(trib(n-i) for i in (1,2,3))
@cache
def new_uses(m):
    return (int(m==0)-int(m==1))+sum(a*new_uses(m-i) for i,a in enumerate((1,3,3,2),1) if m>=i)
def verify():
    for n in range(3,17):
        a=edge(n);h=2**(n-2)+1
        b=1+(2**((n-2)//2) if n%2==0 else 0)
        seen=blue=closed=closed_blue=0
        for p in paths(a):
            seen+=1
            c=bool(a[p[-1]]>>p[0]&1)
            B=all(x+y==n-1 for x,y in zip(p,p[::-1]))
            blue+=B;closed+=c;closed_blue+=c and B
        ck((seen,blue,closed,closed_blue),(h,b,n,1+(n%2==0)))
        ck(H(full(n)),trib(n))
        p=trib(n)-(n-1)+new_uses(n-2)
        ck(H(flipped(full(n),0,n-1)),p)
        if n<=11:print(n,h,trib(n),p,b,n)
    print('Verified literal Hamilton paths, closure, mirror frames and partner counts through 16.')
if __name__=='__main__':verify()
