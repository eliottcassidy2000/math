"""Literal checks for exact complementary-reader returns."""
from itertools import product


def g(x):
    ans=1
    while x:
        x,b=divmod(x,2)
        if not b:
            ans+=x+x%2
    return ans


def rev(x,n):
    return int(f"{x:0{n}b}"[::-1],2)


def masks(w):
    n=len(w)
    return tuple((((1<<n)-1)^((2<<i)-1))|((b=='0')<<n) for i,b in enumerate(w)) + (sum((b=='1')<<i for i,b in enumerate(w)),)


def H(rows):
    n=len(rows)
    pre=[sum(bool(v>>i&1)<<j for j,v in enumerate(rows)) for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        dp[1<<i][i]=1
    for S in range(1,1<<n):
        s=S
        while s:
            b=s&-s
            i=b.bit_length()-1
            s^=b
            prev=pre[i]&(S^b)
            while prev:
                q=prev&-prev
                dp[S][i]+=dp[S^b][q.bit_length()-1]
                prev^=q
    return sum(dp[-1])


def word(x,y):
    return f'{x:b}'+'0'+''.join('1' if b=='0' else '0' for b in f'{y:b}'[::-1])


def pair(k,m):
    n=(m+1)*(k+1)
    x=(2**(k-1)-1)*(2**(n-k-1)-1)//(2**(k+1)-1)
    y,z=2**(n-1)+2**k-1,2**(n-1)-2**k
    return n,x,y,z


if __name__=='__main__':
    for n in range(1,15):
        for x in range(2**n):
            y=2**n-1-x
            assert g(x)-g(y)==x-rev(x,n)+(x>>(n-1))-(x&1)
    for k,m in product(range(2,13),range(1,9)):
        n,x,y,z=pair(k,m)
        a,b=word(x,y),word(x,z)
        h=g(x)+2**n*x+2**(n-2)+2-2**(k-1)
        assert a[0]==b[0]=='1' and a[-1]==b[-1]=='0'
        assert g(int(a,2))==g(int(b,2))==h
        assert 1+int(a,2)==2**(n+1)*x+2**(n-k)-1
        assert 1+int(b,2)==2**n*x+2**(n-1)-2**(n-1-k)+1
    for k,m in [(2,1),(2,2),(3,1),(4,1)]:
        _,x,y,z=pair(k,m)
        a,b=word(x,y),word(x,z)
        A,B=masks(a),masks(b)
        h=g(int(a,2))
        assert H(A)==H(B)==h
        print(k,m,h,(1+int(a,2),1+int(b,2)),A,B)
    print('PASS: all reflection identities through width 14; 88 family pairs; four literal DAG-DP pairs')
