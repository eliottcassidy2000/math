"""Literal DP controls for the signed root wrap."""
from itertools import product

def p(y):
    s=0
    while y:
        y,r=divmod(y,2)
        if not r:s+=y+(y&1)
    return s

def paths(rows,root=None):
    n=len(rows);m=1<<n
    dp=[[0]*n for _ in range(m)]
    for i in range(n):
        if root is None or i==root:dp[1<<i][i]=1
    for s in range(1,m):
        for i,q in enumerate(dp[s]):
            z=rows[i]&~s
            while z and q:
                b=z&-z;dp[s|b][b.bit_length()-1]+=q;z^=b
    return dp

def read(rows,root):
    dp=paths(rows)
    s=(1<<len(rows))-1
    H=[sum(r) for r in dp];H[0]=1
    rootdp=paths(rows,root)
    c=sum(sum(r)*H[s^v] for v,r in enumerate(rootdp))
    assert sum(rootdp[-1])==1
    return H[-1],c,dp[-1][root]

def wrap(rows,root,length,y):
    n=len(rows)
    mask=(1<<(n+length))-(1<<n)
    out=[r|mask for r in rows]
    out.extend(mask&~((1<<(i+1))-1) for i in range(n,n+length))
    for i in range(n):
        if i!=root:out[i]|=1<<(n+length)
    new=1<<root
    for d in range(length):
        if y>>d&1:out[n+d]|=1<<(n+length)
        else:new|=1<<(n+d)
    out.append(new)
    return tuple(out),n+length

def transfer(state,length,y):
    h,c,e=state
    if not length:return h+c-e,h+c,h-e
    B=1<<length
    rev=int(f"{y:0{length}b}"[::-1],2)
    return (1+p(y))*h+(y+1)*c-(not y%2)*e,(B-rev)*h+B*c,int(y>=B//2)*h

S=[((0,),0)]
for n in range(3):
    T=[]
    for rows,r in S:
        H=read(rows,r)
        for ell in range(3):
            for y in range(1<<ell):
                z,q=wrap(rows,r,ell,y)
                assert read(z,q)==transfer(H,ell,y),(rows,r,ell,y,read(z,q),transfer(H,ell,y))
                if (ell,y) in {(0,0),(1,1),(2,1),(2,3)}:T.append((z,q))
    S=T
print("signed root-wrap recurrence controls PASS")
