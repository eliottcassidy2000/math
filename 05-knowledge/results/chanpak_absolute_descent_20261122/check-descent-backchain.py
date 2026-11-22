from itertools import product

def transfer(widths):
    D,U,V=0,0,1
    for m in widths:
        if m:
            D,U,V=(1<<m)*(D+2*U+V),0,D+U+V
        else:
            D,U,V=D+2*U,V,D+U
    return D+U+V,V

def rows(widths):
    boundaries=[0]
    for m in widths:boundaries.append(boundaries[-1]+m+1)
    n=boundaries[-1]+1
    r=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    for a,b in zip(boundaries,boundaries[1:]):
        r[a]^=1<<b
        r[b]^=1<<a
    return r

def DP(r):
    n=len(r); M=1<<n
    dp=[[0]*n for _ in range(M)]
    for i in range(n):dp[1<<i][i]=1
    for S in range(1,M):
        for i,h in enumerate(dp[S]):
            e=r[i]&~S
            while h and e:
                b=e&-e
                dp[S|b][b.bit_length()-1]+=h
                e-=b
    return sum(dp[-1])

witnesses=0
for k in range(6):
    for w in product(range(3),repeat=k):
        r=rows(w)
        if len(r)<=11:
            assert transfer(w)[0]==DP(r),w
            witnesses+=1
for x,y in product(range(1,7),repeat=2):
    assert transfer((x,0,y))[0]==2*((1<<x)+1)*((1<<y)+1)-1
for w in product(range(3),repeat=6):
    if not w[-1]:continue
    F,G=transfer(w)
    for m in range(1,9):
        f,g=transfer(w+(0,m))
        d=f-g;v=(d&-d).bit_length()-1
        old=d>>v
        assert (v-1,old,2*old-g)==(m,F,G)
print('backchain transfer, coupled inverse and exact tournament DP PASS;',witnesses,'tournaments')

def census(X):
    from array import array
    A=set();f=array("I",[0])*(2*X+1)
    for y in range(len(f)):
        if y:
            u=y//2
            f[y]=f[u]+(u+(u&1) if not y&1 else 0)
        if f[y]<X:A.add(f[y]+1)
    stack=[(0,0,1)];B=set()
    while stack:
        D,U,V=stack.pop()
        H=D+U+V
        if H>X:continue
        B.add(H)
        if D or U:stack.append((D+2*U,V,D+U))
        m=2*(D+2*U+V)
        while m+H<=X:
            stack.append((m,0,H))
            m*=2
    return A,B

A,B=census(100001)
assert (len(A),len(B),len(B-A))==(25117,13263,6437)
p=(2,0,0,2,0,0,0,2,0,1,0,1)
assert len(rows(p))==21 and transfer(p)[0]==99991 and 99991 in B-A
print('absolute cutoff census PASS')
