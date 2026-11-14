"""Replay three-slot fallback with explicit Hamilton subset DP."""
def paths(rows):
    n=len(rows); F=(1<<n)-1
    pred=[sum(((r>>i)&1)<<j for j,r in enumerate(rows)) for i in range(n)]
    dp=[[0]*n for _ in range(F+1)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,F+1):
        z=s
        while z:
            b=z&-z; z^=b; i=b.bit_length()-1
            w=pred[i] & (s^b)
            while w:
                c=w&-w; w^=c; j=c.bit_length()-1
                dp[s][i]+=dp[s^b][j]
    return [1]+list(map(sum,dp[1:]))
def conv(a,b):
    c=[]
    for s in range(len(a)):
        v=0; t=s
        while True:
            v+=a[t]*b[s^t]
            if not t: break
            t=(t-1)&s
        c.append(v)
    return c
def signature(r):
    h=paths(r); f=conv(h,h); q=conv(h,f)
    return h[-1],f[-1],q[-1]
def tt(n):
    return tuple((1<<n)-(2<<i) for i in range(n))
def tri(*blocks):
    ret=[]; offset=0; n=sum(map(len,blocks))
    starts=(0,len(blocks[0]),len(blocks[0])+len(blocks[1]))
    for i,r in enumerate(blocks):
        j=(i+1)%3
        bits=((1<<len(blocks[j]))-1)<<starts[j]
        ret.extend((v<<offset)|bits for v in r)
        offset+=len(r)
    return tuple(ret)
if __name__=='__main__':
    case=[
        (tt(2),(1,4,9),45),
        (tt(3),(1,8,27),137),
        (tt(4),(1,16,81),393),
        ((44,49,50,38,41,192,31,95),(345,5640,44505),202425),
        ((92,201,42,112,6,19,52,125),(345,5640,45225),205305)]
    for A,T,H in case:
        h,f,g=signature(A)
        assert (h,f,g)==T
        assert 4*g+5*f-11*h==paths(tri(A,tt(2),tt(2)))[-1]==H
        for k in range(1,5):
            assert paths(tri(A,tt(k),tt(1)))[-1]==(2**k-1)*f+(3-2**k)*h
        print(len(A),T,H)
