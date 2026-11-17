def need(condition):
    if not condition:
        raise AssertionError('failed literal seam replay')

A=(4,5,72,19,39,79,155,63)
T35=(110,108,120,48,99,64,8)
T49=(2,20,17,7,105,15,47)
C3=(2,4,1)

def table(a):
    n=len(a); f=(1<<n)-1
    need(all(0<=r<=f and not(r>>i&1) for i,r in enumerate(a)))
    need(all((r>>j&1)+(a[j]>>i&1)==1 for i,r in enumerate(a) for j in range(i)))
    d=[[0]*n for _ in range(f+1)];h=[1]+[0]*f
    for i in range(n):
        d[1<<i][i]=1
    for s in range(1,f+1):
        for i,r in enumerate(a):
            v=d[s][i]
            h[s]+=v
            if not v:
                continue
            z=r&~s
            while z:
                b=z&-z;z^=b
                d[s|b][b.bit_length()-1]+=v
    return d,h

def inflate(m):
    n=7+m; r=(1<<n)-(1<<8)
    out=[x|((x&1)*r) for x in A]
    out[0]|=r
    out += [A[0]|((1<<n)-(1<<(i+1))) for i in range(8,n)]
    return out

def join(a,v):
    n=len(a); f=(1<<n)-1
    out=list(a)+[f|sum(1<<(n+j) for j in range(i)) for i in range(3)]
    out[v]|=1<<(n+2)
    out[n+2]^=1<<v
    return out

def strong(a):
    f=(1<<len(a))-1
    for i in range(len(a)):
        s=1<<i
        for _ in a:
            z=s
            for j,r in enumerate(a):
                if s>>j&1:
                    z|=r
            s=z
        if s!=f:
            return False
    return True

def ordinal(a,b):
    n=len(a)
    return tuple(x|(((1<<len(b))-1)<<n) for x in a)+tuple(x<<n for x in b)

def apex(w):
    n=len(w)
    return tuple(sum(1<<j for j in range(i+1,n))+((b=='0')<<n)
                 for i,b in enumerate(w))+(sum((b=='1')<<i for i,b in enumerate(w)),)

def run():
    for m in (1,2,3,4):
        a=inflate(m); d,h=table(a)
        need(d[-1][5]>0)
        b=sum(d[s][5]*h[(len(h)-1)^s] for s in range(1,len(h)) if s>>5&1)
        need((h[-1],b)==(64*(1<<m)+1,96*(1<<m)-2))
        out=join(a,5); N=448*(1<<m)-7
        need(strong(out) and table(out)[1][-1]==N)
        print(m, h[-1],b,N)
    for k,t in ((3,T49),(4,ordinal(C3,T35)),(5,apex('100110010')),
                (6,ordinal(T49,ordinal(C3,C3)))):
        need(table(t)[1][-1]==7*((1<<k)-1))
    print('Literal positive family and all four initial counts: PASS')
if __name__=='__main__':
    run()
