def tournament(w):
    n=len(w)+1
    return [[False if i==j else (i<j) if j<n-1 and i<n-1 else (bool(w[j]=='1') if i==n-1 else bool(w[i]=='0')) for j in range(n)] for i in range(n)]

def glue(U,V,reverse=True):
    a,b=len(U),len(V)
    G=[[False]*(a+b) for _ in range(a+b)]
    for i in range(a):
        for j in range(a):
            G[i][j]=U[i][j] if i!=j else False
        for j in range(b):
            G[i][a+j]=True
    for i in range(b):
        for j in range(b):
            G[a+i][a+j]=V[i][j] if i!=j else False
    if reverse:
        G[a-1][a+b-1]=False
        G[a+b-1][a-1]=True
    return G

def count(M):
    n=len(M)
    masks=[sum((1<<j) for j,b in enumerate(r) if i!=j and b) for i,r in enumerate(M)]
    for i in range(n):
        for j in range(i):
            if bool(masks[i]>>j&1)==bool(masks[j]>>i&1):
                raise ValueError('not a tournament')
    d=[[0]*n for i in range(1<<n)]
    for i in range(n):
        d[1<<i][i]=1
    for s in range(1,1<<n):
        for i,t in enumerate(d[s]):
            if t:
                w=masks[i]&~s
                while w:
                    b=w&-w
                    j=b.bit_length()-1
                    d[s|b][j]+=t
                    w^=b
    return sum(d[-1])

for s in range(2,6):
    u='1001'+'0'*s
    S=glue(tournament('10'),glue(tournament(u),tournament('100')),False)
    T=glue(tournament(u[:-1]),tournament('100000'))
    d,f=count(S),count(T)
    if (d,f)!=(324*2**s-27,432*2**s-37) or 3*f!=4*d-3 or d%6!=3:
        raise RuntimeError(('plus',s,d,f))
for w,a in [('11',3),('101',5)]:
    for s in range(1,7):
        u=w+'0'*s
        S=glue(tournament(u),tournament('1110'))
        T=glue(tournament(u+'00'),tournament('10'))
        d,f=count(S),count(T)
        if (d,f)!=(18*a*2**s-1,24*a*2**s-1) or 3*f!=4*d+1 or d%6!=5:
            raise RuntimeError(('minus',s,d,f))
print('Literal-graph controls for all three pencil pairs pass.')
