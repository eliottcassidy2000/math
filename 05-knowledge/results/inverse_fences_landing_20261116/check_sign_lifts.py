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

def rows(a):
    return [[bool(a[i]>>j&1) for j in range(len(a))] for i in range(len(a))]
def endports(a, y):
    from itertools import permutations
    n=len(a); S0=set(range(n)); P=S=E=0
    for m in range(1,1<<n):
        if not (m>>y&1): continue
        prefix=[i for i in range(n) if m>>i&1 and i!=y]
        tail=sorted(S0-set(prefix)-{y})
        def allowed(word):
            return all(a[i]>>j&1 for i,j in zip(word,word[1:]))
        E+=sum(allowed((y,)+w) for w in permutations(prefix)) if not tail else 0
        v=sum(allowed(w+(y,)) for w in permutations(prefix))*sum(allowed(w) for w in permutations(tail))
        P+=v
        S+=v*(-1 if m.bit_count()%2 else 1)
    return P,S,E
for right in [(6,8,2,5),(0,), (2,4,1), (2,4,9,3)]:
    # Reindex it to place the root at the last slot as expected by glue.
    n=len(right)
    order=list(range(1,n))+[0]
    right=tuple(sum(1<<j for j in range(n) if a>>order[j]&1) for a in [right[i] for i in order])
    R=rows(right)
    P,S,E=endports(right,n-1)
    if S!=-E: raise RuntimeError(('sign',right,P,S,E))
    for w in ['110','10100','100','110010']:
        L=tournament(w)
        p=1+int(w,2)
        v=count(L)*count(R)+p*P-E
        if count(glue(L,R))!=v: raise RuntimeError(('bridge',w,right,v))
print('Generalized bridge and signed port controls pass.')

def strong(a):
    n=len(a)
    for i in range(n):
        v=1<<i
        for _ in range(n):
            for j in range(n):
                if v>>j&1: v|=a[j]
        if v!=(1<<n)-1: return False
    return True
V=[(56,37,1,6,14,28),(44,41,2,4,15,28),(16,9,3,5,14),(50,40,3,5,14,28)]
pro=[(37,44,4),(23,31,5),(15,21,3),(43,53,7)]
for a,(h,p,e) in zip(V,pro):
    vv,sv,ev=endports(a,0)
    if not strong(a) or (count(rows(a)),vv,ev,sv)!=(h,p,e,-e):
        raise RuntimeError(('plus chain',a,h,p,e))
for s in [1,2]:
    for W,A in [('11',3),('101',5)]:
        counts=[]
        for a,dt in zip(V,(0,1,2,1)):
            n=len(a); order=list(range(1,n))+[0]
            b=tuple(sum(1<<j for j in range(n) if row>>order[j]&1) for row in [a[i] for i in order])
            counts.append(count(glue(tournament(W+'0'*(s+dt)),rows(b))))
        expected=[A*p*2**s+3 for p in (81,108,144,192)]
        if counts!=expected or any(x%6!=3 or 3*y!=4*x-3 for x,y in zip(counts,counts[1:])):
            raise RuntimeError(('plus chain realized',counts,expected))
print('Both literal strong three-step plus chains pass.')

V=[(92,81,34,6,12,27,60),(52,49,2,7,12,28),(44,17,2,6,13,30),(116,97,34,7,14,24,60),(22,16,2,7,12)]
pro=[(109,134,26),(37,44,8),(25,29,5),(131,157,27),(9,11,3)]
counts=[]
for a,(h,p,e),dt,w in zip(V,pro,(0,2,3,1,6),('101','101','101','101','11')):
    v,sv,ev=endports(a,0)
    if not strong(a) or (count(rows(a)),v,ev,sv)!=(h,p,e,-e):
        raise RuntimeError(('minus chain',a,h,p,e))
    n=len(a); order=list(range(1,n))+[0]
    b=tuple(sum(1<<j for j in range(n) if row>>order[j]&1) for row in [a[i] for i in order])
    counts.append(count(glue(tournament(w+'0'*(1+dt)),rows(b))))
if counts!=[2429,3239,4319,5759,7679] or any(x%6!=5 or 3*y!=4*x+1 for x,y in zip(counts,counts[1:])):
    raise RuntimeError(('minus chain realized',counts))
print('Literal strong four-step minus chain passes.')
