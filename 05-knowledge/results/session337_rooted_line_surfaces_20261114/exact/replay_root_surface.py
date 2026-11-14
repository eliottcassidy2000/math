from itertools import product

def paths(n,flip):
    adj=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    for i,j in flip:adj[i]^=1<<j;adj[j]^=1<<i
    D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s,row in enumerate(D):
        for i,c in enumerate(row):
            if c:
                moves=adj[i]&~s
                while moves:
                    b=moves&-moves;moves^=b;D[s|b][b.bit_length()-1]+=c
    return sum(D[-1])

def marked(U):
    return set(tuple(sorted((i,j))) for i,j in zip(U,U[1:]))

def board(a,q,c):
    p,r=a+2,c+2;n=p+q+r
    Q=list(range(p+q,n))+list(range(p,p+q))+list(range(p))
    assert all(j==i+1 or i-j>=2 for i,j in zip(Q,Q[1:]))
    direct={(p-1,p),(0,p+q-1),(p+q-1,p+q),(p,n-1)}
    assert marked(list(range(n)))^marked(Q)==direct
    assert (Q==[n-1-j for j in reversed(Q)])==(a==c)
    def arc(i,j):return (i<j)^((min(i,j),max(i,j)) in direct)
    spine=[i for i in range(n) if i not in (p,p+q-1)]
    assert all(arc(u,v) for u in spine for v in spine if u<v)
    assert arc(p,p+q-1)
    w=''.join(str(2*arc(p,v)+arc(p+q-1,v)) for v in spine)
    assert w=='1'+'0'*a+'2'*q+'3'*c+'1'
    return paths(n,direct)

def H(a,q,c): return 1+2**(q+1)*(2**a+2**c)+2*3**q*2**(a+c)-(a==0)*2**c-(c==0)*2**a
for a,q,c in product(range(4),range(2,5),range(4)):
    assert board(a,q,c)==H(a,q,c),(a,q,c)
assert H(1,3,3)==1025 and H(2,3,2)==993
print('All 48 independent coefficient / literal line / word checks passed')

from itertools import permutations
def comps(n):
    if n==0:
        yield ()
    else:
        for a in range(1,n+1):
            for s in comps(n-a):
                if a!=1 or not s or s[0]!=1:
                    yield (a,)+s
def lit_paths(t):
    n=len(t)
    D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s,row in enumerate(D):
        for i,c in enumerate(row):
            if c:
                moves=t[i]&~s
                while moves:
                    b=moves&-moves;moves^=b;D[s|b][b.bit_length()-1]+=c
    return sum(D[-1])
for n in range(3,9):
    def allowedQ(P):return all(j==i+1 or i-j>=2 for i,j in zip(P,P[1:]))
    literal={P for P in permutations(range(n)) if allowedQ(P)}
    model=set()
    for k in comps(n):
        borders=[0]
        for s in k:borders.append(borders[-1]+s)
        Q=tuple(u for i in reversed(range(len(k))) for u in range(borders[i],borders[i+1]))
        model.add(Q)
    assert literal==model
print('Composition language verified against literal partner paths through order eight')

def apex(w):
    m=len(w)
    M=[sum(1<<j for j in range(i+1,m))for i in range(m)]+[0]
    for i,d in enumerate(w):
        if d=='1':M[-1]|=1<<i
        else:M[i]|=1<<m
    return M
for k in range(1,6):
    seen=set()
    for v in product('01',repeat=k-1):
        U='1'+''.join(v);A=int(U,2);t=k+1
        h=lit_paths(apex(U+'0'*t))
        assert h==lit_paths(apex(U))+(2**t-1)*A+int(U[-1])
        seen.add(h)
    assert len(seen)==2**(k-1)
print('Injective strong one-hub padding checked for k=1..5')
