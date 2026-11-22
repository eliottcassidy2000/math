from itertools import product, combinations, permutations
from functools import lru_cache

def literal_weights(T):
    n=len(T)
    @lru_cache(None)
    def end(A,i):
        if A==(1<<i):return 1
        return sum(end(A^(1<<i),j) for j in range(n)
                   if j!=i and A>>j&1 and T[j]>>i&1)
    @lru_cache(None)
    def start(A,i):
        if A==(1<<i):return 1
        return sum(start(A^(1<<i),j) for j in range(n)
                   if j!=i and A>>j&1 and T[i]>>j&1)
    full=(1<<n)-1
    def C(u,v):
        return sum(end(A,u)*start(full^A,v) for A in range(1,full)
                   if A>>u&1 and not A>>v&1)
    W=[[0]*n for _ in T]
    for a,b in combinations(range(n),2):
        k=C(a,b)+C(b,a)
        assert not k%2
        W[a][b]=W[b][a]=k//2
        assert k>0
    return W

def matrix(T,W,x,i):
    return sum(w for j,w in enumerate(W[i]) if T[i]>>j&1)-sum(w for j,w in enumerate(W[i]) if x>>j&1)

def check(T,U,W,V,a,b):
    n=len(T); w=W[a][b]; assert w==V[a][b]
    js=[]
    det=lambda i,j: i[0]*j[1]-i[1]*j[0]
    for x in range(1<<n):
        J=det(*((matrix(T,W,x,i),matrix(U,V,x,i)) for i in (a,b)))
        js.append(J)
    z=T[a]|(1<<b)
    assert not (z>>a&1)
    assert js[z|(1<<a)]-js[z]==-w*w
    for p in (2,3,5,7,11,13):
        if w%p:
            assert 4*sum(j%p!=0 for j in js)>=1<<n

pairs=0
for n in range(2,6):
    E=list(combinations(range(n),2));Ts=[];Ws=[]
    for bits in range(1<<len(E)):
        T=[0]*n
        for k,(u,v) in enumerate(E):
            if bits>>k&1:u,v=v,u
            T[u]|=1<<v
        Ts.append(T);Ws.append(literal_weights(T))
    for i,T in enumerate(Ts):
        for k,(a,b) in enumerate(E):
            if i>>k&1:continue
            j=i^(1<<k)
            check(T,Ts[j],Ws[i],Ws[j],a,b)
            pairs+=1
print(pairs,"literal arc reversals: marked determinant identity and modular quarter-density pass")
