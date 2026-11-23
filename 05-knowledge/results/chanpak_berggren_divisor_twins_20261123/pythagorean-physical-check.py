from itertools import combinations,permutations
from functools import cache
@cache
def chart(n):
    E=tuple(combinations(range(n),2)); return E
def rows(n,m):
    A=[0]*n
    for i,(u,v) in enumerate(chart(n)):
        if m>>i&1:A[u]|=1<<v
        else:A[v]|=1<<u
    return tuple(A)
def h(A):
    n=len(A); d=[[0]*n for _ in range(1<<n)]
    for j in range(n):d[1<<j][j]=1
    for s in range(1,1<<n):
        for j,c in enumerate(d[s]):
            if not c:continue
            b=A[j]&~s
            while b:
                e=b&-b;b-=e
                d[s|e][e.bit_length()-1]+=c
    return sum(d[-1])
def strongly(A):
    n=len(A)
    for j in range(n):
        s=1<<j
        for _ in range(n):
            t=s
            for a in range(n):
                if s>>a&1:t|=A[a]
            if t==s:break
            s=t
        if s!=(1<<n)-1:return False
    return True
def marked_canonical(A,u,v):
    return min(sum(((A[P[a]]>>P[b])&1)<<j for j,(a,b) in enumerate(chart(len(A))))
               for W in permutations(set(range(len(A)))-{u,v}) for P in ((u,v)+W,(v,u)+W))
def tagged_mask(A,e,f):
    rest=sorted(set(range(6))-set(e)-set(f))
    assert len(rest)==2
    return min(sum(((A[P[a]]>>P[b])&1)<<j for j,(a,b) in enumerate(chart(6)))
               for X in (e,e[::-1]) for Y in (f,f[::-1])
               for W in (tuple(rest),tuple(rest[::-1])) for P in (X+Y+W,))

def flip(A,e):
    u,v=e; B=list(A)
    B[u]^=1<<v; B[v]^=1<<u
    return tuple(B)

def join(A,B):
    n,m=len(A),len(B)
    return tuple(a|(((1<<m)-1)<<n) for a in A)+tuple(b<<n for b in B)

def shell(B,flag):
    n=len(B);a,b,z=n,n+1,n+2
    A=list(B)+[0]*3
    def edge(i,j):A[i]|=1<<j
    edge(b,a);edge(z,b)
    if flag:edge(z,a)
    else:edge(a,z)
    for i in range(n):edge(a,i);edge(i,b);edge(i,z)
    return tuple(A)

def classify():
    n=6;E=chart(n)
    A=[rows(n,m) for m in range(1<<len(E))]
    H=[h(a) for a in A]
    N=0; names=set()
    for m,t in enumerate(H):
        if t!=5:continue
        for i,e in enumerate(E):
            if H[m^(1<<i)]!=9:continue
            for j,f in enumerate(E):
                if j==i or H[m^(1<<j)]!=19 or H[m^(1<<i)^(1<<j)]!=23:continue
                N+=1
                names.add(tagged_mask(A[m],e,f))
    assert N==1440 and names=={3655,3997}
    assert tagged_mask(tuple(a^63^(1<<i) for i,a in enumerate(A[3655])),(0,1),(2,3))==3997
    print("hole square: 1440 labeled, masks 3655 and 3997, paired by converse")
    a=A[4]; e=(3,4); f=(0,5)
    assert [h(p) for p in [a,flip(a,e),flip(a,f),flip(flip(a,e),f)]]==[5,9,19,23]
    assert len({v for v in e+f})==4
    def canon(B):
        return min(sum(((B[P[a]]>>P[b])&1)<<j for j,(a,b) in enumerate(chart(len(B))))
                   for P in permutations(range(len(B))))
    def sub(B,I):
        return tuple(sum((B[i]>>k&1)<<j for j,k in enumerate(I)) for i in I)
    for m in names:
        B=flip(A[m],(0,1))
        cores=[sub(B,I) for I in combinations(range(6),5) if strongly(sub(B,I))]
        assert len(cores)==1 and h(cores[0])==9 and canon(cores[0])==8
    for word,mask in [('OOOI',8),('OIOI',10)]:
        B=[14,12,8,0,0]
        for i,l in enumerate(word):
            if l=='O':B[4]|=1<<i
            else:B[i]|=16
        assert h(B)==9 and strongly(B) and canon(B)==mask
    print("both H=9 corners are exactly OOOI plus one singleton")

def strong_cubes():
    for b in range(1,5):
        T=tuple(((1<<b)-1)^((1<<(i+1))-1) for i in range(b))
        assert h(shell(T,0))==1+2*2**b and h(shell(T,1))==-1+3*2**b
    B=rows(6,81)
    for e,g in zip(((2,5),(2,4),(1,5)),(6,8,10)):
        B0,B1=B,flip(B,e)
        assert strongly(B0) and strongly(B1)
        assert (h(B0),h(B1))==(23,23+g)
        d=[]
        for fl in (0,1):
            t=[]
            for u,v in ((B0,B0),(B1,B0),(B0,B1),(B1,B1)):
                T=shell(join(u,v),fl)
                assert strongly(T)
                t.append(h(T))
            d.append(t[3]-t[2]-t[1]+t[0])
        assert 3*d[0]-2*d[1]==5*g*g
        print("strong cube:",e,d,3*d[0]-2*d[1],flush=True)

if __name__=="__main__":
    classify()
    strong_cubes()
