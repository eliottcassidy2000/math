from itertools import combinations,permutations,product

def tr(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
def flip(A,a,b):
    A=list(A);A[a]^=1<<b;A[b]^=1<<a;return tuple(A)
def rev(A):
    return tuple(((1<<len(A))-1)^r^(1<<i) for i,r in enumerate(A))
def isos(A,B):
    return any(all((A[i]>>j&1)==(B[P[i]]>>P[j]&1)
               for i,j in combinations(range(len(A)),2)) for P in permutations(range(len(A))))
def support(A):
    return [P for P in permutations(range(len(A))) if all(A[i]>>j&1 for i,j in zip(P,P[1:]))]
for m in range(1,9):
    E={tuple(range(m+2))}
    for S in range(1<<m):
        a=tuple(i+1 for i in range(m) if S>>i&1)
        b=tuple(i+1 for i in range(m) if not (S>>i&1))
        E.add(a+(m+1,0)+b)
    if m<=5:assert set(support(flip(tr(m+2),0,m+1)))==E
    D={}
    for w in E:
        for j in range(m+3):D.setdefault(w[:j],set()).add(w[j:])
    assert len({frozenset(v) for v in D.values()})+1==3*2**m+m+2

def atom(mask):
    A=tr(5)
    for k,(i,j) in enumerate(combinations(range(5),2)):
        if not (mask>>k&1):A=flip(A,i,j)
    return A
A,B=atom(8),atom(10);C=(2,4,1)
assert all(isos(X,rev(X)) for X in (A,B)) and not isos(A,B)
for X in (A,B):
    P=(0,4,3,2,1);assert all(X[i]>>j&1 for i,j in zip(P,P[1:]))
def join(parts):
    A=();s=0
    for B in parts:
        cross=((1<<(s+len(B)))-1)^((1<<s)-1)
        A=tuple(r|cross for r in A)+tuple(r<<s for r in B);s+=len(B)
    return A
def outer(U):
    m=len(U);return tuple(r|(1<<m) for r in U)+(1<<(m+1),(1<<m)-1)
def cond(A,S):
    return all((A[i]&S) in (0,S) for i in range(len(A)) if not (S>>i&1))
def modules(A,m):
    n=len(A)
    return [S for e in combinations(range(n),n-m)
            if cond(A,S:=((1<<n)-1)^sum(1<<i for i in e))]
def Z(A):
    n=len(A);F=(1<<n)-1;D=[[0]*n for _ in range(F+1)]
    for i in range(n):D[1<<i][i]=1
    for S in range(1,F):
        for i,v in enumerate(D[S]):
            if not v:continue
            e=A[i]&~S
            while e:
                b=e&-e;e-=b;D[S|b][b.bit_length()-1]+=v
    h=[1]+list(map(sum,D[1:]))
    return h[F],sum(h[S]*h[F^S] for S in range(F+1))
for K in (1,2):
    for w in product((A,B),repeat=K):
        U=join((C,)+w);T=outer(U);m=len(U)
        assert modules(T,m)==[(1<<m)-1]
        for i in range(m):
            for b in (m,m+1):assert not modules(flip(T,i,b),m)
        if K==1:
            assert Z(U)==(3*9,12*84)
            assert Z(T)[0]==3*9+12*84

def short(T):
    n=len(T);D={}
    for a,b,c in permutations(range(n),3):
        if T[a]>>b&1 and T[b]>>c&1:D.setdefault(tuple(sorted((a,b,c))),set()).add(b)
    return D
def decode(D,n):
    R={(0,1):1};Q=[(0,1)]
    E=set(combinations(range(n),2))
    for e in Q:
        for f in E-dict.keys(R):
            meet=set(e)&set(f)
            if not meet:continue
            v=meet.pop()
            mid=D[tuple(sorted(set(e)|set(f)))]
            s=int(v in mid)
            # 1 means first vertex points to second; outward sense from v changes
            # by one iff v is the second endpoint.
            R[f]=R[e]^(e[1]==v)^(f[1]==v)^s
            Q.append(f)
    assert len(R)==len(E)
    A=[0]*n
    for (i,j),v in R.items():
        if v:A[i]|=1<<j
        else:A[j]|=1<<i
    return tuple(A)
for n in range(3,6):
    es=list(combinations(range(n),2))
    for M in range(1<<len(es)):
        T=tr(n)
        for j,(a,b) in enumerate(es):
            if M>>j&1:T=flip(T,a,b)
        assert decode(short(T),n) in (T,rev(T))
print("DFA residual census, folded short-path parity, strong tape module controls, profiles pass")
