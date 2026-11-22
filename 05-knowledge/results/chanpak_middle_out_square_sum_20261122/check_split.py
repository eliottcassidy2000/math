from collections import Counter
from itertools import permutations, product, combinations
from math import comb
from random import Random
R=Random(53879)

def arc(T,u,v): return T[u]>>v&1
def tournament(n):
    T=[0]*n
    for u,v in combinations(range(n),2):
        if R.randrange(2): u,v=v,u
        T[u]|=1<<v
    return T

def pc(T):
    "Literal permutation/cut enumeration, including one-empty-slot covers."
    C=Counter()
    for p in permutations(range(len(T))):
        for s in range(len(T)+1):
            A,B=p[:s],p[s:]
            if all(arc(T,u,v) for a in (A,B) for u,v in zip(a,a[1:])):
                C[(A[-1] if A else -1,B[0] if B else -1)]+=1
    return C

def cut(T,W,x):
    return sum(W[u,v] for u in range(len(T)) for v in range(len(T))
               if arc(T,u,v) and x>>u&1 and not x>>v&1)

def spectra(T,U,pins=()):
    n=len(T); C,D=pc(T),pc(U)
    W={(u,v):(C[u,v]+C[v,u])//2 for u in range(n) for v in range(n) if u!=v}
    V={(u,v):(D[u,v]+D[v,u])//2 for u in range(n) for v in range(n) if u!=v}
    assert all((C[u,v]+C[v,u])%2==0 for u,v in W)
    free=[i for i in range(n) if i not in dict(pins)]
    pin=sum(b<<i for i,b in pins)
    def lift(x): return pin|sum(((x>>j)&1)<<i for j,i in enumerate(free))
    obs=lambda x:(cut(T,W,lift(x)),cut(U,V,lift(x)))
    target=Counter((a%2,b%2) for a,b in map(obs,range(1<<len(free))))
    predicted=Counter()
    d=len(free)
    for l in (1,2,3):
        evalf=lambda x:sum(z for i,z in enumerate(obs(x)) if l>>i&1)%2
        c=evalf(0)
        f=lambda x:evalf(x)^c
        bil=lambda x,y:f(x^y)^f(x)^f(y)
        basis=[1<<i for i in range(d)]
        summ=(-1)**c
        while True:
            pair=next(((e,g) for e,g in combinations(basis,2) if bil(e,g)),None)
            if pair is None:break
            e,g=pair
            summ*=2*(-1)**(f(e)*f(g))
            basis=[x^(e if bil(g,x) else 0)^(g if bil(e,x) else 0) for x in basis if x not in pair]
        summ*=0 if any(f(x) for x in basis) else 1<<len(basis)
        for a,b in product(range(2), repeat=2):
            predicted[a,b]+=summ*(-1)**((l&1)*a+(l>>1)*b)
    assert all((n0+2**d)//4==target[z] for z,n0 in predicted.items())
    for u,v in combinations(range(n),2):
        if arc(T,u,v)==arc(U,u,v):continue
        assert W[u,v]==V[u,v]
        a,b=free.index(u),free.index(v)
        for x in range(1<<d):
            if (x>>a|x>>b)&1:continue
            vals=[obs(x^y) for y in (0,1<<a,1<<b,1<<a|1<<b)]
            diffs=[s-r for r,s in vals]
            assert diffs[0]+diffs[-1]==diffs[1]+diffs[2]
    return W,V

def module_test(T,m,W):
    # Enumerate ordered distinct-outside skeletons with nonadjacent M runs.
    n=len(T)
    outside=tuple(range(m,n))
    O=Counter()
    for r in range(1,m+1):
        for U in set(permutations(outside+(-2,)*r)):
            for j in range(1,len(U)):
                A,B=U[:j],U[j:]
                ok=lambda p:all(x!=y and arc(T,0 if x==-2 else x,0 if y==-2 else y) for x,y in zip(p,p[1:]))
                if not ok(A) or not ok(B):continue
                if A[-1]==-2 and B[0]>=0: O[r,0,B[0]]+=1
                if A[-1]>=0 and B[0]==-2: O[r,1,A[-1]]+=1
    def end(r,i):
        return sum((-1)**(r-k)*comb(r-1,k-1)*k**i*(k-1)**(m-1-i) for k in range(1,r+1))
    for i in range(m):
        for v in outside:
            assert 2*W[i,v]==sum(end(r,i)*O[r,0,v]+end(r,m-1-i)*O[r,1,v] for r in range(1,m+1))

for n in range(4,8):
    for i in range(6):
        T=tournament(n);m=2+(i%min(3,n-3))
        for u,v in combinations(range(m),2):
            T[u]|=1<<v;T[v]&=~(1<<u)
        for v in range(m,n):
            for u in range(1,m):
                if arc(T,0,v):T[u]|=1<<v;T[v]&=~(1<<u)
                else:T[u]&=~(1<<v);T[v]|=1<<u
        U=T[:];u,v=m,m+1
        U[u]^=1<<v;U[v]^=1<<u
        W,V=spectra(T,U,[(0,i%2)] if i%3 else ())
        module_test(T,m,W);module_test(U,m,V)
print("24 physical seam pairs: pointed module identities, pinned quadratic parity and affine middle-mark census passed")

def symhost(n):
    T=[0]*n; handled=set(); inv=lambda i:n-1-i
    for a,b in combinations(range(n),2):
        if frozenset((a,b)) in handled:continue
        u,v=(a,b) if R.randrange(2) else (b,a)
        T[u]|=1<<v;T[inv(v)]|=1<<inv(u)
        handled.update((frozenset((a,b)),frozenset((inv(a),inv(b)))))
    return T

def split(T,U,W,V,d):
    n=len(T); pair=lambda a,b:(cut(T,W,a|b<<d),cut(U,V,a|b<<d))
    A={};B={}
    for a in range(1<<d):
        key=tuple(sum(M[i,j]*(a>>i&1) for i in range(d)) for M in (W,V) for j in range(d,n))
        A.setdefault(key,{})[pair(a,0)]=a
    for b in range(1<<(n-d)):
        key=tuple(sum(M[i,j]*(b>>(j-d)&1) for j in range(d,n)) for M in (W,V) for i in range(d))
        B.setdefault(key,{})[pair(0,b)]=b
    out=set()
    for D in A.values():
        for E in B.values():
            a0=next(iter(D.values()));b0=next(iter(E.values()))
            corr=tuple(pair(a0,0)[i]+pair(0,b0)[i]-pair(a0,b0)[i] for i in (0,1))
            out.update(tuple(u[i]+v[i]-corr[i] for i in (0,1)) for u in D for v in E)
    return out

for n in range(4,8):
    for i in range(3):
        T=symhost(n); U=T[:]; U[0]^=1<<(n-1); U[-1]^=1
        W,V=spectra(T,U)
        orbs=[tuple(sorted({i,n-1-i})) for i in range((n+1)//2)]
        for mask in range(1<<len(orbs)):
            x=sum(1<<v for i,S in enumerate(orbs) if mask>>i&1 for v in S)
            for M,T0 in ((W,T),(V,U)):
                total=0
                for i,j in combinations(range(len(orbs)),2):
                    a=sum(M[u,v] for u in orbs[i] for v in orbs[j])
                    assert not a%2
                    if (mask>>i^mask>>j)&1:total+=a//2
                assert total==cut(T0,M,x)
        assert split(T,U,W,V,n//2)=={(cut(T,W,x),cut(U,V,x)) for x in range(1<<n)}
print("12 physical common-skew seams: exact orbit folds and exact paired split value matching passed")
