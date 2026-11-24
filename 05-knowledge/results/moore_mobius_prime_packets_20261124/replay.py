from collections import Counter
from itertools import combinations
from math import gcd, prod, comb

def graph(d):
    if d==2:
        return [{(i-1)%5,(i+1)%5} for i in range(5)]
    if d==3:
        pts=list(map(set,combinations(range(5),2)))
        return [{j for j,b in enumerate(pts) if not a&b} for a in pts]
    G=[set() for i in range(50)]
    def edge(u,v):
        G[u].add(v);G[v].add(u)
    for i in range(5):
        for x in range(5):
            edge(5*i+x,5*i+(x+1)%5)
            edge(25+5*i+x,25+5*i+(x+2)%5)
            for j in range(5):edge(5*i+x,25+5*j+(x+i*j)%5)
    return G

def components(G,S):
    left=set(S)
    out=[]
    while left:
        v=min(left);left.remove(v);todo=[v];part={v}
        while todo:
            for w in G[todo.pop()]&left:
                left.remove(w);todo.append(w);part.add(w)
        out.append(part)
    return out

def pentagons(G):
    def grow(p):
        if len(p)==5:
            if p[0] in G[p[-1]] and p[1]<p[-1]:yield p
            return
        for v in G[p[-1]]:
            if v>p[0] and v not in p:yield from grow(p+(v,))
    for u in range(len(G)):yield from grow((u,))

def chart(d):
    G=graph(d); n=len(G); T=Counter(); F=Counter()
    assert all(len(g)==d for g in G)
    assert all(len(G[i]&G[j])==int(j not in G[i]) for i,j in combinations(range(n),2))
    td=Counter();pd=Counter()
    for r in range(n):
        parent=sorted(G[r]); L=[G[i]-{r} for i in parent]
        for i,j,k in combinations(range(d),3):
            parts=components(G,L[i]|L[j]|L[k])
            z=tuple(sorted(len(S)//3 for S in parts))
            assert all(len(S)%3==0 for S in parts)
            T[z]+=1
            for a,b in combinations((parent[i],parent[j],parent[k]),2):
                td[r,*sorted((a,b))]+=z.count(2)
    for C in pentagons(G):
        X=[G[c]-set(C) for c in C]; Xall=set.union(*X)
        parts=components(G,Xall); z=tuple(sorted(len(p)//5 for p in parts))
        assert all(len(p)%5==0 for p in parts)
        F[z]+=1
        for i,r in enumerate(C):
            pd[r,*sorted((C[(i-1)%5],C[(i+1)%5]))]+=z.count(1)
        Y=set(range(n))-Xall-set(C)
        if not all(l==1 for l in z):continue
        tag={x:j for j,p in enumerate(parts) for x in p};t=len(parts)
        word={y:tuple(tag[next(iter(G[y]&A))] for A in X) for y in Y}
        assert all(len(set(w))==5 for w in word.values())
        for i,j in combinations(range(5),2):
            assert Counter((w[i],w[j]) for w in word.values())=={(a,b):1 for a in range(t) for b in range(t) if a!=b}
        for y,w in word.items():
            for i in range(5):
                assert Counter(word[z][i] for z in G[y]&Y)=={a:1 for a in range(t) if a not in (w[i],w[(i-2)%5],w[(i+2)%5])}
    assert all(pd[k]==2*td[k] for k in set(td)|set(pd))
    print('degree',d,'triple',T,'pentagon',F,'root-branch overlap checked')

def f(l,L):return (1<<(L*l-2))+1

def g(v):
    if not v:return 1
    k=v//2
    return g(k) if v%2 else g(k)+k+k%2

def parts(n,j=2):
    if not n:yield ();return
    for l in range(j,n+1):
        for r in parts(n-l,l):yield (l,*r)

def primitives(L,lo,hi):
    # Removing all gcds certifies a primitive prime exists without factoring.
    old=[]
    for j in range(lo,hi+1):
        z=f(j,L)
        for a in old:
            while (b:=gcd(z,a))>1:z//=b
        assert z>1
        old.append(f(j,L))
    return len(old)

def H(z,L): return prod(f(l,L) for l in z)

def spectrum57():
    Z=sorted((H(z,3),z) for z in parts(56))
    assert len({h for h,z in Z})==len(Z)
    assert Z[0]==(17**28,(2,)*28)
    assert Z[1]==(129**2*17**25,(2,)*25+(3,3))
    V=Counter()
    for h,z in Z:
        assert h%16==1
        n17=sum(l%8==2 for l in z)
        A=h
        for _ in range(n17):
            assert A%17==0;A//=17
        assert A%17!=0
        if n17>=25:V[n17,tuple(l for l in z if l!=2)]+=1
        # A logarithm-free recovery of the count using the signless modulo 7 phase.
        C=len(z)
        p=1<<(168-2*C)
        assert p<h and h*16**C<=p*17**C and h%7==pow(3,C,7)
    print('partition deck size',len(Z),'minimum',Z[0][0],'next',Z[1][0],'top 17 bands',V)

if __name__=='__main__':
    for d in (2,3,7):chart(d)
    assert [g(a) for a in (8,16,5024,83922)]==[9,17,4913,59049]
    print('primitive packets',primitives(3,2,56),primitives(5,1,55))
    spectrum57()
