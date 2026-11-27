"""Python 3 replay. Pass a local copy of SpenceMat.html as the first argument."""
from collections import Counter
from itertools import combinations, product
from math import comb
from pathlib import Path
import re
import sys

def mat(rows):
    n = len(rows)
    A = [[(2*int(rows[i][j])-1 if i != j else 0) for j in range(n)]+[-1]
         for i in range(n)]
    return A+[[1]*n+[0]]

def root(A, r):
    e = [a or 1 for a in A[r]]
    return [[a*e[i]*e[j] for j,a in enumerate(row)] for i,row in enumerate(A)]

def shell(A):
    n = len(A)
    counts = Counter(); parity = Counter()
    for r in range(n):
        B = root(A,r)
        for a,b,c in combinations((i for i in range(n) if i != r),3):
            K = B[a][c]*B[b][c]+B[a][b]*B[c][b]+B[b][a]*B[c][a]
            Q = sum(B[r][x]*B[a][x]*B[b][x]*B[c][x] for x in range(n))
            for k,D in enumerate((n-1+K+Q,n-1+K-Q)):
                counts[D] += 1
                if D==0: parity[1-k] += 1
    assert sum(counts.values()) == 2*n*comb(n-1,3)
    return counts,parity

def codes(A):
    n = len(A)
    h = [sum(1<<j for j,a in enumerate(row) if a<0) for row in A]
    universe = (1<<n)-1
    pair = []
    for flip in (False,True):
        C={0,universe}
        for i,word in enumerate(h):
            if flip: word ^= (1<<i) ^ (1<<(n-1))
            C |= {c^word for c in C}
        assert len(C) == 1<<(n//2)
        assert all(c.bit_count()%4==0 for c in C)
        pair.append(C)
    assert pair[0]&pair[1]=={0,universe}
    return pair

def maps(A,B):
    n = len(A)
    r = n-1
    X = root(A,r)
    allv = (1<<n)-1
    for s in range(n):
        Y=root(B,s)
        nb = [sum(1<<j for j,a in enumerate(row) if a==1) for row in Y]
        to = [-1]*n
        to[r]=s
        def search(D):
            if not D:
                yield tuple(to)
                return
            i = min(D,key=lambda i:D[i].bit_count())
            bits=D[i]
            while bits:
                bit=bits&-bits
                bits-=bit
                j=bit.bit_length()-1
                to[i]=j
                d={k:v & (nb[j] if X[i][k]==1 else (allv^nb[j]^bit))
                   for k,v in D.items() if k!=i}
                if all(d.values()):
                    yield from search(d)
        yield from search({i:allv^(1<<s) for i in range(r)})

def check(A):
    n=len(A)
    assert all(A[i][j]==-A[j][i] for i in range(n) for j in range(n))
    assert all(sum(a*b for a,b in zip(A[i],A[j]))==(n-1)*(i==j)
               for i in range(n) for j in range(n))
    C=codes(A)
    f,p=shell(A)
    W=[Counter(c.bit_count() for c in d) for d in C]
    assert all(p[k]==4*W[k][4] for k in (0,1))
    return f,C,W

def test_shell_samples(A):
    n=len(A)
    L=n*(n-2)
    s0=L*L
    s1=((n-2)*(n-4))**2+32*(n-2)
    s2=((n-4)*(n-6))**2+144*(n-4)+16*n
    for r in range(n):
        X=root(A,r)
        for a,b,c in list(combinations((i for i in range(n) if i!=r),3))[::137]:
            K=X[a][c]*X[b][c]+X[a][b]*X[c][b]+X[b][a]*X[c][a]
            Q=sum(X[r][j]*X[a][j]*X[b][j]*X[c][j] for j in range(n))
            for sink in (False,True):
                flips={a,b,c}|({r} if sink else set())
                sums=[sum(-v if j in flips else v for j,v in enumerate(row)) for row in X]
                S=sum((z*z-1)**2 for z in sums)
                assert S==s0-3*s1+3*s2-192*(n-1+K+(-Q if sink else Q))

def unrooted(A):
    n=len(A)
    f=Counter()
    for a,b,c,d in combinations(range(n),4):
        pf=A[a][b]*A[c][d]-A[a][c]*A[b][d]+A[a][d]*A[b][c]
        K=-3 if abs(pf)==3 else 1
        Q=sum(w*x*y*z for w,x,y,z in zip(A[a],A[b],A[c],A[d]))
        f[n-1+K+Q]+=4
        f[n-1+K-Q]+=4
    return f

def order28(html):
    rows=re.findall(r'^\s*([01]{27})\s*$',html,re.M)
    assert len(rows)==65*27
    mat28=[mat(rows[i:i+27]) for i in range(0,len(rows),27)]
    opp=lambda X:[[-a for a in row] for row in X]
    groups={}
    for i,A in enumerate(mat28):
        assert all(A[u][v]==-A[v][u] for u in range(28) for v in range(28))
        assert all(sum(a*b for a,b in zip(A[u],A[v]))==27*(u==v)
                   for u in range(28) for v in range(28))
        f=unrooted(A)
        groups.setdefault(tuple(f[d] for d in range(0,49,8)),[]).append(i)
    bins=[len({k[i] for k in groups}) for i in range(7)]
    print("order 28 distinct one-bin fingerprints:",bins)
    print("order 28 complete shell fingerprints:",len(groups))
    total=0
    collisions=[]
    for g in groups.values():
        assert all(next(maps(mat28[i],mat28[j]),None) is None for i,j in combinations(g,2))
        reps=[]
        for i in g:
            for j in reps:
                if any(next(maps(mat28[i],B),None) is not None
                       for B in (mat28[j],opp(mat28[j]))):
                    break
            else:
                reps.append(i)
        total+=len(reps)
        if len(reps)>1:
            collisions.append([i+1 for i in reps])
    print("order 28 classes up to converse:",total)
    print("order 28 third-shell collision reps:",collisions)
    assert collisions==[[44,48]]
    h = [fourth(mat28[i-1]) for i in collisions[0]]
    dif=[(S,h[0][S],h[1][S]) for S in sorted(h[0].keys()|h[1].keys())
         if h[0][S]!=h[1][S]]
    print("order 28 fourth-shell differing frequencies (S, block44, block48):",dif)
    assert dif
    S,a,b=dif[0]
    n=28
    import math
    N=(math.prod(range(n-5,n+1))+n*(3*(n-1)**2-2*(n-1))-S-n*(2*n-3))//32
    print("order 28 tie breaker (paths, block44, block48):",N,a,b)

def fourth(A):
    n=len(A)
    rowmask=[sum(1<<j for j,x in enumerate(row) if x<0) for row in A]
    diag=[1<<j for j in range(n)]
    L=n*(n-2)
    powers=[((n-1-2*d)**2-1)**2 for d in range(n)]
    h=Counter()
    for r in range(n):
        for t in combinations((i for i in range(n) if i!=r),4):
            mask=rowmask[r]^sum(1<<i for i in t)
            for e in (mask,mask^diag[r]):
                S=sum(powers[((e^h)&~d).bit_count()] for h,d in zip(rowmask,diag))
                assert 135744<S<193600
                h[S]+=1
    assert sum(h.values())==2*n*comb(n-1,4)
    return h

def count_sample(A):
    n=len(A)
    X=root(A,n-1)
    for b in combinations(range(n-1),3):
        e=[-1 if j in b else 1 for j in range(n)]
        T=sum(sum(v*f for v,f in zip(row,e))**4 for row in X)
        if T==94008:
            break
    else:
        raise AssertionError("no score 94008 in this root")
    nb=[sum(1<<j for j in range(n) if X[i][j]*e[i]*e[j]==1)
        for i in range(n)]
    def walk(i,left,used):
        b=nb[i]&~used
        if left==1:
            return b.bit_count()
        s=0
        while b:
            x=b&-b
            b-=x
            s+=walk(x.bit_length()-1,left-1,used|x)
        return s
    h=sum(walk(i,5,1<<i) for i in range(n))
    assert h==3026628
    print("checked actual path count",h,"at flipped triple",b)

def main(path):
    html=Path(path).read_text()
    rows=re.findall(r'^\s*([01]{23})\s*$',html,re.M)
    assert len(rows)==16*23
    A=[mat(rows[i:i+23]) for i in range(0,len(rows),23)]
    expected=[800,684,864,684,864,800,544,696,696,528,528,504,504,1320,1320,0]
    groups={}
    seen=[]
    print("block | frequencies (D=0,8,16,24,32,40) | tetrads (source,sink)")
    for i,X in enumerate(A,1):
        f,C,W=check(X)
        test_shell_samples(X)
        assert f[8]==expected[i-1]
        print(i,[f[d] for d in range(0,41,8)],[w[4] for w in W])
        groups.setdefault(f[8],[]).append(i)
        if any(w[8]==759 and not w[4] for w in W): seen.append(i)
    assert seen==[14,15,16]
    for group in groups.values():
        i,j=group[0]-1,group[-1]-1
        assert next(maps(A[i],[[-a for a in row] for row in A[j]]),None) is not None

    X=A[14]
    C=codes(X)
    tetrads={c for c in C[1] if c.bit_count()==4}
    assert len(tetrads)==66 and len(C[0] & C[1])==2
    pairs=[(0,13),(1,3),(2,8),(4,10),(5,11),(6,7),(9,23),(12,22),
           (14,18),(15,19),(16,21),(17,20)]
    pm=[sum(1<<x for x in t) for t in pairs]
    assert tetrads=={a|b for a,b in combinations(pm,2)}
    ordered=[frozenset(pairs[i]) for i in [0,1,5,11,7,4,10,9,2,3,6,8]]
    perms=list(maps(X,X))
    im={tuple(ordered.index(frozenset(p[i] for i in pair)) for pair in ordered)
        for p in perms}
    def mobius(a,b,c,d,x):
        u,v=(a,c) if x==11 else ((a*x+b)%11,(c*x+d)%11)
        return (u*pow(v,-1,11))%11 if v else 11
    pgl={tuple(mobius(a,b,c,d,x) for x in range(12))
         for a,b,c,d in product(range(11),repeat=4) if (a*d-b*c)%11}
    assert len(im)==len(perms)==1320 and im==pgl
    print("mixed action",len(im))
    print("Golay partner orbits",244823040//6072,244823040//len(im))

    count_sample(A[0])
    for n in (24,32):
        p=n-1
        squares={i*i%p for i in range(1,p)}
        X=mat([''.join(str(int((j-i)%p in squares)) for j in range(p)) for i in range(p)])
        f,_,W=check(X)
        test_shell_samples(X)
        print("Paley",n,"tetrads",[w[4] for w in W],"shell",sorted(f.items()))

    order28(html)

if __name__=='__main__':
    main(sys.argv[1])
