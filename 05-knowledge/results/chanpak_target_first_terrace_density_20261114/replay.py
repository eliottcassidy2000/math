from fractions import Fraction as F
from functools import reduce
from itertools import product
from operator import mul

def kern(q,a,b):
    total=0
    for mask in range(1<<q):
        bits=[(mask>>j)&1 for j in range(q)]
        j=0
        weight=1
        while j<q:
            if not bits[j]:
                j+=1
                continue
            i=j
            while j<q and bits[j]:
                j+=1
            if not ((a and i==0) or (b and j==q)):
                weight*=2 if (j-i)%2==0 else 0
        total+=weight
    return total

def table(q):
    t=[1,1,2]
    while len(t)<=q:
        t.append(sum(t[-3:]))
    def a(i): return t[i] if i>=0 else 0
    return [[a(q)+a(q-2),a(q)+a(q-1)],[a(q)+a(q-1),a(q)+a(q-1)+a(q-2)]]

def monoid_value(gaps):
    p=[i for i,r in enumerate(gaps) if r]
    vec=[1,1<<gaps[0]]
    for i,j in zip(p,p[1:]):
        k=table(j-i-1)
        vec=[sum(vec[a]*k[a][b] for a in (0,1))*((1<<gaps[j]) if b else 1)
             for b in (0,1)]
    return sum(vec)

def direct(gaps):
    anchors=[0]
    for r in gaps:
        anchors.append(anchors[-1]+r+1)
    n=anchors[-1]+1
    rev=set(zip(anchors,anchors[1:]))
    arcs=[[i!=j and ((i<j)!=((min(i,j),max(i,j)) in rev)) for j in range(n)] for i in range(n)]
    return held_karp(arcs)

def held_karp(arcs):
    n=len(arcs)
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for m,u in enumerate(d):
        for i,num in enumerate(u):
            if num:
                for j,ok in enumerate(arcs[i]):
                    if ok and not (m>>j)&1: d[m|(1<<j)][j]+=num
    return sum(d[-1])

def dyadic_lift(t,k):
    assert t&1 and k>=2
    g=[1] if t%4==3 else []
    m=3 if g else 1
    for r in range(2,k):
        if (t-m)%(1<<(r+1)):
            g.append(r)
            m*=1+(1<<r)
    return g or [k]

def upper(x,den=10**6):
    # Strict arithmetic certificate for a ceiling of x^(-49/50).
    x=F(x)
    a,b=x.numerator,x.denominator
    l,r=0,den
    while r-l>1:
        m=(l+r)//2
        if m**50*a**49>=den**50*b**49: r=m
        else: l=m
    assert r**50*a**49>=den**50*b**49
    return F(r,den)

if __name__=='__main__':
    for q,a,b in product(range(10),(0,1),(0,1)):
        assert kern(q,a,b)==table(q)[a][b]
    R=sum(upper(1+(1<<j)) for j in range(1,9))+upper(512)/(1-upper(2))
    Q=sum(upper(F(v,9)) for v in [17,31,57,105,193])+upper(F(193,9))*upper(F(9,5))/(1-upper(F(9,5)))
    assert R<F(793,1000) and Q<F(1202,1000) and R*Q<F(954,1000)
    for k in range(2,12):
        for t in range(1,1<<k,2):
            g=dyadic_lift(t,k)
            assert reduce(mul,(1+(1<<r) for r in g),1)%(1<<k)==t
    checked=0
    for k in range(1,7):
        for g in product(range(4),repeat=k):
            if g[0] and g[-1] and sum(g)+k+1<=10:
                assert direct(g)==monoid_value(g)
                checked+=1
    # A literal 63 witness, two hubs after an ascending spine; first hub -> second.
    def arcs(w):
        m=len(w)
        def a(i,j):
            if i==j: return False
            if j<i: return not a(j,i)
            if i==m: return True
            if j<m: return True
            return not (int(w[i])>>(1-(j-m))&1)
        return [[a(i,j) for j in range(m+2)] for i in range(m+2)]
    witness="100110"
    assert held_karp(arcs(witness))==63
    print(f'kernel checked, {checked} physical words; R<{float(R):.6f}, Q<{float(Q):.6f}; H=63 witness {witness}')

    checked=0
    for n in range(1,6):
        edge=[(i,j) for i in range(n) for j in range(i+1,n)]
        for code in range(1<<len(edge)):
            B=[[False]*n for _ in range(n)]
            for t,(i,j) in enumerate(edge):
                B[i][j]=bool((code>>t)&1)
                B[j][i]=not B[i][j]
            counts=[]
            for mask in range(1<<n):
                S=[i for i in range(n) if (mask>>i)&1]
                counts.append(held_karp([[B[i][j] for j in S] for i in S]) if S else 1)
            signed=0
            L=0
            for m,q in enumerate(counts):
                x=q*counts[((1<<n)-1)^m]
                L+=x
                signed+=(-1)**m.bit_count()*x
            W=[[False]*(n+2) for _ in range(n+2)]
            for i in range(n):
                W[0][i+1]=W[i+1][n+1]=True
                for j in range(n): W[i+1][j+1]=B[i][j]
            W[n+1][0]=True
            assert signed==0 and held_karp(W)==counts[-1]+L
            if n>=2: assert L%4==0
            if n>=4 and n%2==0: assert L%8==0
            # Strongness and sorted-score decomposition.
            def reaches(B,i):
                S={i}; todo=[i]
                for v in todo:
                    for u,a in enumerate(B[v]):
                        if a and u not in S: S.add(u);todo.append(u)
                return len(S)
            if n>=5 and all(reaches(B,i)==n for i in range(n)):
                ranks=sorted(range(n),key=lambda i:-sum(B[i]))
                backward=[(i,j) for i,j in edge if B[ranks[j]][ranks[i]]]
                slack=[sum(i<k<=j for i,j in backward)-1 for k in range(1,n)]
                tri=sum(B[i][j] and B[j][k] and B[k][i] or B[i][k] and B[k][j] and B[j][i]
                        for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
                drop=[(i,j) for i,j in backward if any(slack[i:j])]
                keep=[(i,j) for i,j in backward if (i,j) not in drop]
                assert sum(slack)<=3*(tri-(n-2))
                assert len(drop)<=6*(tri-(n-2))
                assert all(sum(i<k<=j for i,j in keep)<=1 for k in range(1,n))
            checked+=1
    print(f'wrap and edit checks: {checked} labelled reservoirs')
