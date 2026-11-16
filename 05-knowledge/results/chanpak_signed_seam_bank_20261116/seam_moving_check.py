"""Independent literal head/path-cover checks for the two seam banks."""
from itertools import combinations


def tor(n,k):
    g=[[0]*n for _ in range(n)]
    for i,(u,v) in enumerate(combinations(range(n),2)):
        g[u][v]=(k>>i)&1
        g[v][u]=1-g[u][v]
    return g


def ends(g):
    n=len(g)
    words=[sum(1<<v for v in range(n) if g[u][v]) for u in range(n)]
    d=[[0]*n for _ in range(1<<n)]
    for u in range(n):
        d[1<<u][u]=1
    R=(1<<n)-1
    for s,row in enumerate(d):
        outside=R^s
        for u,k in enumerate(row):
            if k:
                bits=words[u]&outside
                while bits:
                    b=bits&-bits
                    d[s|b][b.bit_length()-1]+=k
                    bits^=b
    return d


def profile(g,m):
    d=ends(g); r=ends(tuple(zip(*g))); R=len(d)-1
    p=[1]+list(map(sum,d[1:]))
    E=[sum(k for u,k in enumerate(row) if (m>>u)&1) for row in d]
    L=[sum(k for u,k in enumerate(row) if (m>>u)&1) for row in r]
    return (p[R],sum(p[s]*p[R^s] for s in range(R+1)),
            sum(E[s]*p[R^s] for s in range(R+1)),
            sum(L[s]*p[R^s] for s in range(R+1)),
            sum(E[s]*L[R^s] for s in range(R+1)),E[R],L[R])


def literal(g,m,a,b,c):
    n=len(g); U=a+2*n+b+c
    t,x=U,U+1
    P=[[0]*(U+2) for _ in range(U+2)]
    D=a
    F=a+n+c
    V={j for i in range(n) if (m>>i)&1 for j in (D+i,F+i)}
    for u,v in combinations(range(U+2),2):
        if v==x:
            uv=u==t or u in V
        elif D<=u<v<D+n:
            uv=g[u-D][v-D]
        elif F<=u<v<F+n:
            uv=g[v-F][u-F]
        else:
            uv=True
        P[u][v]=int(uv)
        P[v][u]=1-int(uv)
    if len(P)!=U+2:
        raise AssertionError
    return sum(ends(P)[-1])


def require(cond):
    if not cond:
        raise AssertionError('physical comparison failed')


def verify(n,code,m,c,budget,s,expect,en,h,q):
    g=tor(n,code)
    H,Q,E,L,J,T,S=profile(g,m)
    require((H,Q,E,L,J,T,S)==(h,q,*expect))
    energy=2*h*J+2**c*L*L+(T*T if c==0 else 0)
    require(energy==en)
    a=s+(budget-c)//2
    require(literal(g,m,a,a,c)==q*q*2**(2*a+c)+h*h-energy)


if __name__=='__main__':
    for m,e,E in [(50,(403,364,114,28,23),146732),
                  (11,(356,365,110,23,24),146734),
                  (38,(313,370,80,20,23),146740)]:
        verify(7,1149404,m,0,0,1,e,E,59,940)
    for m,c,e,E in [(9,0,(12,12,4,2,2),188),
                    (8,2,(5,7,0,1,1),196),
                    (10,0,(14,13,2,3,2),198)]:
        verify(4,4,m,c,2,1,e,E,5,32)
    print('Six exact profiles and six literal tournament reads passed')
