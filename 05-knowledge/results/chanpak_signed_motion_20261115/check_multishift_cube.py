"""Independent exact fixtures for projected moving-register head windows."""
from itertools import product
def tournament(n,c):
    a=[[0]*n for _ in range(n)]
    k=0
    for u in range(n):
        for v in range(u+1,n):
            a[u][v]=(c>>k)&1;a[v][u]=1-a[u][v];k+=1
    return a
def ends(a):
    n=len(a); f=1<<n; out=[[0]*n for _ in range(f)]
    for u in range(n):out[1<<u][u]=1
    for s in range(1,f):
        for u,w in enumerate(out[s]):
            if w:
                for v in range(n):
                    if not (s>>v&1) and a[u][v]:
                        out[s|1<<v][v]+=w
    return out
def coordinates(a,V):
    d=ends(a);r=ends(list(zip(*a)));f=len(d)-1
    cnt=[1]+[sum(row) for row in d[1:]]
    V=[i for i in range(len(a)) if (V>>i)&1]
    t=[sum(row[i] for i in V) for row in d]
    s=[sum(row[i] for i in V) for row in r]
    h=cnt[-1]
    q=sum(cnt[i]*cnt[f^i] for i in range(f+1))
    E=sum(t[i]*cnt[f^i] for i in range(f+1))
    L=sum(cnt[i]*s[f^i] for i in range(f+1))
    J=sum(t[i]*s[f^i] for i in range(f+1))
    return (h,q,E,L,t[-1],s[-1],J)
def lift(a,V,x,y):
    n=len(a);m=x+n+y
    out=[[0]*(m+2) for _ in range(m+2)]
    for u in range(m):
        for v in range(m):
            out[u][v]=a[u-x][v-x] if x<=u<x+n and x<=v<x+n else int(u<v)
        out[u][m]=1
        b=int(x<=u<x+n and (V>>(u-x))&1)
        out[u][m+1]=b;out[m+1][u]=1-b
    out[m][m+1]=1
    return out
def strong(a):
    def cover(b):
        s={0}
        for _ in a:s|={v for u in s for v in range(len(a)) if b[u][v]}
        return len(s)==len(a)
    return cover(a) and cover(list(zip(*a)))
def check(c,V,f,a,b):
    t=tournament(8 if c>=1<<21 else 7,c)
    p=coordinates(t,V)
    assert strong(t) and strong(lift(t,V,a,b))
    assert sum(ends(lift(t,V,a,b))[-1])==f(a,b)
    return p
if __name__=='__main__':
    for j,V in enumerate((14,38,42)):
        f=lambda a,b:1276*4**a-42*2**a-29-2*j
        p=check(1265214,V,f,1,1)
        check(1265214,V,f,2,2)
        assert p[2]-p[3]==-42 and p[-1]==124+2*j
        print('fixed-core',1265214,V,p)
    for j,(c,V) in enumerate(((1564245,15),(1888791,82),(101621,83),(758958,28))):
        f=lambda a,b:2360*2**b+1271+2*j
        p=check(c,V,f,1,1)
        check(c,V,f,1,3)
        assert 2*p[1]-p[3]==2360 and p[0]+2*p[2]-p[-1]==1271+2*j
        print('moving-core',c,V,p)

    n8=((75823804,45),(47535549,98),(207924374,106),(257190115,116),(253178809,51),(114827917,58))
    for j,(c,V) in enumerate(n8):
        f=lambda a,b:7617*2**b+3617+2*j
        p=check(c,V,f,1,1)
        check(c,V,f,1,3)
        assert 2*p[1]-p[3]==7617 and p[0]+2*p[2]-p[-1]==3617+2*j
        print('longer route',c,V,p)
    n7=((1564245,15),(1888791,82),(101621,83),(758958,28))
    routes=[
      (7,n7,((5,3,4,6,0,2,1),(5,2,0,1,6,4,3),(5,4,1,0,3,2,6)), ((2,3),)*3),
      (8,n8,((0,2,6,5,3,1,4,7),(7,5,3,0,2,1,4,6),(7,6,2,0,1,4,3,5),
              (2,4,7,1,6,0,3,5),(4,1,2,0,6,5,3,7)),((5,1),(4,3),(4,0),(4,2),(2,4)))]
    for n,sequence,ps,profile in routes:
        a=tournament(n,sequence[0][0]);V=sequence[0][1]
        for (code,W),p,(exp_internal,exp_contacts) in zip(sequence[1:],ps,profile):
            b=tournament(n,code)
            c=[[b[p[u]][p[v]] for v in range(n)] for u in range(n)]
            U=sum(((W>>p[u])&1)<<u for u in range(n))
            internal=sum(a[u][v]!=c[u][v] for u in range(n) for v in range(u+1,n))
            contacts=(U^V).bit_count()
            assert (internal,contacts)==(exp_internal,exp_contacts)
            a,V=c,U
        print('coherent route flips',n,profile)
