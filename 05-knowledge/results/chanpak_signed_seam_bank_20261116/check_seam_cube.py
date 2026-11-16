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
def catalog(n,c,V):
    return coordinates(tournament(n,c),V)

def value(A,B,x,y,g,z):
    h,q,E,L,T,S,J=A; H,Q,EE,LL,TT,SS,JJ=B
    p=h*H+g*q*Q*x*y
    e=g*q*EE+H*E
    l=h*LL+g*Q*L
    j=H*J+h*JJ+g*EE*L+z*T*SS
    return p+x*e-y*l-j,(e,l,j)

def actual(ds,rs,guards):
    a,c,b=guards
    tapes=[tournament(a,(1<<(a*(a-1)//2))-1),ds[0],tournament(c,(1<<(c*(c-1)//2))-1),ds[1],tournament(b,(1<<(b*(b-1)//2))-1)]
    x=sum(map(len,tapes))+1;t=x-1
    A=[[0]*(x+1) for _ in range(x+1)];ix=[]
    for k,D in enumerate(tapes):
        for u in range(len(D)):ix.append((k,u))
    for a0,(k,u) in enumerate(ix):
        for b0,(j,v) in enumerate(ix):
            A[a0][b0]=tapes[j][u][v] if k==j else int(k<j)
        m=rs[(k-1)//2] if k in (1,3) else 0
        A[a0][x]=(m>>u)&1;A[x][a0]=1-A[a0][x];A[a0][t]=1
    A[t][x]=1
    return sum(ends(A)[-1])
def check():
    examples=[
      ((4,17),(3,2),((9,1),(10,1)),[(b+1,0,b) for b in (1,2,3)]),
      ((5,41),(5,41),((3,24),(9,9)),[(a,0,b) for a,b in ((1,1),(1,3),(3,2))]),
      ((5,855),(5,221),((23,13),(13,25)),[(2,0,b) for b in (1,2,3)]),
    ]
    for a,b,ps,guards in examples:
        counts=[]
        for V,W in ps:
            p=catalog(*a,V);q=catalog(*b,W)
            print('cores',a,b,'masks',V,W,'profiles',p,q)
            for (i,j,k) in guards:
                h,coeffs=value(p,q,1<<i,1<<k,1<<j,int(j==0))
                actual_h=actual([tournament(*a),tournament(*b)],[V,W],(i,j,k))
                assert h==actual_h
                print('guard lengths',i,j,k,'projection',coeffs,'H',h)
            counts.append(h)
        assert abs(counts[0]-counts[1])==2
    D=tournament(9,53826135213)
    move=(53,77,141)
    for a,b in [(1,1),(1,2),(2,1),(2,3)]:
        for i,V in enumerate(move):
            assert coordinates(D,V)==(1521,26412,10422,11997,631,750,3908-2*i)
            p=actual([D,[]],[V,0],(a,0,b))
            assert p==1521+26412*2**(a+b)+10422*2**a-11997*2**b-(3908-2*i)
    assert tuple((move[i]^move[i+1]).bit_count() for i in (0,1))==(4,2)
    print('three-mask fixed-core route checked at four independently varying guard pairs')
    print('fixture replay passes')
if __name__=='__main__':check()
