"""Exact storage bounds and independent physical/synthetic port controls."""
from fractions import Fraction as F
from math import factorial, comb, ceil
from itertools import combinations, permutations

def nbits(n):
    return (n-1).bit_length()

def bounds():
    out=[]
    for c in (25,100,175):
        t=c//5
        u=(2065+t)//15-25
        M=7020-c
        U=M+2*(675-c)-2*(135-t)
        a,b=F(U,225),F(40,3)
        err=next(i for i in range(200) if i>=a+b and (i-a-b)**2>=4*a*b)
        Q=F((270-c)*49+268*t,225)-8
        L=F(11*u*u,100)-4*Q-F(268*t,225)-err+2*F(36*t-c,15)
        J=F(36*t,15)-F(c,27)
        R=ceil(52*u-L+F(5,12)*J)
        assert 5**12<=50**5
        port=nbits(120*11**5*factorial(54)**15)
        basis=nbits(comb(270,u))
        diag=u
        rows=nbits(50**(R-2*u)*54**(2*u))
        out.append((c,u,err,R,port,basis,diag,rows,2+port+basis+diag+rows))
    print('capacity:',out)
    assert max(row[-1] for row in out)<=comb(249,2)
    assert max(row[-1] for row in out)==30868==comb(249,2)-8
    lo,hi,scale=27182818,27182819,10**7
    e0=sum(F(1,factorial(i)) for i in range(15))
    assert F(lo,scale)<e0<e0+F(16,15*factorial(15))<F(hi,scale)
    assert all(factorial(j)*hi**j<=((j+5)*scale)**j for j in range(1,251))
    for n,bits in [(249,1397),(250,1404)]:
        a=(n*(n+3)+10*(n+1))*scale
        b=2*(n+1)*lo
        assert a**(n+1)<(b**(n+1))<<bits
    return out

def rankperm(p):
    s=list(range(len(p)))
    x=0
    for v in p:
        i=s.index(v)
        x=x*len(s)+i
        s.pop(i)
    return x

def unrankperm(n,k):
    digits=[]
    for r in range(1,n+1):
        k,i=divmod(k,r)
        digits.append(i)
    assert not k
    s=list(range(n))
    return [s.pop(i) for i in reversed(digits)]

def fp(a,b,n,M=(0,1,2,3,4)):
    return tuple((a+b*m)%n for m in M)

def affine_atlas(n):
    # Example of an equivariant five-colour incidence table; not a Moore graph
    # on n=55. n=5 below additionally gets a full graph audit.
    return [fp(a,b,n) for a in range(n) for b in range(1,n)]

def allowed(x):
    return [set((x[i],x[(i+2)%5],x[(i+3)%5])) for i in range(5)]

def possible(x,y):
    return all(z not in t for z,t in zip(y,allowed(x)))

def check_55():
    n=55
    table=affine_atlas(n)
    assert len(table)==n*(n-1)
    for i,j in combinations(range(5),2):
        cells={(x[i],x[j]) for x in table}
        assert len(cells)==n*(n-1) and all(a!=b for a,b in cells)
    # Normalize the first contact to each of the 5 orbit representatives.
    # Adding 5 to all five coordinates is the order-11 action.
    packet=[]
    for a in range(5):
        rows=[x for x in table if x[0]==a]
        rows.sort(key=lambda x:x[1])
        rem=[i for i in range(n) if i!=a]
        for i in range(2,5):
            p=[rem.index(x[i]) for x in rows]
            q=rankperm(p)
            assert unrankperm(54,q)==p
            packet.append(q)
    assert len(packet)==15
    # Exhaust the permitted port-free candidate lists in one physical row
    # for each slope. The affine translation then covers all other rows.
    for b in range(1,n):
        x=fp(0,b,n)
        f=allowed(x)
        for y0 in range(n):
            if y0 in f[0]: continue
            # Same-remote-orbit candidates have intercept 0 mod 5 and
            # identical slope. This test is for the stored off-diagonal rows.
            no_given_port=[y for d in range(1,n) if not(d==b and y0%5==0)
                           for y in [fp(y0,d,n)] if all(a!=z for a,z in zip(x,y))]
            exc=y0 in {x[1],x[4]} # holes to the two distance-one colours
            assert len(no_given_port)<=(54 if exc else 50)
            # At least one entry per cell survives the stronger port-edge
            # exclusions, so the simplest zero-domain test does not reject.
            assert any(possible(x,y) for y in no_given_port)
    print('55-symbol synthetic near-array: all ten projections and 15 rank round-trips')

def audit_graph(adj,d):
    assert len(adj)==d*d+1 and all(len(r)==d for r in adj)
    for i,row in enumerate(adj):
        for j in range(i):
            assert len(row&adj[j])==int(j not in row)

def full_hs(table):
    # 0..4 fixed roots; 5+5c+p the port of colour c and phase p.
    a=[set() for _ in range(50)]
    def edge(u,v):a[u].add(v);a[v].add(u)
    for c in range(5):
        edge(c,(c+1)%5)
        for p in range(5):
            edge(c,5+5*c+p)
            edge(5+5*c+p,5+5*((c+2)%5)+p)
    for i,x in enumerate(table):
        for c,p in enumerate(x):
            edge(30+i,5+5*c+p)
        for j in range(i):
            if possible(x,table[j]):edge(30+i,30+j)
    audit_graph(a,7)

def hs_atlas():
    id=tuple(range(5))
    def compatible(a,b):return len({(x-y)%5 for x,y in zip(a,b)})==5
    opts=[a for p in permutations(range(1,5)) if compatible(a:=(0,*p),id)]
    assert opts==[fp(0,b,5) for b in (2,3,4)]
    assert all(compatible(a,b) for a,b in combinations(opts,2))
    ret=[]
    for order in permutations((2,3,4)):
        table= [fp(a,b,5,(0,1,*order)) for a in range(5) for b in range(1,5)]
        degrees={sum(possible(x,y) for y in table) for x in table}
        ret.append(degrees)
        if degrees=={2}:full_hs(table)
    assert ret==[{2}]+[{0}]*5
    print('HS: only first of six exact pairwise arrays extends to degree 7')

def petersen():
    V=list(combinations(range(5),2))
    A=[[int(not(set(x)&set(y))) for y in V] for x in V]
    adj=[{j for j,x in enumerate(a) if x} for a in A]
    audit_graph(adj,3)
    root=V.index((3,4))
    port=sorted(adj[root])
    rem=[i for i in range(10) if i not in port and i!=root]
    # The order-3 cycle (0 1 2): unique phase fibre for each remote orbit.
    def move(i):
        return V.index(tuple(sorted((j+1)%3 if j<3 else j for j in V[i])))
    orb=[]
    while rem:
        seed=rem[0]
        row={seed,move(seed),move(move(seed))}
        orb.append(row)
        rem=[i for i in rem if i not in row]
    assert len(orb)==2
    for r in orb:
        for v in r:
            own,=adj[v]&set(port)
            nbrs=adj[v]&orb[1-orb.index(r)]
            assert len(nbrs)==2
            assert {next(iter(adj[w]&set(port))) for w in nbrs}==set(port)-{own}
    print('Petersen: one port colour forces the two crossed phases')

def metagraph():
    out=[]
    for m in (2,5,3,4,0,7):
        A=[[int(i<j) for j in range(4)] for i in range(4)]
        for k,(i,j) in enumerate(((0,2),(0,3),(1,3))):
            if m>>k&1:A[i][j],A[j][i]=0,1
        P=list(permutations(range(4)))
        def canon(converse):
            return min(sum((A[p[i]][p[j]]^converse)<<b
                           for b,(i,j) in enumerate(combinations(range(4),2))) for p in P)
        h=sum(all(A[i][j] for i,j in zip(p,p[1:])) for p in P)
        out.append((m,h,canon(0),min(canon(0),canon(1)),(m&1)==(m>>2&1)))
    assert out[0][:4]==(2,5,4,4) and out[1][:4]==(5,5,4,4)
    print('4-site tiling atlas (mask,H,class,merged,blue):',out)

if __name__=='__main__':
    bounds()
    petersen()
    hs_atlas()
    metagraph()
    check_55()
