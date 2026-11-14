"""Literal bilateral splice and the (11,13,26) transverse transfer."""
from bridge import bridge,rooted_profile
from physical_seed import C,adjacency,transfer

W=(47,56); V=(56,-47); Q=5345
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def apply(q,z):
    e,f,g=C(q);a,b=z
    return e*a+f*b,f*a+g*b
def value(q,x,y):return dot(x,apply(q,y))
def ready_pair(R,root):
    H,A,B=rooted_profile(R,root)
    assert H==A and -B>=1 and H+B>=2
    return -B,H+B

def splice(A,u,B,v,q):
    assert q>=0
    A=list(A)
    for _ in range(max(0,q-1)):
        A=bridge(A,u,0);u=len(A)-1
    n,m=len(A),len(B)
    right=((1<<m)-1)<<n
    rows=[a|right for a in A]+[(((1<<m)-1)^(1<<j)^b)<<n for j,b in enumerate(B)]
    if q:
        rows[u]^=1<<(n+v)
        rows[n+v]|=1<<u
    return rows

def encode(z):
    u=dot(W,z); v=dot(V,z)
    j,s=divmod(v-2958*u,2*Q)
    assert s==0
    return u,j

def lookup_target(N,q,left,right):
    """Return a literal graph from bank entries (rows, root, ready_pair)."""
    if q not in (26,29,30): raise ValueError(q)
    ell={26:9,29:56,30:103}[q]
    e,f,g=C(q)
    D=(e+g+2*ell*2958)//Q
    L={}; R={}
    for bank,out in ((left,L),(right,R)):
        for rows,root,z in bank:
            u,j=encode(z)
            if u not in out:out[u]={}
            out[u][j]=(rows,root,z)
    for u,As in L.items():
        for v,Bs in R.items():
            if e*u*v>56*56*N:continue
            c,r=divmod(N-D*u*v,2*ell)
            if r:continue
            if len(As)<=len(Bs):
                keys=[(j,(c-v*j)//u) for j in As if (c-v*j)%u==0]
            else:
                keys=[((c-u*k)//v,k) for k in Bs if (c-u*k)%v==0]
            for j,k in keys:
                if j in As and k in Bs:
                    a,x,z=As[j]; b,y,w=Bs[k]
                    assert value(q,z,w)==N
                    return splice(a,x,b,y,q),(z,w)
    return None

def path_count(rows):
    from functools import cache
    n=len(rows)
    @cache
    def last(S,i):
        T=S^(1<<i)
        if not T:return 1
        return sum(last(T,j) for j in range(n) if T>>j&1 and rows[j]>>i&1)
    return sum(last((1<<n)-1,i) for i in range(n))

if __name__=='__main__':
    t={30:0,31:1,32:1}
    for r in range(29,-1,-1): t[r]=t[r+3]-t[r+2]-t[r+1]
    for r in range(33,65): t[r]=t[r-1]+t[r-2]+t[r-3]
    assert all(value(q,V,V)==t[q] for q in t)
    assert [r for r in t if not t[r]]==[30,29,26,13]
    assert apply(11,(1,2))==(2209,2632)==tuple(47*x for x in W)
    assert apply(13,V)==(0,0) and apply(26,V)==tuple(9*x for x in W)
    for q,l in ((26,9),(29,56),(30,103)):
        assert apply(q,V)==tuple(l*x for x in W)
        trace=C(q)[0]+C(q)[2]
        assert C(q)[0]*C(q)[2]-C(q)[1]**2==-l*l
        D=(trace+2*l*2958)//Q
        assert D*Q==trace+2*l*2958
        print(q,l,D)
    for a in (1,2):
        for b in (1,2):
            P=bridge((0,),0,a);R=bridge((0,),0,b)
            for q in range(4):
                H=path_count(splice(P,len(P)-1,R,len(R)-1,q))
                assert H==value(q,(1,1<<a),(1,1<<b))
    O=(2,20,9,3,13)
    assert rooted_profile(O,1)[:2]==(11,13)
    O=bridge(O,1,1)
    assert ready_pair(O,6)==(11,26)
    R=bridge((0,),0,1)
    for q in range(3):
        assert path_count(splice(O,len(O)-1,R,len(R)-1,q))==value(q,(11,26),(1,2))
    # Two independent literal switches; one pair on either side.
    A=((1,0,1,0,2,1),(1,0,0,2,3))
    B=((1,0,1,0,4,0,1,0,1),(4,0,0,2,0,1,0,0,1))
    X=tuple(map(transfer,A));Y=tuple(map(transfer,B))
    assert tuple(a-b for a,b in zip(*X))==tuple(2*x for x in V)
    assert tuple(a-b for a,b in zip(*Y))==tuple(2*x for x in V)
    for q in (13,23,26,27,29,30):
        checks=[value(q,x,y) for x in X for y in Y]
        assert checks[0]+checks[3]-checks[1]-checks[2]==4*t[q]
        if q==26:
            assert checks[1]-checks[3]==18*dot(W,Y[1])
            assert checks[2]-checks[3]==18*dot(W,X[1])
        # Independent interval packing verifies all four actual mixed graphs.
        from physical_seed import interval
        assert checks==[interval(a+(0,)*q+tuple(reversed(b))) for a in A for b in B]
    # Make the 11/13 source interact with the reset-free source pair at 26.
    H=[value(26,(11,26),y) for y in Y]
    assert H[0]-H[1]==18*1973
    P=[splice(O,len(O)-1,adjacency(w),len(adjacency(w))-1,26) for w in B]
    assert all(len(x)==50 for x in P)
    print('literal nonterrace q26 pair',H)
    from math import isqrt
    p=H[1]
    assert p&1 and all(p%d for d in range(3,isqrt(p)+1,2))
    left=((O,len(O)-1,(11,26)),)
    right=tuple((adjacency(w),len(adjacency(w))-1,z) for w,z in zip(B,Y))
    assert lookup_target(p,26,left,right)==(P[1],((11,26),Y[1]))
    assert lookup_target(p+2,26,left,right) is None
