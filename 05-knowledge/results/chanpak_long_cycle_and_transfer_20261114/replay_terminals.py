"""Finite local Jacobian certificates and independent literal path checks."""
from itertools import product
from functools import lru_cache

C0 = [(1,1,1),(1,2,2),(3,3,4)]
def seq(n, modulus=None):
    C = C0.copy()
    for q in range(3,n+1):
        row=tuple(sum(v[i] for v in C[-3:]) for i in range(3))
        C.append(tuple(x % modulus for x in row) if modulus else row)
    return C
C=seq(160)
def gate(q,r):
    e,f,g=C[q]
    y=1<<r
    return e+y*f,f+y*g

def jac(p,L,R,q,r):
    u,v=gate(q,r); U,V=gate(q+L,r); W,Z=gate(q,r+R)
    assert all(x%p==0 for x in (U-u,V-v,W-u,Z-v))
    return tuple((x//p)%p for x in (U-u,V-v,W-u,Z-v))

def literal_tables(count, arcs):
    m=1<<count
    ends=[[0]*count for _ in range(m)]
    for v in range(count):ends[1<<v][v]=1
    for S in range(1,m):
        for v,t in enumerate(ends[S]):
            if t:
                mask=arcs[v]&~S
                while mask:
                    bit=mask&-mask
                    mask^=bit
                    ends[S|bit][bit.bit_length()-1]+=t
    return [1]+[sum(row) for row in ends[1:]]

def profile(arcs,root):
    n=len(arcs); h=literal_tables(n,arcs); full=(1<<n)-1
    started=[0]*len(h)
    @lru_cache(None)
    def paths(S,v):
        if S==1<<v:return 1
        return sum(paths(S^(1<<v),i) for i in range(n) if (arcs[v]&S)>>i&1)
    for S in range(1,full+1):
        if S>>root&1:started[S]=paths(S,root)
    R=sum(started[S]*h[full^S] for S in range(full+1))
    return h[-1],R

def strong(arcs):
    n=len(arcs); full=(1<<n)-1
    assert all(not (arcs[i]>>i&1) for i in range(n))
    assert all(((arcs[i]>>j&1)+(arcs[j]>>i&1))==1
               for i in range(n) for j in range(i))
    for v in range(n):
        reached=1<<v
        for _ in range(n):
            for i in range(n):
                if reached>>i&1:reached|=arcs[i]
        if reached!=full:return False
    return True

def bridge(arcs,root,r):
    n=len(arcs); k=n+r; full=(1<<(k+1))-1
    added=full^((1<<n)-1)
    return (tuple(a|added for a in arcs[:root])
        +((arcs[root]|added)^(1<<k),)
        +tuple(a|added for a in arcs[root+1:])
        +tuple(full^((1<<(v+1))-1) for v in range(n,k))+(1<<root,))

def terminal(arcs,root,t,q,r):
    arcs=bridge(arcs,root,t)
    for _ in range(q):arcs=bridge(arcs,len(arcs)-1,0)
    return bridge(arcs,len(arcs)-1,r)

# (a, core bitmasks, marked root, R_+)
SOURCES=[
    (1, (0,),0,1),
    (19,(2,4,41,19,7,27),0,23),
    (5,(2,4,9,3),2,7),
    (23,(2,4,25,19,35,15),1,23),
    (25,(2,12,9,49,7,23),0,28),
    (11,(2,20,9,3,13),1,13),
    (13,(2,12,9,17,7),0,14),
    (31,(2,12,33,21,7,27),0,34),
]
if __name__=='__main__':
    for p,L in ((3,13),(5,31),(7,48)):
        assert all((C[q+L][i]-C[q][i])%p==0 for q in range(3) for i in range(3))
    zeros7=[(q,r) for q,r in product(range(48),range(3)) if not any(u%7 for u in gate(q,r))]
    assert zeros7==[(7,2),(12,0),(23,2),(28,0),(39,2),(44,0)]
    def determinant(J):
        x,y,u,v=J;return x*v-y*u
    assert [determinant(jac(7,48,3,q,r))%7 for q,r in zeros7]==[5,4,3,1,6,2]
    nonsingular5={}
    for q,r in product(range(31),(1,3)):
        if determinant(jac(5,31,4,q,r))%5:
            nonsingular5.setdefault(tuple(x%5 for x in gate(q,r)),(q,r))
    assert set(nonsingular5)==set(product(range(5),repeat=2))
    print("Nonsingular odd-r mod-5 coefficient clocks:", sorted(nonsingular5.items()))
    # Show the full source set supplies all mod-2 base classes and all mod-3 cases.
    assert {a%16 for a,_,_,_ in SOURCES}==set(range(1,16,2))
    for a,T,u,P in SOURCES:
        assert strong(T) and profile(T,u)==(a,P)
        assert a%7 and a%3 and P%3 and (a%5 or P%5)
        for t in range(4,10):
            b=(1<<t)*P
            unit={((a*U+b*V)%3) for q in range(13)
                 for U,V in [gate(q,3)]
                 for dU,dV,_,_ in [jac(3,13,2,q,3)]
                 if (a*dU+b*dV)%3}
            assert {1,2}<=unit
        for residue in (0,3,6):
            assert any(
                (a*U+b*V)%9==residue and (a*dU+b*dV)%3
                for q,t in product(range(13),range(4,10))
                for b in [(1<<t)*P] for U,V in [gate(q,3)]
                for _,_,dU,dV in [jac(3,13,2,q,3)])
    for q in (12,28,44):
        assert C[q][0]%16==1
        assert tuple((u-v)%32 for u,v in zip(C[q+48],C[q]))==(16,8,24)
    # Verify M^48 == I mod 8, not only on one state.
    for j in range(3):
        s=[int(j==i) for i in range(3)]
        for _ in range(48):s=s[1:]+[sum(s)%8]
        assert s==[int(j==i) for i in range(3)]
    # Literal bypass of a provably impossible length-seven ready target.
    A=(6,28,104,113,165,195,147,15); root=5
    assert strong(A) and profile(A,root)==(563,1066)
    a,b=563,8*1066
    assert gate(3,1)==(17,20)
    B=terminal(A,root,3,3,1)
    assert len(B)==17 and strong(B) and literal_tables(len(B),B)[-1]==17*a+20*b==180131
    print("Literal terminal-three bypass:",B)
