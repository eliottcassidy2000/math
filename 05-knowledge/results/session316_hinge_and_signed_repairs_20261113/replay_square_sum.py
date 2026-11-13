from itertools import combinations
from math import isqrt
from collections import Counter
from functools import lru_cache

def graph(n):
    return {i:{j for j in range(1,n+1) if i!=j and isqrt(i+j)**2==i+j} for i in range(1,n+1)}

def component_count(A,S=()):
    U=set(A)-set(S)
    c=0
    while U:
        V=[U.pop()]
        c+=1
        while V:
            I=A[V.pop()] & U
            U-=I
            V.extend(I)
    return c

def forced_graph(n):
    A=graph(n)
    for x,S in ((1,{3,8}),(3,{1,6})):
        for y in A[x]-S:
            A[y].remove(x)
        A[x]=S
    return A

def propagate(A,forced,ends):
    E={tuple(sorted((x,y))) for x in A for y in A[x]}
    f=set(forced)
    for x,y in f:assert (x,y) in E
    N={x:{e for e in E if x in e} for x in A}
    while True:
        change=False
        for x in sorted(A):
            d=1 if x in ends else 2
            avail=N[x]&E
            used=N[x]&f
            if len(used)>d or len(avail)<d:return "degree"
            if len(used)==d and avail-used:
                E-=avail-used
                change=True
            if len(avail)==d and avail-used:
                f|=avail
                change=True
        if change:continue
        P={x:x for x in A}
        def root(x):
            while P[x]!=x:x=P[x]
            return x
        for x,y in f:
            x,y=root(x),root(y)
            if x==y:return "cycle"
            P[x]=y
        return "open"

def segment_witness(n,seg):
    A=graph(n)
    U=sum(1<<x for x in seg)
    allv=sum(1<<x for x in A)
    @lru_cache(None)
    def walk(x,M):
        if M==allv:return ()
        for y in sorted(A[x]):
            if M>>y&1:continue
            if y not in seg:
                V=walk(y,M|1<<y)
                if V is not None:return (y,)+V
            elif y in (seg[0],seg[-1]):
                Q=seg if y==seg[0] else seg[::-1]
                V=walk(Q[-1],M|U)
                if V is not None:return Q+V
    for x in A:
        if x in seg:
            if x not in (seg[0],seg[-1]):continue
            Q=seg if x==seg[0] else seg[::-1]
            V=walk(Q[-1],U)
        else:
            Q=(x,)
            V=walk(x,1<<x)
        if V is not None:return Q+V

for n in range(8,21):
    assert component_count(forced_graph(n))>1
cuts=[(4,),(4,),(4,7,14),(2,5,12),(5,7,12,14),(7,10,14,23)]
for n,S in enumerate(cuts,21):
    assert component_count(forced_graph(n),S)>len(S)+1
wanted=[{"degree":18,"cycle":8},{"degree":27},{"degree":28},{"degree":29},{"degree":465}]
for n,w in enumerate(wanted,27):
    A=forced_graph(n)
    leaves={x for x in A if len(A[x])==1}
    R={(1,8),(1,3),(3,6)}
    assert Counter(propagate(A,R,set(E)) for E in combinations(A,2) if leaves<=set(E))==w

for s,n in ((3,32),(4,27),(5,25)):
    seg=(s*s-1,1,3,s*s-3)
    assert all(segment_witness(m,seg) is None for m in range(max(seg),n))
    P=segment_witness(n,seg)
    assert P is not None and set(P)==set(range(1,n+1))
    assert all(y in graph(n)[x] for x,y in zip(P,P[1:]))
    assert any(P[i:i+4] in (seg,seg[::-1]) for i in range(n-3))
    print(s,n,P)

def rooted_paths(n,root):
    A=graph(n)
    mask=sum(1<<x for x in A)
    def visit(P,U):
        if U==mask:
            yield P
            return
        for x in sorted(A[P[-1]]):
            if not U>>x&1:
                yield from visit(P+(x,),U|1<<x)
    yield from visit((root,),1<<root)

for n,W in [(23,(3,3)),(25,(10,7)),(26,(12,0))]:
    P=list(rooted_paths(n,18))
    assert (len(P),sum(abs(p.index(1)-p.index(3))==1 for p in P))==W

pairs=[(s,t) for s in range(3,7) for t in range(3,8) if isqrt(s*s+t*t-4)**2==s*s+t*t-4]
assert pairs==[(6,7)]
print("Hinge cut/forcing certificates, first segment orders, and edge support passed")



PORT=(1,8,28,78,66,16,84,34,2,7,9,72,22,27,77,4,5,31,81,19,69,
75,25,11,14,64,80,30,70,20,29,79,21,71,73,23,13,63,18,68,
76,26,74,24,12,37,87,82,32,17,67,33,83,86,36,85,35,65,15,10,6,3)

def tok(P):
    seen=set()
    i=0
    while i<len(P):
        a=P[i]
        if a<=12:
            yield (False,a,a)
            i+=1
            continue
        b=P[i+1]
        c=min(a,b)-25
        assert abs(b-a)==50 and -12<=c<=12 and c not in seen
        seen.add(c)
        yield (True,a,b)
        i+=2
    assert seen==set(range(-12,13))

def squares(w):
    return all(isqrt(x+y)**2==x+y for x,y in zip(w,w[1:]))

T=tuple(tok(PORT))
assert sorted(a for isblock,a,b in T if not isblock)==list(range(1,13))
assert len(T)==37 and all(isqrt(x[2]+y[1])**2==x[2]+y[1] for x,y in zip(T,T[1:]))

def lift(P):
    assert len(P)%2 and P[0]==1 and P[-1]==3 and squares(P)
    D=[]
    for isblock,a,b in T:
        if not isblock:
            D.append(a)
            continue
        c=min(a,b)-25
        B=[25*x+(c if i%2==0 else -c) for i,x in enumerate(P)]
        D.extend(B if a<b else reversed(B))
    assert len(D)==25*len(P)+12 and set(D)==set(range(1,len(D)+1))
    assert D[:2]==[1,8] and D[-2:]==[6,3] and squares(D)
    E={frozenset((a,b)) for a,b in zip(D,D[1:]+D[:1])}
    assert all(frozenset(e) in E for e in ((8,1),(1,3),(3,6),(2,7),(4,5)))
    return D

P=[1,8,28,21,4,32,17,19,6,30,34,15,10,26,23,13,12,24,25,11,5,20,29,35,14,2,7,18,31,33,16,9,27,22,3]
for m in range(1,4):
    P=lift(P)
    assert len(P)==(71*25**m-1)//2
    print("marked 25-fold lift",m,len(P))

