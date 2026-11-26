"""Integer Wilson-action tests on exact small Moore charts and the marked board."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import sys

from nested_core_check import hs, graph as fixture, orbit_chart, legal
from slot_interface_check import SOURCE, leaf, read_board, BRANCHES

def inv(p):
    return tuple(sorted(range(len(p)),key=p.__getitem__))

def mul(p,q):
    return tuple(q[x] for x in p)

def petersen():
    c=list(map(set,combinations(range(5),2)))
    return [{j for j in range(10) if not x.intersection(c[j])} for x in c]

def chart(A,r=0):
    heads=sorted(A[r])
    S=[sorted(A[h]-{r}) for h in heads]
    k=len(S);m=k-1
    assert len(set.union(*map(set,S)))==k*m
    P={}
    for i,j in combinations(range(k),2):
        p=[]
        for v in S[i]:
            Y=[y for y,w in enumerate(S[j]) if w in A[v]]
            assert len(Y)==1
            p.append(Y[0])
        P[i,j]=tuple(p)
        P[j,i]=inv(P[i,j])
    return P,k

def energy(P,k):
    e=0
    a,b=0,0
    for i,j in combinations(range(k),2):
        d=P[i,j]
        x=[mul(P[i,l],P[l,j]) for l in range(k) if l not in (i,j)]
        a+=sum(sum(u==v for u,v in zip(t,d)) for t in x)
        b+=sum(sum(u==v for u,v in zip(s,t)) for s,t in combinations(x,2))
        e+=sum(sum((n-1)**2 for n in Counter(row).values())+(k-1-len(set(row)))
               for row in zip(d,*x))
    assert a%3==b%2==0 and e==2*(a+b)
    return e,a//3,b//2

def leafgraph(P,k):
    m=k-1
    G=[set() for _ in range(k*m)]
    for i,j in combinations(range(k),2):
        for s,t in enumerate(P[i,j]):
            a,b=i*m+s,j*m+t
            G[a].add(b);G[b].add(a)
    return G

def counts(G):
    t=a=0
    for x,y in combinations(range(len(G)),2):
        c=len(G[x]&G[y])
        t+=c*(y in G[x])
        a+=c*(c-1)//2
    assert t%3==a%2==0
    return t//3,a//2

def closure(G,m):
    ID=tuple(range(m))
    K={ID}; Q=[ID]
    for v in Q:
        for g in G:
            w=mul(v,g)
            if w not in K:K.add(w);Q.append(w)
    return K

def controls():
    for A,expected in [(petersen(),2),(hs(),720)]:
        P,k=chart(A);m=k-1
        G=[mul(mul(P[0,i],P[i,j]),P[j,0]) for i,j in combinations(range(1,k),2)]
        H=closure(G,m)
        print(k,"rooted sheet holonomy group order",len(H))
        assert len(H)==expected
        assert energy(P,k)==(0,0,0)
        edit_hist=Counter()
        for g in H:
            gg=[inv(g)]+[mul(mul(P[i,0],inv(g)),P[0,i]) for i in range(1,k)]
            gi=list(map(inv,gg))
            Q={(i,j):mul(mul(gg[i],p),gi[j]) for (i,j),p in P.items()}
            assert all(Q[i,0]==P[i,0] for i in range(1,k))
            E=[(i,j) for i,j in combinations(range(k),2) if Q[i,j]!=P[i,j]]
            if not E:continue
            degrees=Counter(v for ed in E for v in ed)
            m0=k-len(degrees)
            assert all(d>=m0+1 for d in degrees.values())
            edit_hist[m0,len(E)]+=1
        print("Isomorphic controls fixing one full branch star:",sorted(edit_hist.items()))
        for i,j in combinations(range(k),2):
            for r,s in combinations(range(m),2):
                Q=P.copy()
                p=list(P[i,j])
                p[r],p[s]=p[s],p[r]
                Q[i,j]=tuple(p);Q[j,i]=inv(p)
                assert energy(Q,k)==(4*(2*k-3),2,2*(k-3))
                assert counts(leafgraph(Q,k))==(2,2*(k-3))
    print("All 318 single-link transpositions checked via both routes and leaf-graph cycles.")

def marked(board, log):
    ref,D=read_board(board)
    C,S=orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
    A=fixture(C,S)
    for row in Path(log).read_text().splitlines():
        budget,claimed,*ys=map(int,row.split())
        assert 0<=sum(y!=r for y,r in zip(ys,ref) if r>=0)<=budget and ys[0]==7
        W=[leaf(h,y) for h,y in zip(BRANCHES,ys)]
        assert all(legal(A,SOURCE,w) for w in W)
        collisions=Counter()
        for g in range(3):
            for t,k in Counter(R[y][g] for R,y in zip(D,ys)).items():
                if k>1:collisions[g,t]=k*(k-1)//2
        assert sum(collisions.values())==claimed
        # Independent physical check: for each pair the number of 4-cycles
        # through both new edges is their number of old common neighbors.
        assert not any(y in A[x] for x,y in combinations(W,2))
        assert sum(len(A[x]&A[y]) for x,y in combinations(W,2))==claimed
        print("At most",budget,"changes, collision witness",claimed,
              " ; changes",[(h,r,y) for h,r,y in zip(BRANCHES,ref,ys) if r!=y])

if __name__=='__main__':
    controls()
    if len(sys.argv)>1:marked(*sys.argv[1:])
