from collections import Counter
from itertools import combinations, permutations
import json

def graph(n, edges):
    A=[set() for _ in range(n)]
    for a,b in edges:
        assert a!=b
        A[a].add(b);A[b].add(a)
    return A

def verify(A,k):
    assert len(A)==k*k+1 and all(len(s)==k for s in A)
    for x,y in combinations(range(len(A)),2):
        assert int(y in A[x])+len(A[x]&A[y])==1

def petersen():
    E=[(i,(i+1)%5) for i in range(5)]+[(5+i,5+(i+2)%5) for i in range(5)]
    E += [(i,5+i) for i in range(5)]
    return graph(10,E)

def hs():
    E=[]
    for i in range(5):
        for x in range(5):
            E.extend([(5*i+x,5*i+(x+1)%5),(25+5*i+x,25+5*i+(x+2)%5)])
            E.extend((5*i+x,25+5*j+(x+i*j)%5) for j in range(5))
    return graph(50,E)

def cycles5(A):
    def visit(p):
        if len(p)==5:
            if p[1]<p[-1] and p[0] in A[p[-1]]:
                yield tuple(sorted(p))
        else:
            for y in A[p[-1]]:
                if y>p[0] and y not in p:
                    yield from visit(p+[y])
    return {t for i in range(len(A)) for t in visit([i])}

def comps(A,R):
    unseen=set(R); sizes=[]
    while unseen:
        x=unseen.pop(); C={x}; L=[x]
        for y in L:
            for z in A[y]&unseen:
                unseen.remove(z);C.add(z);L.append(z)
        sizes.append(len(C))
    return tuple(sorted(sizes))

def shell(A,H,k,d):
    assert len(H)==d*d+1
    assert all(len(A[u]&H)==d for u in H)
    F={u:A[u]-H for u in H}
    Ftot=set().union(*F.values())
    assert sum(map(len,F.values()))==len(Ftot)==len(H)*(k-d)
    assert all(len(U)==k-d for U in F.values())
    Z=set(range(len(A)))-H-Ftot
    assert len(Z)==(k-d)*(k-1-d*(d-1))
    assert all(len(A[z]&Z)==k-len(H) for z in Z)
    codes={x:A[x]&Z for x in Ftot}
    for u in H:
        for v in H:
            for x in F[u]:
                assert len(A[x]&F[v])==int(u!=v and v not in A[u])
    for z in Z:
        assert all(len(A[z]&U)==1 for U in F.values())
    for x,S in codes.items():
        covered=set()
        for z in S:
            ball={z}|(A[z]&Z)
            assert not (covered&ball)
            covered |= ball
        remaining=Z-covered
        covered2=set()
        for y in A[x]&Ftot:
            B=codes[y]
            assert not (covered2&B)
            covered2|=B
        assert remaining==covered2
    if d==2 and k==7:
        assert comps(A,Ftot)==(5,5,5,5,5)
        C5=[]
        for u in sorted(Z):
            if any(u in S for S in C5): continue
            S={u}; Q=[u]
            for v in Q:
                for w in A[v]&Z:
                    if w not in S:
                        S.add(w);Q.append(w)
            C5.append(S)
        anchors=[H]+C5; assert len(anchors)==5 and all(len(S)==5 for S in anchors)
        labels=[]
        for x in Ftot:
            name=[]
            for H0 in anchors:
                E=A[x]&H0; assert len(E)==1
                name.append(next(iter(E)))
            labels.append(name)
        for i,j in combinations(range(5),2):
            assert len({(x[i],x[j]) for x in labels})==25
    return {"fibers":len(F),"width":k-d,"residual":len(Z),"components":comps(A,Z)}

def edge_grid(A,k,x,y):
    assert y in A[x]
    L=A[x]-{y}; R=A[y]-{x}
    U=set(range(len(A)))-{x,y}-L-R
    labels={};n=k-1
    assert len(U)==n*n
    for a in L:
        for b in R:
            C=A[a]&A[b]
            assert len(C)==1
            z=next(iter(C));assert z in U and z not in labels
            labels[z]=(a,b)
    assert set(labels)==U
    assert all(len(A[z]&U)==k-2 for z in U)
    for z in U:
        a,b=labels[z]
        nbrs=[labels[t] for t in A[z]&U]
        assert all(c!=a and d!=b for c,d in nbrs)
        assert len({c for c,d in nbrs})==len({d for c,d in nbrs})==k-2
    for v,w in combinations(U,2):
        a,b=labels[v];c,d=labels[w]
        assert int(w in A[v])+len(A[v]&A[w]&U)==int(a!=c and b!=d)
    return comps(A,U)


def best_cyclic_rows(A,k):
    n=k-1
    groups=set()
    for p in permutations(range(1,n)):
        q=(0,)+p
        f=tuple(q[(q.index(j)+1)%n] for j in range(n))
        g=tuple(range(n));S=[]
        for _ in range(n):
            S.append(g);g=tuple(f[t] for t in g)
        groups.add(frozenset(S))
    hist=Counter()
    for u,N in enumerate(A):
        for v in N:
            if v<u: continue
            left=sorted(N-{v}); right=sorted(A[v]-{u})
            grid=[[next(iter(A[a]&A[b])) for b in right] for a in left]
            for i,line in enumerate(grid):
                row=set()
                for j,other in enumerate(grid):
                    if j==i:continue
                    row.add(tuple(next(b for b,z in enumerate(other) if z in A[x]) for x in line))
                if k==7:
                    assert all(all(t[x]!=x and t[t[x]]==x for x in range(n)) for t in row)
                hit=max(len(row&G) for G in groups)
                hist[hit]+=1
    return dict(sorted(hist.items()))

def main():
    P=petersen();H=hs()
    verify(P,3);verify(H,7)
    pcycles=cycles5(P);hcycles=cycles5(H)
    assert len(pcycles)==12 and len(hcycles)==1260
    assert all(shell(P,set(C),3,2)["residual"]==0 for C in pcycles)
    base=set(range(5))|set(range(25,30))
    nest=shell(H,base,7,3)
    htypes=Counter(shell(H,set(C),7,2)["components"] for C in hcycles)
    gridsP=Counter(edge_grid(P,3,x,y) for x,N in enumerate(P) for y in N if y>x)
    gridsH=Counter(edge_grid(H,7,x,y) for x,N in enumerate(H) for y in N if y>x)
    codec=lambda D:{",".join(map(str,k)):v for k,v in sorted(D.items())}
    print(json.dumps({"pentagons_in_P":len(pcycles),"pentagon_shell_types_in_HS":codec(htypes),
         "P_in_HS":nest,"best_cyclic_rows_P":best_cyclic_rows(P,3),"best_cyclic_rows_HS":best_cyclic_rows(H,7),"P_edge_grids":codec(gridsP),"HS_edge_grids":codec(gridsH)},indent=2))

if __name__=="__main__":
    main()
