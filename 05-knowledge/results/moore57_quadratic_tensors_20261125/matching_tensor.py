"""Exact matched-pair quotient and midpoint tests for Petersen and Hoffman--Singleton."""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import sys
import sympy as sp

sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"moore57_thirteen_20261125"))
from moore_chart import hoffman_edge


def quotient(edges, pairs):
    n=len(pairs)
    ids={x:(i,1-2*s) for i,pair in enumerate(pairs) for s,x in enumerate(pair)}
    assert len(ids)==2*n and all(frozenset(p) in edges for p in pairs)
    R,D=sp.zeros(n),sp.zeros(n)
    for edge in edges:
        (a,sa),(b,sb)=(ids[v] for v in edge)
        if a!=b:
            assert not R[a,b]
            R[a,b]=R[b,a]=1
            D[a,b],D[b,a]=sa,sb
    return R,D


def check(edges, pairs, d, matrices=True):
    R,D=quotient(edges,pairs)
    n=R.rows
    neighbors=[{j for j in range(n) if R[i,j]} for i in range(n)]
    signs=[[int(D[i,j]) for j in range(n)] for i in range(n)]
    assert all(len(N)==2*d-2 for N in neighbors) and all(sum(row)==0 for row in signs)
    triangles=Counter()
    for i,j in combinations(range(n),2):
        common=neighbors[i]&neighbors[j]
        routes=Counter((signs[i][k],signs[j][k]) for k in common
                       if signs[k][i]==signs[k][j])
        if R[i,j]:
            assert routes=={(-signs[i][j],-signs[j][i]):1}
            for k in common:
                if k>j:
                    c=sum([signs[i][j]==signs[i][k],signs[j][i]==signs[j][k],signs[k][i]==signs[k][j]])
                    triangles[c]+=1
        else:
            assert routes=={(-1,-1):1,(-1,1):1,(1,-1):1,(1,1):1}
    assert triangles[1]==n*(d-1) and not (set(triangles)-{0,1})
    S=R.multiply_elementwise(D.multiply_elementwise(D.T))
    if matrices:
        I,J=sp.eye(n),sp.ones(n)
        C=R+3*I
        V=S-I
        assert C*C+D.T*D == (4*d-3)*I+8*J
        assert D*D.T+V*V == (4*d-3)*I
        assert D*C+V*D == sp.zeros(n)
    return R,S,D,triangles


def petersen():
    vertices=list(combinations(range(5),2))
    E={frozenset((i,j)) for i,j in combinations(range(10),2) if set(vertices[i]).isdisjoint(vertices[j])}
    def factors(todo):
        if not todo:
            yield ()
            return
        x=min(todo)
        for y in todo:
            if frozenset((x,y)) in E:
                for tail in factors(todo-{x,y}):
                    yield ((x,y),)+tail
    cases=list(factors(set(range(10))))
    assert len(cases)==6
    for pairs in cases:
        R,S,D,tri=check(E,pairs,3)
        assert R==sp.ones(5)-sp.eye(5) and tri=={1:10}
    print("Petersen: all six contractions give K5 and ten single-coherent triangles.")


def hoffman():
    E={frozenset((i,j)) for i,j in combinations(range(50),2) if hoffman_edge(i,j)}
    hist=Counter()
    char_types={}
    for sigma in permutations(range(5)):
        pairs=[(5*i+x,25+5*sigma[i]+(x+i*sigma[i])%5) for i in range(5) for x in range(5)]
        R,S,D,tri=check(E,pairs,7,matrices=False)
        # A characteristic polynomial is invariant under relabelling a contraction.
        char=tuple(R.charpoly().all_coeffs())
        hist[(tri[0],char)] +=1
        if char not in char_types:
            check(E,pairs,7)
            u=5
            eye=sp.eye(25)
            fplus=25-(R-(u-3)*eye).to_DM().rank()
            fminus=25-(R+(u+3)*eye).to_DM().rank()
            splus=25-(S-(u+1)*eye).to_DM().rank()
            sminus=25-(S+(u-1)*eye).to_DM().rank()
            assert (fplus-sminus,splus-fminus)==(3,4)
            # Symbolically same degree-57 calculation, calibrated in degree seven.
            char_types[char]=(tri[0],fplus,fminus,splus,sminus,D.to_DM().rank())
    print("HS row-permutation matchings, by zero-coherent triangles and quotient spectrum:")
    for (nzero,char),v in sorted(hist.items(),key=lambda x:(x[0][0],x[1])):
        print(nzero,v,char_types[char][1:])
    print(len(hist),"classes by statistic,",len(char_types),"distinct characteristic polynomials")


if __name__=="__main__":
    petersen()
    hoffman()
