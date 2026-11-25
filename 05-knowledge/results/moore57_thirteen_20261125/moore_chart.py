"""Exact small certificates for induced square-sum prefixes and mod-13 excision."""
from itertools import combinations
from math import comb, isqrt, prod
import sympy as sp

from sympy.polys.domains import GF


def square_sum(i, j):
    return i != j and isqrt(i+j)**2 == i+j


def square_rows(n):
    return {i: {j for j in range(1, n+1) if square_sum(i,j)} for i in range(1, n+1)}


def hoffman_edge(x, y):
    hx, j = divmod(x, 5)
    hy, k = divmod(y, 5)
    if hx == hy:
        steps = (1, 4) if hx < 5 else (2, 3)
        return (j-k) % 5 in steps
    if (hx < 5) == (hy < 5):
        return False
    if hx >= 5:
        return hoffman_edge(y, x)
    return k == (hx*(hy-5)+j) % 5


def small_thresholds():
    witnesses = (
        ([0,3,1], lambda a,b: (a-b) % 5 in (1,4)),
        ([(0,1),(0,2),(2,3),(0,3),(1,2)], lambda a,b: set(a).isdisjoint(b)),
        ([0,24,1,5,43,2,33,35,7,32,21,9,36,20], hoffman_edge),
    )
    for (vertices, edge), d, determinant in zip(witnesses,(2,3,7),(-1,-32,-931724928)):
        for i,j in combinations(range(len(vertices)),2):
            assert edge(vertices[i],vertices[j]) == square_sum(i+1,j+1)
        n = len(vertices)+1
        A = sp.Matrix(n,n,lambda i,j: int(square_sum(i+1,j+1)))
        gram = (d-1)*sp.eye(n)+sp.ones(n)-A-A*A
        assert gram.det(method='domain-ge') == determinant


# One line removes the centre and these degree-two middle vertices.
# Reverse order builds the full graph out of an empty set.
SCHEDULE = (
    (17,(19,)), (9,(16,27)), (10,(6,26)), (21,(15,28)), (5,(4,20)),
    (13,(12,23)), (14,(2,11,22)), (18,(7,)), (25,(24,)), (3,(1,)),
    (8,()), (29,()),
)


def triple_port_certificate():
    vertices = list(combinations(range(5),2))
    rows = [{j for j,v in enumerate(vertices) if set(u).isdisjoint(v)} for u in vertices]
    all_sites = set(range(10))
    hexagons = [set(c) for c in combinations(range(10),6)
                if all(len(rows[u] & set(c)) == 2 for u in c)]
    assert len(hexagons) == 10
    for a,b in combinations(all_sites,2):
        if b not in rows[a]:
            continue
        broken = [set(neighbors) for neighbors in rows]
        broken[a].remove(b)
        broken[b].remove(a)
        recovered = []
        for z,other in ((a,b),(b,a)):
            y1,y2 = broken[z]
            ring = all_sites - {z,other,y1,y2}
            assert ring in hexagons
            assert all(len(broken[v]&ring) == 2 for v in (other,y1,y2))
            assert len((broken[y1]&broken[y2])) == 1
            recovered.append(ring)
        assert recovered[0] != recovered[1]
    n,d = 3250,57
    nonstar = comb(n,3) - n*d//2*(n-2) + n*comb(d,2) - n*comb(d,3)
    assert nonstar == 20*266266000


def embedding_certificate():
    d = 57
    rows = square_rows(29)
    old = set(rows)
    floors = []
    for centre, mids in SCHEDULE:
        assert all(len(rows[m]&old) == 2 and centre in rows[m] for m in mids)
        old -= {centre,*mids}
        anchors = rows[centre] & old
        assert len(anchors) <= 1
        terminals = []
        for m in mids:
            terminal, = rows[m] & old
            terminals.append(terminal)
        assert len(set(terminals)) == len(terminals)
        assert anchors.isdisjoint(terminals)
        assert all((rows[u] & anchors) == set() for u in terminals)

        def uncovered(u):
            return {v for v in old if v != u and v not in rows[u]
                    and not rows[u] & rows[v] & old}

        forbidden_pairs = {tuple(sorted((u,v))) for u in terminals for v in uncovered(u)}
        if anchors:
            anchor, = anchors
            bound = d-len(rows[anchor]&old)-len(uncovered(anchor))-len(forbidden_pairs)
        else:
            twice_edges = sum(len(rows[u] & old) for u in old)
            bound = d*d+1-(d+1)*len(old)+twice_edges-(d-2)*len(forbidden_pairs)
        assert bound > 0
        floors.append(bound)
    assert not old
    assert floors == [16,4,17,20,28,39,2193,2635,2967,3079,3192,3250]
    assert prod(floors) == 13012905007307981605939200000
    assert all(square_sum(i,j) for i,j in combinations([6,19,30],2))
    seeds = {1,2,4,10,11,13,16,17,18}
    assert all(not rows[u] & seeds for u in seeds)
    for v in [15,14,12,7,9,8,3,6,5, *range(19,30)]:
        assert len(rows[v] & seeds) == 2
        seeds.add(v)
    assert seeds == set(rows)
    for n in range(18,30):
        assert sum(len(row) for row in square_rows(n).values()) == 4*(n-9)
    classes = {u: 0 for u in rows}
    for expected in [4,12,28,29]:
        signatures = {u:(classes[u], tuple(sorted(classes[v] for v in rows[u]))) for u in rows}
        codes = {key:i for i,key in enumerate(sorted(set(signatures.values())))}
        classes = {u:codes[key] for u,key in signatures.items()}
        assert len(codes) == expected
    print('Q29 lower bounds', floors, 'distinct induced copies',prod(floors))


def exterior_certificate(n, cycle):
    A = sp.Matrix(n,n, lambda i,j: int((i-j)%n in (1,n-1) if cycle else square_sum(i+1,j+1)))
    K = sp.Matrix([[57]]).col_join(-sp.ones(n,1)).row_join(sp.ones(1,n).col_join(-sp.eye(n)-A))
    x = sp.Symbol('x')
    char = sp.Poly(K.charpoly(x).as_expr(),x,modulus=13)
    factors = char.factor_list()[1]
    defects = []
    for poly, mult in factors:
        mat = sp.zeros(n+1)
        for coeff in poly.all_coeffs():
            mat = mat*K+int(coeff)*sp.eye(n+1)
        reduced = mat.to_DM().convert_to(GF(13))
        nullities = [(n+1-(reduced**j).rank())//poly.degree() for j in range(1,mult+1)]
        # Over the algebraic closure each root has these block lengths.
        tails = [nullities[0]]+[b-a for a,b in zip(nullities,nullities[1:])]+[0]
        blocks = [k for k in range(1,len(tails)) for _ in range(tails[k-1]-tails[k])]
        assert sum(blocks) == mult
        if poly.degree() == 1:
            root = -int(poly.TC()) % 13
            if root in (5,7):
                blocks = [k+1 for k in blocks]
            defects.extend((root,k) for k in blocks if k > 1)
        else:
            assert all(k == 1 for k in blocks)
    expected = {(13,True):[(5,2),(10,6),(10,7)], (25,True):[], (32,True):[(7,2)],
                (13,False):[(7,2),(8,2)], (25,False):[], (29,False):[]}
    assert sorted(defects) == expected[n,cycle]
    if not cycle and n in (25,29):
        assert sorted(f.degree() for f,_ in factors)==([4,4,18] if n==25 else [1,3,26])
    print('C' if cycle else 'Q', n, 'nontrivial exterior Jordan blocks', sorted(defects))


if __name__ == '__main__':
    small_thresholds()
    triple_port_certificate()
    embedding_certificate()
    for size in (13,25,32):
        exterior_certificate(size,True)
    for size in (13,25,29):
        exterior_certificate(size,False)
