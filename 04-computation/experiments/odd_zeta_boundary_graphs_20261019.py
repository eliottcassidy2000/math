#!/usr/bin/env python3
"""Exact finite tests for odd-zeta boundary graph geometry; Python 3 + SymPy.
Universe: all 2-subsets of five labels; every pair of r-fold corners, r<=3;
all labelled internal skeletons for the four displayed stable-tree profiles.
No floating-point arithmetic, no network, no numerical zeta inference.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import factorial, prod
import sympy as sp

V = list(combinations(range(1, 6), 2))
E = {frozenset((a,b)) for a,b in combinations(V,2) if set(a).isdisjoint(b)}
def adjacent(a,b): return frozenset((a,b)) in E
A = sp.Matrix(10,10,lambda i,j: int(adjacent(V[i],V[j])))
B = sp.Matrix(5,10,lambda i,j: int(i+1 in V[j]))
Q = A-sp.eye(10)
H = sp.ones(5)/4-sp.eye(5)
assert len(E)==15 and all(sum(A[i,j] for j in range(10))==3 for i in range(10))
assert Q==sp.ones(10)-B.T*B==B.T*H*B
assert A.eigenvals()=={3:1,1:5,-2:4}
assert Q.eigenvals()=={2:1,0:5,-3:4}
assert B.rank()==Q.rank()==5
basis=[(1,2),(1,3),(2,3),(1,4),(1,5)]
cols=[V.index(x) for x in basis]
T=B[:,cols]
assert abs(T.det())==2
coordinates=T.inv()*B
assert all(x.q==1 for x in coordinates)
Q0=Q.extract(cols,cols)
assert Q0.det()==1 and Q==coordinates.T*Q0*coordinates
relations=[]
for a,b,c,d in combinations(range(1,6),4):
    for plus,minus in [([(a,b),(c,d)],[(a,c),(b,d)]),
                       ([(a,b),(c,d)],[(a,d),(b,c)])]:
        row=[0]*10
        for x in plus: row[V.index(x)]+=1
        for x in minus: row[V.index(x)]-=1
        relations.append(row)
R=sp.Matrix(relations)
assert R.rank()==5 and R*B.T==sp.zeros(10,5) and R*Q==sp.zeros(10,10)
print('Petersen/K5 intersection identity: Q=A-I=J-B^T B=B^T((J5/4)-I)B')
print('spectrum Q:',sorted(Q.eigenvals().items()))
print('integral basis:',basis,'incidence determinant:',T.det(),'Gram determinant:',Q0.det())
print('Keel relations rank:',R.rank(),'all annihilate B^T and Q')
print('integral Gram matrix:',Q0.tolist())
F=sp.ones(5)-sp.eye(5)
kappa=2*sp.ones(5,1)
assert F.T*H*B==B and F.T*H*F==sp.ones(5)-sp.eye(5)
assert (kappa.T*H*kappa)[0]==5
simple=sp.Matrix.hstack(*[sp.eye(5)[:,i]-sp.eye(5)[:,i+1] for i in range(4)])
assert simple.T*H*simple==sp.Matrix(4,4,lambda i,j: -2 if i==j else int(abs(i-j)==1))
roots=[u for u in product((-1,0,1),repeat=5) if sum(u)==0 and sum(x*x for x in u)==2]
assert len(roots)==20
fiber_colors=Counter()
for edge in E:
    u,v=tuple(edge); omitted=next(iter(set(range(1,6))-set(u)-set(v)))
    assert B[:,V.index(u)]+B[:,V.index(v)]==F[:,omitted-1]
    fiber_colors[omitted]+=1
assert fiber_colors==Counter({i:3 for i in range(1,6)})
print('Geometric K5 vertices: five forgetful fibre classes F_i=1-e_i; intersections F^T H B=B.')
print('Petersen edge colors/reducible fibres:',dict(sorted(fiber_colors.items())))
print('E4=A4: anticanonical square5; perpendicular sum-zero lattice;20 roots; negative A4 Cartan verified.')


# Explicit standard Petersen model and Kuratowski witnesses.
a=[(1,2),(3,4),(1,5),(2,3),(4,5)]
b=[(3,5),(2,5),(2,4),(1,4),(1,3)]
assert set(a+b)==set(V)
assert all(adjacent(a[i],a[(i+1)%5]) and adjacent(a[i],b[i]) and
           adjacent(b[i],b[(i+2)%5]) for i in range(5))
branch5=[{a[i],b[i]} for i in range(5)]
assert all(any(adjacent(u,v) for u in branch5[i] for v in branch5[j])
           for i,j in combinations(range(5),2))
branch33=[{a[2],a[1]},{b[3],b[0]},{b[4],a[4]}, {a[3]},{b[1]},{b[2]}]
assert len(set.union(*branch33))==9 and a[0] not in set.union(*branch33)
q33={(i,j) for i,j in combinations(range(6),2)
     if any(adjacent(u,v) for u in branch33[i] for v in branch33[j])}
assert q33=={(i,j) for i in range(3) for j in range(3,6)}
print('K5 minor: contract five disjoint spokes a_i-b_i.')
print('K3,3 minor: delete a0; branch sets',[[a.index(v) if v in a else 5+b.index(v) for v in s] for s in branch33])
print('Petersen contains no K5 or K3,3 subgraph: degree3 and triangle-free; a K3,3 subgraph would use all degrees and be a proper connected component.')

for r in range(1,4):
    corners=list(product(range(10),repeat=r))
    intersection_dimensions=Counter()
    quotient_edges=set()
    pmap={V.index(v):i for i,S in enumerate(branch5) for v in S}
    for x,y in combinations(corners,2):
        changes=sum(xi!=yi for xi,yi in zip(x,y))
        if all(xi==yi or A[xi,yi] for xi,yi in zip(x,y)):
            intersection_dimensions[r-changes]+=1
            if changes==1:
                qx=tuple(pmap[i] for i in x); qy=tuple(pmap[i] for i in y)
                if qx!=qy: quotient_edges.add(frozenset((qx,qy)))
    cart=intersection_dimensions[r-1]
    full=sum(intersection_dimensions.values())
    assert cart==3*r*10**r//2
    assert full==10**r*(4**r-1)//2
    assert len(quotient_edges)==2*r*5**r
    positive=(5**r+(-3)**r)//2; negative=(5**r-(-3)**r)//2
    t=sp.Symbol('t')
    middle=sp.expand((1+5*t+t*t)**r).coeff(t,r)
    fpoly=sp.Poly((1+10*t+15*t*t)**r,t)
    red_euler=-fpoly.eval(-1)
    assert red_euler==-6**r
    print(f'r={r}: corners={10**r}; Cartesian edges={cart}; all nonempty pair intersections={full}; dimension histogram={dict(sorted(intersection_dimensions.items()))}')
    print(f'  actual boundary dual complex: {10*r} vertices, top dimension={2*r-1}, reduced Euler={red_euler}, top Betti={6**r}')
    print(f'  corner Gram: rank={5**r}, nullity={10**r-5**r}, signature=({positive},{negative}); full middle rank={middle}')
    print(f'  K5 Cartesian-power minor: vertices={5**r}, edges={len(quotient_edges)}')
    if r<=2:
        Qt=Q if r==1 else sp.kronecker_product(Q,Q)
        assert Qt.rank()==5**r
        ev=Qt.eigenvals()
        assert sum(m for x,m in ev.items() if x>0)==positive
        assert sum(m for x,m in ev.items() if x<0)==negative
        print('  explicit exact tensor spectrum:',sorted(ev.items()))

# Stable tree profile formula, independently counted via internal skeletons.
def decode_prufer(word,k):
    degree=[1]*k
    for u in word: degree[u]+=1
    edges=[]
    for u in word:
        v=next(i for i,d in enumerate(degree) if d==1)
        edges.append((u,v)); degree[u]-=1; degree[v]-=1
    last=[i for i,d in enumerate(degree) if d==1]
    edges.append(tuple(last))
    return edges

def profile_count(vals):
    k=len(vals); n=2+sum(d-2 for d in vals)
    denominator=prod(factorial(m) for m in Counter(vals).values())
    closed=Fraction(factorial(n+k-2),denominator*prod(factorial(d-1) for d in vals))
    total=Fraction(0)
    skeletons=[[]] if k==1 else (decode_prufer(w,k) for w in product(range(k),repeat=k-2))
    for edges in skeletons:
        deg=Counter(v for e in edges for v in e)
        leaves=[vals[i]-deg[i] for i in range(k)]
        if min(leaves)<0: continue
        total+=Fraction(factorial(n),denominator*prod(factorial(m) for m in leaves))
    assert total==closed and total.denominator==1
    return n,int(total)

for vals,expected in [([5,5],35),([5,4,5],5775),([5,3,3,5],51975),([5,5,5,3],525525)]:
    n,count=profile_count(vals); assert count==expected
    print(f'stable profile {sorted(vals)}: N={n}, codim={len(vals)-1}, dimension={sum(d-3 for d in vals)}, labelled strata={count}')

def verify_tree(n,marks,edges,expected):
    assert set.union(*map(set,marks))==set(range(1,n+1))
    assert sum(map(len,marks))==n and len(edges)==len(marks)-1
    deg=Counter(v for e in edges for v in e)
    vals=[len(marks[i])+deg[i] for i in range(len(marks))]
    assert vals==expected and min(vals)>=3
    assert n==2+sum(d-2 for d in vals)
    assert n-3-(len(edges))==sum(d-3 for d in vals)
    return vals
print('explicit N8:',verify_tree(8,[range(1,5),range(5,9)],[(0,1)],[5,5]))
print('explicit N10:',verify_tree(10,[range(1,5),[5],[6],range(7,11)],[(0,1),(1,2),(2,3)],[5,3,3,5]))
print('explicit N12:',verify_tree(12,[[],range(1,5),range(5,9),range(9,13)],[(0,1),(0,2),(0,3)],[3,5,5,5]))
assert factorial(12)//(factorial(4)**3*factorial(3))==5775
print('N12 central-trivalent star subfamily: 5775 strata (of 525525 total profile strata)')
for n in (8,10,12):
    r=(n-2)//3
    assert 3*r+2<=n<3*(r+1)+2
    print(f'N={n}, weight cap={n-3}: maximum number of Mbar0,5 stable-tree factors={r}')
# Valuation sidecar for Brown's actual odd-zeta cellular configurations.
# It is precisely this data that a bare Petersen graph forgets.
def runs(order,S):
    return sum(order[i] in S and order[i-1] not in S for i in range(len(order)))
def odd_order(m):
    out=[2*m]
    for k in range(2,m+1):
        out.append(k)
        if k<m: out.append(2*m-k+1)
    return out+[1,m+1]
def stable_splits(n):
    # Unique representative containing label 1; both sides have size at least 2.
    for mask in range(1<<(n-1)):
        S={1}|{i+2 for i in range(n-1) if mask>>i&1}
        if 2<=len(S)<=n-2: yield S
for n,selected in [(8,[set(range(1,5))]),
                   (10,[set(range(1,5)),set(range(1,6)),set(range(1,7))]),
                   (12,[set(range(1,5)),set(range(5,9)),set(range(9,13))])]:
    delta=list(range(1,n+1)); sigma=odd_order(n//2)
    assert sorted(sigma)==delta
    allsplits=list(stable_splits(n)); assert len(allsplits)==2**(n-1)-n-1
    hist=Counter();facets=0
    for S in allsplits:
        r0,r1=runs(delta,S),runs(sigma,S)
        # Independent form of Brown's ID definition: half of same-side edges.
        ID=lambda order: sum((order[i] in S)==(order[i-1] in S) for i in range(n))
        assert ID(delta)==n-2*r0 and ID(sigma)==n-2*r1
        slope,constant=r1-r0,r1-2
        assert 2*slope==ID(delta)-ID(sigma)
        assert 2*constant==n-4-ID(sigma)
        if r0==1:
            facets+=1; assert r1>=2 and slope>=1 and constant>=0
        if n==8 and len(S)==4: hist[(slope,constant)]+=1
    assert facets==n*(n-3)//2
    node_orders=[(runs(sigma,S)-runs(delta,S),runs(sigma,S)-2) for S in selected]
    assert all(a>=1 and b>=0 for a,b in node_orders)
    print(f'Brown N={n}, odd order={sigma}: checked {len(allsplits)} stable splits; {facets} cell facets all convergent')
    print('  chosen Petersen-product stratum node orders (slope,constant):',node_orders)
    print('  therefore every ordinary logarithmic residue along these nodes is zero for integral exponent k>=0')
    if hist: print('  N8 all 35 4|4 Petersen-square divisors, order histogram:',dict(sorted(hist.items())))
print('ALL CHECKS PASSED')
