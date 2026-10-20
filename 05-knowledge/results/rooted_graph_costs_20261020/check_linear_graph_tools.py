#!/usr/bin/env python3
"""Exact finite controls for rooted star tomography and finite probability patching.
Standard library only. Every check remains active under python -O.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb
import json

checks = 0

def check(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)

def subsets(mask):
    s = mask
    while True:
        yield s
        if not s:
            break
        s = (s - 1) & mask

def chi(x, mask):
    return -1 if (x & mask).bit_count() & 1 else 1

def star(x, mask):
    k = (x & mask).bit_count()
    return min(k, mask.bit_count() - k)

def star_fourier(n, k):
    return Q(sum((-1)**a * comb(k,a) * comb(n-k,b) * min(a+b,n-a-b)
                 for a in range(k+1) for b in range(n-k+1)), 2**n)

def root_cost(n, r, edges):
    return tuple(min(sum(((x>>u)&1)==((x>>v)&1) for u,v in edges)
                     for x in range(t,1<<n,1<<r)) for t in range(1<<r))

def fourier(table, mask):
    return Q(sum(v*chi(x,mask) for x,v in enumerate(table)),len(table))

def components(n,r,edges):
    unseen=set(range(r,n)); result=[]
    while unseen:
        component={min(unseen)}; unseen-=component
        stack=list(component)
        while stack:
            u=stack.pop()
            for a,b in edges:
                v=b if a==u else a if b==u else None
                if v in unseen:
                    unseen.remove(v); component.add(v); stack.append(v)
        roots={v for a,b in edges for u,v in [(a,b),(b,a)] if u in component and v<r}
        result.append((component,roots))
    return result

# Exact nonzero top coefficients, independently by direct cube and binomial sums.
coefficient_rows=[]
for n in range(2,17,2):
    m=n//2
    top=Q((-1)**m*comb(2*m-2,m-1),2**(2*m-1))
    direct=Q(sum(star(x,(1<<n)-1)*chi(x,(1<<n)-1) for x in range(1<<n)),1<<n)
    check(top==direct==star_fourier(n,n),f'top coefficient n={n}')
    coefficient_rows.append([n,str(top)])

# Arbitrary nonsymmetric rational laws: recover every even moment and the
# symmetrized atom law using ONLY the even-star expectations.
tomography_rows=[]
for n in range(1,9):
    weights=[(17*x*x+31*x+23)%97+1 for x in range(1<<n)]
    den=sum(weights); law=[Q(w,den) for w in weights]
    even=sorted((s for s in range(1<<n) if s.bit_count()%2==0),key=lambda s:(s.bit_count(),s))
    moments={0:Q(1)}
    for s in even[1:]:
        val=sum(p*star(x,s) for x,p in enumerate(law))
        lower=sum(star_fourier(s.bit_count(),t.bit_count())*moments[t]
                  for t in subsets(s) if t!=s and t.bit_count()%2==0)
        moments[s]=(val-lower)/star_fourier(s.bit_count(),s.bit_count())
        check(moments[s]==sum(p*chi(x,s) for x,p in enumerate(law)),f'moment {n},{s}')
    for x in range(1<<n):
        recovered=sum(moments[s]*chi(x,s) for s in even)/2**n
        check(recovered==(law[x]+law[x^((1<<n)-1)])/2,f'atom {n},{x}')
    tomography_rows.append([n,len(even)-1])

# All labelled edge subsets touching the interior, on prefix roots in 0..5
# vertices. Prefix roots lose no isomorphism type of labelled pattern.
pattern_count=0; low_boundary_count=0
for n in range(6):
    for r in range(n+1):
        allowed=[(u,v) for u in range(n) for v in range(u+1,n) if v>=r]
        for flags in range(1<<len(allowed)):
            edges=tuple(e for j,e in enumerate(allowed) if (flags>>j)&1)
            table=root_cost(n,r,edges)
            pattern_count+=1
            for x in range(1<<r):
                check(table[x]==table[x^((1<<r)-1)],'complement')
                # independent proper-colouring existence, not derived from min cost
                proper=any(all(((z>>u)&1)!=((z>>v)&1) for u,v in edges)
                           for z in range(x,1<<n,1<<r))
                check((table[x]==0)==proper,'zero-cost support')
            coeff={s:fourier(table,s) for s in range(1<<r)}
            for s in coeff:
                if s.bit_count()%2:
                    check(coeff[s]==0,'odd Fourier coefficient')
            if all(len(B)<=3 for _,B in components(n,r,edges)):
                low_boundary_count+=1
                for s in coeff:
                    if s.bit_count()>2:
                        check(coeff[s]==0,'three-port component compression')
            # Independently reduce every graph cost to a linear star-basis expansion.
            residual=dict(coeff); alpha={}
            even=sorted((s for s in coeff if s and s.bit_count()%2==0),
                        key=lambda s:(-s.bit_count(),s))
            for s in even:
                a=residual[s]/star_fourier(s.bit_count(),s.bit_count())
                alpha[s]=a
                for t in subsets(s):
                    residual[t]-=a*star_fourier(s.bit_count(),t.bit_count())
            for s,v in residual.items():
                if s:
                    check(v==0,'star-basis residual')
            for x,v in enumerate(table):
                check(v==residual[0]+sum(a*star(x,s) for s,a in alpha.items()),'star-basis value')

# Minimal four-root witness, including all proper marginals, and unit geometry.
R=4; allmask=15
laws={eps:[Q(1,8) if chi(x,15)==eps else Q(0) for x in range(16)] for eps in [1,-1]}
for s in range(15):
    check(sum(p*chi(x,s) for x,p in enumerate(laws[1]))==
          sum(p*chi(x,s) for x,p in enumerate(laws[-1])),'all proper marginals')
parity_costs={str(eps):sum(p*star(x,15) for x,p in enumerate(law)) for eps,law in laws.items()}
check(parity_costs=={'1':Q(3,2),'-1':Q(1)},'four-star gap')
for x in range(16):
    check(star(x,15)==Q(5,4)-Q(1,4)*sum(chi(x,(1<<i)|(1<<j)) for i,j in combinations(range(4),2))+Q(1,4)*chi(x,15),'four-star polynomial')
points=[(0,0),(1,0),(0,1),(-1,0),(0,-1)]
unit_edges=[(i,j) for i,j in combinations(range(5),2) if sum((points[i][k]-points[j][k])**2 for k in range(2))==1]
check(unit_edges==[(0,1),(0,2),(0,3),(0,4)],'induced unit-distance four-star')

# Gluing cost additivity for disjoint interiors/owned edges: all fragments
# on roots 0,1,2, with private centres 3 and 4, no boundary edges.
glue_count=0
for left in range(8):
    E=tuple((i,3) for i in range(3) if (left>>i)&1)
    for right in range(8):
        F=tuple((i,4) for i in range(3) if (right>>i)&1)
        a=root_cost(5,3,E); b=root_cost(5,3,F); c=root_cost(5,3,E+F)
        check(c==tuple(x+y for x,y in zip(a,b)),'glue disjoint interiors')
        glue_count+=1
# Shared-interior hostile: u-y and y-v separately cost zero; their union may not.
a=root_cost(3,2,((0,2),)); b=root_cost(3,2,((1,2),)); c=root_cost(3,2,((0,2),(1,2)))
check(a[1]+b[1]==0 and c[1]==1,'shared-interior hostile')
# Shared edge is counted once by union, twice by naive sum.
check(root_cost(2,2,((0,1),))[0]==1,'shared-edge union')
check(2*root_cost(2,2,((0,1),))[0]==2,'shared-edge double count')

# Tomography does not make one-sided star bounds sufficient for all gadgets.
short_path=root_cost(3,2,((0,2),(2,1)))
long_path=root_cost(4,2,((0,2),(2,3),(3,1)))
check(all(a+b==1 for a,b in zip(short_path,long_path)),'opposite path parity costs')
check(short_path[0]==0 and long_path[0]==1,'one-sided star bound hostile')

# Exact finite 2x2 minimax, all 81 payoffs in {0,1,2} and 25 rational targets.
def minmax_two(a):
    # a[constraint][strategy], choose probability p of strategy zero
    candidates={Q(0),Q(1)}
    den=(a[0][0]-a[0][1])-(a[1][0]-a[1][1])
    if den:
        p=Q(a[1][1]-a[0][1],den)
        if 0<=p<=1: candidates.add(p)
    return min(max(p*a[i][0]+(1-p)*a[i][1] for i in range(2)) for p in candidates)

def maxmin_dual(a):
    candidates={Q(0),Q(1)}
    den=(a[0][0]-a[1][0])-(a[0][1]-a[1][1])
    if den:
        p=Q(a[1][1]-a[1][0],den)
        if 0<=p<=1: candidates.add(p)
    return max(min(p*a[0][j]+(1-p)*a[1][j] for j in range(2)) for p in candidates)

minimax_cases=0
for raw in product(range(3),repeat=4):
    f=[raw[:2],raw[2:]]
    for c in product([Q(i,2) for i in range(5)],repeat=2):
        a=[[Q(f[i][j])-c[i] for j in range(2)] for i in range(2)]
        check(minmax_two(a)==maxmin_dual(a),'finite minimax')
        minimax_cases+=1

# Five-vertex path exposes bounded replication. Roots a,b have opposite signs;
# W={w}; private interiors s,t. Every pattern over {a,b} has cost <= |S|/3.
# Enumerate every allowed S and every touching edge subpattern, not merely full path.
path_edges=((0,3),(3,2),(2,4),(4,1))
path_J_patterns=0
for S in subsets((1<<2)|(1<<3)|(1<<4)):
    sverts=[v for v in range(2,5) if (S>>v)&1]
    allowed=[e for e in path_edges if all(v<2 or v in sverts for v in e) and any(v in sverts for v in e)]
    for F in range(1<<len(allowed)):
        E=[e for i,e in enumerate(allowed) if (F>>i)&1]
        cost=min(sum(((x>>u)&1)==((x>>v)&1) for u,v in E)
                 for x in range(2,32,4)) # bit0=0, bit1=1
        check(cost<=Q(len(sverts),3),'path J at C=1/3')
        path_J_patterns+=1
# Over N={a,b,w}, costs for each private two-edge path are complementary.
P=root_cost(5,3,((0,3),(3,2)))
Qtab=root_cost(5,3,((1,4),(4,2)))
check(P[2]+Qtab[2]==P[6]+Qtab[6]==1,'complementary extension costs')
check(min(P[2],P[6])==min(Qtab[2],Qtab[6])==0,'individually extendible')
check(all(min(a,b)<=Q(1+a+b,3) for a,b in product(range(2),repeat=2)),'one-copy tests pass')
check(min(2,2)>Q(1+2+2,3),'two-copy test fails')

# Round-up tariff bound used in finite successor theorem, all simplex grid
# vectors with denominator <=12, 1..4 constraints, 1..4 copy scales,
# and all size vectors in {0,1,2}. All arithmetic exact.
rounding_cases=0
for r in range(1,5):
    def compositions(total,length,prefix=()):
        if length==1:
            yield prefix+(total,)
        else:
            for h in range(total+1):
                yield from compositions(total-h,length-1,prefix+(h,))
    for den in range(1,13):
        for nums in compositions(den,r):
            for k in range(1,5):
                delta=[Q((k*a+den-1)//den)-Q(k*a,den) for a in nums]
                for sizes in product(range(3),repeat=r):
                    check(sum(d*s for d,s in zip(delta,sizes))<=sum(sizes)-min(sizes),'rounding tariff')
                    rounding_cases+=1

result={
 'status':'FINITE-EXACT controls; mathematical proofs are in README.md',
 'checks':checks,
 'top_star_coefficients':coefficient_rows,
 'tomography_root_and_probe_counts':tomography_rows,
 'all_patterns_vertices_0_to_5':pattern_count,
 'three_port_component_patterns':low_boundary_count,
 'four_star_expected_costs':{k:str(v) for k,v in parity_costs.items()},
 'disjoint_interior_gluing_pairs':glue_count,
 'minimax_2x2_cases':minimax_cases,
 'five_vertex_path_rooted_J_patterns':path_J_patterns,
 'rounding_tariff_cases':rounding_cases,
 'hostiles':['same proper marginals, different four-star cost','shared interior','shared boundary edge','one-sided star bounds do not control every gadget','one copy is insufficient for successor J'],
 'not_claimed':['a new finite chromatic bound','source Lean validation','the countable-colouring theorem proved here','geometric realizability of arbitrary copied gadgets']
}
print(json.dumps(result,indent=2,sort_keys=True))
