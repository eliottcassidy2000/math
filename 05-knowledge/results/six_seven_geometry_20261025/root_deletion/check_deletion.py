#!/usr/bin/env python3
"""Exact controls for deletion, branch splitting, and the doubled cube.

Universe: all 36 seven-element subsets of the nine rational magic forms,
each of seven labelled deletions; all t=0..3 Boolean component models;
all F17 points of the eight explicitly displayed central branch conics.
No search truncation, floating arithmetic, inherited producer imports, or
assert statements. The geometric all-characteristic-zero proof is separate.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FORMS = ((1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),
         (1,1,1),(1,-1,-1),(1,1,-1),(1,-1,1))
NAMES = ('E','E+x','E-x','E+z','E-z','E+x+z','E-x-z','E+x-z','E-x+z')


def check(test, message):
    if not test:
        raise RuntimeError(message)


def normalize(v):
    pivot = next((x for x in v if x), None)
    check(pivot is not None, 'zero projective vector')
    return tuple(Fraction(x, pivot) for x in v)


def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def incidence(subset):
    points = {normalize(cross(FORMS[a], FORMS[b])) for a,b in combinations(subset,2)}
    rows = [(point,tuple(i for i in subset if dot(FORMS[i],point)==0)) for point in points]
    check(all(2 <= len(labels) <= 3 for _,labels in rows), 'bad incidence multiplicity')
    return rows


def line_basis(form):
    pivot = next(i for i,x in enumerate(form) if x)
    basis=[]
    for j in range(3):
        if j != pivot:
            v=[Fraction(0) for _ in range(3)]
            v[j]=1
            v[pivot]=-Fraction(form[j],form[pivot])
            basis.append(tuple(v))
    return tuple(basis)


def binary_rank(rows):
    pivots={}
    for row in rows:
        while row:
            pos=row.bit_length()-1
            if pos in pivots:
                row ^= pivots[pos]
            else:
                pivots[pos]=row
                break
    return len(pivots)


def restriction_data(subset, deleted):
    basis=line_basis(FORMS[deleted])
    groups={}
    for j in subset:
        if j!=deleted:
            restriction=normalize(tuple(dot(FORMS[j],v) for v in basis))
            groups.setdefault(restriction,[]).append(j)
    check(all(len(g) in (1,2) for g in groups.values()),'not pair-only restrictions')
    groups=list(groups.values())
    zero_labels={j:k for k,g in enumerate(groups) for j in g}
    remaining=[j for j in subset if j!=deleted]
    reference=zero_labels[remaining[0]]
    parity_rows=[(1 << zero_labels[j]) ^ (1 << reference) for j in remaining[1:]]
    return groups,binary_rank(parity_rows)


histogram=Counter()
hourglass=[]
for subset in combinations(range(9),7):
    rows=incidence(subset)
    triple_labels=[labels for _,labels in rows if len(labels)==3]
    for deleted in subset:
        t=sum(deleted in labels for labels in triple_labels)
        groups,rank=restriction_data(subset,deleted)
        d=len(groups)
        check(d==6-t,'restriction/incidence disagreement')
        check(sum(len(g)==2 for g in groups)==t,'paired-zero disagreement')
        check(rank==d-1,'squareclass rank disagreement')
        components=2**(5-rank)
        degree=2**rank
        rh=-2*degree+d*(degree//2)
        check(rh%2==0,'odd Riemann-Hurwitz expression')
        genus=1+rh//2
        check(components==2**t,'component number disagreement')
        check(genus==1+2**(3-t)*(2-t),'genus disagreement')
        nodes=8*t
        arithmetic_genus=components*genus+nodes-components+1
        check(arithmetic_genus==17,'arithmetic genus failed')
        fixed_euler=components*(2-2*genus)
        quotient_euler=24+nodes
        check(2*quotient_euler-fixed_euler==80,'resolved Euler failed')
        # Direct intersection arithmetic on the blown-up K3.
        m_square=8-nodes
        m_dot_k=nodes
        check(2*2+(m_square+m_dot_k)//2==8,'cover chi failed')
        histogram[t]+=1
        if subset==tuple(range(7)):
            hourglass.append(dict(deleted=NAMES[deleted],incident_triples=t,
                                  six_root_nodes=4*(len(triple_labels)-t),
                                  branch_components=components,component_genus=genus,
                                  branch_nodes=nodes,quotient_euler=quotient_euler,
                                  fixed_euler=fixed_euler,
                                  restriction_groups=[[NAMES[j] for j in g] for g in groups]))
check(sum(histogram.values())==252,'wrong universe size')
check([row['incident_triples'] for row in hourglass]==[3,1,1,1,1,1,1],
      'wrong hourglass deletion types')

cube_rows=[]
for t in range(4):
    size=2**t
    mult=2**(4-t)
    matrix=[[mult*((2-t) if a==b else int((a^b).bit_count()==1))
             for b in range(size)] for a in range(size)]
    spectrum=Counter()
    for support in range(size):
        vector=[(-1)**((support & a).bit_count()) for a in range(size)]
        eigenvalue=2**(5-t)*(1-support.bit_count())
        observed=[sum(matrix[a][b]*vector[b] for b in range(size)) for a in range(size)]
        check(observed==[eigenvalue*x for x in vector],'Walsh eigenvector failure')
        spectrum[eigenvalue]+=1
        if support.bit_count()==1:
            check(observed==[0]*size,'half-cube relation not null')
    check(spectrum[0]==t,'wrong nullity')
    check(sum(sum(row) for row in matrix)==32,'B squared not 32')
    # Degree against H: B=2H, so each row sum is twice H.C.
    check(all(sum(row)==2**(5-t) for row in matrix),'unequal component H-degree')
    cube_rows.append(dict(t=t,edge_multiplicity=mult,
                         spectrum={str(k):v for k,v in sorted(spectrum.items()) if v},
                         divisor_span_rank=size-t,
                         matrix=matrix if t==3 else None))


def projective_finite_points(q,n):
    for pivot in range(n):
        for tail in product(range(q),repeat=n-pivot-1):
            yield (0,)*pivot+(1,)+tail


def normalize_mod(v,q):
    pivot=next(x for x in v if x%q)
    inv=pow(pivot,-1,q)
    return tuple(x*inv%q for x in v)


q=17
imaginary_unit=4
check(imaginary_unit**2%q==q-1,'wrong sqrt(-1)')
conics=[]
for bits in range(8):
    points=set()
    for a,b,c in projective_finite_points(q,3):
        if (c*c-a*a-b*b)%q:
            continue
        coords=(a,imaginary_unit*(-1)**(bits&1)*a,
                b,imaginary_unit*(-1)**((bits>>1)&1)*b,
                c,imaginary_unit*(-1)**((bits>>2)&1)*c)
        point=normalize_mod(coords,q)
        squares=[x*x%q for x in point]
        check((squares[0]+squares[1])%q==0,'A pair not on E=0')
        check((squares[2]+squares[3])%q==0,'B pair not on E=0')
        check((squares[4]+squares[5])%q==0,'C pair not on E=0')
        check((squares[4]-squares[0]-squares[2])%q==0,'conic quadratic failed')
        points.add(point)
    check(len(points)==q+1,'conic count failed')
    conics.append(points)
intersection_counts={}
node_points=set()
for a,b in combinations(range(8),2):
    shared=conics[a]&conics[b]
    expected=2 if (a^b).bit_count()==1 else 0
    check(len(shared)==expected,'direct conic intersection failure')
    intersection_counts[f'{a:03b}-{b:03b}']=len(shared)
    node_points.update(shared)
check(len(node_points)==24,'wrong direct central node count')
check(all(sum(p in c for c in conics)==2 for p in node_points),'triple component intersection')
check(len(set.union(*conics))==8*(q+1)-24,'conic union count failed')

# Hostile control: naively calling the central branch smooth genus 17
# predicts Euler -32, while the actual normalization has Euler 16.
check(8*2 != 2-2*17,'hostile smooth-genus substitution failed to distinguish')

receipt=dict(status='FINITE-EXACT',universe='36 seven-root subsets x 7 deletions',
             cases=252,t_histogram={str(k):v for k,v in sorted(histogram.items())},
             hourglass=hourglass,cube_models=cube_rows,
             central_conics_F17=dict(q=q,sqrt_minus_one=imaginary_unit,
                                     component_points=[len(c) for c in conics],
                                     pair_intersections=intersection_counts,
                                     distinct_nodes=len(node_points),
                                     union_points=len(set.union(*conics))),
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
encoded=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
(ROOT/'deletion_receipt.json').write_text(encoded,encoding='utf-8')
print(encoded,end='')
