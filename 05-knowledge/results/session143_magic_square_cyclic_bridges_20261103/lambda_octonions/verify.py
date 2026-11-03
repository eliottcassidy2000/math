#!/usr/bin/env python3
"""Portable exact checker: stdlib only; all checks survive python -O.
Run: python3 verify.py [--mutate-sign]
The mutation is an intentionally failing hostile control, not a valid output.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
import json, sys

def require(ok, why):
    if not ok:
        raise RuntimeError(why)

@lru_cache(None)
def cd(n, a, b):
    """Recursive (a,b)(c,d)=(ac-db*,a*d+cb); signed basis output."""
    if n == 0:
        return 1, 0
    h = 1 << (n-1)
    i, j = a % h, b % h
    if a < h and b < h:
        return cd(n-1, i, j)
    if a < h <= b:
        s, k = cd(n-1, i, j)
        return s * (-1 if i else 1), k+h
    if b < h <= a:
        s, k = cd(n-1, j, i)
        return s, k+h
    s, k = cd(n-1, j, i)
    return -s * (-1 if i else 1), k

def bit(x, i):
    return (x >> i) & 1

def exponent(a, b):
    x = [bit(a, i) for i in range(3)]
    y = [bit(b, i) for i in range(3)]
    ans = (sum(x[i]*y[j] for i in range(3) for j in range(i,3))
           + x[1]*x[2]*y[0] + x[0]*x[2]*y[1] + x[0]*x[1]*y[2]) % 2
    if '--mutate-sign' in sys.argv and (a,b) == (1,2):
        ans ^= 1
    return ans

def determinant(a,b,c):
    return sum(bit(a,i)*bit(b,j)*bit(c,k)
               for i,j,k in permutations(range(3))) % 2

def assoc(a,b,c):
    s1, ab = cd(3,a,b)
    s2, abc1 = cd(3,ab,c)
    s3, bc = cd(3,b,c)
    s4, abc2 = cd(3,a,bc)
    require(abc1 == abc2 == a^b^c, 'associator output label')
    return s1*s2-s3*s4

# Independent Boolean-polynomial certificate of delta f = det.
def add(*ps):
    z = set()
    for p in ps:
        z.symmetric_difference_update(p)
    return z

def mul(p,q):
    z = set()
    for a in p:
        for b in q:
            m = a | b
            if m in z: z.remove(m)
            else: z.add(m)
    return z

def fp(x,y):
    return add(*(mul(x[i],y[j]) for i in range(3) for j in range(i,3)),
               mul(mul(x[1],x[2]),y[0]), mul(mul(x[0],x[2]),y[1]),
               mul(mul(x[0],x[1]),y[2]))

X = [{1 << i} for i in range(3)]
Y = [{1 << (i+3)} for i in range(3)]
Z = [{1 << (i+6)} for i in range(3)]
XY = [add(x,y) for x,y in zip(X,Y)]
YZ = [add(y,z) for y,z in zip(Y,Z)]
delta_poly = add(fp(X,Y),fp(XY,Z),fp(Y,Z),fp(X,YZ))
det_poly = add(*(mul(mul(X[i],Y[j]),Z[k]) for i,j,k in permutations(range(3))))
require(delta_poly == det_poly and len(det_poly) == 6, 'symbolic coboundary identity')

for a,b in product(range(8),repeat=2):
    require(cd(3,a,b) == ((-1)**exponent(a,b), a^b), 'cochain vs Cayley-Dickson')
for a,b,c in product(range(8),repeat=3):
    df = exponent(a,b)^exponent(a^b,c)^exponent(b,c)^exponent(a,b^c)
    require(df == determinant(a,b,c), 'basis coboundary determinant')
    require((assoc(a,b,c) != 0) == bool(determinant(a,b,c)), 'associator support')
    require(assoc(a,b,c) == -assoc(b,a,c) == -assoc(a,c,b), 'alternating associator')
require(assoc(1,2,4) == 2, 'nonassociative hostile witness')

subspaces = []
for mask in range(256):
    W = frozenset(x for x in range(8) if mask >> x & 1)
    if 0 in W and all(a^b in W for a in W for b in W):
        subspaces.append(W)
require(Counter(len(W) for W in subspaces) == {1:1,2:7,4:7,8:1}, 'subspace universe')
for W in subspaces:
    require(all(cd(3,a,b)[1] in W for a in W for b in W), 'subalgebra closure')
    is_associative = all(assoc(a,b,c) == 0 for a,b,c in product(W,repeat=3))
    require(is_associative == (len(W)<=4), 'rank <= 2 associative iff')

def span(gens):
    W={0}
    for g in gens:
        W |= {x^g for x in list(W)}
    return frozenset(W)

cubes=set()
for basis in combinations(range(1,8),3):
    if determinant(*basis):
        cube=frozenset(span([basis[i] for i in range(3) if mask>>i&1]) for mask in range(8))
        require(len(cube)==8,'coordinate cube distinct vertices')
        require(sum(W<V and len(V)==2*len(W) for W in cube for V in cube)==12,'cube edges')
        cubes.add(cube)
require(len(cubes)==28,'28 basis-independent coordinate cubes')

# Exact diagonal automorphisms and fixed-subalgebra dictionary.
def chi(u,x):
    return (-1)**((u&x).bit_count()%2)
for u,a,b in product(range(8),repeat=3):
    require(chi(u,a)*chi(u,b)==chi(u,a^b),'character automorphism')
for H in subspaces:
    fixed=frozenset(x for x in range(8) if all(chi(u,x)==1 for u in H))
    require(fixed in subspaces and len(fixed)*len(H)==8,'fixed support annihilator')

# Q(sqrt(2),sqrt(3),sqrt(5)) has the same label/action but different product.
def q(a,b):
    ans=1
    for i,d in enumerate((2,3,5)):
        if bit(a,i) and bit(b,i):ans*=d
    return ans
for a,b,c in product(range(8),repeat=3):
    require(q(a,b)==q(b,a),'field product symmetry')
    require(q(a,b)*q(a^b,c)==q(b,c)*q(a,b^c),'field cocycle associativity')
require(q(1,1)==2 and cd(3,1,1)==(-1,0),'field vs octonion square hostile')

# All 128 sign gauges leave the associator unchanged; precisely eight fix table.
gauge_tables=set();kernel=[]
for mask in range(128):
    g=lambda a: ((mask>>(a-1))&1) if a else 0
    fg=lambda a,b: exponent(a,b)^g(a)^g(b)^g(a^b)
    table=tuple(fg(a,b) for a,b in product(range(8),repeat=2))
    gauge_tables.add(table)
    if table==tuple(exponent(a,b) for a,b in product(range(8),repeat=2)):kernel.append(mask)
    for a,b,c in product(range(8),repeat=3):
        require(fg(a,b)^fg(a^b,c)^fg(b,c)^fg(a,b^c)==determinant(a,b,c),'gauge cannot remove associator')
require(len(gauge_tables)==16 and len(kernel)==8,'gauge counts')

# The given nine-cycle is translation +4 in Z/9, with representative set 2..10.
cycle=[4,8,3,7,2,6,10,5,9]
require(len(set(x%9 for x in cycle))==9,'nine distinct residues')
require(all((cycle[(i+1)%9]-cycle[i])%9==4 for i in range(9)),'translation +4 mod9')
require([x for x in range(8) if (x if 9%2 else 0)==0]==[0],'C9 to F2^3 homomorphism trivial')

# Lucas Lambda_5 is a different object: 11 vertices, not 8.
words=[x for x in range(32) if all(not(bit(x,i) and bit(x,(i+1)%5)) for i in range(5))]
edges=[(x,y) for x,y in combinations(words,2) if (x^y).bit_count()==1]
require(len(words)==11 and len(edges)==15,'Lucas cube order-five universe')

# Exact perimeter-root lift: clockwise cells, with opposite labels XOR 4.
perimeter = [0,1,2,5,8,7,6,3]
cell_to_label = {cell:a for a,cell in enumerate(perimeter)}
lines = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
r=4
left_table=[cd(3,r,a) for a in range(8)]
for a in range(8):
    s,b=left_table[a]
    t,c=left_table[b]
    require(c==a and s*t==-1,'L_r squared is minus identity')
    require(perimeter[b]==8-perimeter[a],'left permutation is halfturn')
require(sorted(b for s,b in left_table)==list(range(8)),'signed permutation preserves norm')
rotated_lines=[tuple(sorted(8-cell for cell in line)) for line in lines]
require(sorted(rotated_lines)==sorted(tuple(sorted(line)) for line in lines),'all eight magic lines permuted')
require(cd(3,r,0)==(1,r),'left multiplication is not unital algebra automorphism')
# Symbolic square-root lift forgets signs: squared coordinate y_b^2 = x_(b XOR r)^2.
# Test distinct integer squares, positive root convention, and nonmagic hostile.
roots=list(range(1,10))
lifted=[0]*9
lifted[4]=roots[4]
for a,cell in enumerate(perimeter):
    sign,b=left_table[a]
    lifted[perimeter[b]]=sign*roots[cell]
require([x*x for x in lifted]==[roots[8-cell]**2 for cell in range(9)],'square map intertwines halfturn')
require(sum(x*x for x in lifted)==sum(x*x for x in roots),'full nine-root norm')
require(any(x<0 for x in lifted),'positive-root convention is not preserved')

def magic(values):
    sums=[sum(values[cell] for cell in line) for line in lines]
    return len(set(sums))==1
require(not magic([x*x for x in roots]) and not magic([x*x for x in lifted]),'lift cannot turn hostile into magic')
require(magic([1]*9),'nonvacuity degenerate square-magic control')
loshu=[8,1,6,3,5,7,4,9,2]
translation_magic=[]
translation_hostiles={}
for v in range(8):
    moved=loshu[:]
    for a,cell in enumerate(perimeter):
        moved[perimeter[a^v]]=loshu[cell]
    if magic(moved):translation_magic.append(v)
    else:translation_hostiles[v]={'values':moved,'line_sums':[sum(moved[cell] for cell in line) for line in lines]}
require(translation_magic==[0,4],'not all octonion label translations preserve magic')

# Norm and opposite pairs are necessary but not the full magic constraints.
paired_roots=[13,23,35,47,65,79,85,89,91]
paired_squares=[x*x for x in paired_roots]
center=paired_squares[4]
require(all(paired_squares[c]+paired_squares[8-c]==2*center for c in (0,1,2,3)), 'opposite-pair hostile assumptions')
require(sum(paired_squares[c] for c in perimeter)==8*center, 'perimeter norm hostile assumption')
paired_sums=[sum(paired_squares[c] for c in line) for line in lines]
require(paired_sums==[1923,12675,23427,9603,12675,15747,12675,12675], 'paired hostile exact line sums')
require(not magic(paired_squares), 'norm plus opposite pairs is insufficient')

certificate={
 'status':'VERIFIED exact finite identities; accompanying proof states uniform claims',
 'cochain_pairs':64,'associator_triples':512,
 'nonzero_associator_ordered_triples':sum(assoc(*t)!=0 for t in product(range(8),repeat=3)),
 'nonzero_associator_unordered_imaginary_triples':sum(assoc(*t)!=0 for t in combinations(range(1,8),3)),
 'associator_coefficient_census':dict(sorted(Counter(assoc(*t) for t in product(range(8),repeat=3)).items())),
 'symbolic_coboundary_monomial_masks':sorted(delta_poly),
 'subspace_support_cardinality_census':dict(sorted(Counter(len(W) for W in subspaces).items())),
 'coordinate_boolean_cubes':len(cubes),
 'diagonal_character_automorphisms':8,
 'gauge_count':128,'gauge_orientations':len(gauge_tables),'gauge_kernel_count':len(kernel),
 'field_cocycle_triples':512,
 'hostile_nonassociativity':{'labels':[1,2,4],'left':[1,7],'right':[-1,7],'difference':[2,7]},
 'nine_cycle':cycle+[cycle[0]],'nine_cycle_step_mod9':4,
 'Lucas_Lambda5_vertices':len(words),'Lucas_Lambda5_edges':len(edges),
 'magic_perimeter_clockwise_cells':perimeter,'magic_left_unit_label':r,
 'magic_left_signed_table':left_table,'magic_lines_permuted':8,'magic_lift_order':4,
 'magic_quotient_halfturn_order':2,
 'magic_XOR_translations_preserving_all_lines':[0,4],
 'magic_other_translation_hostiles':translation_hostiles,
 'norm_opposite_pairs_hostile_roots':paired_roots,
 'norm_opposite_pairs_hostile_line_sums':paired_sums,
 'norm_opposite_pairs_hostile_perimeter_norm':8*center,
 'scope':'No interpretation of user lambda-cube term is assumed; no LRC, magic-square, Baby Monster, or n5 spectral theorem is inferred.'}
print(json.dumps(certificate,indent=2,sort_keys=True))
