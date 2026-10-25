"""Exact P1(F7) -> affine eight -> marked conic-component comparison.

Universe: all 70 four-subsets; all 336 projective transformations; all 168
square-determinant transformations. Integer arithmetic only. No asserts:
normal and optimized Python execute the same validation gates.
"""
import itertools as it
import json

P = tuple(range(8))                 # 7 denotes infinity
LABEL = (1, 2, 4, 3, 6, 7, 5, 0)  # chosen U design; not an additive F7 map
D = (7, 3, 5)                     # inherited conic-action direction basis

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def dot(a, b):
    return (a & b).bit_count() % 2

def hp(t):
    return (1, 0) if t == 7 else (t, 1)

def det(v, w):
    return (v[0]*w[1] - v[1]*w[0]) % 7

def cross(q):
    a, b, c, d = map(hp, q)
    den = det(d, a)*det(c, b) % 7
    return det(c, a)*det(d, b)*pow(den, -1, 7) % 7

def transform(m, t):
    a,b,c,d = m
    x,y = hp(t)
    X,Y = (a*x+b*y)%7, (c*x+d*y)%7
    return X*pow(Y,-1,7)%7 if Y else 7

def normalized(m):
    lead = next(x for x in m if x)
    return tuple(x*pow(lead,-1,7)%7 for x in m)

def block_image(p, block):
    return frozenset(p[x] for x in block)

def cube_edges():
    return {frozenset((x,x^d)) for x in P for d in (1,2,4)}

def linear(cols, x):
    out = 0
    for j in range(3):
        if x & (1<<j):
            out ^= cols[j]
    return out

def affine_data(p):
    """Permutation p acts on the residue/P1 labels; return binary affine data."""
    f = {LABEL[t]: LABEL[p[t]] for t in P}
    shift = f[0]
    cols = tuple(f[1<<j]^shift for j in range(3))
    need(all(f[x] == (linear(cols,x)^shift) for x in P), 'not affine')
    return cols, shift, f

def main():
    matrices = sorted({normalized(m) for m in it.product(range(7), repeat=4)
                       if (m[0]*m[3]-m[1]*m[2])%7})
    maps = [tuple(transform(m,t) for t in P) for m in matrices]
    ps = [p for m,p in zip(matrices,maps)
          if (m[0]*m[3]-m[1]*m[2])%7 in (1,2,4)]
    need(len(maps)==len(set(maps))==336 and len(ps)==168,'group orders')
    fours = list(map(frozenset,it.combinations(P,4)))
    equi = {q for q in fours if cross(sorted(q)) in (3,5)}
    A = {block_image(p,(7,0,1,3)) for p in ps}
    B = {block_image(p,(7,0,1,5)) for p in ps}
    need(len(equi)==28 and len(A)==len(B)==14 and not A&B and A|B==equi,
         'equianharmonic orbit partition')
    for blocks in (A,B):
        need(all(Pset-q in blocks for q in blocks for Pset in [frozenset(P)]),
             'complement closure')
        need(all(sum(set(t)<=q for q in blocks)==1 for t in it.combinations(P,3)),
             'Steiner triple coverage')
    U = {frozenset((t,(t+1)%7,(t+3)%7)) for t in range(7)}
    need({q-{7} for q in A if 7 in q}==U,'derived Fano design')
    affine_planes = {frozenset(x for x in P if dot(s,x)==b)
                     for s in range(1,8) for b in (0,1)}
    need({frozenset(LABEL[t] for t in q) for q in A}==affine_planes,
         'affine identification')
    affine = [affine_data(p) for p in ps]
    need({p[7] for p in ps}==set(P),'projective transitivity')
    need(not any(all(p[t]==t for p in ps) for t in P),'common fixed point hostile')
    need(len({cols for cols,shift,f in affine})==168,'faithful linear-part map')
    inversion = tuple(transform((0,6,1,0),t) for t in P)  # -1/z
    ci,ai,fi = affine_data(inversion)
    need(not any(fi[x]==x for x in P),'fixed-point-free affine involution')
    need(ai not in {linear(ci,x)^x for x in P},'nontrivial affine gain')
    # Specific component marking: choose dot-duality between inertia and
    # the simplex subgroup, then use its already proved component action.
    component = tuple(sum(dot(LABEL[t],D[j])<<j for j in range(3)) for t in P)
    need(len(set(component))==8 and component[7]==0 and component[0]==7,
         'specific conic marking')
    need({frozenset(component[t] for t in q) for q in A}==affine_planes,
         'component block transport')
    # Raw equianharmonic blocks are TWO designs: forgetting the orbit makes
    # every triple have two completions, not a Steiner system.
    need(all(sum(set(t)<=q for q in equi)==2 for t in it.combinations(P,3)),
         'lost-chirality negative control')
    # The other chirality is not the same affine design under this marking.
    need({frozenset(LABEL[t] for t in q) for q in B} != affine_planes,
         'marking/chirality negative control')
    # Exact matrix/Walsh check, and the additional cube frame it requires.
    edges=cube_edges()
    M=[[(-2 if x==y else 2 if frozenset((x,y)) in edges else 0)
        for y in P] for x in P]
    eigen=[]
    for s in P:
        v=[(-1)**dot(s,x) for x in P]
        lam=4*(1-s.bit_count())
        need(all(sum(M[x][y]*v[y] for y in P)==lam*v[x] for x in P),'Walsh spectrum')
        eigen.append(lam)
    # Preserve all affine planes, destroy one genuine conic intersection.
    shear=tuple(x ^ (((x>>1)&1)<<0) for x in P)
    need({block_image(shear,q) for q in affine_planes}==affine_planes,'shear planes')
    need(M[0][2]==2 and M[shear[0]][shear[2]]==0,'shear intersection hostile')
    cube_preservers=sum({frozenset((f[x],f[y])) for x,y in map(tuple,edges)}==edges
                        for cols,shift,f in affine)
    # A projective group preserving the design still need not preserve the
    # marked cube graph; do not silently promote its action to geometry.
    need(cube_preservers < 168,'full projective group is not cube symmetry')
    order3=tuple(transform((2,0,0,4),t) for t in P)  # z -> 4z
    fixed_pairs=[q for q in it.combinations(P,2) if block_image(order3,q)==frozenset(q)]
    fixed_equi=[q for q in equi if block_image(order3,q)==q]
    need(len(fixed_pairs)==1 and len(fixed_equi)==4,'inequivalent 28-set actions')
    report={
        'status':'PASS', 'universe':{'four_subsets':70,'PGL':336,'PSL':168},
        'equianharmonic_blocks':28,'chiral_SQS_blocks':[14,14],
        'derived_Fano_triples':sorted(map(sorted,U)),
        'label_order':'0,1,2,3,4,5,6,infinity',
        'binary_labels':LABEL,'conic_component_labels':component,
        'conic_basis':D,'walsh_eigenvalues_mask_order':eigen,
        'minus_inverse_affine':{'columns':ci,'translation':ai,
                               'fixed_points':[x for x in P if fi[x]==x]},
        'PSL_elements_preserving_standard_cube':cube_preservers,
        'order3_fixed_counts_28sets':{'two_subsets':len(fixed_pairs),
                                    'equianharmonic_quadruples':len(fixed_equi)},
        'negative_controls':['28 blocks have two completions per triple',
                             'opposite chirality requires a new affine marking',
                             'shear preserves design but destroys intersection',
                             'transitive affine PSL has no common origin',
                             '28 sextic labels and 28 quadrangles have different characters']}
    print(json.dumps(report,sort_keys=True,indent=2))

if __name__=='__main__':
    main()
