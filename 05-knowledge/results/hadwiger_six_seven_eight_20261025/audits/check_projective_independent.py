#!/usr/bin/env python3
"""Independent generator/orbit audit of the projective eight; no producer imports."""
from itertools import combinations,permutations,product
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
P=tuple(range(8)); infinity=7
labels=(1,2,4,3,6,7,5,0)
def need(ok,message):
    if not ok:
        raise RuntimeError(message)
def compose(a,b):
    return tuple(a[b[x]] for x in P)
def generated(generators):
    seen={P};todo=[P]
    while todo:
        a=todo.pop()
        for b in generators:
            c=compose(a,b)
            if c not in seen:
                seen.add(c);todo.append(c)
    return seen
def act(g,A):
    return frozenset(g[x] for x in A)
def dot(a,b):
    return (a&b).bit_count()%2
def lin(columns,x):
    ans=0
    for j in range(3):
        if x>>j&1:
            ans^=columns[j]
    return ans
T=tuple((x+1)%7 if x!=7 else 7 for x in P)
S=tuple(7 if x==0 else 0 if x==7 else -pow(x,-1,7)%7 for x in P)
N=tuple(3*x%7 if x!=7 else 7 for x in P)
G=generated((T,S));wide=generated((T,S,N))
need(len(G)==168 and len(wide)==336,'generated group orders')
seed=frozenset((7,0,1,3))
A={act(g,seed) for g in G}
B={act(N,b) for b in A}
need(len(A)==len(B)==14 and not A&B,'two chiral systems')
need({act(g,seed) for g in wide}==A|B,'PGL orbit')

# Independent homogeneous polynomial test for equianharmonic ratio.
def hp(x):
    return (1,0) if x==7 else (x,1)
def determinant(a,b):
    return (a[0]*b[1]-a[1]*b[0])%7
def equi(block):
    a,b,c,d=[hp(x) for x in sorted(block)]
    X=determinant(a,c)*determinant(b,d)%7
    Y=determinant(a,d)*determinant(b,c)%7
    return (X*X-X*Y+Y*Y)%7==0
all_blocks={frozenset(b) for b in combinations(P,4)}
need({b for b in all_blocks if equi(b)}==A|B,'crossratio polynomial orbit identity')
for family in (A,B):
    need(all(frozenset(P)-b in family for b in family),'complement closure')
    need(all(sum(set(t)<=b for b in family)==1 for t in combinations(P,3)),
         'Steiner triple multiplicity')
need(all(sum(set(t)<=b for b in A|B)==2 for t in combinations(P,3)),
     'unchiral hostile')
planes={frozenset(x for x in P if dot(a,x)==b) for a in range(1,8) for b in (0,1)}
need({frozenset(labels[x] for x in b) for b in A}==planes,'chosen marking')
need({frozenset(labels[x] for x in b) for b in B}!=planes,'opposite marking hostile')
conic=tuple(sum(dot(labels[x],a)<<j for j,a in enumerate((7,3,5))) for x in P)
need(conic==(7,3,5,4,6,1,2,0),'conic correspondence table')
need({frozenset(conic[x] for x in b) for b in A}==planes,'conic planes')

linear_parts=set();affine={}
for g in G:
    f={labels[x]:labels[g[x]] for x in P}
    shift=f[0]
    columns=tuple(f[1<<j]^shift for j in range(3))
    need(all(f[x]==(lin(columns,x)^shift) for x in P),'affine structure')
    linear_parts.add(columns);affine[g]=(columns,shift,f)
need(len(linear_parts)==168,'linear part faithfulness')
need(affine[T][0]==(2,4,3) and affine[T][1]==0,'translation affine formula')
need(affine[S][0]==(1,4,2) and affine[S][1]==1,'inversion affine formula')
need(not any(affine[S][2][x]==x for x in P),'affine inversion fixed points')
need(1 not in {x^lin((1,4,2),x) for x in P},'crossed cocycle obstruction')
need(sum(lin((1,4,2),x)==x for x in P)==4,'linear/affine action hostile')
edges={frozenset((x,x^(1<<j))) for x in P for j in range(3)}
preservers=[g for g,(_,_,f) in affine.items()
            if {frozenset(f[x] for x in edge) for edge in edges}==edges]
need(len(preservers)==6,'marked cube preserving subgroup')
M=[[-2 if x==y else 2 if frozenset((x,y)) in edges else 0 for y in P] for x in P]
for mask in P:
    v=[(-1)**dot(mask,x) for x in P]
    need([sum(M[x][y]*v[y] for y in P) for x in P]
         ==[4*(1-mask.bit_count())*a for a in v],'Walsh divisor modes')
shear=tuple(x^(((x>>1)&1)<<0) for x in P)
need({act(shear,b) for b in planes}==planes and M[0][2]==2 and M[shear[0]][shear[2]]==0,
     'affine incidence loss')
order3=tuple(4*x%7 if x!=7 else 7 for x in P)
fixed_pairs=sum(act(order3,pair)==frozenset(pair) for pair in combinations(P,2))
fixed_quad=sum(act(order3,b)==b for b in A|B)
need((fixed_pairs,fixed_quad)==(1,4),'two28sets character hostile')

result=dict(status='PASS / FINITE-EXACT independent review',
            generated_PSL=len(G),generated_PGL=len(wide),blocks=[len(A),len(B)],
            affine_linear_parts=len(linear_parts),cube_preservers=len(preservers),
            conic_marking=conic,linear_inversion_fixed=4,affine_inversion_fixed=0,
            order3_fixed_pairs=fixed_pairs,order3_fixed_quadrangles=fixed_quad,
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
text=json.dumps(result,indent=2,sort_keys=True)+'\n'
(HERE/'projective_independent.json').write_text(text,encoding='utf-8')
print(text,end='')
