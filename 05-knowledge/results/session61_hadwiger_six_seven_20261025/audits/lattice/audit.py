#!/usr/bin/env python3
"""Independent Session61 exact auditor. Does not import producer scripts."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from collections import Counter, defaultdict
import json
from pathlib import Path

OUT=Path(__file__).parent

def need(v,s):
    if not v: raise RuntimeError(s)

def q(v):
    a,b=v
    return a*a-a*b+b*b

def inner(u,v):
    a,b=u;c,d=v
    return a*c+b*d-F(a*d+b*c,2)

def sub(u,v):return tuple(a-b for a,b in zip(u,v))
def add(u,v):return tuple(a+b for a,b in zip(u,v))
def scale(t,u):return tuple(t*a for a in u)
def det(u,v):return u[0]*v[1]-u[1]*v[0]

def det3(M):
    a,b,c=M
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])

# Veronese identity, evaluated on exhaustive small independent triples;
# symbolic coefficient verification is in the accompanying proof, not inferred from this census.
vecs=list(product(range(-2,3),repeat=2))
count=0
for u,v,w in combinations(vecs,3):
    lhs=det3([(a*a,2*a*b,b*b) for a,b in (u,v,w)])
    rhs=2*det(u,v)*det(u,w)*det(v,w)
    need(lhs==rhs,'Veronese factorization')
    count+=1
# Hostile to old claim: exactly two directions can occur in nonarithmetic Gram.
# Q=x²+sqrt(2)xy+y² is positive definite and Q(1,0)=Q(0,1)=1.
# Rational/irrational comparison forces xy=0 and x²+y²=1, yielding four vectors.

B=((0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1))
c=lambda v,k=2:(v[0]+k*v[1])%7
valid=[]
for A,C in [(1,j) for j in range(7)]+[(0,1)]:
    if len({(A*a+C*b)%7 for a,b in B})==7:valid.append((A,C))
need(valid==[(1,2),(1,4)],'two index7 kernels')
# Every nonzero ideal vector of norm below21 occurs in this complete box:
# q(a,b)>=3a²/4,3b²/4, hence q<21 implies |a|,|b|<=5.
small=[v for v in product(range(-6,7),repeat=2) if v!=(0,0) and c(v)==0 and q(v)<=21]
need(sorted({q(v) for v in small})==[7,21],'first ideal shells')
short=[v for v in small if q(v)==7]
need(len(short)==6,'six shortest ideal vectors')
# Difference body 2H, cyclically ordered rational vertices.
V=[(F(4,3),F(2,3)),(F(2,3),F(4,3)),(F(-2,3),F(2,3)),
   (F(-4,3),F(-2,3)),(F(-2,3),F(-4,3)),(F(2,3),F(-2,3))]
def to_segment(v,p,r):
    d=sub(r,p)
    t=max(F(0),min(F(1),inner(sub(v,p),d)/q(d)))
    z=add(p,scale(t,d))
    return q(sub(v,z)),z
nearest={str(v):min(to_segment(v,V[i],V[(i+1)%6]) for i in range(6)) for v in short}
need({d for d,z in nearest.values()}=={F(7,3)},'exact cell distance')
need(nearest[str((2,-1))][1]==(F(2,3),F(-2,3)),'nearest vertex')

# Construct geometric faces first, then quotient; no preloaded Fano design.
triangles=(((0,0),(1,0),(1,1)),((0,0),(1,1),(0,1)))
faces_by_chirality=[]
for k in (2,4):
    systems=[]
    for tri in triangles:
        systems.append(frozenset(frozenset((t+c(v,k))%7 for v in tri) for t in range(7)))
    faces_by_chirality.append(tuple(systems))
up,down=faces_by_chirality[0]
faces=up|down
need((len(up),len(down),len(faces))==(7,7,14),'14 distinct faces')
edges=Counter(e for f in faces for e in combinations(sorted(f),2))
need(len(edges)==21 and set(edges.values())=={2},'closed pseudomanifold')
links={}
for v in range(7):
    E={tuple(sorted(f-{v})) for f in faces if v in f}
    adj=defaultdict(set)
    for a,b in E:adj[a].add(b);adj[b].add(a)
    need(len(adj)==6 and all(len(z)==2 for z in adj.values()),'link degree2')
    reached={next(iter(adj))};todo=list(reached)
    while todo:
        w=todo.pop()
        for z in adj[w]-reached:reached.add(z);todo.append(z)
    need(len(reached)==6,'connected link')
    links[str(v)]=sorted(E)
# Explicit orientation inherited from both lattice triangles.
orient_edges=Counter()
for tri in triangles:
    for t in range(7):
        f=[(t+c(v))%7 for v in tri]
        for i in range(3):
            a,b=f[i],f[(i+1)%3]
            orient_edges[tuple(sorted((a,b)))]+=1 if a<b else -1
need(set(orient_edges.values())=={0},'coherent orientation')
need(7-21+14==0,'Euler zero')
for design in (up,down):
    pairs=Counter(e for f in design for e in combinations(sorted(f),2))
    need(len(pairs)==21 and set(pairs.values())=={1},'Steiner triple system')

# Work in original branch ordering, deliberately unlike residue ordering.
branch_coords=B
inv_c={c(v):i for i,v in enumerate(branch_coords)}
def pull(design,k=2):
    inv={c(v,k):i for i,v in enumerate(branch_coords)}
    return frozenset(frozenset(inv[t] for t in f) for f in design)
branches=[pull(d,k) for k,sys in zip((2,4),faces_by_chirality) for d in sys]
need(len(set(branches))==4,'four chosen branch designs')
actual=frozenset(frozenset(t) for t in ((0,1,2),(0,3,4),(0,5,6)))
need(all(not d&actual for d in branches),'all four avoid actual crossings')

# All thirty labelled Fano designs, generated via F2 incidence from scratch.
xor_lines=[t for t in combinations(range(1,8),3) if t[0]^t[1]^t[2]==0]
all_designs=set()
markings={}
for lab in permutations(range(1,8)):
    inv={v:i for i,v in enumerate(lab)}
    d=frozenset(frozenset(inv[t] for t in f) for f in xor_lines)
    all_designs.add(d)
    markings.setdefault(d,lab)
need(len(all_designs)==30,'30 labelled designs')
need(Counter(len(d&actual) for d in all_designs)=={0:16,1:12,3:2},'16/12/2')

# Independently construct complete metric dihedral group by unit multiplication
# and conjugation, transport to branch coordinates, and compute all orbits.
R=lambda v:(v[0]-v[1],v[0])
C=lambda v:(v[0]-v[1],-v[1])
idx={v:i for i,v in enumerate(branch_coords)}
G=set()
for flip in (0,1):
    for turns in range(6):
        p=[]
        for v in branch_coords:
            w=C(v) if flip else v
            for _ in range(turns):w=R(w)
            p.append(idx[w])
        G.add(tuple(p))
need(len(G)==12,'D12')
def move(d,p):return frozenset(frozenset(p[i] for i in f) for f in d)
remaining=set(all_designs);orbits=[]
while remaining:
    d=min(remaining,key=lambda z:sorted(map(sorted,z)))
    orbit={move(d,p) for p in G}
    remaining-=orbit
    orbits.append((len(d&actual),len(orbit)))
need(sorted(orbits)==[(0,4),(0,12),(1,6),(1,6),(3,2)],'all D12 orbits')
orb={move(branches[0],p) for p in G}
need(orb==set(branches),'exact special orbit4 identification')
# Verify binary kernels, not merely face systems.
def kernel(design):
    lab=markings[design]
    def parity(mask):
        z=0
        for i in range(7):
            if (mask>>i)&1:z^=lab[i]
        return z
    return frozenset(mask for mask in range(128) if mask.bit_count()%2==0 and parity(mask)==0)
for d in orb:
    ker=kernel(d)
    need(Counter(m.bit_count() for m in ker)=={0:1,4:7},'simplex kernel')
    need(all((127^sum(1<<i for i in f)) not in ker for f in actual),'fixed support absence')
    stab=[p for p in G if move(d,p)==d]
    need(len(stab)==3,'C3 stabilizer')
    lab=markings[d]
    nontrivial=next(p for p in stab if p!=tuple(range(7)))
    act={0:0}|{lab[i]:lab[nontrivial[i]] for i in range(7)}
    need(all(act[x^y]==act[x]^act[y] for x,y in product(range(8),repeat=2)),'quotient linear action')
    need(sum(act[x]==x for x in range(8))==2,'fixed subspace dim1')

# Hostile controls: branch crossings are NOT Fano lines, and whole torus faces
# are NOT one Fano plane (two occurrences per pair); T² is not unit distance.
need(up.isdisjoint(down),'two Fano systems disjoint')
need(set(q(v) for a,b in product(B[1:],repeat=2) for v in [add(a,b)] if v!=(0,0))=={1,3,4},'T-square nonunit lengths')

summary={
 'veronese_integer_triples_checked':count,
 'perfect_index7_kernels':valid,
 'cell_distance_squared':'7/3',
 'strict_scale_interval':['sqrt(3/7)','sqrt(3)/2'],
 'torus':{'V':7,'E':21,'F':14,'links':links,'orientable':True,'Euler':0},
 'fano_face_systems':[[sorted(map(sorted,d)) for d in sys] for sys in faces_by_chirality],
 'branch_designs': [sorted(map(sorted,d)) for d in branches],
 'D12_orbits':sorted(orbits),
 'special_orbit4_verified':True,
 'chosen_simplex_kernels':[sorted(kernel(d)) for d in branches],
 'all_four_free_inherited_surface_criterion':True,
 'stabilizer_C3_binary_fixed_dimension':1,
 'status':'PASS; exact finite tests; lattice general proof and geometric/surface deductions in AUDIT.md'
}
(OUT/'certificate.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ('torus','fano_face_systems','branch_designs')},sort_keys=True))
