"""Exact metric cube-fold and recovered Moser controls; standard library only."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json

checks = 0
def need(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)

def q(p):
    a,b = p
    return a*a-a*b+b*b

def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))

cube_edges = {(a,b) for a,b in combinations(range(8),2) if (a^b).bit_count()==1}
faces = []
for i,j in combinations(range(3),2):
    k = next(k for k in range(3) if k not in (i,j))
    for bit in (0,1):
        a = bit << k
        faces.append((a,a^(1<<i),a^(1<<i)^(1<<j),a^(1<<j)))

def rref(rows):
    a = [list(map(F,row)) for row in rows]
    pivot = 0
    for col in range(8):
        nz = next((i for i in range(pivot,len(a)) if a[i][col]),None)
        if nz is None:
            continue
        a[pivot],a[nz] = a[nz],a[pivot]
        d = a[pivot][col]
        a[pivot] = [x/d for x in a[pivot]]
        for i in range(len(a)):
            if i!=pivot and a[i][col]:
                d = a[i][col]
                a[i] = [x-d*y for x,y in zip(a[i],a[pivot])]
        pivot += 1
        if pivot==len(a):
            break
    return a[:pivot]

def difference(a,b):
    row=[0]*8
    row[a]=1
    row[b]=-1
    return row

collision_census = Counter()
forced_witnesses = []
for a,b in combinations(range(8),2):
    d=(a^b).bit_count()
    if d==1:
        collision_census['edge_forbidden'] += 1
        continue
    rows=[difference(a,b)]
    for face in faces:
        if a in face and b in face:
            continue  # Such a face is not injective; do not impose its rhombus law.
        row=[0]*8
        for sign,v in zip((1,-1,1,-1),face):
            row[v]=sign
        rows.append(row)
    rank=len(rref(rows))
    forced=[(u,v) for u,v in combinations(range(8),2)
            if len(rref(rows+[difference(u,v)]))==rank]
    if d==2:
        need(len(forced)>1,'distance-two collision must propagate')
        collision_census['second_collision_forced'] += 1
        forced_witnesses.append({'pair':[a,b],'forced':forced})
    else:
        need(forced==[(a,b)],'opposite collision has exactly one forced pair')
        need(8-rank==3,'opposite scalar space is constant plus two dimensions')
        collision_census['opposite_fold'] += 1
need(dict(collision_census)=={'edge_forbidden':12,'second_collision_forced':12,
                             'opposite_fold':4},'full 28-pair census')

# u=(1,0),v=(0,1),w=(-1,-1) in the Eisenstein norm a²-ab+b².
vectors=[(1,0),(0,1),(-1,-1)]
fold=[tuple(sum(((m>>i)&1)*vectors[i][j] for i in range(3))
            for j in range(2)) for m in range(8)]
need(len(set(fold))==7 and fold[0]==fold[7]==(0,0),'single opposite fold')
for a,b in cube_edges:
    need(q(sub(fold[a],fold[b]))==1,'all twelve cube edges stay unit')
points=sorted(set(fold))
index={p:i for i,p in enumerate(points)}
quotient={tuple(sorted((index[fold[a]],index[fold[b]]))) for a,b in cube_edges}
unit={(a,b) for a,b in combinations(range(7),2) if q(sub(points[a],points[b]))==1}
need(quotient==unit and len(unit)==12,'no extra image unit edge')
distance_counts=Counter(q(sub(points[a],points[b])) for a,b in combinations(range(7),2))
need(distance_counts=={1:12,3:6,4:3},'full seven-point metric')
pullback_counts=Counter(q(sub(fold[a],fold[b])) for a,b in combinations(range(8),2))
need(pullback_counts=={0:1,1:18,3:6,4:3},'labelled pullback hostile')
for m in range(8):
    need(fold[m^7]==tuple(-x for x in fold[m]),'complement becomes central inversion')
    for axis in range(2):
        need(-2*fold[m][axis]+2*sum(fold[m^(1<<i)][axis] for i in range(3))==0,
             'coordinate functions lie in K3 conic-matrix kernel')

# All cube graph automorphisms and the stabilizer of the collision fibre.
stabilizers=[]
for shift in range(8):
    for perm in permutations(range(3)):
        image=[shift ^ sum(((m>>i)&1)<<perm[i] for i in range(3)) for m in range(8)]
        need({tuple(sorted((image[a],image[b]))) for a,b in cube_edges}==cube_edges,
             'cube automorphism')
        descends={image[0],image[7]}=={0,7}
        need(descends==(shift in (0,7)),'exact descent stabilizer')
        if descends:
            p1,p2=fold[image[1]],fold[image[2]]
            A=(p1,p2)
            for m,p in enumerate(fold):
                expected=tuple(p[0]*p1[k]+p[1]*p2[k] for k in range(2))
                need(expected==fold[image[m]],'linear metric action')
            need(q(p1)==q(p2)==q(tuple(p1[k]+p2[k] for k in range(2)))==1,
                 'Eisenstein metric preserved')
            stabilizers.append(image)
need(len(stabilizers)==12,'D12 stabilizer')

# A faithful, injective eight-point control: every nonedge is checked.
unit_vectors=[(F(1),F(0)),(F(3,5),F(4,5)),(F(5,13),F(12,13))]
rational=[tuple(sum(((m>>i)&1)*unit_vectors[i][j] for i in range(3))
                for j in range(2)) for m in range(8)]
need(len(set(rational))==8,'faithful cube is injective')
rational_edges={(a,b) for a,b in combinations(range(8),2)
                if sum(t*t for t in sub(rational[a],rational[b]))==1}
need(rational_edges==cube_edges,'faithful cube has exactly intended edges')

# Independent arithmetic in Q(s,t), s²=-3,t²=-11, conjugation s,t -> -s,-t.
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def scale(a,c): return tuple(x*F(c) for x in a)
def mul(a,b):
    out=[F(0)]*4
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            common=i&j
            factor=(-3 if common&1 else 1)*(-11 if common&2 else 1)
            out[i^j] += x*y*factor
    return tuple(out)
def conj(a): return tuple((-1 if i in (1,2) else 1)*x for i,x in enumerate(a))
zero=(F(0),)*4
one=(F(1),F(0),F(0),F(0))
alpha=(F(1,2),F(1,2),F(0),F(0))
beta=(F(5,6),F(0),F(1,6),F(0))
spindle=[zero,one,alpha,add(one,alpha),beta,mul(alpha,beta),mul(beta,add(one,alpha))]
need(len(set(spindle))==7,'Moser seven distinct points')
spindle_edges=set()
for a,b in combinations(range(7),2):
    delta=sub(spindle[a],spindle[b])
    dist=mul(delta,conj(delta))
    need(dist!=zero,'no Moser collision')
    if dist==one:
        spindle_edges.add((a,b))
expected={(0,1),(0,2),(1,2),(1,3),(2,3),
          (0,4),(0,5),(4,5),(4,6),(5,6),(3,6)}
need(spindle_edges==expected,'full induced Moser graph')

def coloring_count(n,edges,k):
    count=0
    witness=None
    for col in product(range(k),repeat=n):
        if all(col[a]!=col[b] for a,b in edges):
            count+=1
            if witness is None: witness=col
    return count,witness

cube_chi2=coloring_count(8,cube_edges,2)
flower_chi2=coloring_count(7,unit,2)
flower_chi3=coloring_count(7,unit,3)
moser_chi3=coloring_count(7,spindle_edges,3)
moser_chi4=coloring_count(7,spindle_edges,4)
need(cube_chi2[0]==2,'cube is bipartite connected')
need(flower_chi2[0]==0 and flower_chi3[0]==6,'flower chromatic three')
need(moser_chi3[0]==0 and moser_chi4[0]>0,'Moser chromatic four')
need(moser_chi4[0]==384,'recovered exact Moser four-color count')
need(beta[2]!=0 and beta[3]==0,'spindle rotation outside Q(s)')

# Independent coordinate basis rank: 1,alpha,beta,alpha*beta span degree four.
basis=[one,alpha,beta,mul(alpha,beta)]
det=F(0)
for perm in permutations(range(4)):
    inversions=sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
    term=F((-1)**inversions)
    for i in range(4): term*=basis[i][perm[i]]
    det+=term
need(det!=0,'spindle coordinate additive rank four')

receipt={
    'status':'FINITE-EXACT; all-image classification is proved in RESULT.md',
    'checks':checks,
    'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'collision_pair_census':dict(collision_census),
    'distance_two_witness':forced_witnesses[0],
    'fold_coordinates':fold,
    'image_squared_distance_counts':dict(distance_counts),
    'labelled_squared_distance_counts':dict(pullback_counts),
    'descended_cube_automorphisms':len(stabilizers),
    'faithful_cube_vectors':[[str(x) for x in p] for p in unit_vectors],
    'unit_graphs':{'cube':[8,12,2],'flower':[7,12,3],'Moser':[7,11,4]},
    'Moser_four_colorings':moser_chi4[0],
    'Moser_four_color_witness':moser_chi4[1],
    'Moser_basis_determinant':str(det)
}
print(json.dumps(receipt,indent=2,sort_keys=True))
