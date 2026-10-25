#!/usr/bin/env python3
"""Exact standalone certificate for the Session61 A2 / seven-color bridge.

Standard library only. Explicit failures survive `python3 -O`.
The unbounded claims are proved in RESULTS.md; finite checks here certify
matrices, complete finite classifications, local graph quotients, and metric
algebra. This file never writes outside its explicitly requested output.
"""
from fractions import Fraction as F
from itertools import product, combinations
import json

checks = 0

def require(value, label):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(label)

def add(u,v): return (u[0]+v[0],u[1]+v[1])
def sub(u,v): return (u[0]-v[0],u[1]-v[1])
def scale(k,u): return (k*u[0],k*u[1])
def mul(u,v):
    a,b=u; c,d=v
    return (a*c-b*d, a*d+b*c-b*d)
def norm(v):
    a,b=v
    return a*a-a*b+b*b

def dot(u,v): return (norm(add(u,v))-norm(u)-norm(v))/2

def qdot(u,v):
    a,b=u; c,d=v
    return a*c+b*d-F(a*d+b*c,2)

def mm(A,B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def mv(A,v): return tuple(sum(A[i][j]*v[j] for j in range(2)) for i in range(2))
def transpose(A): return tuple(zip(*A))
def det(A): return A[0][0]*A[1][1]-A[0][1]*A[1][0]
def matpow(A,n):
    out=((1,0),(0,1))
    for _ in range(n): out=mm(out,A)
    return out

def c(k,v): return (v[0]+k*v[1])%7

O=(0,0)
U={(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)}
B={O}|U
G=((F(1),F(-1,2)),(F(-1,2),F(1)))
R=((1,-1),(1,0))
C=((1,-1),(0,-1))
I=((1,0),(0,1))
require(matpow(R,6)==I and matpow(R,3)==((-1,0),(0,-1)), 'rotation order6')
require(mm(C,C)==I and mm(mm(C,R),C)==matpow(R,5),'dihedral relations')
D12={matpow(R,i) for i in range(6)}|{mm(matpow(R,i),C) for i in range(6)}
require(len(D12)==12, 'twelve matrices')
for A in D12:
    require({mv(A,v) for v in U}==U,'root permutation')
    require(mm(mm(transpose(A),G),A)==G,'Gram preservation')

alpha=(2,-1); alphabar=(3,1)
M2=((2,1),(-1,3)); M4=((3,-1),(1,2))
require(mul(alpha,alphabar)==(7,0), 'split 7 product')
for M in (M2,M4):
    require(det(M)==7,'index7 determinant')
    require(mm(mm(transpose(M),G),M)==tuple(tuple(7*t for t in row) for row in G),'norm scales7')
require(mm(M2,M4)==((7,0),(0,7)), 'multiplication product matrices')
for v in B:
    require(c(2,mv(R,v))==3*c(2,v)%7,'c2 rotation')
    require(c(4,mv(R,v))==5*c(4,v)%7,'c4 rotation')
    require(c(2,mv(C,v))==c(4,v),'conjugation exchange')
require(c(2,alpha)==0 and c(4,alpha)!=0,'distinct kernels hostile')

# Split residue ring and all norm-seven elements; conjugate shells differ.
for v in product(range(7),repeat=2):
    require(norm(v)%7==c(2,v)*c(4,v)%7,'split norm identity')
require(len({(c(2,v),c(4,v)) for v in product(range(7),repeat=2)})==49,'CRT bijective')
require({(c(2,v),c(4,v)) for v in U}=={(u,pow(u,-1,7)) for u in range(1,7)},'six roots exhaust split norm1 torus')
norm7={v for v in product(range(-3,4),repeat=2) if norm(v)==7}
require(norm7=={mul(alpha,u) for u in U}|{mul(alphabar,u) for u in U},'all twelve norm7 vectors')
require(len(norm7)==12,'norm7 count twelve not six')

# Complete classification of all 48 nonzero linear maps F7^2 -> F7.
valid=[]
for a,b in product(range(7),repeat=2):
    if (a,b)==(0,0): continue
    if len({(a*x+b*y)%7 for x,y in B})==7: valid.append((a,b))
expected={(a,(a*k)%7) for a in range(1,7) for k in (2,4)}
require(set(valid)==expected and len(valid)==12,'all linear maps')
# Each of the 8 lines in the projective dual represents one index7 subgroup.
projective=[(1,k) for k in range(7)]+[(0,1)]
perfect=[(a,b) for a,b in projective if len({(a*x+b*y)%7 for x,y in B})==7]
require(perfect==[(1,2),(1,4)],'exactly two index7 perfect sublattices')
for k in (2,4):
    require({c(k,v) for v in B}==set(range(7)),'complete flower reps')

# Explicit coset decomposition in a diagnostic box, in addition to algebraic proof.
reps={c(2,v):v for v in B}
for v in product(range(-12,13),repeat=2):
    r=reps[c(2,v)]; w=sub(v,r)
    n1=F(3*w[0]-w[1],7); n2=F(w[0]+2*w[1],7)
    require(n1.denominator==n2.denominator==1 and mv(M2,(n1,n2))==w,'perfect ball decomposition')

# T itself versus T^2. Exact graph-radius-two difference set has three norms.
D2=U|{add(u,v) for u in U for v in U}; D2.discard(O)
require(len(D2)==18 and sorted({norm(v) for v in D2})==[1,3,4],'T2 metric support')
for u,v in combinations(B,2):
    require(sub(v,u) in D2,'flower clique in square')
for k in (2,4):
    for v in D2: require(c(k,v)!=0,'proper T2 coloring')
for v in U: require((v[0]+v[1])%3!=0,'proper T coloring')
triangle=[O,(1,0),(1,1)]
require(all(sub(v,u) in U for u,v in combinations(triangle,2)),'unit triangle')
flower_edges=[(u,v) for u,v in combinations(sorted(B),2) if sub(v,u) in U]
require(len(flower_edges)==12,'wheel has twelve edges')

# Voronoi hexagon H, all rational coordinates in the Eisenstein basis.
H=[(F(2,3),F(1,3)),(F(1,3),F(2,3)),(F(-1,3),F(1,3)),
   (F(-2,3),F(-1,3)),(F(-1,3),F(-2,3)),(F(1,3),F(-1,3))]
for v in H:
    require(norm(v)==F(1,3),'circumradius squared')
    x,z=2*v[0]-v[1],-v[0]+2*v[1]
    require(max(abs(x),abs(z),abs(x+z))==1,'hourglass boundary')
    require(norm(v)==(x*x+x*z+z*z)/3,'primal metric coordinate map')
require(max(norm(sub(u,v)) for u,v in combinations(H,2))==F(4,3),'diameter squared')
p=(F(2,3),F(-2,3)); delta=sub(alpha,p)
require(p in [scale(2,v) for v in H], 'closest 2H vertex')
require(delta==add(scale(F(4,3),(1,0)),scale(F(1,3),(0,-1))),'normal cone certificate')
require(all(qdot(delta,sub(scale(2,v),p))<=0 for v in H),'support minimizer certificate')
require(norm(delta)==F(7,3),'exact separation squared')
# Rational segment minimization independently certifies every short ideal vector.
def point_segment_sq(p,a,b):
    w=sub(b,a); t=qdot(sub(p,a),w)/norm(w)
    t=max(F(0),min(F(1),t))
    return norm(sub(p,add(a,scale(t,w))))
H2=[scale(2,v) for v in H]
short={mul(alpha,u) for u in U}
require(len(short)==6 and all(norm(v)==7 for v in short),'all short ideal directions')
for v in short:
    require(min(point_segment_sq(v,H2[i],H2[(i+1)%6]) for i in range(6))==F(7,3),'short metric distances')
# At h=3/4 the squared diameter and separation straddle 1 strictly.
h=F(3,4)
require(h*h*F(4,3)==F(3,4)<1,'diameter margin')
require(h*h*F(7,3)==F(21,16)>1,'separation margin')
require(F(3,7)<h*h<F(3,4),'strict scale interval')

# Actual triangular quotient face designs (combinatorial incidence, not collinearity).
up={frozenset(((t+x)%7 for x in (0,1,3))) for t in range(7)}
down={frozenset(((t+x)%7 for x in (0,2,3))) for t in range(7)}
for faces in (up,down):
    require(len(faces)==7,'seven face orbits')
    pairs=[frozenset(e) for f in faces for e in combinations(f,2)]
    require(len(pairs)==len(set(pairs))==21,'Fano pair multiplicity one')
require(up.isdisjoint(down),'opposite face designs disjoint')
require(7-21+14==0,'torus Euler characteristic')
for t in range(7):
    require({(t+c(2,u))%7 for u in U}==set(range(7))-{t},'graph covering K7 neighborhood')
perm_faces=lambda F,m:{frozenset(m*x%7 for x in f) for f in F}
require(perm_faces(up,3)==down and perm_faces(down,3)==up,'rotation switches faces')
require(perm_faces(up,2)==up and perm_faces(down,2)==down,'120 rotation preserves faces')
require(perm_faces(up,6)==down,'central inversion exchanges Fano structures')

# Hostiles: nonideal seven-color map is not even T^2 proper; metric data cannot drop.
require(c(1,(1,-1))==0 and (1,-1) in D2,'bad linear quotient hostile')
require(c(2,(0,0))==c(2,alpha) and norm(alpha)==7,'same color has nonunit distance')
require(norm((1,1))==1 and (1+1)%3!=0,'three-color positive control')
require(F(1,4)*F(7,3)<1, 'too-small h fails strict separation condition')
require(F(4,3)>1,'unscaled hexagon fails within-cell diameter condition')

out={
 'status':'PASS', 'checks':checks,
 'scope':'finite exact algebra/classification; universal implications proved in RESULTS.md',
 'gram':[[str(v) for v in row] for row in G],
 'flower_labels_c2':[[*v,c(2,v)] for v in sorted(B)],
 'index_seven_subgroups':len(projective), 'valid_kernels':perfect,
 'valid_linear_maps':len(valid), 'dihedral_order':len(D12),
 'total_norm7_vectors':len(norm7), 'split_residue_ring_size':49,
 'T_flower_edges':len(flower_edges),'T_chromatic':3,'T_square_chromatic':7,
 'T_square_squared_edge_lengths':sorted({norm(v) for v in D2}),
 'unit_spacing_hex_diameter_squared':'4/3',
 'unit_spacing_same_color_separation_squared':'7/3',
 'h':'3/4','scaled_diameter_squared':'3/4','scaled_separation_squared':'21/16',
 'torus_VEF':[7,21,14],
 'up_Fano_faces':sorted([sorted(f) for f in up]),
 'down_Fano_faces':sorted([sorted(f) for f in down]),
 'negative_controls':['c1 kills a graph-distance-two vector','h=1/2 violates strict separation','h=1 violates cell diameter','conjugation sends kernel I2 to distinct I4'],
}
print(json.dumps(out,sort_keys=True,indent=2))
