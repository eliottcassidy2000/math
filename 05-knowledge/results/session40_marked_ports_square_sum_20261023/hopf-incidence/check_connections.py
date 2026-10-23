#!/usr/bin/env python3
"""Exact hostile controls, standard-library only. No asserts: also works under -O."""
import json, math, itertools, hashlib
from fractions import Fraction as F
from pathlib import Path
counts={}
def check(ok, label):
    if not ok: raise RuntimeError(label)
    counts[label]=counts.get(label,0)+1
# Pascal support and five-tap carry repair.
for r in range(2,11):
    n=2**r
    for j in range(n+1):
        expected={0:1,n//4:4,n//2:6,3*n//4:4,n:1}.get(j,0)
        check(math.comb(n,j)%8==expected,'mod8_five_tap_coefficients')
for r in range(1,10):
    n=2**r
    for s in range(1,r+2):
        support=[j for j in range(n+1) if math.comb(n,j)%(2**s)]
        step=2**max(0,r-s+1)
        check(support==list(range(0,n+1,step)),'prime_power_support')
for r in range(2,9):
    n=2**r
    for j in (0,n//4,n//2,3*n//4,n):
        check(math.comb(n,j)%8!=0,'each_sample_indispensable')
    check(math.comb(n,n//4)%8==4,'endpoint_only_hostile')
    for j in range(n+1):
        check(math.comb(n,j)%8==math.comb(n,n-j)%8,'reversal_blindness')
# Every additive C3 -> Z/(2^s) map is zero.
for s in range(1,13):
    check([a for a in range(2**s) if 3*a%(2**s)==0]==[0],'coprime_phase_firewall')
# Exact Gaussian rational arithmetic for projectors and the Bargmann identity.
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
def mm(a,b):
    return [[add(mul(a[i][0],b[0][j]),mul(a[i][1],b[1][j])) for j in range(2)] for i in range(2)]
def tr(a): return add(a[0][0],a[1][1])
def proj(u):
    x,y,z=u
    return [[(F(1+z,2),0),(F(x,2),F(-y,2))],[(F(x,2),F(y,2)),(F(1-z,2),0)]]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
axes=[tuple(sign if i==j else 0 for i in range(3)) for j in range(3) for sign in (-1,1)]
for u,v,w in itertools.product(axes,repeat=3):
    value=tr(mm(mm(proj(u),proj(v)),proj(w)))
    expected=(F(1+dot(u,v)+dot(v,w)+dot(w,u),4),F(dot(u,cross(v,w)),4))
    check(value==expected,'bargmann_pauli_identity')
u,v,w=(0,0,1),(1,0,0),(0,1,0)
wm=(0,-1,0)
check(tr(mm(mm(proj(u),proj(v)),proj(w)))==(F(1,4),F(1,4)),'oriented_phase_positive')
check(tr(mm(mm(proj(u),proj(v)),proj(wm)))==(F(1,4),F(-1,4)),'oriented_phase_negative')
for a,b in ((u,v),(v,w),(w,u),(v,wm),(wm,u)):
    check(tr(mm(proj(a),proj(b)))==(F(1,2),0),'equal_pair_magnitudes')
# Fixed-root Hopf/Bargmann sign tournaments: exact integer planar controls.
def det2(a,b): return a[0]*b[1]-a[1]*b[0]
for points in [[(1,0),(-1,1),(-1,-1)],[(10,0),(3,10),(-8,6),(-8,-6),(3,-10)]]:
    n=len(points)
    for i in range(n):
        for j in range(n):
            if i!=j:
                d=det2(points[i],points[j])
                check(d!=0 and (d>0)==(0<(j-i)%n<=(n-1)//2),'hopf_cyclic_realizations')
qr={1,2,4}
for i in range(7):
    out=[j for j in range(7) if (j-i)%7 in qr]
    check(all(sum((v-u)%7 in qr for v in out if v!=u)==1 for u in out),'paley7_local_cycle_obstruction')
local=0
edges=list(itertools.combinations(range(4),2))
for bits in range(64):
    adj=[[False]*4 for _ in range(4)]
    for k,(i,j) in enumerate(edges):
        adj[i][j]=bool(bits>>k&1);adj[j][i]=not adj[i][j]
    good=True
    for i in range(4):
        for incoming in (False,True):
            neigh=[j for j in range(4) if adj[j][i]] if incoming else [j for j in range(4) if adj[i][j]]
            if len(neigh)==3 and all(sum(adj[u][v] for v in neigh if v!=u)==1 for u in neigh): good=False
    local+=good
check(local==48,'four_vertex_local_transitivity_count')
# Square-shell reflection composition, including invalid path controls.
shells=[4,9,16,25,36,49]
valid=0
for length in range(1,6):
    for word in itertools.product(shells,repeat=length):
        shift=sum((-1)**(length-j-1)*s for j,s in enumerate(word))
        for start in (1,2,4,6,13,25):
            x=start; seq=[x]
            for s in word: x=s-x; seq.append(x)
            check(x==(-1)**length*start+shift,'square_shell_affine_formula')
            if all(1<=z<=25 for z in seq) and all(a!=b for a,b in zip(seq,seq[1:])):
                valid+=1
                if x==start:
                    check(shift==(2*start if length%2 else 0),'closed_walk_square_relation')
# A balanced even shell word can still leave the interval / hit forbidden loops.
word=(4,9,9,4); seq=[25]
for s in word: seq.append(s-seq[-1])
check(seq[-1]==seq[0] and min(seq)<1,'balanced_word_not_realizable')
# One vertex repeated: closed-walk relation does not imply a simple cycle.
check([1,3,1]==[1,4-1,4-(4-1)],'backtracking_not_cycle')
# Chordal actual transport: two support triangles generate all closed walks.
perms=list(itertools.permutations(range(3)))
def compose(p,q): return tuple(p[q[i]] for i in range(3))
def inv(p): return tuple(p.index(i) for i in range(3))
identity=(0,1,2)
flat_count=0
for p01,p02,p03,p12,p23 in itertools.product(perms,repeat=5):
    triangle_flat=(compose(p12,p01)==p02 and compose(p23,p02)==p03)
    frame_flat=(p12==compose(p02,inv(p01)) and p23==compose(p03,inv(p02)))
    check(triangle_flat==frame_flat,'chordal_triangle_gluing')
    flat_count+=triangle_flat
check(flat_count==6**3,'chordal_flat_count')
nonflat=0
for p01,p12,p23,p30 in itertools.product(perms,repeat=4):
    hol=compose(p30,compose(p23,compose(p12,p01)))
    nonflat+=hol!=identity
    # Every triangle condition is vacuous on this support.
    check(sorted(hol)==[0,1,2],'chordless_cycle_holonomy')
check(nonflat==6**4-6**3,'vacuous_triangle_hostile_count')
receipt={'status':'PASS','counts':counts,'total_checks':sum(counts.values()),'valid_square_walk_instances':valid,
 'scope':'Exact algebraic controls; not a native tournament wiring census, Hopf realization, or square-sum Hamiltonian enumeration.',
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(receipt,indent=2,sort_keys=True))
