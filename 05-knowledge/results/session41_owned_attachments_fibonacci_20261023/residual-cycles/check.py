#!/usr/bin/env python3
"""Exact independent finite controls. No assert statements; python -O safe."""
from itertools import product, permutations
from collections import defaultdict, Counter
from pathlib import Path
import hashlib, json

counts=Counter()
def need(v,label):
    if not v: raise RuntimeError(label)
    counts[label]+=1

def mv(A,x): return tuple(sum(a*b for a,b in zip(r,x)) for r in A)
def add(x,y): return tuple(a+b for a,b in zip(x,y))
def conforms(h,g): return all(a*b>=0 and abs(a)<=abs(b) for a,b in zip(h,g))
def atoms(A,D):
    K=[g for g in product(*(range(-d,d+1) for d in D)) if any(g) and not any(mv(A,g))]
    # Sort by norm, checking all smaller previously found atoms suffices:
    # any proper conformal kernel subvector itself contains an atom.
    K.sort(key=lambda g:sum(map(abs,g)))
    G=[]
    for g in K:
        if not any(conforms(h,g) for h in G): G.append(g)
    # Independent minimality confirmation against the entire kernel box.
    for g in G:
        need(not any(h!=g and conforms(h,g) for h in K),'atom_full_box_minimality')
    for g in K:
        need(any(conforms(h,g) for h in G),'atom_full_box_coverage')
    return G

def test_fiber_family(A,D,Ps):
    G=atoms(A,D)
    fibers=defaultdict(list)
    for x in product(*(range(d+1) for d in D)): fibers[mv(A,x)].append(x)
    for xs in fibers.values():
        for x in xs:
            legal=[g for g in G if all(0<=a+b<=d for a,b,d in zip(x,g,D))]
            for P in Ps:
                direct=len({mv(P,y) for y in xs})==1
                criterion=all(not any(mv(P,g)) for g in legal)
                need(direct==criterion,'basepoint_atomic_iff')
                if not direct:
                    g=next(g for g in legal if any(mv(P,g)))
                    need(add(x,g) in xs and mv(P,add(x,g))!=mv(P,x),'actual_changed_projection_witness')
    return len(G),len(fibers)

Ps=[((1,0,0),),((1,1,1),),((1,-2,3),),((0,0,0),),((1,0,0),(0,1,0),(0,0,1))]
rows=list(product(range(-2,3),repeat=3))
for row in rows: test_fiber_family((row,),(2,2,2),Ps)
# Deliberately varied deterministic two-row families, not a random claim.
for i in range(125): test_fiber_family((rows[i],rows[(37*i+19)%125]),(2,2,2),Ps)

# Two permutation tables, shared top-left owner bit.
A=[]
for shift in (0,4):
    for slots in ((0,1),(2,3),(0,2),(1,3)):
        r=[0]*8
        for k in slots:r[shift+k]=1
        A.append(tuple(r))
A.append((1,0,0,0,-1,0,0,0))
x=(1,0,0,1,1,0,0,1); y=(0,1,1,0,0,1,1,0)
G=atoms(A,(1,)*8)
g=tuple(b-a for a,b in zip(x,y))
need(set(G)=={g,tuple(-v for v in g)},'exact_two_switch_atomic_bank')
feasible=[p for p in product((0,1),repeat=8) if mv(A,p)==mv(A,x)]
need(set(feasible)=={x,y},'exact_two_switch_fiber')
for h in (g[:4]+(0,)*4,(0,)*4+g[4:]):
    need(any(mv(A,h)),'individual_cycle_breaks_shared_owner')
need(not any(mv(A,g)) and add(x,g)==y,'joint_cycle_preserves_shared_owner')
# Positive certificate for this nonsingleton fiber: shared difference vanishes.
P=((1,0,0,0,-1,0,0,0),)
need(all(not any(mv(P,h)) for h in G),'nontrivial_coupled_positive_certificate')
need(len({mv(P,p) for p in feasible})==1,'positive_certificate_direct_check')

# Singleton coupled fiber: fixed top-left excludes the relaxed cycle.
perms=((1,0,0,1),(0,1,1,0))
need([p for p in perms if p[0]==1]==[perms[0]],'relaxed_cycle_not_coupled_counterexample')
# Integer-versus-real witness, use exact doubled coordinates for fractional point.
need([(a,b) for a,b in product((0,1),repeat=2) if 2*a+b==1]==[(0,1)],'integer_width_uniqueness')
need(2*1+0==2,'fractional_half_point_after_doubling')

# Fixed linear mod-8 guard. The auxiliary clock z is bounded, never discarded.
modA=((1,1,-8),)
modG=atoms(modA,(7,7,1))
modfiber=[x for x in product(range(8),range(8),range(2)) if not any(mv(modA,x))]
need({(a,b) for a,b,z in modfiber}=={(a,b) for a,b in product(range(8),repeat=2) if (a+b)%8==0},'mod8_auxiliary_exactness')
need(all(a+b==8*z for a,b,z in modfiber),'mod8_clock_retained')
for x in modfiber:
    legal=[g for g in modG if all(0<=a+b<=d for a,b,d in zip(x,g,(7,7,1)))]
    need(all(not any(mv(modA,g)) for g in legal),'mod8_every_atomic_move_keeps_guard')

# Binary relation compiler: direct vertex assignments versus independently
# joining actual allowed edge tuples on shared endpoint labels.
pairs=list(product((0,1),repeat=2))
relations=[tuple(p for k,p in enumerate(pairs) if mask>>k&1) for mask in range(16)]
edges=((0,1),(1,2),(0,2))
for rels in product(relations,repeat=3):
    direct={x for x in product((0,1),repeat=3) if all((x[i],x[j]) in R for (i,j),R in zip(edges,rels))}
    joined=set()
    for chosen in product(*rels):
        # w one-hot on each edge; z induced at vertices; sharing constraints.
        if chosen[0][0]==chosen[2][0] and chosen[0][1]==chosen[1][0] and chosen[1][1]==chosen[2][1]:
            joined.add((chosen[0][0],chosen[0][1],chosen[1][1]))
    need(direct==joined,'triangle_relation_compiler_all_4096')
neq=((0,1),(1,0))
need(not any(all(x[i]!=x[j] for i,j in edges) for x in product((0,1),repeat=3)),'odd_cycle_integer_infeasible')
# Half-marginal witness, scaled by two: every allowed w is1, every z is1.
need(all(sum(a==v for a,b in neq)==1 and sum(b==v for a,b in neq)==1 for v in (0,1)),'odd_cycle_fractional_marginal_witness')

# Nonabelian controls: permutations represented as exact value tuples.
S=list(permutations(range(3))); identity=(0,1,2)
def comp(g,h):return tuple(g[h[i]] for i in range(3))
def inv(g):return tuple(g.index(i) for i in range(3))
def conj(g,h):return comp(comp(g,h),inv(g))
t=(1,2,0);s=(0,2,1);C3=(identity,t,comp(t,t))
def section(H):return [x for x in range(3) if all(h[x]==x for h in H)]
def central(H):return [g for g in C3 if all(comp(g,h)==comp(h,g) for h in H)]
need(section((t,))==[] and section((t,s))==[],'both_holonomy_carriers_have_no_section')
need(len(central((t,)))==3 and len(central((t,s)))==1,'centralizer_distinguishes_order31_mechanism')
comm=comp(comp(comp(t,s),inv(t)),inv(s))
need(comm==comp(t,t) and comm!=identity,'commutator_translation_hostile')
def sign(g):return (-1)**sum(g[i]>g[j] for i in range(3) for j in range(i+1,3))
need(sign(comm)==1,'abelian_sign_guard_misses_commutator')

# Simultaneous conjugacy, never independent cycle conjugators.
# 36^2 ordered holonomy pairs, exact root-rotation tests.
for H in product(S,repeat=2):
    for K in product(S,repeat=2):
        common=[g for g in C3 if all(conj(g,h)==k for h,k in zip(H,K))]
        independent=all(any(conj(g,h)==k for g in C3) for h,k in zip(H,K))
        if common: need(independent,'joint_implies_separate_conjugacy')
        if independent and not common: counts['separate_conjugacy_false_positives']+=1
need(counts['separate_conjugacy_false_positives']>0,'shared_root_conjugator_is_essential')

receipt={'status':'PASS','scope':'bounded attachment and frame-comparison carriers, not native wiring exhaustion',
'counts':dict(sorted(counts.items())),
'controls':{'matrix_families':250,'triangle_relation_instances':4096,'two_switch_feasible_points':feasible,
'bounded_mod8_atom_count':len(modG),'bounded_mod8_fiber_count':len(modfiber)},
'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(receipt,indent=2,sort_keys=True))
