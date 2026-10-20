#!/usr/bin/env python3
"""Independent bounded audit; imports no producer modules and writes only stdout."""
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CHECKS = Counter()

def check(ok, label):
    CHECKS[label] += 1
    if not ok:
        raise RuntimeError(label)

def odd(g):
    return sum(a > b for i, a in enumerate(g) for b in g[i+1:]) % 2

def edge_mask(word, edges):
    return sum(1 << edges.index(tuple(sorted(e))) for e in zip(word, word[1:]))

# Six-bit graph carrier, independent of the producer's frozenset carrier.
E4 = tuple(combinations(range(4), 2))
def graph_action(m, g):
    return sum(1 << E4.index(tuple(sorted((g[a], g[b]))))
               for i, (a, b) in enumerate(E4) if m >> i & 1)
def cp(m):
    return min(m, 63 ^ m)
def opposition(m):
    return sum(1 << E4.index(tuple(v for v in range(4) if v not in e))
               for i, e in enumerate(E4) if m >> i & 1)
def letter_op(e):
    return e if e[0] == e[1] else tuple(v for v in range(4) if v not in e)
def phi(m):
    values = set()
    for z in (m, 63 ^ m):
        deg = [sum(v in e for i, e in enumerate(E4) if z >> i & 1) for v in range(4)]
        if 3 in deg:
            values.add((deg.index(3),) * 2)
        elif 0 in deg:
            values.add((deg.index(0),) * 2)
        else:
            for q in permutations(range(4)):
                if edge_mask(q, E4) == z:
                    values.add(tuple(sorted((q[1], q[2]) if odd(q) else (q[0], q[3]))))
    check(len(values) == 1, 'all_four_representatives_give_one_phi')
    return values.pop()

X4 = sorted({cp(m) for m in range(64) if m.bit_count() == 3})
PHI = {m: phi(m) for m in X4}
check(len(set(PHI.values())) == 10, 'tetrahedral_bijection')
for g in permutations(range(4)):
    for m in X4:
        y = tuple(sorted(g[a] for a in PHI[m]))
        check(PHI[cp(graph_action(m, g))] == (letter_op(y) if odd(g) else y), 'full_S4_twist')
for m in X4:
    check(PHI[cp(opposition(m))] == letter_op(PHI[m]), 'tetrahedral_opposition')
fixed = sum(cp(opposition(m)) == m for m in X4)
check(fixed == 4, 'four_singletons_three_doubletons')
tau = {2:4,4:2,3:12,12:3,5:20,20:5,6:6,10:10,15:15,30:30}
B = {2:(0,1),3:(0,2),5:(0,3),4:(2,3),12:(1,3),20:(1,2),30:(0,0),15:(1,1),10:(2,2),6:(3,3)}
invphi = {v:k for k,v in PHI.items()}
arithmetic_graph_rows = []
for d in sorted(B):
    g = invphi[B[d]]
    check(invphi[B[tau[d]]] == cp(opposition(g)), 'whole_arithmetic_graph_conjugacy')
    arithmetic_graph_rows.append({'divisor':d,'letter':B[d],'complement_pair_masks':[g,63^g]})
check(6 % 2 == 0 and tau[6] % tau[2] != 0, 'divisibility_hostile_2_divides_6')
check(all(60//d != d for d in B), 'divisor_complement_no_fixed_points')
trans = (1,0,2,3)
check(sum(cp(graph_action(g,trans)) == g for g in X4) == 2, 'S4_source_fixed_count')
check(sum(tuple(sorted(trans[a] for a in e)) == e for e in PHI.values()) == 4, 'S4_target_fixed_count')
graph_types = {min(graph_action(m,g) for g in permutations(range(4))) for m in range(64)}
check(len(graph_types) == 11, 'eleven_unlabelled_four_graphs')

def canon(v):
    return min(tuple(x % 5 for x in v), tuple(-x % 5 for x in v))
def determinant(v,w):
    return (v[0]*w[1]-v[1]*w[0]) % 5
def adjacent(v,w):
    return determinant(v,w) in (1,4)

V = sorted({canon(v) for v in product(range(5),repeat=2) if any(v)})
N,S = canon((1,0)),canon((2,0))
U = [canon((i,1)) for i in range(5)]
L = [canon((2*i-1,2)) for i in range(5)]
CODE = {N:947,S:76,**dict(zip(U,[613,217,710,341,668])),**dict(zip(L,[313,682,355,410,806]))}
check(len(CODE) == 12 and set(CODE) == set(V), 'literal_code_table_complete')
distance_rows=[]
for v in V:
    dist={v:0}; todo=deque([v])
    while todo:
        a=todo.popleft()
        for b in V:
            if b not in dist and adjacent(a,b):
                dist[b]=dist[a]+1;todo.append(b)
    check(Counter(dist.values()) == {0:1,1:5,2:5,3:1}, 'icosahedral_distance_profile')
    for w in V:
        check((dist[w]==3)==(w==canon((2*v[0],2*v[1]))), 'antipode_exactly_distance_three')
        check((dist[w]==2)==(determinant(v,w) in (2,3)), 'determinant_four_exactly_distance_two')
    distance_rows.append({'vertex':v,'distances':[[w,dist[w]] for w in V]})

# Tournament matrices and physical path-complement returns, not just XOR distance.
E5=tuple(combinations(range(5),2))
def beats(c,a,b,edges=E5):
    bit=(c>>edges.index(tuple(sorted((a,b)))))&1
    return bit==int(a<b)
def relabel(c,g,edges=E5):
    inv=[g.index(i) for i in range(len(g))]
    return sum(int(beats(c,inv[a],inv[b],edges))<<i for i,(a,b) in enumerate(edges))
def paths(c):
    return [q for q in permutations(range(5)) if all(beats(c,a,b) for a,b in zip(q,q[1:]))]
def regular(c):
    return all(sum(beats(c,a,b) for b in range(5) if b!=a)==2 for a in range(5))
def complement_path(c,q):
    return c ^ (1023 ^ edge_mask(q,E5))

check(len({min(relabel(c,g,E4) for g in permutations(range(4))) for c in range(64)})==4,'four_unlabelled_four_tournaments')
legal_return_rows=[]
for v,c in CODE.items():
    check(regular(c),'literal_code_regular')
    check(c.bit_count()%2==CODE[N].bit_count()%2,'literal_code_common_parity_sheet')
    weights=Counter()
    for p in paths(c):
        middle=complement_path(c,p)
        for q in paths(middle):
            target=complement_path(middle,q)
            same_endpoints=(p[0],p[-1])==(q[0],q[-1])
            check(regular(target)==same_endpoints,'physical_return_endpoint_criterion')
            if regular(target):weights[target]+=1
    check(weights[c]==15 and sum(weights.values())==25,'return_presentation_weights')
    neighbours={x for x in weights if x!=c}
    check(neighbours=={CODE[w] for w in V if adjacent(v,w)},'physical_return_equals_field_adjacency')
    check(all(weights[x]==2 for x in neighbours),'each_real_edge_has_two_presentations')
    legal_return_rows.append({'code':c,'weights':dict(sorted(weights.items()))})

# Row-major matrices and explicit generator closure.
def mv(m,v):
    a,b,c,d=m
    return ((a*v[0]+b*v[1])%5,(c*v[0]+d*v[1])%5)
def mm(m,n):
    a,b,c,d=m;e,f,g,h=n
    return ((a*e+b*g)%5,(a*f+b*h)%5,(c*e+d*g)%5,(c*f+d*h)%5)
matrices=list(product(range(5),repeat=4))
bulk=[m for m in matrices if (m[0],m[2])!=(0,0) and (m[1],m[3])!=(0,0)]
dc=Counter((a*d-b*c)%5 for a,b,c,d in bulk)
check(len(matrices)==625 and len(bulk)==576 and dc=={0:96,1:120,2:120,3:120,4:120},'independent_625_matrix_census')
SL=[m for m in matrices if (m[0]*m[3]-m[1]*m[2])%5==1]
actions={tuple(CODE[canon(mv(m,v))] for v in V) for m in SL}
even_actions={tuple(relabel(CODE[v],g) for v in V) for g in permutations(range(5)) if not odd(g)}
check(len(actions)==60 and actions==even_actions,'whole_SL2_action_equals_A5')
generators=[((1,1,0,1),(2,3,4,0,1)),((0,4,1,0),(3,4,2,0,1))]
for m,g in generators:
    for v in V:check(CODE[canon(mv(m,v))]==relabel(CODE[v],g),'literal_generator_correspondence')
generated={(1,0,0,1)};todo=list(generated)
while todo:
    m=todo.pop()
    for g,_ in generators:
        z=mm(g,m)
        if z not in generated:generated.add(z);todo.append(z)
check(generated==set(SL),'generators_reach_all_120_matrices')
for v in V:
    check(CODE[canon((2*v[0],2*v[1]))]==(1023^CODE[v]),'literal_antipode_converse')
    check(CODE[canon((v[1],v[0]))]==(1023^relabel(CODE[v],(1,0,2,4,3))),'literal_reflection_converse_permutation')

# Start with literal seven-vertex paths, then recover the runs rather than encode them.
def sb(a,b):
    return (a<b) ^ (tuple(sorted((a,b))) in [(0,3),(3,6)])
HP=[q for q in permutations(range(7)) if all(sb(a,b) for a,b in zip(q,q[1:]))]
label={frozenset():1,frozenset([1]):2,frozenset([2]):3,frozenset([1,2]):4,frozenset([0,1,2]):0}
decoded={}
for q in HP:
    before=set(q[:q.index(3)])
    I=frozenset(before&{0,1,2})
    K=frozenset({0,1,2}-{6-x for x in before&{4,5,6}})
    check(I in label and K in label,'direct_path_has_admissible_source_states')
    v=(label[I],label[K]);decoded[v]=q
check(len(HP)==25 and len(decoded)==25,'independent_direct_25_path_decoder')
for v,q in decoded.items():
    check(tuple(6-x for x in reversed(q))==decoded[(v[1],v[0])],'actual_source_R_is_swap')
    check(tuple(-x%5 for x in v) in decoded,'source_J_has_actual_path_target')
blue=sum(tuple(6-x for x in reversed(q))==q for q in HP)
check(blue==5,'five_actual_blue_paths')
check(decoded[(0,0)]==tuple(range(7)),'standard_path_is_zero')
automorphisms=[g for g in permutations(range(7)) if all(sb(a,b)==sb(g[a],g[b]) for a,b in combinations(range(7),2))]
check(automorphisms==[tuple(range(7))],'source_tournament_is_rigid')
balanced=json.loads((ROOT/'session17-balanced-fibres/marked_f5_paths.json').read_text())
check({tuple(r['vector']):tuple(r['path']) for r in balanced}==decoded,'separate_balanced_path_file_exact_match')
field_receipt=json.loads((HERE/'frozen/field_receipt.json').read_text())
check({tuple(r['field_vector']):tuple(r['path']) for r in field_receipt['literal_source_paths']}==decoded,'producer_path_file_exact_match')
markings=[]
for q in permutations(range(1,5)):
    h=(0,)+q
    if all(h[-x%5]==-h[x]%5 for x in range(5)):
        markings.append(h)
check(len(markings)==8,'eight_allowed_markings')
cube=lambda v:tuple(x**3%5 for x in v)
check(adjacent((1,1),(2,1)) and not adjacent(cube((1,1)),cube((2,1))),'literal_cubing_hostile')

# The three small named Collatz block claims, independently recomputed exactly.
def block_value(word,x):
    for a in sorted(word,reverse=True):x=(3*x+1)/2**a
    return x
block_rows=[]
for n,expected in [(27,4),(63,5),(75,1)]:
    x=n;word=[]
    while x>=n:
        x=3*x+1;a=0
        while x%2==0:x//=2;a+=1
        word.append(a)
    best={(0,0):Fraction(n)};record=[]
    for b in range(1,expected+1):
        for j in range(b,len(word)+1):
            best[b,j]=min(block_value(word[i:j],best[b-1,i]) for i in range(b-1,j) if (b-1,i) in best)
        record.append(str(best[b,len(word)]))
    check(all(best[b,len(word)]>=n for b in range(1,expected)) and best[expected,len(word)]<n,'named_minimum_block_claim')
    block_rows.append({'source':n,'first_descent_depth':len(word),'endpoint':x,'minimum_blocks':expected,'best_envelopes':record})

manifest=json.loads((HERE/'frozen/publication_manifest.json').read_text())
for item in manifest['files']:
    data=(HERE/'frozen'/Path(item['source']).name).read_bytes()
    check(len(data)==item['bytes'] and hashlib.sha256(data).hexdigest()==item['sha256'],'exact_frozen_manifest')
print(json.dumps({'status':'PASS','checks':dict(CHECKS),'total_checks':sum(CHECKS.values()),'arithmetic_graph_rows':arithmetic_graph_rows,'distance_rows':distance_rows,'physical_return_rows':legal_return_rows,'collatz_controls':block_rows,'eight_markings':markings},sort_keys=True,indent=2))
