#!/usr/bin/env python3
"""Literal finite replay of U21's 128 internal triangle twists; no HP DP."""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import hashlib
import json

checks = 0

def require(test, label):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(label)


V = [(a,b) for a in range(7) for b in range(3)]
INDEX = {x:i for i,x in enumerate(V)}
S = ({1,2,3},{0,1,2,4},{1,2,4})
D = {1,2,4}
FULL = (1<<21)-1
PAIRS = list(combinations(range(21),2))


def base_arc(x,y):
    a,b=x
    c,d=y
    # x^{-1} y = (2^{-b}(c-a), d-b) in C7 semidirect C3.
    return ((pow(pow(2,b,7),-1,7)*(c-a))%7) in S[(d-b)%3]


def rows(mask):
    result=[]
    for i,x in enumerate(V):
        value=0
        for j,y in enumerate(V):
            if i==j:
                continue
            arc=base_arc(x,y)
            if x[0]==y[0] and mask>>x[0]&1:
                arc=not arc
            value |= int(arc)<<j
        result.append(value)
    return result


def bit(r,i,j):
    return r[i]>>j&1


def reachable(r,start):
    seen=1<<start
    while True:
        nxt=seen
        for i in range(21):
            if seen>>i&1:
                nxt |= r[i]
        if nxt==seen:
            return seen
        seen=nxt


def paley(a,c):
    return (c-a)%7 in D


def image_mask(mask,p):
    return sum(1<<p[a] for a in range(7) if mask>>a&1)


ordinary_maps=[tuple((t+pow(2,k,7)*a)%7 for a in range(7)) for t in range(7) for k in range(3)]
anti_maps=[tuple((t-pow(2,k,7)*a)%7 for a in range(7)) for t in range(7) for k in range(3)]
actual_auto=[]
actual_anti=[]
for p in permutations(range(7)):
    if all(paley(a,c)==paley(p[a],p[c]) for a,c in combinations(range(7),2)):
        actual_auto.append(p)
    if all(paley(a,c)!=paley(p[a],p[c]) for a,c in combinations(range(7),2)):
        actual_anti.append(p)
require(set(actual_auto)==set(ordinary_maps),'full Paley quotient automorphisms')
require(set(actual_anti)==set(anti_maps),'full Paley quotient anti-maps')

allrows=[rows(mask) for mask in range(128)]
cross_max=0
common_values=Counter()
for mask,r in enumerate(allrows):
    require(all(r[i]>>i&1==0 for i in range(21)),'no loops')
    require(all(bit(r,i,j)+bit(r,j,i)==1 for i,j in PAIRS),'tournament relation')
    require([x.bit_count() for x in r]==[10]*21,'all twists regular')
    reverse=[FULL ^ (1<<i) ^ r[i] for i in range(21)]
    require(reachable(r,0)==FULL and reachable(reverse,0)==FULL,'all twists strongly connected')
    for i,j in PAIRS:
        co=(r[i]&r[j]).bit_count()
        ci=(reverse[i]&reverse[j]).bit_count()
        require(co==ci,'common in/out equality in regular tournament')
        require((co==7)==(V[i][0]==V[j][0]),'weight7 recovers exactly triangle fibres')
        if V[i][0]!=V[j][0]:
            cross_max=max(cross_max,co)
        common_values[co]+=1
    for a in range(7):
        triangle=[INDEX[(a,b)] for b in range(3)]
        require(all(sum(bit(r,i,j) for j in triangle)==1 for i in triangle),'each internal triangle directed')
        ports=[]
        for c in range(7):
            if not paley(a,c):
                continue
            between=[(b,e) for b in range(3) for e in range(3)
                     if bit(r,INDEX[(a,b)],INDEX[(c,e)])]
            require(len(between)==8,'Paley majority is exactly eight to one')
            dissent=set((b,e) for b in range(3) for e in range(3))-set(between)
            expected_b=next(b for b in range(3) if pow(2,b+2,7)==(c-a)%7)
            require(dissent=={(expected_b,expected_b)},'literal dissent-port formula')
            ports.append(expected_b)
        require(sorted(ports)==[0,1,2],'outgoing majority arcs distinguish all three ports')
    for t in range(7):
        for k in range(3):
            a_map=tuple((t+pow(2,k,7)*a)%7 for a in range(7))
            lift=[INDEX[(a_map[a],(b+k)%3)] for a,b in V]
            target=allrows[image_mask(mask,a_map)]
            require(all(bit(r,i,j)==bit(target,lift[i],lift[j]) for i,j in PAIRS),'every affine mask map literally lifts')
            anti=tuple((t-pow(2,k,7)*a)%7 for a in range(7))
            anti_lift=[INDEX[(anti[a],(b+k)%3)] for a,b in V]
            anti_target_mask=127^image_mask(mask,anti)
            target=allrows[anti_target_mask]
            require(all(bit(r,i,j)!=bit(target,anti_lift[i],anti_lift[j]) for i,j in PAIRS),'inversion-complement literally anti-lifts')
            require(anti_target_mask.bit_count()==7-mask.bit_count(),'anti-map complements odd-size mask universe')
            require(anti_target_mask.bit_count()!=mask.bit_count(),'no mask can be self-converse')

ordinary={min(image_mask(mask,p) for p in ordinary_maps) for mask in range(128)}
merged={min([image_mask(mask,p) for p in ordinary_maps]+[127^image_mask(mask,p) for p in anti_maps]) for mask in range(128)}
require(len(ordinary)==12,'twelve ordinary mask orbits')
require(merged=={0,1,3,7,11,13},'six merged mask orbits and representatives')
require(sum(sum(image_mask(mask,p)==mask for mask in range(128)) for p in ordinary_maps)==21*12,'Burnside direct fixed-mask sum')
cycle_types=Counter()
for p in ordinary_maps:
    left=set(range(7)); lengths=[]
    while left:
        start=min(left); here=start; length=0
        while here in left:
            left.remove(here); length+=1; here=p[here]
        lengths.append(length)
    cycle_types[tuple(sorted(lengths))]+=1
require(cycle_types==Counter({(1,1,1,1,1,1,1):1,(7,):6,(1,3,3):14}),'Burnside cycle types')

# The seven Fano lines from the named F8 polynomial, independently of orbit types.
def f8mul(a,b):
    z=0
    while b:
        if b&1: z ^= a
        a <<= 1
        if a&8: a ^= 11
        b >>= 1
    return z
powers=[1]
for k in range(6): powers.append(f8mul(powers[-1],2))
require(len(set(powers))==7 and f8mul(powers[-1],2)==1,'F8 alpha order7')
lines={frozenset((a+d)%7 for d in D) for a in range(7)}
require(len(lines)==7,'seven translated Fano lines')
for line in lines:
    values=[powers[a] for a in line]
    require(values[0]^values[1]^values[2]==0,'field-vector Fano dependence')
for a,c in combinations(range(7),2):
    require(sum(a in L and c in L for L in lines)==1,'unique Fano line through pair')
triple_types=[]
for mask in (7,11,13):
    vertices=[a for a in range(7) if mask>>a&1]
    scores=sorted(sum(paley(a,c) for c in vertices if c!=a) for a in vertices)
    orbit={image_mask(mask,p) for p in ordinary_maps}
    triple_types.append({'mask':mask,'type':('transitive' if scores==[0,1,2] else 'Fano-line cycle' if frozenset(vertices) in lines else 'nonline cycle'),'orbit_size':len(orbit)})
require([r['type'] for r in triple_types]==['transitive','Fano-line cycle','nonline cycle'],'three mask-size3 types')
require([r['orbit_size'] for r in triple_types]==[21,7,7],'all35 triples accounted for')

# Arithmetic checks on incoming H values, explicitly NOT a new H computation.
inherited=[(0,21,123376391895435,5875066280735),(1,3,123381234928431,41127078309477),
 (3,1,123384468287435,123384468287435),(7,1,123386086601285,123386086601285),
 (11,3,123385527546237,41128509182079),(13,3,123386644838787,41128881612929)]
weights=[]
for mask,aut,H,degree in inherited:
    stabilizer=sum(image_mask(mask,p)==mask for p in ordinary_maps)
    require(stabilizer==aut,'literal mask stabilizer matches claimed full automorphism order')
    require(H==aut*degree,'inherited H divided by exact automorphism order')
    weights.append({'mask':mask,'automorphism_order':aut,'inherited_H_not_recomputed':H,
                    'fundamental_black_degree_arithmetic':degree,'full_merged_black_degree_arithmetic':2*degree})
arc_counts=[585723718963,592255137742,594981636439,571650274131,584204812473,
585280890671,590182955148,585182864820,595482271370,590121718978]
require(sum(arc_counts)==inherited[0][3],'incoming arc incidences sum to H/21 arithmetic only')

print(json.dumps({'status':'FINITE-EXACT literal128-twist replay; inherited H only arithmetic-checked',
 'checks':checks,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'vertices':'(a,b), a mod7, b mod3; index3a+b', 'twists':128,
 'all_regular_outdegree':10,'all_strong':True,'canonical_pair_weight':7,'maximum_cross_pair_common_out':cross_max,
 'common_out_histogram_across_all_twists':dict(sorted(common_values.items())),
 'quotient_auto_order':len(actual_auto),'quotient_anti_maps':len(actual_anti),
 'ordinary_orbit_count':len(ordinary),'ordinary_representatives':sorted(ordinary),'merged_representatives':sorted(merged),
 'self_converse_count':0,'three_point_types':triple_types,'weights':weights,
 'H_scope':'No Hamilton-path subset DP or large path-count verification was performed.',
 'family_boundary':'128 internal triangle-orientation masks on one fixed U21 port lift; not the30 labelled Fano-plane/simplex-code family.'},indent=2,sort_keys=True))
