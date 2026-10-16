#!/usr/bin/env python3
"""Exact palette graph, component flips, and combinatorial hostile controls.

All 2^15 admissible-pair banks are covered by matching/vertex-cover witnesses.
The Mycielski control is an abstract graph, NOT an embedding in Q(zeta55).
"""
from collections import Counter
from itertools import combinations, product
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
PAIRS=list(combinations(range(6),2))
BITS=[sum(1<<a for a in p) for p in PAIRS]
EDGES=[]
for i,a in enumerate(BITS):
    for j,b in enumerate(BITS):
        if a&b==0: EDGES.append((i,j,BITS.index(63^a^b)))


def need(p,m):
    if not p: raise AssertionError(m)


def matching_cover(bank):
    adj=[[j for a,j,q in EDGES if a==i and (bank>>q)&1] for i in range(15)]
    right={}
    def augment(i,seen):
        for j in adj[i]:
            if j in seen: continue
            seen.add(j)
            if j not in right or augment(right[j],seen):
                right[j]=i;return True
        return False
    for i in range(15): augment(i,set())
    left={i:j for j,i in right.items()}
    zl=set(range(15))-left.keys();zr=set();todo=list(zl)
    while todo:
        i=todo.pop()
        for j in adj[i]:
            if left.get(i)==j or j in zr: continue
            zr.add(j)
            if j in right and right[j] not in zl:
                zl.add(right[j]);todo.append(right[j])
    cl=set(range(15))-zl;cr=zr
    need(len(cl)+len(cr)==len(right),'matching/cover equality')
    need(all(i in cl or j in cr for i in range(15) for j in adj[i]),'actual vertex cover')
    need(len(set(left.values()))==len(left) and all(j in adj[i] for i,j in left.items()),'actual matching')
    return left,cl,cr


def permitted(components):
    bank=0
    for q,(a,b) in enumerate(PAIRS):
        if all((not(S>>a&1) and not(T>>b&1)) or
               (not(S>>b&1) and not(T>>a&1)) for S,T in components): bank|=1<<q
    return bank


def graph_control():
    # Core: the triangle-free, four-chromatic Mycielskian of C5.
    edges=set()
    def edge(a,b): edges.add(tuple(sorted((a,b))))
    for i in range(5):
        j=(i+1)%5;edge(i,j);edge(i,5+j);edge(j,5+i);edge(10,5+i)
    colours=[0,1,0,1,2]*2+[3]
    need(all(colours[a]!=colours[b] for a,b in edges),'supplied proper four-colouring')
    adj=[set() for _ in range(11)]
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    need(not any(b in adj[c] for a in range(11) for b,c in combinations(adj[a],2)), 'triangle-free core')
    order=sorted(range(11),key=lambda i:-len(adj[i]));assignment={}
    def three(j):
        if j==11:return True
        i=order[j]
        for c in range(3):
            if all(assignment.get(a)!=c for a in adj[i]):
                assignment[i]=c
                if three(j+1):return True
                del assignment[i]
        return False
    need(not three(0),'independent exhaustive non-three-colourability')
    # Six single-port leaves in each exterior colour class; no exterior edges
    # are added to them. Their boundary labels are all six colours per class.
    forbidden=[set() for _ in range(11)]
    for cls in range(4):
        for c in range(6):colours.append(cls);forbidden.append({c})
    adj.extend(set() for _ in range(24))
    unseen={i for i,c in enumerate(colours) if c in [2,3]};comps=[];masks=[]
    while unseen:
        start=min(unseen);todo=[start];unseen.remove(start);component=[]
        while todo:
            i=todo.pop();component.append(i)
            for j in adj[i]:
                if j in unseen:unseen.remove(j);todo.append(j)
        comps.append(component)
        masks.append(tuple(sum(1<<c for c in set().union(*(forbidden[i] for i in component if colours[i]==side)))
                           for side in [2,3]))
    bank=permitted(masks);need(bank==(1<<15)-1,'component flips retain every pair')
    for cls in range(4):
        need(set().union(*(forbidden[i] for i,c in enumerate(colours) if c==cls))==set(range(6)),
             'no class is port-free or has a globally unused colour')
    # Literal six-colour lift through one palette edge.
    a,b,q=next(e for e in EDGES if bank>>e[2]&1)
    result=[None]*len(colours)
    for i,c in enumerate(colours):
        if c<2:result[i]=next(x for x in PAIRS[[a,b][c]] if x not in forbidden[i])
    c,d=PAIRS[q]
    for component in comps:
        for x,y in [(c,d),(d,c)]:
            if all((x if colours[i]==2 else y) not in forbidden[i] for i in component):
                for i in component:result[i]=x if colours[i]==2 else y
                break
        else:raise AssertionError('missing component lift')
    need(all(result[a]!=result[b] for a,b in edges) and
         all(result[i] not in forbidden[i] for i in range(len(colours))),'literal full exterior/port lift')
    # Same actual eight-port edge relation as the inherited geometry.
    ports=[0,1,2,3,4,5,1,2]
    pe=[(i,(i+1)%5) for i in range(5)]+[(0,5),(5,6),(6,7),(7,0)]
    need(all(ports[a]!=ports[b] for a,b in pe),'proper six-colouring of the unchanged port graph')
    return dict(status='COMBINATORIAL ONLY; no unit-distance realization asserted',
                exterior_vertices=len(colours),core_edges=len(edges),admissible_pairs=PAIRS,
                supplied_four_colours=colours,lift=result,
                component_count=len(comps),port_colours=ports)


def geometric_pair_bank():
    source=ROOT/'05-knowledge/results/six_colour_core_doubling_20261012.json'
    inherited=json.loads(source.read_text())
    D=inherited['apex_port_pairs']
    expected=sorted([sorted((i,(i+2)%5)) for i in range(5)]+
                    [[1,5],[1,7],[4,5],[4,7]])
    need(sorted(D)==expected,'exact inherited nine apex port-pair incidences')
    B=[(i,(i+1)%5) for i in range(5)]+[(0,5),(5,6),(6,7),(7,0)]
    patterns=[]
    def extend(c):
        i=len(c)
        if i==8:patterns.append(c);return
        for x in range(min(6,max(c,default=-1)+2)):
            if all(c[j]!=x for a,b in B for j in
                   ([b] if a==i else [a] if b==i else []) if j<i):extend(c+[x])
    extend([])
    counts=Counter(max(c)+1 for c in patterns)
    need(dict(counts)=={3:30,4:180,5:247,6:113},'complete proper boundary-pattern universe')
    exceptions=[];hist=Counter();assignments=0;distinct_assignments=0
    for c in patterns:
        types=[(i,PAIRS.index(tuple(sorted((c[a],c[b])))))
               for i,(a,b) in enumerate(D) if c[a]!=c[b]]
        mult=Counter(p for i,p in types)
        eligible=[]
        # Independent capacity criterion for blocking all six options of
        # a SINGLE complementary pair. Its three complementary-pair groups
        # are disjoint, and each needs two physical blocker copies.
        for q in range(15):
            groups=[(i,j) for i,j,k in EDGES if k==q and i<j]
            if all(mult[i]+mult[j]>=2 for i,j in groups):eligible.append(q)
        need(len(eligible)<=1,'at most one geometrically blockable pair per boundary pattern')
        seen=set();possible=0;witness=None
        for mask in range(1<<len(types)):
            assignments+=1;left=right=0;classes=[None]*9
            for j,(i,p) in enumerate(types):
                if mask>>j&1:left|=1<<p;classes[i]=1
                else:right|=1<<p;classes[i]=2
            if (left,right) in seen:continue
            seen.add((left,right));distinct_assignments+=1
            live=0
            for i,j,q in EDGES:
                if not(left>>i&1 or right>>j&1):live|=1<<q
            dead=((1<<15)-1)^live
            need(dead.bit_count()<=1,'no assignment of the actual nine points blocks two pairs')
            possible|=dead
            if dead:witness=dict(classes=classes,left=[PAIRS[i] for i in range(15) if left>>i&1],
                                 right=[PAIRS[i] for i in range(15) if right>>i&1])
        need(possible==sum(1<<q for q in eligible),'capacity test agrees with every literal blocker assignment')
        hist[(max(c)+1,len(eligible))]+=1
        if eligible:exceptions.append(dict(port_colours=c,exceptional_pair=PAIRS[eligible[0]],
                                           blocker_witness=witness))
    independent=[]
    for mask in range(256):
        if all(not(mask>>a&1 and mask>>b&1) for a,b in B):
            cut=sum(bool(mask>>a&1)!=bool(mask>>b&1) for a,b in D)
            independent.append((mask,cut))
    need(len(independent)==46 and max(v for m,v in independent)==6,
         'one boundary colour is incident with at most six distinct-colour blockers')
    return dict(status='FINITE-EXACT; full nine-point incidence bank, arbitrary C1/C2 assignments',
                inherited_source=dict(path=str(source.relative_to(ROOT)),
                                      sha256=hashlib.sha256(source.read_bytes()).hexdigest()),
                apex_port_pairs=D,patterns=len(patterns),raw_assignments=assignments,
                distinct_class_pair_banks=distinct_assignments,exceptional_patterns=exceptions,
                histogram=[dict(colours=a,exceptional_pairs=b,patterns=n) for (a,b),n in sorted(hist.items())],
                independent_port_sets=len(independent),max_one_colour_blockers=6)


def small_boundary_controls():
    count=0
    for a,b,c,d in product(range(12),repeat=4):
        if a+b+c+d<=11:
            need(sum(sorted([a,b,c,d])[:2])<=5,'two least loaded classes')
            count+=1
    q=PAIRS.index((4,5));options=[(i,j) for i,j,k in EDGES if k==q]
    slots=sorted({(0,i) for i,j in options}|{(1,j) for i,j in options})
    banks=0
    for n in range(6):
        for bank in combinations(slots,n):
            forbidden=set(bank)
            need(any((0,i) not in forbidden and (1,j) not in forbidden for i,j in options),
                 'five blockers cannot destroy six complementary options')
            banks+=1
    # Complete four-partite hostile: every part has one single-port vertex
    # forbidding each used boundary colour 0..4. Assigning all six colours
    # to disjoint part palettes loses no possible colouring.
    surviving=0
    for owners in product(range(4),repeat=6):
        palettes=[{c for c,p in enumerate(owners) if p==i} for i in range(4)]
        surviving+=all(all(p-{c} for c in range(5)) for p in palettes)
    need(surviving==0,'five-colour boundary hostile with no double-port vertex')
    palettes=[{4},{5},{0,1},{2,3}]
    need(all(all(p-{c} for c in range(4)) for p in palettes),'four-colour positive control')
    return dict(load_compositions=count,blocker_banks=banks,
                five_colour_hostile_palette_assignments=4**6,
                five_colour_hostile_surviving=surviving)


if __name__=='__main__':
    hist=Counter();passing=0
    for bank in range(1<<15):
        matching,cl,cr=matching_cover(bank)
        qs=[q for q in range(15) if bank>>q&1]
        disjoint=any(BITS[a]&BITS[b]==0 for a,b in combinations(qs,2))
        star=any(all(bank>>q&1 for q,p in enumerate(PAIRS) if c in p) for c in range(6))
        need((len(matching)>=10)==(disjoint or star),'complete nine-blocker classification')
        hist[(len(qs),len(matching))]+=1;passing+=len(matching)>=10
    examples={}
    for name,qs in [('disjoint',[(0,1),(2,3)]),
                    ('star4',[(0,1),(0,2),(0,3),(0,4)]),
                    ('star5',[(0,i) for i in range(1,6)]),
                    ('triangle',[(0,1),(0,2),(1,2)])]:
        bank=sum(1<<PAIRS.index(q) for q in qs);m,cl,cr=matching_cover(bank)
        examples[name]=dict(bank=qs,matching=[[PAIRS[i],PAIRS[j]] for i,j in sorted(m.items())],
                            cover_left=[PAIRS[i] for i in sorted(cl)],cover_right=[PAIRS[i] for i in sorted(cr)])
    # Three-vertex path: both endpoints are in the same colour class.
    need(not(permitted([(3,0)])&1),'two locally nonempty lists can disagree on a component')
    need(permitted([(1,2),(2,1)])&1,'separate components can choose opposite orientations')
    out=ROOT/'05-knowledge/results/six_colour_component_pairs_20261016a.json'
    result=dict(status='FINITE-EXACT palette and abstract-graph controls',banks=1<<15,
                nine_blocker_safe_banks=passing,
                histogram=[dict(pairs=a,matching=b,count=n) for (a,b),n in sorted(hist.items())],
                examples=examples,graph_control=graph_control(),geometry=geometric_pair_bank(),
                small_boundary=small_boundary_controls())
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('banks',1<<15,'safe',passing,'examples',examples)
    print('geometry',result['geometry']['histogram'],'raw assignments',result['geometry']['raw_assignments'])
    print('sha256',hashlib.sha256(out.read_bytes()).hexdigest())
