#!/usr/bin/env python3
"""Independent frozen-packet audit. Standard library only; no producer imports.

Inputs: --source frozen packet, --repo checkout containing its dependency pins.
This script authenticates, reconstructs, and parses literal witnesses. It never
writes its inputs. JSON stdout is deterministic and unchanged by python -O.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, permutations, product
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import re

FROZEN = 'a5490b17a42062db222676a44a1eb2ef1c3f0b55499b2389cc70691886f9112c'

def require(condition, message):
    if not condition:
        raise ValueError(message)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def pairings(vertices):
    if not vertices:
        yield ()
    else:
        x = vertices[-1]
        for j,y in enumerate(vertices[:-1]):
            for tail in pairings(vertices[:j]+vertices[j+1:-1]):
                yield tuple(sorted(tail+((y,x),)))

def stars(k):
    vertices = tuple(range(2*k-1))
    for middle in vertices:
        for leaves in combinations([v for v in vertices if v != middle], 2):
            used = {middle,*leaves}
            star = tuple(tuple(sorted((middle,v))) for v in leaves)
            for matching in pairings(tuple(v for v in vertices if v not in used)):
                yield tuple(sorted(star+matching))

def cut_test(edges,n):
    # Independent cut-coverage difference array.
    difference = [0]*(n+1)
    for a,b in edges:
        difference[a+1] += 1
        difference[b+1] -= 1
    current = 0
    for cut in range(1,n):
        current += difference[cut]
        if current == 0:
            return False
    return True

def words(total,parts,cap):
    if total < 0 or total > parts*cap:
        return
    if parts == 0:
        if total == 0:
            yield ()
        return
    for value in range(min(total,cap)+1):
        for tail in words(total-value,parts-1,cap):
            yield (value,)+tail

def tournament(edges,n,gaps,strong):
    gs = ((0,)+gaps+(0,)) if strong else gaps
    labels = []; weights = []; endpoint = {}
    for j in range(n+1):
        if gs[j]:
            labels.append(('gap',j)); weights.append(gs[j])
        if j<n:
            endpoint[j] = len(labels)
            labels.append(('endpoint',j)); weights.append(1)
    rows = [sum(1<<j for j in range(i+1,len(labels))) for i in range(len(labels))]
    for a,b in edges:
        x,y = endpoint[a],endpoint[b]
        rows[x] &= ~(1<<y)
        rows[y] |= 1<<x
    return tuple(rows),tuple(weights)

def homogeneous_pairs(rows):
    n=len(rows)
    return [(a,b) for a,b in combinations(range(n),2)
            if all(((rows[a]>>v)&1)==((rows[b]>>v)&1)
                   for v in range(n) if v not in (a,b))]

def strong_actual(rows):
    # Warshall reachability, independently of the endpoint-cut test.
    reach = [r|(1<<i) for i,r in enumerate(rows)]
    for v in range(len(rows)):
        for u in range(len(rows)):
            if reach[u] & (1<<v):
                reach[u] |= reach[v]
    return all(r==(1<<len(rows))-1 for r in reach)

def canonical(rows,weights):
    # Generic weighted directed colour refinement followed by exhaustive
    # permutations of every unresolved invariant cell. No large-gap decoder
    # or endpoint/gap classification is assumed.
    palette = sorted(set(weights)); colour = [palette.index(w) for w in weights]
    while True:
        groups = [[v for v,c in enumerate(colour) if c==j] for j in range(max(colour)+1)]
        masks = [sum(1<<v for v in group) for group in groups]
        signatures = [(colour[v],tuple((rows[v]&mask).bit_count() for mask in masks))
                      for v in range(len(rows))]
        unique=sorted(set(signatures)); lookup={sig:j for j,sig in enumerate(unique)}
        refined=[lookup[sig] for sig in signatures]
        if refined==colour:
            break
        colour=refined
    count=1
    for group in groups:
        count*=factorial(len(group))
    require(count<=200000,'unexpected canonical search size')
    best=None
    for block_order in product(*(permutations(group) for group in groups)):
        order=tuple(v for group in block_order for v in group)
        candidate=(tuple(weights[v] for v in order),tuple(sum(((rows[v]>>w)&1)<<j
                    for j,w in enumerate(order)) for v in order))
        if best is None or candidate<best:
            best=candidate
    return best,count

def feedback_dp(rows,weights):
    # Add first vertex, unlike the producer's append-last recurrence.
    n=len(rows); costs=[10**9]*(1<<n); costs[0]=0
    for mask in range(1,1<<n):
        rest=mask
        while rest:
            bit=rest&-rest; v=bit.bit_length()-1; tail=mask^bit
            penalty=sum(weights[u]*weights[v] for u in range(n)
                        if tail&(1<<u) and rows[u]&(1<<v))
            costs[mask]=min(costs[mask],costs[tail]+penalty)
            rest^=bit
    return costs[-1]

def binomial(n,r):
    if r<0:
        return 0
    require(n>=r>=0,'invalid present binomial')
    return comb(n,r)

def formula(k,delta,strong,F,I,P):
    m=2*k-1 if strong else 2*k+1
    E=I[k] if strong else F[k]
    J=(2*k-3)*I[k-1] if strong else F[k]
    p=P[k] if strong else (k-1)*F[k]
    return (E*(binomial(m+delta-1,delta)-m*binomial(m+delta-k-1,delta-k))
       +(m*E-J)*binomial(m+delta-k-2,delta-k)
       +p*binomial(m+delta-k-3,delta-k-1))

def coefficient(total,parts,cap):
    if total<0:
        return 0
    polynomial=[1]+[0]*total
    for _ in range(parts):
        polynomial=[sum(polynomial[i-v] for v in range(min(cap,i)+1)) for i in range(total+1)]
    return polynomial[total]

def audit_zero(source,repo):
    F=[1]; I=[0]; P=[0,0]
    for k in range(1,33):
        F.append(F[-1]*(2*k-1))
        I.append(F[k]-sum(I[a]*F[k-a] for a in range(1,k)))
        if k>=2:
            P.append((k-1)*F[k]-sum(P[j]*sum(F[a]*F[k-j-a] for a in range(k-j+1)) for j in range(2,k)))
    graphs=[]
    for k in range(2,7):
        ms=list(pairings(tuple(range(2*k))))
        ps=list(stars(k))
        irreducible=[e for e in ms if cut_test(e,2*k)]
        jp=sum(sum(b-a==1 for a,b in e) for e in irreducible)
        pc=sum(cut_test(e,2*k-1) for e in ps)
        require(len(ms)==F[k] and len(irreducible)==I[k], 'matching graph recurrence')
        require(len(ps)==(k-1)*F[k] and pc==P[k], 'P3 graph recurrence')
        require(jp==(2*k-3)*I[k-1], 'adjacent-pair deletion/insertion')
        graphs.append(dict(k=k,matching=len(ms),irreducible=len(irreducible),p3=len(ps),irreducible_p3=pc,adjacent=jp))
    formula_rows=0
    for k in range(3,33):
        for delta in range(2*k-2):
            permitted=[(d,h) for d in range(2*k+1) for h in range(2*k+2)
                       if (k+1)*d+k*h<=delta]
            require(set(permitted)<={(0,0),(0,1),(1,0)},'corridor record exhaustion')
            for strong in (False,True):
                m=2*k-1 if strong else 2*k+1; E=I[k] if strong else F[k]
                J=(2*k-3)*I[k-1] if strong else F[k]; p=P[k] if strong else (k-1)*F[k]
                independent=(E*coefficient(delta,m,k-1)
                  +(m*E-J)*coefficient(delta-k,m-1,k-1)
                  +p*coefficient(delta-k-1,m-1,k-1))
                require(independent==formula(k,delta,strong,F,I,P),'corridor formula vs bounded polynomial coefficients')
                formula_rows+=1
    expected=json.loads((source/'zero-heavy-check.json').read_text())
    expected_rows={(r['k'],r['delta'],r['strong']):r for r in expected['boundary_layers']}
    rows=[]; checked=0; dp=0; relabels=0; max_permutations=0
    # Includes the ENTIRE stated corridor for k=3,4, plus scoped small exceptions.
    cases=[(k,d) for k in (3,4) for d in range(2*k-2)]+[(2,2),(2,3),(3,4)]
    for k,delta in cases:
        for strong in (False,True):
            seen=set(); counts=Counter(); reject=Counter(); excluded_nonstrong=0
            n=2*k; m=n-1 if strong else n+1
            def records():
                for edges in pairings(tuple(range(n))):
                    for deficit in words(delta,m,k-1):
                        yield 'full',edges,n,tuple(k-v for v in deficit)
                    for hole in range(m):
                        for deficit in words(delta-k,m-1,k-1):
                            gaps=list(k-v for v in deficit);gaps.insert(hole,0)
                            yield 'hole',edges,n,tuple(gaps)
                if delta>=k+1:
                    for edges in stars(k):
                        for deficit in words(delta-k-1,m-1,k-1):
                            yield 'p3',edges,n-1,tuple(k-v for v in deficit)
            for kind,edges,n_end,gaps in records():
                literal,weights=tournament(edges,n_end,gaps,strong)
                pairs=homogeneous_pairs(literal)
                forbidden=False
                if kind=='hole':
                    hole=gaps.index(0)+(1 if strong else 0)
                    forbidden=0<hole<n_end and (hole-1,hole) in edges
                require(bool(pairs)==forbidden,'literal clone filter')
                if forbidden:
                    if not strong or cut_test(edges,n_end):
                        reject[kind]+=1
                    continue
                if strong:
                    actual=strong_actual(literal)
                    require(actual==cut_test(edges,n_end),'SCC and proper-cut equivalence')
                    if not actual:
                        excluded_nonstrong+=1;continue
                require(sum(weights)==(2*k*k+(k if strong else 3*k))-delta,'literal mass')
                key,perms=canonical(literal,weights)
                max_permutations=max(max_permutations,perms)
                require(key not in seen,'distinct records have identical weighted tournament class')
                seen.add(key);counts[kind]+=1;checked+=1
                if sum(counts.values())<=2:
                    order=tuple(reversed(range(len(literal))))
                    changed=tuple(sum(((literal[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order)
                    require(canonical(changed,tuple(weights[v] for v in order))[0]==key,'generic canonical relabel invariance')
                    relabels+=1
                if k==2 and counts[kind]<=2:
                    require(feedback_dp(literal,weights)==k,'literal weighted feedback optimum')
                    dp+=1
            target=formula(k,delta,strong,F,I,P)
            require(len(seen)==target,'literal class count vs formula')
            ref=expected_rows.get((k,delta,strong))
            if ref:
                require(ref['classes']==len(seen),'producer receipt class count')
                require(ref['types']=={dict(full='A',hole='B',p3='C')[a]:b for a,b in counts.items()},'producer types')
                require(ref['clone_rejections']==({'B':reject['hole']} if reject['hole'] else {}),'producer clone reject count')
            rows.append(dict(k=k,delta=delta,strong=strong,classes=len(seen),types=dict(counts),clone_rejects=dict(reject),nonstrong_excluded=excluded_nonstrong))
    old=json.loads((repo/'05-knowledge/results/charge_palettes_boundary_20261020g/zero-heavy-mass.json').read_text())
    old2=next(row for row in old['complete_budget_compiler'] if row['k']==2)
    for delta in (2,3):
        for strong in (False,True):
            mass=10-delta if strong else 14-delta
            mapping=old2['strong_zero_heavy_by_mass' if strong else 'all_zero_heavy_by_mass']
            require(mapping[str(mass)]==formula(2,delta,strong,F,I,P),'independent complete k2 compiler source')
    return dict(graphs=graphs,formula_coefficient_checks=formula_rows,
      literal_class_checks=checked,canonical_max_permutations=max_permutations,
      weighted_feedback_dp_checks=dp,relabel_checks=relabels,layers=rows)

ZERO=(Q(0),)*4

def subtract(a,b):
    return tuple(x-y for x,y in zip(a,b))

def plus(a,b):
    return tuple(x+y for x,y in zip(a,b))

def distance_squared(a,b):
    # Direct real-coordinate quadratic form, independent of generic field multiplication.
    x,y,z,w=subtract(a,b)
    return (x*x+3*y*y+11*z*z+33*w*w, 2*x*w-2*y*z)

def rotate(p,t):
    u=Q(1-3*t*t,1+3*t*t);v=Q(2*t,1+3*t*t)
    a,b,c,d=p
    return (u*a-3*v*b,u*b+v*a,u*c-3*v*d,u*d+v*c)

def edges_for(points):
    return [(u,v) for u,v in combinations(range(len(points)),2)
            if distance_squared(points[u],points[v])==(1,0)]

def colouring_count(n,edges,boundary,root,q=6):
    fixed=dict(boundary);fixed[0]=root
    if any(a in fixed and b in fixed and fixed[a]==fixed[b] for a,b in edges):
        return 0,None
    earlier=[[] for _ in range(n)]
    adjacency=[[] for _ in range(n)]
    for a,b in edges:
        adjacency[a].append(b);adjacency[b].append(a)
    variables=[v for v in range(n) if v not in fixed]
    total=0;first=None
    def visit(j):
        nonlocal total,first
        if j==len(variables):
            total+=1
            if first is None:
                first=[fixed[v] for v in range(n)]
            return
        v=variables[j];forbidden={fixed[u] for u in adjacency[v] if u in fixed}
        for c in range(q):
            if c not in forbidden:
                fixed[v]=c;visit(j+1)
        fixed.pop(v,None)
    visit(0)
    return total,first

def proper(colours,edges,q):
    return all(isinstance(c,int) and 0<=c<q for c in colours) and all(colours[a]!=colours[b] for a,b in edges)

def audit_palette(source,repo):
    inherited=json.loads((repo/'05-knowledge/results/first_exit_complement_games_20261020/palette-glue.json').read_text())
    receipt=json.loads((source/'palette-check.json').read_text())
    # Literal coordinate reconstruction from pinned forcing.md seven-vertex table
    # and its nine reflection triples, rather than adopting receipt coordinates.
    points=[ZERO,(Q(1),Q(0),Q(0),Q(0)),(Q(1,2),Q(1,2),Q(0),Q(0)),
      (Q(3,2),Q(1,2),Q(0),Q(0)),(Q(5,6),Q(0),Q(1,6),Q(0)),
      (Q(5,12),Q(5,12),Q(1,12),Q(1,12)),(Q(15,12),Q(5,12),Q(3,12),Q(1,12))]
    triples=[(1,3,2),(2,4,0),(5,6,4),(1,4,0),(2,6,3),(3,5,6),(1,5,0),(2,3,1),(4,6,5)]
    for a,b,c in triples:
        points.append(subtract(plus(points[a],points[b]),points[c]))
    require(points==[tuple(map(Q,p)) for p in inherited['points'][:16]],'source coordinate table and reflections')
    edges=edges_for(points)
    require(len(set(points))==16 and len(edges)==31,'distinct source points and unit incidence')
    require(edges==[tuple(e) for e in receipt['base_edges']],'source edge witness')
    require(edges==[tuple(e) for e in inherited['edges'] if max(e)<16],'inherited source edge witness')
    boundary={v:(v-7)//3 for v in range(7,16)}
    repair=dict(boundary);repair[7]=1
    normal=[];repaired=[]
    for c in range(6):
        normal.append(colouring_count(16,edges,boundary,c)[0])
        repaired.append(colouring_count(16,edges,repair,c)[0])
    require(normal==[24,24,24,0,0,0]==receipt['source_counts'],'exact original root fibres')
    require(repaired==[60,84,84,8,8,8]==receipt['repair_counts'],'exact repaired root fibres')
    witness_text=(source/'palettes.md').read_text()
    matches=re.findall(r'\(3,0,4,5,4,5,3,1,0,0,1,1,1,2,2,2\)',witness_text)
    require(len(matches)==1,'locate exact note literal witness')
    note_literal=json.loads('['+matches[0][1:-1]+']')
    require(proper(note_literal,edges,6) and note_literal[0]==3,'parse note literal colouring')
    require(all(note_literal[v]==repair[v] for v in repair),'note literal input colours')
    for target,witness in enumerate(receipt['repair_witnesses']):
        require(len(witness)==16 and proper(witness,edges,6),'literal repaired root witness')
        require(witness[0]==target and all(witness[v]==repair[v] for v in repair),'literal root and port identities')
    placed=[tuple(map(Q,p)) for p in receipt['points']]
    require(len(placed)==len(set(placed))==91 and placed[:16]==points,'actual 91-point retained geometry')
    rotations=receipt['rotation_parameters']
    require(rotations==[None,2,3,4,5,6],'declared finite placement controls')
    expected=set(edges)
    physical_ports=[]
    for i,t in enumerate(rotations):
        embedding=[0]+list(range(15*i+1,15*i+16))
        physical_ports.extend(embedding[v] for v in range(7,16))
        if i:
            require(1<=t<=900*i+1,'finite avoidance bound')
            require([placed[v] for v in embedding]==[rotate(p,t) for p in points],'literal rotation reconstruction')
            # Every earlier integer must fail, and the chosen integer must avoid
            # ALL nonroot coincidences and unit contacts.
            old=placed[:15*i+1]
            for trial in range(1,t+1):
                copied=[rotate(p,trial) for p in points[1:]]
                bad=any(a==b or distance_squared(a,b)==(1,0) for a in old[1:] for b in copied)
                require(bad==(trial<t),'first surviving exact rotation')
        expected.update(tuple(sorted((embedding[a],embedding[b]))) for a,b in edges)
    actual=edges_for(placed)
    require(set(actual)==expected=={tuple(e) for e in receipt['actual_edges']},'full actual induced union')
    require(len(actual)==186 and len(physical_ports)==len(set(physical_ports))==54,'edge and distinct physical-port count')
    require(receipt['unordered_distances_checked']==comb(91,2)==4095,'all pair distances')
    require(len(receipt['four_colouring'])==91 and proper(receipt['four_colouring'],actual,4),'literal global four-colour hostile')
    direct_counts=[]
    for item in receipt['nearest_repairs']:
        n=item['copies'];colours=item['colouring'];palettes=item['palettes'];c=item['root']
        require(len(colours)==15*n+1==item['vertices'],'nearest literal vertex count')
        local_edges=[e for e in actual if e[1]<len(colours)]
        require(proper(colours,local_edges,6) and colours[0]==c,'nearest literal on ACTUAL induced graph')
        loads=[sum(x in a for a in palettes) for x in range(6)]
        require(loads==item['loads'],'literal palette loads')
        changes=0;nearest_boundary={}
        for i,palette in enumerate(palettes):
            embedding=[0]+list(range(15*i+1,15*i+16))
            changed=[]
            for v in range(7,16):
                physical=embedding[v];nearest_boundary[physical]=colours[physical]
                if colours[physical]!=palette[(v-7)//3]:
                    changed.append(v);changes+=1
            require(changed==([] if c in palette else [7]),'copy-specific physical-port one-change identity')
        require(changes==n-max(loads)==n//2==item['radius'],'physical Hamming radius and balanced optimum')
        require(item['target_extensions']==24**loads[c]*8**(n-loads[c]),'literal nearest repair product count')
        if n<=3:
            count,_=colouring_count(len(colours),local_edges,nearest_boundary,c)
            require(count==item['target_extensions'],'direct whole-graph enumeration of nearest fibre')
            direct_counts.append(dict(copies=n,root=c,extensions=count))
    # All palette multisets in the source-declared universe; independent
    # subfamily hitting-set optimization, not assumed repair/load identity.
    bank=list(combinations(range(6),3));multisets=[]
    for n in range(1,5):
        histogram=Counter();checked=0
        for palettes in combinations_with_replacement(bank,n):
            loads=[sum(c in p for p in palettes) for c in range(6)]
            best=n
            for size in range(n+1):
                if any(set.intersection(set(range(6)),*(set(palettes[i]) for i in range(n) if i not in chosen))
                       for chosen in combinations(range(n),size)):
                    best=size;break
            require(best==n-max(loads)<=n//2,'minimum repaired-copy subset')
            require(18*best<=9*n,'normalized radius bound')
            histogram[best]+=1;checked+=1
        literal=receipt['palette_multiset_controls'][n-1]
        require(literal['copies']==n and literal['palette_multisets']==checked,'declared palette multiset universe')
        require(literal['radius_histogram']=={str(k):v for k,v in histogram.items()},'palette histogram')
        multisets.append(dict(copies=n,count=checked,histogram=dict(histogram)))
    return dict(source_vertices=16,source_edges=31,root_fibres=normal,repair_fibres=repaired,
      literal_repair_witnesses=7,actual_vertices=91,actual_edges=186,physical_ports=54,
      exact_unordered_distances=4095,rotations=rotations,direct_global_fibres=direct_counts,
      nearest_literal_replays=6,palette_multisets=multisets)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parent.parent/'session18-zero-heavy')
    parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parent.parent/'math')
    args=parser.parse_args();source=args.source.resolve();repo=args.repo.resolve()
    require(digest(source/'manifest.json')==FROZEN,'FROZEN manifest identity')
    manifest=json.loads((source/'manifest.json').read_text())
    for name,sha in manifest['files'].items():
        require(digest(source/name)==sha,'frozen file '+name)
    for name,sha in manifest['dependencies'].items():
        require(digest(repo/name)==sha,'source dependency '+name)
    result=dict(status='PASS independent finite controls and frozen source pins; all-order proofs reviewed separately',
      frozen_manifest_sha256=FROZEN,source_files_authenticated=len(manifest['files']),
      dependency_files_authenticated=len(manifest['dependencies']),
      zero_heavy=audit_zero(source,repo),palette=audit_palette(source,repo))
    # Reauthenticate after the entire audit to detect concurrent changes.
    require(digest(source/'manifest.json')==FROZEN,'frozen manifest unchanged after audit')
    for name,sha in manifest['files'].items():
        require(digest(source/name)==sha,'frozen source unchanged after audit '+name)
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
