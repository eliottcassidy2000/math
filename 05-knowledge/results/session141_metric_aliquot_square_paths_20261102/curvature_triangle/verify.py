#!/usr/bin/env python3
"""Portable finite checks; proofs of all-level claims are in PROOF.md."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import acosh, asinh, cosh, log
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def source_checks():
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
    for p in pins['git_sources']:
        data = (HERE / p['snapshot']).read_bytes()
        need(len(data) == p['bytes'], 'source byte count: ' + p['snapshot'])
        need(sha256(data).hexdigest() == p['sha256'], 'source digest: ' + p['snapshot'])
    return len(pins['git_sources'])


def C(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def predecessors(n):
    out = [2 * n]
    if n % 6 == 4:
        out.append((n - 1) // 3)
    return out


def native_checks():
    depth = {5: 0}
    front = {5}
    counts, layers = [1], [1]
    for r in range(1, 27):
        new = {u for n in front for u in predecessors(n)}
        need(not (new & depth.keys()), 'repeated root5 depth')
        for u in new:
            need(C(u) in front, 'inverse edge')
            x = u
            for t in range(r):
                need(x != 5, 'first-hit prematurely reached')
                x = C(x)
            need(x == 5, 'whole first-hit path')
            depth[u] = r
        front = new
        counts.append(len(depth))
        layers.append(len(new))
    rows = ['integer\tfirst_hit_steps'] + [f'{n}\t{d}' for n, d in sorted(depth.items())]
    literal = ('\n'.join(rows) + '\n').encode()
    need(literal == (HERE / 'sources/inverse_root5_depth26.tsv').read_bytes(), 'pinned state/clock TSV')
    need(counts[21] == 704 and counts[26] == 2263, 'finite Collatz counts')
    need(all(counts[r] <= 1 + 3*r*(r+1)//2 for r in range(21)), 'earlier honeycomb capacity')
    need(all(counts[r] <= 1 + 3*r*(r+1) for r in range(26)), 'earlier triangular capacity')
    need(704 - 694 == 10 and 2263 - 2107 == 156, 'finite excess arithmetic')
    edges = {(n, C(n)) for n in depth if n != 5}
    need(len(edges) == len(depth)-1 and all(v in depth for _, v in edges), 'stopped tree edge count')
    return {'counts': counts, 'levels': layers, 'tsv_sha256': sha256(literal).hexdigest(),
            'honeycomb': {'radius': 21, 'sources': 704, 'capacity': 694, 'excess': 10},
            'triangular': {'radius': 26, 'sources': 2263, 'capacity': 2107, 'excess': 156}}


def honey_neighbors(v):
    a, b, s = v
    return [(a,b,1),(a-1,b,1),(a,b-1,1)] if not s else [(a,b,0),(a+1,b,0),(a,b+1,0)]


def tri_neighbors(v):
    a, b = v
    return [(a+x,b+y) for x,y in [(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]]


def lattice_checks():
    summaries = []
    for name, root, neighbors, factor in [('honeycomb',(0,0,0),honey_neighbors,Q(3,2)),
                                          ('triangular',(0,0),tri_neighbors,Q(3))]:
        seen = {root}
        front = {root}
        for r in range(106):
            if r:
                front = {w for v in front for w in neighbors(v)} - seen
                seen.update(front)
            need(len(seen) == 1+factor*r*(r+1), name + ' exact BFS ball')
        summaries.append({'graph':name,'radii_checked':106,'last_radius':105,'last_count':len(seen)})
    # Independent coordinate region descriptions, both honeycomb parity cases.
    for k in range(15):
        tri = {(a,b) for a in range(-k,k+1) for b in range(-k,k+1) if abs(a+b)<=k}
        white = tri | {(a-1,b) for a,b in tri} | {(a,b-1) for a,b in tri}
        region = {(a,b) for a in range(-k-1,k+1) for b in range(-k-1,k+1) if -k-1<=a+b<=k}
        need(white == region and len(white)==3*(k+1)**2, 'white union/region')
        need(len(tri)==1+3*k*(k+1), 'triangular coordinate count')
    return summaries


def odd_C(n):
    m = 3*n+1
    a = 0
    while m % 2 == 0:
        m //= 2
        a += 1
    return m, a


def sf5(n):
    exponents = [2,4,6] if n%3==1 else [3,5,7]
    out = [((2**a*n-1)//3,a) for a in exponents if (2**a*n-1)%9]
    need(len(out)==2 and all(u>n and u%2 and u%3 and odd_C(u)==(n,a) for u,a in out), 'SF5 legality')
    return out


def union_checks():
    nodes = {(): 5}
    exponents = {}
    front = [()]
    for _ in range(8):
        new = []
        for w in front:
            for bit,(n,a) in enumerate(sf5(nodes[w])):
                v = w+(bit,)
                nodes[v] = n
                exponents[v] = a
                new.append(v)
        front = new
    need(len(set(nodes.values()))==len(nodes), 'distinct odd nodes')
    paths = {}
    for w,n in nodes.items():
        p = [n]
        limit = sum(exponents[w[:i]]+1 for i in range(1,len(w)+1))
        while p[-1]!=5:
            need(len(p)<=limit, 'stored clock path bound')
            p.append(C(p[-1]))
        need(len(p)-1==limit, 'stored first-hit clock')
        paths[w] = p
    table = []
    for h in range(9):
        selected = [w for w in nodes if len(w)<=h]
        source_union = {n for w in selected for n in paths[w]}
        states = {}
        edges = set()
        for w in selected:
            p = paths[w]
            for i,n in enumerate(p):
                d = len(p)-1-i
                need(n not in states or states[n]==d, 'shared union first-hit clocks')
                states[n] = d
            edges.update(zip(p,p[1:]))
        parents = [w for w in selected if len(w)<h]
        formula = len(selected)+sum(max(a for _,a in sf5(nodes[w])) for w in parents)
        naive = 1+sum(exponents[w]+1 for w in selected if w)
        excess = sum(min(a for _,a in sf5(nodes[w])) for w in parents)
        need(len(source_union)==formula==naive-excess, 'native union formula')
        need(len(edges)==len(source_union)-1, 'native union is stopped tree')
        need({n for w in selected if len(w)==h for n in paths[w]}==source_union, 'leaf paths span same union')
        table.append({'odd_depth':h,'odd_vertices':len(selected),'native_vertices':len(source_union),
                      'naive_subdivision_vertices':naive,'duplicate_even_states':excess,
                      'native_edges':len(edges),'max_first_hit_clock':max(states.values())})
    need((table[-1]['native_vertices'],table[-1]['naive_subdivision_vertices'],table[-1]['duplicate_even_states'])==(2001,2811,810), 'depth8 union witness')
    p = paths[(0,)]
    q = paths[(1,)]
    meet = next(n for n in p if n in q)
    need(nodes[(0,)]==13 and nodes[(1,)]==53 and meet==40 and (p.index(meet),q.index(meet))==(1,3), 'shared-tail first meeting')
    return {'rows':table,'sibling_meeting':{'sources':[13,53],'meeting':40,'first_hit_clocks':[1,3],'native_distance':4,'independent_subdivision_distance':10}}


def five_step_checks():
    k_by_residue = {1:2,2:1,4:2,5:3,7:4,8:1}
    def children(n):
        k = k_by_residue[n%9]
        need(2**k*n%9 in (4,7), 'five-step residue selector')
        nums = [32*n,(16*n-2**(4-k))//3]
        need((16*n-2**(4-k))%3==0 and len(set(nums))==2, 'five-step integral distinct')
        for u in nums:
            need(u>n and u%3, 'five-step positive monotone unit')
            v = u
            for _ in range(5):
                v = C(v)
            need(v==n, 'five-step forward replay')
        return nums
    for n in range(1,2001):
        if n%3:
            children(n)
    front, all_nodes = {5}, {5}
    rows = []
    for h in range(1,11):
        new = {u for n in front for u in children(n)}
        need(len(new)==2**h and not new & all_nodes, 'five-step tree collision')
        for u in new:
            v = u
            for _ in range(5*h):
                need(v!=5, 'five-step root5 first hit')
                v = C(v)
            need(v==5, 'five-step exact clock to root5')
        all_nodes.update(new)
        front = new
        rows.append({'macro_depth':h,'endpoints':len(new),'cumulative':len(all_nodes),'native_first_hit_steps':5*h})
    return {'selector_parent_checks':1334,'tree_rows':rows}


def placement(q,W,depth):
    pos = {(): (W/2,Q(1))}
    front = [()]
    for _ in range(depth):
        new = []
        for w in front:
            x,y = pos[w]
            for bit,sign in [(0,-1),(1,1)]:
                v = w+(bit,)
                pos[v] = (x+sign*W*(1-q)*y/2,q*y)
                new.append(v)
        front = new
    return pos


def cosh_distance(p,q):
    x,y = p
    X,Y = q
    return ((x-X)**2+y*y+Y*Y)/(2*y*Y)


def common_prefix(u,v):
    h=0
    while h<min(len(u),len(v)) and u[h]==v[h]:
        h+=1
    return h


def geometry_checks(q,W,depth,gapped):
    pos = placement(q,W,depth)
    beta = W*(1-q)/2
    target = 1+(1+W*W/4)*(1-q)**2/(2*q)
    edge_count = 0
    intervals = {}
    for w,(x,y) in pos.items():
        if len(w)<depth:
            intervals.setdefault(len(w),[]).append((x-beta*y,x+beta*y))
        if not w:
            continue
        parent = pos[w[:-1]]
        X,Y = parent
        sign = -1 if w[-1]==0 else 1
        center = X+sign*(beta*beta+q*q-1)*Y/(2*beta)
        need(cosh_distance(parent,(x,y))==target, 'equal edge cosh')
        need((X-center)**2+Y*Y==(x-center)**2+y*y, 'edge circle radii')
        need((X-center)*(x-center)>0, 'edge monotone semicircle side')
        need(sign*(X-center)>0 and y==q*Y, 'height-band monotonicity')
        edge_count += 1
    for level,spans in intervals.items():
        spans.sort()
        need(all(left[1]<right[0] for left,right in zip(spans,spans[1:])), 'different-parent horizontal edge intervals')
    need(len(set(pos.values()))==len(pos), 'injectivity')
    pair_count, equality_count = 0,0
    same_level_floor = 1+(W*(1-q)/q)**2/2
    g = W*(1-2*q)
    for u,v in combinations(pos,2):
        h=common_prefix(u,v)
        m=len(u)+len(v)-2*h
        x,y=pos[u]
        X,Y=pos[v]
        cd=cosh_distance(pos[u],pos[v])
        pair_count+=1
        need(cd<=(W*W+2)*q**(-m)/2,'QI upper rational precursor')
        need(cd>=(q**(len(u)-len(v))+q**(len(v)-len(u)))/2,'vertical lower precursor')
        if len(u)==len(v):
            need(abs(x-X)>=W*(1-q)*q**(len(u)-1), 'same-level center separation')
            need(cd>=same_level_floor,'same-level metric separation')
        if gapped and h<min(len(u),len(v)):
            need(abs(x-X)>=g*q**h,'relative sibling gap')
            need(2*cd-2>=g*g*q**(-m),'QI lower rational precursor')
        if q==Q(1,2) and W==1:
            is_edge=(m==1)
            need(cd>=Q(21,16),'dyadic minimum')
            need((cd==Q(21,16))==is_edge,'dyadic minimum exactly edges')
            equality_count+=(cd==Q(21,16))
    row={'q':str(q),'W':str(W),'depth':depth,'vertices':len(pos),'edges':edge_count,
         'pairs':pair_count,'edge_cosh':str(target),'minimum_equality_pairs':equality_count,
         'gapped':gapped}
    return row


def radial_checks():
    pos=placement(Q(1,2),Q(1),9)
    row=[]
    for k in range(10):
        vals=[cosh_distance(pos[()],p) for w,p in pos.items() if len(w)==k]
        exact=(5*Q(2)**k-2+5*Q(2)**(-k))/8
        need(max(vals)==exact,'exact radial maximum')
        need((Q(5,4)*Q(2)**k+Q(4,5)*Q(2)**(-k))/2-exact==Q(1,4)-Q(9,40)*Q(2)**(-k),'radial upper difference')
        if k:
            u=(0,)+(1,)*(k-1)
            v=(1,)+(0,)*(k-1)
            need(common_prefix(u,v)==0 and cosh_distance(pos[u],pos[v])==Q(3,2),'cross-root non-QI hostile')
        row.append({'depth':k,'vertices_through_depth':2**(k+1)-1,'outer_radius_cosh':str(exact),'cross_root_tree_distance':2*k if k else None})
    return row


def scalar_controls():
    # Rational specialization of Euler's incidence identities, not a classification claim.
    V,E,F=12,30,20
    need(V-E+F==2 and 3*F==2*E and 5*V==2*E and (6-5)*V==12,'degree5 sphere curvature')
    # Numerical values are illustrative only. Rational certificates carry finite geometry checks.
    examples=[]
    for q,W in [(Q(1,3),Q(1)),(Q(49,100),Q(1,10)),(Q(499,1000),Q(1,100))]:
        qf,Wf=float(q),float(W)
        ell=acosh(1+(1+Wf*Wf/4)*(1-qf)**2/(2*qf))
        s=log(1/qf)
        delta=min(s,2*asinh(Wf*(1-qf)/(2*qf)))
        need(ell>log(2) and delta>0,'illustrative threshold/separation')
        examples.append({'q':str(q),'W':str(W),'edge_length_float':ell,'separation_lower_float':delta,
                         'qi_additive_lower_float':-2*log(Wf*(1-2*qf)),
                         'siblings_closer_than_edges':2*asinh(Wf*(1-qf)/(2*qf))<ell})
    return {'sphere_degree5':{'vertices':V,'edges':E,'faces':F},'floating_examples_not_proofs':examples}


def main():
    return {'status':'PASS','source_pins_checked':source_checks(),
            'native_collatz':native_checks(),'lattice_bfs':lattice_checks(),
            'native_source_unions':union_checks(),'five_step_supply':five_step_checks(),
            'geometry':[geometry_checks(Q(1,2),Q(1),7,False),
                        geometry_checks(Q(1,3),Q(1),6,True),
                        geometry_checks(Q(2,5),Q(1,2),6,True),
                        geometry_checks(Q(49,100),Q(1,10),6,True)],
            'dyadic_radial':radial_checks(),'controls':scalar_controls(),
            'limits':['Finite exact checks corroborate but do not replace all-level proofs.',
                      'The sharp threshold is an infimum; equality is not settled.',
                      'No native Collatz hyperbolic embedding or convergence theorem is claimed.']}


if __name__=='__main__':
    print(json.dumps(main(),sort_keys=True,indent=2))
