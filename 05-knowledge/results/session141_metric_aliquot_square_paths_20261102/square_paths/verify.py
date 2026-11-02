#!/usr/bin/env python3
"""Standalone exact square-guard census and port-normal-form controls.

No repository, nonstandard dependency, network, or tournament-space enumeration.
All checks remain active under python -O. Writes JSON only to stdout.
"""
from collections import Counter
from itertools import combinations, permutations
from math import comb, isqrt
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def arcset(path):
    return frozenset(zip(path, path[1:]))


def square_paths(n):
    """All-start depth-first search, without degree/leaf pruning."""
    adjacency = {u: tuple(v for v in range(1, n+1)
                         if u != v and isqrt(u+v)**2 == u+v)
                 for u in range(1, n+1)}
    found = []
    def visit(path, used):
        if len(path) == n:
            found.append(path)
            return
        for v in adjacency[path[-1]]:
            if not used & (1 << (v-1)):
                visit(path+(v,), used | (1 << (v-1)))
    for start in range(1, n+1):
        visit((start,), 1 << (start-1))
    found.sort()
    require(len(found) == len(set(found)), 'duplicate path')
    require(set(found) == {p[::-1] for p in found}, 'reversal closure')
    edges = tuple((u,v) for u in adjacency for v in adjacency[u] if u < v)
    return edges, found


def compatible(arcs):
    return all((v,u) not in arcs for u,v in arcs)


def common_blocks(paths):
    """Canonical maximal common-directed-arc components, ordered by first path."""
    require(bool(paths), 'nonempty family')
    shared = set.intersection(*(set(arcset(p)) for p in paths))
    blocks = [[paths[0][0]]]
    for u,v in zip(paths[0], paths[0][1:]):
        if (u,v) in shared:
            blocks[-1].append(v)
        else:
            blocks.append([v])
    which = {v:i for i,b in enumerate(blocks) for v in b}
    orders = []
    for p in paths:
        order = []
        for v in p:
            if not order or order[-1] != which[v]:
                order.append(which[v])
        require(sorted(order) == list(range(len(blocks))), 'blocks contiguous')
        require(tuple(v for i in order for v in blocks[i]) == p, 'block reconstruction')
        orders.append(order)
    return blocks, orders


def singleton_test(blocks, orders):
    joins = {(i,j) for p in orders for i,j in zip(p,p[1:])}
    return not any(len(blocks[i]) == len(blocks[j]) == 1 and (j,i) in joins
                   for i,j in joins)


def clique_moments(paths, edge_count):
    arcs = [arcset(p) for p in paths]
    neighbours = [{j for j in range(i+1,len(paths)) if compatible(arcs[i]|arcs[j])}
                  for i in range(len(paths))]
    moments = Counter({0: 1 << edge_count})
    union_profiles = {}
    cliques = []
    def visit(indices, candidates, union):
        for i in sorted(candidates):
            selected = indices+(i,)
            joint = union|arcs[i]
            k = len(selected)
            moments[k] += 1 << (edge_count-len(joint))
            union_profiles.setdefault(k, Counter())[len(joint)] += 1
            cliques.append(selected)
            blocks, orders = common_blocks([paths[j] for j in selected])
            require(singleton_test(blocks,orders), 'singleton conflict criterion')
            joins = {(u,v) for p in orders for u,v in zip(p,p[1:])}
            require(len(joint) == len(paths[0])-len(blocks)+len(joins), 'port union formula')
            if k == 2:
                d = len(joint)-(len(paths[0])-1)
                require(len(blocks) == d+1, 'pair common block normal form')
            visit(selected, candidates & neighbours[i], joint)
    visit((), set(range(len(paths))), frozenset())
    top = max(moments)
    histogram = {k:sum((-1)**(j-k)*comb(j,k)*moments[j] for j in range(k,top+1))
                 for k in range(top+1)}
    require(min(histogram.values()) >= 0, 'moment positivity')
    return histogram, dict(sorted(moments.items())), union_profiles, cliques


def active_cubes(edges, paths):
    """Distinct algorithm: enumerate ONLY each path's free-edge cube."""
    index = {e:i for i,e in enumerate(edges)}
    source_count = Counter()
    for path in paths:
        used = positive = 0
        for u,v in zip(path,path[1:]):
            bit = 1 << index[min(u,v),max(u,v)]
            used |= bit
            if u < v:
                positive |= bit
        free_bits = [1 << i for i in range(len(edges)) if not (used >> i) & 1]
        sources = [positive]
        for bit in free_bits:
            sources += [x|bit for x in sources]
        require(len(sources) == (1 << (len(edges)-len(path)+1)), 'cube cardinality')
        source_count.update(sources)
    result = dict(Counter(source_count.values()))
    result[0] = (1 << len(edges))-len(source_count)
    return dict(sorted(result.items())),len(source_count)


def controls():
    # Every pair of permutations through n=5; includes many conflicting families.
    tests = 0
    for n in range(2,6):
        ps = list(permutations(range(n)))
        for p,q in combinations(ps,2):
            blocks,orders = common_blocks([p,q])
            require(singleton_test(blocks,orders) == compatible(arcset(p)|arcset(q)),
                    'normal-form equivalence, including negative cases')
            tests += 1
    # Exact block-shuffle obstruction: singleton A,B demand opposite same arc.
    require(not compatible(arcset((0,1,2))|arcset((0,2,1))), 'singleton hostile')
    require(compatible(arcset((0,1,2,3))|arcset((0,3,1,2))), 'one long block suffices')
    # A two-block swap can have d=1, but it necessarily closes a Hamilton cycle.
    p,q=(0,1,2,3),(2,3,0,1)
    require(compatible(arcset(p)|arcset(q)), 'd1 positive cycle control')
    require(len(arcset(p)|arcset(q))==4,'d1 union')
    # A fixed directed source makes subdivision/block blow-up HP-faithful.
    D={(0,1),(0,3),(1,2),(1,4),(2,3),(3,1),(3,4),(4,2)}
    port_orders=[p for p in permutations(range(5)) if arcset(p) <= D]
    require(port_orders==[(0,1,2,3,4),(0,1,4,2,3),(0,3,1,4,2)],'three-path port bank')
    blowups=0
    for lengths in ((1,1,1,1,1),(2,2,2,2,2),(3,1,2,1,2)):
        blocks=[];v=0
        for length in lengths:
            blocks.append(tuple(range(v,v+length)));v+=length
        arcs=set().union(*(arcset(b) for b in blocks))
        arcs |= {(blocks[i][-1],blocks[j][0]) for i,j in D}
        adj={u:[] for u in range(v)}
        for u,w in arcs:adj[u].append(w)
        found=[]
        def visit(p):
            if len(p)==v:found.append(p);return
            for w in adj[p[-1]]:
                if w not in p:visit(p+(w,))
        for u in range(v):visit((u,))
        expected={tuple(x for i in p for x in blocks[i]) for p in port_orders}
        require(set(found)==expected,'source-preserving blow-up correspondence')
        blowups+=1
    # Without a source, expanding a directed cycle adds Hamiltonian cuts.
    cycle3={(0,1),(1,2),(2,0)}
    cycle4={(0,1),(1,2),(2,3),(3,0)}
    require(sum(arcset(p)<=cycle3 for p in permutations(range(3)))==3, 'source-free base')
    require(sum(arcset(p)<=cycle4 for p in permutations(range(4)))==4, 'source-free hostile')
    return {'pair_normal_form_cases':tests,'source_blowup_controls':blowups,
            'source_free_cycle_hostile':'passed',
            'negative_singleton_control':'passed','positive_cycle_control':'passed'}


def main():
    expected_paths={1:1,15:2,16:2,17:2,23:6,25:20,26:24,27:70,28:104,29:38,30:40}
    expected_hist={
      25:{0:4294962432,1:4608,2:256},
      26:{0:17179857952,1:10240,2:928,3:64},
      27:{0:68719413784,1:54976,2:7360,3:488,4:120,5:8},
      28:{0:274877723488,1:158096,2:21840,3:2928,4:528,5:64},
      29:{0:1099511486464,1:129024,2:10240,3:2048},
      30:{0:4398046212096,1:274432,2:20480,3:4096}}
    rows=[];witnesses={}
    for n in range(1,31):
        edges,paths=square_paths(n)
        require(len(paths)==expected_paths.get(n,0),'square path census')
        hist,moments,profiles,cliques=clique_moments(paths,len(edges))
        literal,active=active_cubes(edges,paths)
        require(hist==literal,'clique versus active cube algorithms')
        require(sum(hist.values())==1<<len(edges),'guard total')
        require(sum(k*v for k,v in hist.items())==len(paths)*(1<<(len(edges)-n+1))
                if paths else hist=={0:1<<len(edges)},'first moment')
        if n in expected_hist:require(hist==expected_hist[n],'new spectrum controls')
        if n>1 and paths:
            require(len(paths)%2==0,'directed/reversal count')
        rows.append({'n':n,'guard_edges':len(edges),'total_tournament_bits':comb(n,2),
                     'outside_free_bits':comb(n,2)-len(edges),'directed_paths':len(paths),
                     'guard_histogram':hist,'guard_binomial_moments':moments,
                     'compatible_subset_union_profiles':profiles,'active_orientations':active})
        if n in (25,26):
            nontrivial=[]
            for selected in cliques:
                if len(selected)<2:continue
                blocks,orders=common_blocks([paths[j] for j in selected])
                nontrivial.append({'path_indices':selected,'blocks':blocks,'block_orders':orders,
                                   'union_arcs':len(set.union(*(set(arcset(paths[j])) for j in selected)))})
            witnesses[n]={'path_bank':paths,'compatible_subsets':nontrivial}
    return {'status':'passed','scope':'square guards on 1..n, 1 <= n <= 30; labelled tournaments',
            'methods':['all-start path DFS','compatible-subset union moments','active-cube orientation counting'],
            'controls':controls(),'rows':rows,'witnesses':witnesses}

if __name__=='__main__':
    print(json.dumps(main(),sort_keys=True,indent=2))
