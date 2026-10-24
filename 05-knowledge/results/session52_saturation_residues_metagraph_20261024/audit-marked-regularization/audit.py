#!/usr/bin/env python3
"""Independent bit-row audit. No imports from the candidate; no assert statements.

The lift is constructed by XORing copy parity on every source pair. Exact graph
maps use domain propagation and neighbor-subtournament score invariants, rather
than the candidate's matrix/common-neighbor search.
"""
import argparse
import collections
import itertools
import json
import math
from pathlib import Path


def require(ok, why):
    if not ok:
        raise ValueError(why)


def vertices(bits):
    while bits:
        low = bits & -bits
        yield low.bit_length() - 1
        bits ^= low


def arc(g, u, v):
    return (g[u] >> v) & 1


def from_mask(n, mask):
    g = [0] * n
    for k, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        x, y = (u, v) if (mask >> k) & 1 else (v, u)
        g[x] |= 1 << y
    return tuple(g)


def from_arcs(n, arcs):
    g = [0] * n
    for u, v in arcs:
        g[u] |= 1 << v
    return tuple(g)


def dual(g):
    all_bits = (1 << len(g)) - 1
    return tuple(all_bits ^ (1 << i) ^ r for i, r in enumerate(g))


def lift(g):
    n = len(g)
    edges = []
    for u, v in itertools.combinations(range(2*n+1), 2):
        if v == 2*n:
            direction = u >= n
        elif u % n == v % n:
            direction = True
        else:
            direction = arc(g, u % n, v % n) ^ (u // n != v // n)
        edges.append((u, v) if direction else (v, u))
    return from_arcs(2*n+1, edges)


def valid(g):
    n = len(g)
    return all(not arc(g, i, i) and 0 <= g[i] < 1 << n for i in range(n)) and all(
        arc(g, i, j) != arc(g, j, i) for i, j in itertools.combinations(range(n), 2))


def induced(g, ids):
    return from_arcs(len(ids), ((i,j) for i in range(len(ids)) for j in range(len(ids)) if arc(g,ids[i],ids[j])))


def is_map(g, h, p):
    return len(p) == len(g) == len(h) and sorted(p) == list(range(len(g))) and all(
        arc(g, i, j) == arc(h, p[i], p[j]) for i in range(len(g)) for j in range(len(g)))


def invariant(g, v):
    ins = dual(g)[v]
    outs = g[v]
    return (outs.bit_count(), tuple(sorted((g[w] & outs).bit_count() for w in vertices(outs))),
            tuple(sorted((g[w] & ins).bit_count() for w in vertices(ins))))


def graph_maps(g, h, fixed=(), first_only=False):
    """All bijections, with dynamic arc-consistent bitset candidate domains."""
    if len(g) != len(h):
        return []
    n = len(g)
    gi = [invariant(g,i) for i in range(n)]
    hi = [invariant(h,i) for i in range(n)]
    if collections.Counter(gi) != collections.Counter(hi):
        return []
    domains = [sum(1 << j for j in range(n) if gi[i] == hi[j]) for i in range(n)]
    hd = dual(h)
    output = []
    def add(d, p, i, j):
        if not d[i] >> j & 1:
            return None
        pp = p.copy()
        dd = d.copy()
        pp[i] = j
        dd[i] = 1 << j
        for k in range(n):
            if pp[k] >= 0:
                if k != i and arc(g,i,k) != arc(h,j,pp[k]):
                    return None
            else:
                dd[k] &= (h[j] if arc(g,i,k) else hd[j])
                if not dd[k]:
                    return None
        return dd, pp
    p = [-1] * n
    for i,j in fixed:
        result = add(domains,p,i,j)
        if result is None:
            return []
        domains,p = result
    def walk(d,p):
        choices = [(d[i].bit_count(),i) for i in range(n) if p[i] < 0]
        if not choices:
            require(is_map(g,h,p), 'invalid complete graph map')
            output.append(tuple(p))
            return
        _,i = min(choices)
        for j in vertices(d[i]):
            nxt = add(d,p,i,j)
            if nxt is not None:
                walk(*nxt)
                if first_only and output:
                    return
    walk(domains,p)
    return output


def canon(g):
    return min(tuple(arc(g,i,j) for i in p for j in p) for p in itertools.permutations(range(len(g))))


def brute_maps(g):
    return {p for p in itertools.permutations(range(len(g))) if is_map(g,g,p)}


def reverse(g, pairs):
    h = list(g)
    for a,b in pairs:
        h[a] ^= 1 << b
        h[b] ^= 1 << a
    return tuple(h)


def edges(n):
    return set(itertools.combinations(range(n),2))


def pair(a,b):
    return tuple(sorted((a,b)))


def block(n,a,b):
    return {pair(x,y) for x in (a,n+a) for y in (b,n+b)}


def cycle_ok(g, cyc):
    return len(cyc) == len(set(cyc)) and all(arc(g,cyc[i-1],cyc[i]) for i in range(len(cyc)))


def paths(g):
    # Recursive path extension, not filtering all permutations.
    def walk(p,left):
        if not left:
            yield p
        else:
            for v in vertices(left if not p else left & g[p[-1]]):
                yield from walk(p+(v,),left ^ (1 << v))
    return set(walk((),(1 << len(g))-1))


def root_predicate(g):
    return [i for i in range(len(g)) if all((g[i] & g[j]) != 0 for j in range(len(g)) if j != i)]


def lifted(p):
    n = len(p)
    return (2*n,) + tuple(x for i in p[::-1] for x in (i,n+i))


def exhaustive():
    rows = []
    ordinary_counts = [1]
    self_converse_counts = [1]
    for n in range(1,6):
        hp_total = edge_total = special = root_moves = 0
        sa_hist = collections.Counter()
        da_hist = collections.Counter()
        source_reps = {}
        for mask in range(1 << (n*(n-1)//2)):
            g = from_mask(n,mask)
            d = lift(g)
            require(valid(g) and valid(d),'tournament validation')
            require([x.bit_count() for x in d] == [n] * (2*n+1),'regular degree')
            require(induced(d,list(range(n))) == g == induced(d,list(range(n,2*n))), 'same copies')
            nonextreme = all(0 < r.bit_count() < n-1 for r in g)
            special += nonextreme
            require(root_predicate(d) == ([2*n] if nonextreme else []),'canonical-root predicate')
            for i in range(n):
                require(d[i] & d[n+i] == 0,'matched zero-common-out')
                require((d[2*n] & d[i]).bit_count() == g[i].bit_count(),'z U common-out formula')
                require((d[2*n] & d[n+i]).bit_count() == n-1-g[i].bit_count(),'z V common-out formula')
            source_aut = brute_maps(g)
            require(set(graph_maps(g,g)) == source_aut,'map-search independently checked by permutation brute force')
            target_aut = set(graph_maps(d,d))
            expected = {p + tuple(n+i for i in p) + (2*n,) for p in source_aut}
            require({p for p in target_aut if p[2*n] == 2*n} == expected,'root-stabilizer equality')
            if nonextreme:
                require(target_aut == expected,'full group on nonextreme domain')
            root_moves += bool(target_aut - expected)
            sa_hist[len(source_aut)] += 1
            da_hist[len(target_aut)] += 1
            swap = tuple((i+n) % (2*n) for i in range(2*n)) + (2*n,)
            require(is_map(lift(dual(g)),dual(d),swap),'dual copy swap')
            supports = set()
            for a,b in edges(n):
                support = block(n,a,b)
                require(not support & supports,'non-disjoint source blocks')
                supports |= support
                x,y = (a,b) if arc(g,a,b) else (b,a)
                require(cycle_ok(d,(x,y,n+x,n+y)), 'four-cycle orientation')
                require(lift(reverse(g,{(a,b)})) == reverse(d,support),'single-edge lift equality')
                edge_total += 1
            hp = paths(g)
            hc = set()
            for q in itertools.permutations(range(n)):
                candidate = (2*n,) + tuple(x for i in q for x in (i,n+i))
                if cycle_ok(d,candidate):
                    hc.add(candidate)
            require({lifted(p) for p in hp} == hc,'path to distinguished cycle exactness')
            require(len(hp) % len(source_aut) == 0,'free orbit integral')
            for p in hp:
                require(sum(tuple(a[i] for i in p) == p for a in source_aut) == 1,'free path action')
                c = lifted(p)
                path_edges = {pair(p[i-1],p[i]) for i in range(1,n)}
                changed_source = edges(n)-path_edges
                changed_target = set().union(*(block(n,*e) for e in changed_source))
                lifted_partner = lift(reverse(g,changed_source))
                require(lifted_partner == reverse(d,changed_target),'path complement block transport')
                require(cycle_ok(lifted_partner,c),'cycle retained by block move')
                require(len(changed_target) == 2*(n-1)*(n-2),'block-support cardinality')
                cycle_edges = {pair(c[i-1],c[i]) for i in range(len(c))}
                full_partner = reverse(d,edges(2*n+1)-cycle_edges)
                difference = sum(arc(full_partner,i,j) != arc(lifted_partner,i,j) for i,j in edges(2*n+1))
                require(difference == 5*(n-1),'native full-complement difference ledger')
            hp_total += len(hp)
            source_reps.setdefault(canon(g),g)
        reps = list(source_reps.values())
        self_count = sum(canon(g) == canon(dual(g)) for g in reps)
        core_reps = [g for g in reps if all(0 < r.bit_count() < n-1 for r in g)]
        core_self_count = sum(canon(g) == canon(dual(g)) for g in core_reps)
        core_merged_count = len({min(canon(g),canon(dual(g))) for g in core_reps})
        ordinary_counts.append(len(reps))
        self_converse_counts.append(self_count)
        if n >= 2:
            require(len(core_reps) == ordinary_counts[n]-2*ordinary_counts[n-1]+ordinary_counts[n-2], 'no-extreme ordinary count formula')
            require(core_self_count == self_converse_counts[n]-self_converse_counts[n-2], 'no-extreme self-converse formula')
            require(2*core_merged_count == len(core_reps)+core_self_count, 'no-extreme merged formula')
        target_reps = []
        merged_targets = []
        merged_sources = set()
        rooted_pairs = 0
        for i,g in enumerate(reps):
            d = lift(g)
            merged_sources.add(min(canon(g),canon(dual(g))))
            if not any(graph_maps(d,q,first_only=True) for q in target_reps):
                target_reps.append(d)
            if not any(graph_maps(d,q,first_only=True) or graph_maps(d,dual(q),first_only=True) for q in merged_targets):
                merged_targets.append(d)
            for j,h in enumerate(reps):
                require(bool(graph_maps(d,lift(h),fixed=((2*n,2*n),),first_only=True)) == (i == j), 'rooted class-faithfulness')
                require(bool(graph_maps(d,dual(lift(h)),fixed=((2*n,2*n),),first_only=True)) == (canon(g) == canon(dual(h))), 'rooted converse faithfulness')
                rooted_pairs += 1
        require(hp_total == math.factorial(n) * 2**((n-1)*(n-2)//2),'global path incidence')
        rows.append({'n':n,'labelled_tournaments':2**(n*(n-1)//2),'labelled_nonextreme_tournaments':special,
                     'labelled_hamilton_paths':hp_total,'single_arc_block_checks':edge_total,
                     'source_automorphism_order_histogram':dict(sorted(sa_hist.items())),
                     'target_automorphism_order_histogram':dict(sorted(da_hist.items())),
                     'targets_with_nonfixed_original_root':root_moves,'source_ordinary_classes':len(reps),
                     'target_unmarked_classes':len(target_reps),'source_converse_merged_classes':len(merged_sources),
                     'target_unmarked_converse_merged_classes':len(merged_targets),
                     'target_classes_lost_by_forgetting_root':len(reps)-len(target_reps),
                     'independent_rooted_class_pair_tests':rooted_pairs,
                     'source_self_converse_ordinary_classes':self_count,
                     'source_nonextreme_ordinary_classes':len(core_reps),
                     'source_nonextreme_self_converse_classes':core_self_count,
                     'source_nonextreme_merged_classes':core_merged_count})
    return rows


def join(*blocks):
    g = []
    n = sum(map(len,blocks))
    offset = 0
    for b in blocks:
        later = ((1 << n)-1) ^ ((1 << (offset+len(b)))-1)
        g.extend((r << offset) | later for r in b)
        offset += len(b)
    return tuple(g)


def scc_profile(g):
    # Discover SCCs by forward/reverse reachability from each remaining vertex.
    def reach(h,v):
        seen = 1 << v
        while True:
            nxt = seen
            for u in vertices(seen):
                nxt |= h[u]
            if nxt == seen:
                return seen
            seen = nxt
    rem = (1 << len(g))-1
    comps=[]
    while rem:
        v = next(vertices(rem))
        comp = reach(g,v) & reach(dual(g),v)
        comps.append(comp)
        rem ^= comp
    rank = {c: sum(bool(g[next(vertices(c))] & d) for d in comps if c != d) for c in comps}
    comps.sort(key=lambda c: -rank[c])
    return [c.bit_count() for c in comps]


def matrix(g):
    return [''.join(str(arc(g,i,j)) for j in range(len(g))) for i in range(len(g))]


def collision(q):
    a,b = join((0,),q),join(q,(0,))
    n = len(a)
    p = (2*n,) + tuple(range(2*n))
    da,db = lift(a),lift(b)
    require(is_map(da,db,p),'recenter map fails an adjacency')
    return {'source_1_plus_Q':matrix(a),'source_Q_plus_1':matrix(b),
            'D_1_plus_Q':matrix(da),'D_Q_plus_1':matrix(db),
            'explicit_D_isomorphism_old_to_new':list(p),
            'source_scc_sizes':scc_profile(a),'other_source_scc_sizes':scc_profile(b),
            'source_isomorphic':bool(graph_maps(a,b,first_only=True)),
            'source_isomorphic_to_other_converse':bool(graph_maps(a,dual(b),first_only=True)),
            'source_aut_order':len(graph_maps(a,a)),'target_aut_order':len(graph_maps(da,da)),
            'source_hamilton_paths':len(paths(a))}


def controls():
    c3 = (2,4,1)
    ordinary = collision(c3)
    merged = collision(join((0,),c3))
    larger = collision(join(c3,(6,12,8,1)))
    require(not larger['source_isomorphic'] and not larger['source_isomorphic_to_other_converse'],'larger merged collision scope')
    require(larger['source_scc_sizes'] == [1,3,4] and larger['other_source_scc_sizes'] == [3,4,1],'larger SCC distinction')
    require(not ordinary['source_isomorphic'] and ordinary['source_isomorphic_to_other_converse'],'ordinary collision scope')
    require(not merged['source_isomorphic'] and not merged['source_isomorphic_to_other_converse'],'merged collision scope')
    require(merged['source_scc_sizes'] == [1,1,3] and merged['other_source_scc_sizes'] == [1,3,1],'SCC distinction')
    recenter_checks = 0
    for n in range(6):
        for mask in range(2**(n*(n-1)//2)):
            q = from_mask(n,mask)
            a,b=join((0,),q),join(q,(0,))
            require(is_map(lift(a),lift(b),(2*n+2,)+tuple(range(2*n+2))),'universal small recenter')
            recenter_checks += 1
    t = (2,0)
    d = lift(t)
    rooted_cycles = {(4,)+p for p in itertools.permutations(range(4)) if cycle_ok(d,(4,)+p)}
    require(len(rooted_cycles) == 2 and len(paths(t)) == 1,'all-vs-distinguished cycles')
    require(len(graph_maps(d,d)) == 5 and len(graph_maps(t,t)) == 1,'full-Aut false extension')
    require(root_predicate(d) == [],'unrestricted root predicate false')
    bad = reverse(d,{(2,3)})
    require(any(r.bit_count() != 2 for r in bad),'converse-copy false extension')
    # Important non-necessity control: source extremes do not force a larger group.
    require(merged['source_aut_order'] == merged['target_aut_order'] == 3,'sufficient-not-necessary control')
    return {'ordinary_collision_order_4':ordinary,'merged_collision_order_5':merged,
            'merged_collision_order_8':larger,
            'all_small_recenter_checks_Q_orders_0_through_5':recenter_checks,
            'all_vs_distinguished_cycle_control':{'D_transitive_2_all_cycles_rooted_at_z':len(rooted_cycles),'distinguished_cycles':1},
            'false_extensions_rejected':['replace second copy by converse','drop root for unrestricted full automorphism equality',
                'drop no-source/no-sink condition from the canonical-root predicate','identify all Hamilton cycles with matched alternating cycles',
                'identify native full-target noncycle complement with source complement transport',
                'infer unmarked ordinary injection','infer unmarked converse-merged injection'],
            'no_source_no_sink_is_sufficient_not_necessary':True}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    rows=exhaustive()
    result={'status':'PASS','method':'independent bit rows; copy-parity XOR lift; domain-propagation isomorphisms; DFS paths',
            'no_extreme_class_count_and_regular_lower_bound_controls':[{'source_order':r['n'],'target_regular_order':2*r['n']+1, 'ordinary_lower_bound':r['source_nonextreme_ordinary_classes'],'self_converse_lower_bound':r['source_nonextreme_self_converse_classes'],'merged_lower_bound':r['source_nonextreme_merged_classes']} for r in rows if r['n'] >= 2],
            'exhaustive_rows':rows,'total_labelled_tournaments':sum(r['labelled_tournaments'] for r in rows),**controls()}
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    Path(args.output).write_text(payload)
    print(payload,end='')

if __name__ == '__main__':
    main()
