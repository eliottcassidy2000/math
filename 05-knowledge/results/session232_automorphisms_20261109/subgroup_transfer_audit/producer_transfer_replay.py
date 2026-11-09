#!/usr/bin/env python3
"""Exact finite groupoid transfer for directed path prefixes. No repo imports."""
import itertools as it
import json
from collections import Counter, defaultdict
from math import gcd
from pathlib import Path


def need(x, msg):
    if not x:
        raise RuntimeError(msg)


def tournament(n, mask):
    a = [0] * n
    for bit, (u, v) in enumerate(it.combinations(range(n), 2)):
        if mask >> bit & 1:
            a[u] |= 1 << v
        else:
            a[v] |= 1 << u
    return a


def autos(a):
    n = len(a)
    degree = [x.bit_count() for x in a]
    groups = [tuple(i for i in range(n) if degree[i] == d) for d in sorted(set(degree))]
    answer = []
    for blocks in it.product(*(it.permutations(c) for c in groups)):
        p = list(range(n))
        for c, images in zip(groups, blocks):
            for x, y in zip(c, images):
                p[x] = y
        if all(((a[u] >> v) & 1) == ((a[p[u]] >> p[v]) & 1)
               for u in range(n) for v in range(u + 1, n)):
            answer.append(tuple(p))
    return answer


def image_mask(s, p):
    return sum(1 << p[v] for v in range(len(p)) if s >> v & 1)


def analyze(a, check_literal=True, chosen_group=None):
    n = len(a)
    group = autos(a) if chosen_group is None else chosen_group
    full = (1 << n) - 1
    canon = {}
    reps = {}
    for s in range(1, full + 1):
        for u in range(n):
            if s >> u & 1:
                key = min((image_mask(s, p), p[u]) for p in group)
                canon[s, u] = key
                reps[key] = key
    kernel = {key: sum(all(p[v] == v for v in range(n) if key[0] >> v & 1)
                       for p in group) for key in reps}
    dp = defaultdict(int)
    naive = defaultdict(int)
    transfers = []
    stabilization_checked = 0
    for b in sorted(reps, key=lambda x: (x[0].bit_count(), x)):
        s, u = b
        if s.bit_count() == 1:
            dp[b] += 1
            naive[b] += 1
        if s.bit_count() >= n - 2:
            need(kernel[b] == 1, ('terminal stabilization', a, b, kernel[b]))
            stabilization_checked += 1
        counts = Counter(canon[s | (1 << v), v] for v in range(n)
                         if not s >> v & 1 and a[u] >> v & 1)
        for target, multiplicity in counts.items():
            numerator = multiplicity * kernel[target]
            need(numerator % kernel[b] == 0, ('nonintegral transfer', b, target))
            w = numerator // kernel[b]
            dp[target] += w * dp[b]
            naive[target] += multiplicity * naive[b]
            transfers.append((b, target, multiplicity, kernel[b], kernel[target], w))
    orbit_total = sum(c for (s, u), c in dp.items() if s == full)
    naive_total = sum(c for (s, u), c in naive.items() if s == full)
    direct = defaultdict(set)
    h = 0
    if check_literal:
        for k in range(1, n + 1):
            for path in it.permutations(range(n), k):
                if all(a[path[i]] >> path[i + 1] & 1 for i in range(k - 1)):
                    s = sum(1 << v for v in path)
                    key = canon[s, path[-1]]
                    direct[key].add(min(tuple(p[v] for v in path) for p in group))
                    if k == n:
                        h += 1
        need(all(dp[key] == len(direct[key]) for key in reps), ('prefix orbit mismatch', a))
        need(h == len(group) * orbit_total, ('H mismatch', a))
    return dict(n=n, adj=a, aut=len(group), H=h, path_orbits=orbit_total,
                naive_transfer=naive_total, endpoint_state_orbits=len(reps),
                corrected_transfers=sum(w != m for b, t, m, kb, ka, w in transfers),
                stabilization_states=stabilization_checked), transfers


def ordered_sum(a, b):
    k = len(a)
    right = ((1 << len(b)) - 1) << k
    return [row | right for row in a] + [row << k for row in b]


def seven_states(a):
    """Literal ordered two-path covers; empty paths allowed."""
    n = len(a)
    h = q = 0
    aw = [0]*n; gw = [0]*n; sw = [0]*n; ew = [0]*n
    cross = [[0]*n for _ in range(n)]
    for p in it.permutations(range(n)):
        good = [bool(a[p[i]] >> p[i+1] & 1) for i in range(n-1)]
        if all(good):
            h += 1; sw[p[0]] += 1; ew[p[-1]] += 1
        for k in range(n+1):
            if all(good[:max(k-1, 0)]) and all(good[k:]):
                q += 1
                if k < n: aw[p[k]] += 1
                if k: gw[p[k-1]] += 1
                if 0 < k < n: cross[p[k-1]][p[k]] += 1
    states = []
    for mask in range(1 << n):
        vs = [v for v in range(n) if mask >> v & 1]
        states.append((h, q, sum(aw[v] for v in vs), sum(gw[v] for v in vs),
                       sum(cross[u][v] for u in vs for v in vs),
                       sum(sw[v] for v in vs), sum(ew[v] for v in vs)))
    return states


def main():
    universe = Counter(); total_states = corrected = 0
    for n in range(1, 6):
        for mask in range(1 << (n*(n-1)//2)):
            a = tournament(n, mask)
            result, _ = analyze(a)
            universe[n] += 1
            total_states += result['endpoint_state_orbits']
            corrected += result['corrected_transfers']
    c3 = [2,4,1]
    source_c3 = ordered_sum([0], c3)
    hostile, transfers = analyze(source_c3)
    need((hostile['H'],hostile['aut'],hostile['path_orbits'],hostile['naive_transfer']) == (3,3,1,3), 'hostile')
    changing = [dict(source=b,target=t,N=m,k_source=kb,k_target=ka,W=w)
                for b,t,m,kb,ka,w in transfers if w != m]
    suffix_control, _ = analyze(ordered_sum(c3, [2,0]))
    two_c3 = ordered_sum(c3, c3)
    diagonal = [tuple((v//3)*3 + (v+j)%3 for v in range(6)) for j in range(3)]
    subgroup_control, _ = analyze(two_c3, chosen_group=diagonal)
    need((subgroup_control['H'], subgroup_control['aut'], subgroup_control['path_orbits']) == (9,3,3), 'proper subgroup')
    # For every cut on every tournament <=4, all seven counts divide the cut stabilizer.
    cut_cases = 0
    for n in range(1, 5):
        for mask in range(1 << (n*(n-1)//2)):
            a = tournament(n, mask); group = autos(a)
            for cut, st in enumerate(seven_states(a)):
                stabilizer = sum(image_mask(cut,p) == cut for p in group)
                need(all(x % stabilizer == 0 for x in st), ('seven divisibility',a,cut,st,stabilizer))
                cut_cases += 1
    short_adj = [4,1,18,7,75,31,175,63]
    short_group = autos(short_adj)
    short_state = seven_states(short_adj)[124]
    need(short_state == (45,1440,878,577,332,25,14), 'short state')
    need(len(short_group)==1, 'short core rigidity')
    output = dict(status='PASS', labelled_tournament_counts=dict(universe),
                  total_tournaments=sum(universe.values()), prefix_state_orbits_checked=total_states,
                  corrected_transfer_edges=corrected,
                  hostile_source_triangle=hostile, hostile_edges=changing,
                  suffix_triangle_control=suffix_control, proper_subgroup_control=subgroup_control, cut_state_cases=cut_cases,
                  rigid_cut_nonrigid_core=dict(adj=c3,cut=1,state=seven_states(c3)[1],aut=len(autos(c3)),cut_stabilizer=sum(image_mask(1,p)==1 for p in autos(c3))),
                  session231_short12=dict(adj=short_adj,cut=124,state=short_state,
                                         aut=len(short_group),content=gcd(gcd(short_state[0],short_state[2]),short_state[3])))
    print(json.dumps(output, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
