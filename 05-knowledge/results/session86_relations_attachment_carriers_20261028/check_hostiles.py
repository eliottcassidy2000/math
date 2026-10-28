#!/usr/bin/env python3
"""Independent finite checks for Session86 port-information hostiles.
Only the standard library; all validity gates survive python -O.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
import json

CHECKS = 0

def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def tournament(bits):
    """C3[A_b0,A_b1,A_b2]; A1 is literal converse of A0."""
    n = 4 * len(bits)
    a = [[0] * n for _ in range(n)]
    for u, v in combinations(range(n), 2):
        bu, qu = divmod(u, 4)
        bv, qv = divmod(v, 4)
        if bu == bv:
            forward = qu == 0 or (qu != 0 and qv != 0 and (qv-qu) % 3 == 1)
            forward ^= bool(bits[bu])
        else:
            forward = (bv-bu) % 3 == 1
        a[u][v] = int(forward)
        a[v][u] = 1-int(forward)
    return tuple(map(tuple, a))


def automorphisms(a, colors=None):
    """Exhaustive arc-preserving backtracking; no module assumption."""
    n = len(a)
    colors = tuple([0] * n if colors is None else colors)
    degrees = tuple(map(sum, a))
    domains = [tuple(v for v in range(n) if degrees[u] == degrees[v] and colors[u] == colors[v]) for u in range(n)]
    p = [-1]*n
    used = set()
    result = []
    def visit():
        if len(used) == n:
            result.append(tuple(p))
            return
        options = []
        for u in range(n):
            if p[u] == -1:
                candidates = [v for v in domains[u] if v not in used and all(p[w] == -1 or a[u][w] == a[v][p[w]] for w in range(n))]
                if not candidates:
                    return
                options.append((len(candidates),u,candidates))
        _, u, candidates = min(options)
        for v in candidates:
            p[u] = v
            used.add(v)
            visit()
            used.remove(v)
            p[u] = -1
    visit()
    return result


def cyclic_triples(a):
    return tuple(t for t in combinations(range(len(a)),3) if a[t[0]][t[1]] == a[t[1]][t[2]] == a[t[2]][t[0]])


def path_profiles(a):
    n = len(a)
    dp = [[0]*n for _ in range(1<<n)]
    h = [0]*(1<<n)
    h[0] = 1
    for mask in range(1,1<<n):
        for v in range(n):
            if mask >> v & 1:
                rest = mask ^ (1<<v)
                dp[mask][v] = 1 if rest == 0 else sum(dp[rest][u]*a[u][v] for u in range(n))
        h[mask] = sum(dp[mask])
    return tuple(h), tuple(dp[-1])


def path_cover_profiles(h, n):
    """Unordered nonempty path components, each block weighted by H."""
    pc = [[0]*(n+1) for _ in h]
    pc[0][0] = 1
    for mask in range(1,len(h)):
        anchor = mask & -mask
        sub = mask
        while sub:
            if sub & anchor:
                rest = mask ^ sub
                for d in range(1,mask.bit_count()+1):
                    pc[mask][d] += h[sub]*pc[rest][d-1]
            sub = (sub-1) & mask
    return tuple(map(tuple, pc))


def orbit_partition(group, points):
    unseen = set(points)
    orbits = []
    while unseen:
        x = min(unseen)
        orbit = {tuple(g[v] for v in x) for g in group}
        require(orbit <= unseen, 'orbit overlaps previously removed orbit')
        unseen -= orbit
        orbits.append(tuple(sorted(orbit)))
    return tuple(orbits)


def attach(a, ports):
    """Distinctly colored marker i beats exactly port i among blue vertices."""
    n = len(a)
    b = [list(row)+[0]*len(ports) for row in a]
    b += [[0]*(n+len(ports)) for _ in ports]
    for i,p in enumerate(ports):
        for v in range(n):
            b[n+i][v] = int(v == p)
            b[v][n+i] = int(v != p)
        for j in range(i):
            b[n+j][n+i] = 1
    return tuple(map(tuple,b)), [0]*n + list(range(1,len(ports)+1))


def permutation_order(p):
    q = tuple(range(len(p)))
    for k in range(1,100):
        q = tuple(p[v] for v in q)
        if q == tuple(range(len(p))):
            return k
    raise RuntimeError('unexpected order')


def digest(value):
    return sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()


def cyclic_alignment_checks():
    cube = set(product((0,1), repeat=3))
    rot = lambda x: (x[2],x[0],x[1])
    diagonal = {x for x in cube if rot(x) == x}
    even = {x for x in cube if sum(x) % 2 == 0}
    odd = cube - even
    dot = lambda x,y: sum(a*b for a,b in zip(x,y)) % 2
    require(all(dot(rot(x),rot(y)) == dot(x,y) for x in cube for y in cube), 'cyclic action must preserve common form')
    require({rot(x) for x in odd} == odd and not(even & odd), 'disjoint invariant parity cosets')
    require(len(cube) == 8 and len(diagonal) == 2 and len(even) == 4, 'alignment sizes')
    for relation in (cube,diagonal,even):
        require({rot(x) for x in relation} == relation, 'cyclic invariance')
    pair_scopes = ((0,1),(1,2),(2,0))
    projections = [{tuple(x[i] for i in scope) for x in even} for scope in pair_scopes]
    require(all(len(p) == 4 for p in projections), 'even parity pairwise full')
    for scope in pair_scopes:
        require(Counter(tuple(x[i] for i in scope) for x in even) == Counter(tuple(x[i] for i in scope) for x in odd), 'same exact two-port marginal counts, different joint parity')
    join = {x for x in cube if all(tuple(x[i] for i in scope) in p for scope,p in zip(pair_scopes,projections))}
    require(join == cube, 'marginal join expands even parity')
    require((0,0,1) in cube and rot((0,0,1)) != (0,0,1), 'invariance does not force fixed points')
    return {'rotation_invariant_uniform_full_cube':8,'pointwise_fixed_diagonal':2,'rotation_invariant_even_parity':4,'all_three_pairwise_projection_join':8,'same_form_same_cyclic_action_disjoint_parity_cosets':[sorted(even),sorted(odd)],'cube_rotation_orbit_sizes':sorted(len({x,rot(x),rot(rot(x))}) for x in sorted(cube) if x == min(x,rot(x),rot(rot(x))))}


def main():
    words = list(product((0,1),repeat=3))
    baseline = tournament(words[0])
    c0 = cyclic_triples(baseline)
    h0,e0 = path_profiles(baseline)
    pc0 = path_cover_profiles(h0,12)
    all_cases = []
    saved = {}
    for bits in words:
        a = tournament(bits)
        c = cyclic_triples(a)
        h,e = path_profiles(a)
        pc = path_cover_profiles(h,12)
        require(c == c0, 'literal cyclic-triple mismatch')
        require(h == h0, 'labelled subset Hamiltonian mismatch')
        require(pc == pc0, 'labelled subset path-cover mismatch')
        g = automorphisms(a)
        expected = 81 if len(set(bits)) == 1 else 27
        require(len(g) == expected, 'full automorphism group order')
        require(all(all(a[u][v] == a[p[u]][p[v]] for u,v in combinations(range(12),2)) for p in g), 'automorphism arc verification')
        vertex_orbits = orbit_partition(g, [(i,) for i in range(12)])
        sizes = sorted(map(len,vertex_orbits))
        require(sizes == ([3,9] if expected == 81 else [1,1,1,3,3,3]), 'vertex orbit sizes')
        orders = dict(sorted(Counter(map(permutation_order,g)).items()))
        result = {'leaf_bits':bits,'automorphism_order':len(g),'element_order_histogram':orders,'vertex_orbits':[[x[0] for x in orb] for orb in vertex_orbits],'vertex_orbit_sizes':sizes,'hamiltonian_count':h[-1],'end_vertex_path_counts':e}
        if bits in ((0,0,0),(0,0,1)):
            for r in (1,2,3):
                orbits = orbit_partition(g, permutations(range(12),r))
                result[f'ordered_distinct_{r}_port_orbit_size_histogram'] = dict(sorted(Counter(map(len,orbits)).items()))
                result[f'ordered_distinct_{r}_port_orbit_count'] = len(orbits)
            attachment_aut = []
            attachment_h = []
            for p in range(12):
                b, colors = attach(a,(p,))
                ag = automorphisms(b, colors)
                stabilizer = [q for q in g if q[p] == p]
                require({q[:12] for q in ag} == set(stabilizer), 'role-colored marker attachment stabilizer')
                attachment_aut.append(len(ag))
                attachment_h.append(path_profiles(b)[0][-1])
            result['single_red_marker_automorphism_orders_by_named_port'] = attachment_aut
            result['single_red_marker_hamiltonian_counts_by_named_port'] = attachment_h
            # Full independent colored-graph checks on one representative of every ordered pair orbit.
            two_orbits = orbit_partition(g, permutations(range(12),2))
            for orb in two_orbits:
                ports = orb[0]
                b, colors = attach(a, ports)
                ag = automorphisms(b,colors)
                stabilizer = [q for q in g if all(q[v] == v for v in ports)]
                require({q[:12] for q in ag} == set(stabilizer), 'two-role attachment stabilizer')
            result['two_role_attachment_representatives_checked'] = len(two_orbits)
            saved[bits] = (a,g,h,e,pc)
        all_cases.append(result)
    t, u = saved[(0,0,0)], saved[(0,0,1)]
    group_g, group_h = set(t[1]), set(u[1])
    require(group_h < group_g, 'actual subgroup inclusion')
    outside = group_g - group_h
    require(len(outside) == 54 and all(all(g[v] != v for v in range(12)) for g in outside), 'outside cosets must be fixed-point-free')
    stabilizer_sizes = []
    for mask in range(1,1<<12):
        gs = {g for g in group_g if all(not(mask >> v & 1) or g[v] == v for v in range(12))}
        hs = {g for g in group_h if all(not(mask >> v & 1) or g[v] == v for v in range(12))}
        require(gs == hs, 'actual nonempty pointwise stabilizer subgroup equality')
        stabilizer_sizes.append(len(gs))
    transporter_g = {g for g in group_g if g[0] == 4}
    transporter_h = {g for g in group_h if g[0] == 4}
    require(len(transporter_g) == 27 and not transporter_h, 'missing cross-block transporter')
    require(t[3] != u[3], 'endpoint profile should distinguish')
    require(t[2] == u[2] and t[4] == u[4], 'unmarked profiles must not distinguish')
    bt,_ = attach(t[0],(8,))
    bu,_ = attach(u[0],(8,))
    require((8,9,12) in cyclic_triples(bt) and (8,9,12) not in cyclic_triples(bu), 'explicit marker triple must split')
    require(path_profiles(bt)[0][-1] == 1395279 and path_profiles(bu)[0][-1] == 1061775, 'explicit marker H must split')
    result = {'status':'PASS','universe':'All 8 fixed-label 12-vertex depth-one binary leaf choices; unrestricted automorphism backtracking.', 'source_commit':'5ed1835235b36a57b8cd11eb705aa259b6bc0f08','cyclic_triple_count':len(c0),'literal_cyclic_triples':c0,'all_subset_h_count':len(h0),'all_subset_h_sha256':digest(h0),'all_subset_path_cover_sha256':digest(pc0),'full_path_cover_counts_d0_to_d12':pc0[-1],'tournament_cases':all_cases,'cyclic_and_marginal_hostiles':cyclic_alignment_checks(),'pointwise_stabilizer_hostile':{'literal_subgroup_inclusion':True,'outside_fixed_point_free_permutations':len(outside),'all_nonempty_named_subsets_checked':len(stabilizer_sizes),'pointwise_stabilizer_size_sha256':digest(stabilizer_sizes),'transporter_0_to_4_cardinalities':[len(transporter_g),len(transporter_h)]}}
    result['explicit_checks'] = CHECKS
    print(json.dumps(result,sort_keys=True,indent=2))

if __name__ == '__main__':
    main()
