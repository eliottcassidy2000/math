"""Full graph group of the pinned625 source by exhaustive rooted refinement."""
from collections import Counter
from pathlib import Path
import argparse
import json
import hashlib
import check as c

HERE = Path(__file__).resolve().parent


def source_rows():
    q = 25
    points = list(c.product(range(q), repeat=2))
    left, owner, index = set(points)-{(0, 0)}, {}, 0
    while left:
        v, orb = min(left), []
        while v not in orb:
            orb.append(v)
            v = c.mul((2, 3), v, q, 1)
        for x, y in orb:
            owner[x, y] = index, 1
            owner[-x % q, -y % q] = index, -1
        left -= set(orb) | {(-x % q, -y % q) for x, y in orb}
        index += 1
    eps = json.loads((c.RESULTS/'source_cutoffs_balanced_lifts_20261029/lifts.json').read_text())['level25']['signs']
    conn = {v for v, (i, sign) in owner.items() if eps[i]*sign == 1}
    return [sum(1 << (((x+u) % q)*q+(y+v) % q) for u, v in conn) for x, y in points]


def individualize(parts, cell, vertex):
    return [d for group in parts for d in (([vertex], [x for x in group if x != vertex])
                                          if group == cell else (group,))]


def universal_socle_checks():
    reports = []
    for p, k, a in ((5, 1, 1), (5, 2, 1), (5, 3, 1), (13, 1, 3), (13, 2, 3), (17, 1, 8)):
        q, powers = c.clock(p, k, a)
        z = powers[p**(k-1)]
        eta = ((z[0]-1) % q, z[1])
        residue = lambda v: (v[0]+a*pow(2, -1, p)*v[1]) % p
        c.need(eta != (0, 0) and c.mul(eta, eta, q, a) == (0, 0), 'nonzero square-zero socle')
        c.need(residue(eta) == 0 and c.mul(powers[1], eta, q, a) == eta, 'clock fixes socle pointwise')
        units = 0
        for v in c.product(range(q), repeat=2):
            c.need(c.mul(eta, v, q, a) == tuple(x*residue(v) % q for x in eta), 'socle product factors through residue')
            if residue(v):
                units += 1
                orbit = {c.mul(powers[j*p**(k-1)], v, q, a) for j in range(p)}
                line = {tuple((v[i]+j*eta[i]) % q for i in range(2)) for j in range(p)}
                c.need(orbit == line and len(line) == p, 'all unit differences admit every socle shift')
        reports.append({'p': p, 'k': k, 'a': a, 'eta': eta,
                        'all_ring_elements': q*q, 'all_units': units})
    return reports


def enlarged_cycle_checks():
    reports = []
    # Complete F_2 and all A_2/(A_2 intersect F_2) cosets at p=5;
    # three nonlinear fibre laws per nonidentity coset, not all of H_2.
    p, k, a = 5, 2, 1
    q, powers = c.clock(p, k, a)
    eta = ((powers[p][0]-1) % q, powers[p][1])
    points = list(c.product(range(q), repeat=2))
    residue = lambda v: (v[0]+3*v[1]) % p

    def permutation(t, b, h):
        result = []
        for x in points:
            y = c.mul(powers[t], x, q, a)
            u = h[residue(x)]
            result.append(((y[0]+b[0]+eta[0]*u) % q)*q
                          +(y[1]+b[1]+eta[1]*u) % q)
        return result

    actual_f = Counter()
    for h in c.product(range(p), repeat=p):
        actual_f[c.cycle_type(permutation(0, (0, 0), h))] += 1
    formula_f = Counter()
    from math import comb
    for j in range(p+1):
        typ = tuple((length, count) for length, count in
                    ((1, (p-j)*p**(2*k-1)), (p, j*p**(2*k-2))) if count)
        formula_f[typ] = comb(p, j)*(p-1)**j
    c.need(actual_f == formula_f, 'complete fibre group cycle inventory')
    quotient_reps = sorted({min(tuple((x[i]+eta[i]*u) % q for i in range(2))
                               for u in range(p)) for x in points})
    c.need(len(quotient_reps) == p**(2*k-1), 'all translation quotient representatives')
    laws = ((1, 0, 0, 0, 0), (0, 1, 4, 4, 1), (4, 1, 3, 0, 2))
    checked = 0
    coset_types = Counter()
    for t in range(p**(k-1)):
        for b in quotient_reps:
            if t == 0 and b == (0, 0):
                continue
            expected = c.cycle_type(permutation(t, b, (0,)*p))
            coset_types[expected] += 1
            for h in laws:
                c.need(c.cycle_type(permutation(t, b, h)) == expected, 'nonidentity coset cycle type independent of fibre law')
                checked += 1
    affine = c.formula_inventory(p, k)
    intersection = Counter({((1, 625),): 1, ((1, 125), (5, 100)): 20, ((5, 125),): 4})
    actual_coset_mass = Counter({typ: count*p**p for typ, count in coset_types.items()})+actual_f
    formula_h = Counter({typ: (count-intersection[typ])*p**(p-2) for typ, count in affine.items()})+formula_f
    c.need(actual_coset_mass == formula_h and sum(formula_h.values()) == 5**9, 'enlarged inventory mass and independent coset formula')
    c.need(min(625-sum(n for _, n in typ) for typ in formula_h if typ != ((1, 625),)) == 100, 'enlarged minimum index')
    coloring_orbits = {}
    for Q in (2, 3, 4):
        total = sum(count*Q**sum(n for _, n in typ) for typ, count in formula_h.items())
        c.need(total % 5**9 == 0, 'full625 coloring Burnside divisibility')
        coloring_orbits[str(Q)] = str(total//5**9)
    # The k>=2 boundary is necessary: a top-residue cocycle creates a25-cycle.
    q1, powers1 = c.clock(5, 1, 1)
    eta1 = ((powers1[1][0]-1) % 5, powers1[1][1])
    hostile = [((x+1+eta1[0]*((x+3*y) % 5 == 0)) % 5)*5
               +(y+eta1[1]*((x+3*y) % 5 == 0)) % 5 for x, y in c.product(range(5), repeat=2)]
    c.need(c.cycle_type(hostile) == ((25, 1),), 'depth-one hostile to coset cycle constancy')
    higher = []
    for p2, k2, a2 in ((5, 3, 1), (13, 2, 3)):
        q2, powers2 = c.clock(p2, k2, a2)
        eta2 = ((powers2[p2**(k2-1)][0]-1) % q2, powers2[p2**(k2-1)][1])
        tests = ((0, (1, 0)), (0, (p2**(k2-1), 0)), (1, (0, 0)),
                 (p2, (p2, 0)), (1, (p2, 0)))
        for t, b in tests:
            plain, changed = [], []
            for x in c.product(range(q2), repeat=2):
                y = c.mul(powers2[t], x, q2, a2)
                z = (x[0]+a2*pow(2, -1, p2)*x[1]) % p2
                h = (z*z+2*z+1) % p2
                plain.append(((y[0]+b[0]) % q2)*q2+(y[1]+b[1]) % q2)
                changed.append(((y[0]+b[0]+eta2[0]*h) % q2)*q2
                               +(y[1]+b[1]+eta2[1]*h) % q2)
            typ = c.cycle_type(plain)
            c.need(c.cycle_type(changed) == typ, 'higher nonlinear coset control')
            higher.append({'p': p2, 'k': k2, 'a': a2, 't': t, 'b': b,
                           'h': 'z^2+2z+1', 'cycles': list(typ)})
    return {'all_fibre_permutations': 5**5, 'all_nonidentity_cosets': len(quotient_reps)*5-1,
            'nonlinear_coset_probes': checked, 'full_group_inventory':
            [{'cycles': list(typ), 'multiplicity': count} for typ, count in sorted(formula_h.items())],
            'full625_Q_color_orbits': coloring_orbits,
            'full625_prime_child_classes': str(int(coloring_orbits['2'])-4),
            'higher_nonlinear_controls': higher,
            'depth_one_hostile_cycles': [[25, 1]]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'parent_group.json')
    args = parser.parse_args()
    rows = source_rows()
    old = json.loads((HERE/'certificate.json').read_text())['large_native']
    c.need(c.digest(list(map(hex, rows))) == old['parent_relation_sha256'], 'same pinned native source')
    rooted, _ = c.refinement(rows, [[0], list(range(1, 625))])
    cell = next(group for group in rooted if 1 in group)
    orbit, v = set(), (0, 1)
    while v[0]*25+v[1] not in orbit:
        orbit.add(v[0]*25+v[1])
        v = c.mul((2, 3), v, 25, 1)
    c.need(set(cell) == orbit and len(cell) == 25, 'root-stabilizer orbit exactly25')
    fixed, _ = c.refinement(rows, individualize(rooted, cell, 1))
    leaves, profiles = [], Counter()

    def visit(parts, address):
        c.need(len(leaves) < 10000, 'exhaustive branch safety cap')
        parts, _ = c.refinement(rows, parts)
        profiles[tuple(sorted(Counter(map(len, parts)).items()))] += 1
        if len(parts) == 625:
            leaves.append((address, [group[0] for group in parts]))
            return
        cell = min((group for group in parts if len(group) > 1), key=lambda group: len(group))
        for vertex in cell:
            visit(individualize(parts, cell, vertex), address+[vertex])

    visit(fixed, [])
    baseline = leaves[0][1]
    automorphisms, rejected = [], 0
    for address, order in leaves:
        perm = [0]*625
        for x, y in zip(baseline, order):
            perm[x] = y
        c.need(perm[0] == 0 and perm[1] == 1 and sorted(perm) == list(range(625)), 'two marked vertices and bijection')
        valid = all(((rows[x] >> y) & 1) == ((rows[perm[x]] >> perm[y]) & 1)
                    for x in range(625) for y in range(x+1, 625))
        if valid:
            automorphisms.append(perm)
        else:
            rejected += 1
    c.need(len({tuple(p) for p in automorphisms}) == len(automorphisms), 'distinct certified maps')
    # An independent formula reconstructs every one of the enumerated maps.
    # Here ell(x,y)=2x+y mod5, d=(5,15), and h(0)=h(1)=0.
    formula_maps = set()
    for values in c.product(range(5), repeat=3):
        h = (0, 0)+values
        perm = tuple(((x+5*h[(2*x+y) % 5]) % 25)*25
                     +(y+15*h[(2*x+y) % 5]) % 25
                     for x, y in c.product(range(25), repeat=2))
        formula_maps.add(perm)
    c.need(formula_maps == {tuple(p) for p in automorphisms}, 'complete stabilizer equals independent fibre shifts')
    for residue in range(5):
        perm = [((x+5*((2*x+y) % 5 == residue)) % 25)*25
                +(y+15*((2*x+y) % 5 == residue)) % 25
                for x, y in c.product(range(25), repeat=2)]
        c.need(all(((rows[x] >> y) & 1) == ((rows[perm[x]] >> perm[y]) & 1)
                   for x in range(625) for y in range(x+1, 625)), 'literal single-fibre generator preserves all arcs')
    socle_reports = universal_socle_checks()
    enlarged_report = enlarged_cycle_checks()
    answer = {'status': 'FINITE-EXACT full graph group, exhaustive refinement after intrinsic marked orbit',
              'source_relation_sha256': old['parent_relation_sha256'], 'root_orbit_size': 25,
              'complete_leaves': len(leaves), 'rejected_maps': rejected,
              'two_point_stabilizer_order': len(automorphisms),
              'full_graph_automorphism_order': 625*25*len(automorphisms),
              'rooted_profile': dict(Counter(map(len, rooted))),
              'two_point_profile': dict(Counter(map(len, fixed))),
              'refinement_profile_counts': [(list(key), value) for key, value in sorted(profiles.items())],
              'two_point_group_structure': 'C5^3: x -> x+(5,15)h(2x+y), h(0)=h(1)=0',
              'full_group_structure': 'H_2=F_2 A_2, |F_2|=5^5, |F_2 intersect A_2|=5^2',
              'universal_socle_controls': socle_reports,
              'enlarged_cycle_controls': enlarged_report,
              'complete_stabilizer_maps': automorphisms,
              'complete_leaf_addresses': [address for address, _ in leaves], 'gates': c.GATES}
    text = json.dumps(answer, indent=2, sort_keys=True)+'\n'
    args.output.write_bytes(text.encode())
    print(json.dumps({'leaves': len(leaves), 'two_point_stabilizer': len(automorphisms),
                      'full_graph_group': answer['full_graph_automorphism_order'],
                      'sha256': hashlib.sha256(text.encode()).hexdigest()}))


if __name__ == '__main__':
    main()
