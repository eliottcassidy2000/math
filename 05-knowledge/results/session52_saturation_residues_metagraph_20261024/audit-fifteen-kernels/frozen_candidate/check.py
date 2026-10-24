#!/usr/bin/env python3
"""Complete finite replay of the fifteen-kernel carrier. Standard library only.
All gates use explicit exceptions and remain active under python -O.
"""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0

def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)

def compose(g, h):
    return tuple(g[h[i]] for i in range(len(g)))

def cycle_type(g):
    unseen = set(range(len(g)))
    lengths = []
    while unseen:
        v = min(unseen)
        k = 0
        while v in unseen:
            unseen.remove(v)
            k += 1
            v = g[v]
        lengths.append(k)
    return tuple(sorted(lengths))

def group_closure(gens):
    identity = tuple(range(len(gens[0])))
    seen, frontier = {identity}, [identity]
    while frontier:
        g = frontier.pop()
        for h in gens:
            gh = compose(g, h)
            if gh not in seen:
                seen.add(gh)
                frontier.append(gh)
    return sorted(seen)

def point_orbits(group, size):
    unseen, out = set(range(size)), []
    while unseen:
        v = min(unseen)
        orbit = sorted({g[v] for g in group})
        unseen.difference_update(orbit)
        out.append(orbit)
    return out

def matching_generator(labels):
    if not labels:
        yield ()
        return
    a = labels[0]
    for b in labels[1:]:
        for tail in matching_generator(tuple(x for x in labels if x not in (a, b))):
            yield tuple(sorted(((a, b),) + tail))

def symplectic(x, y):
    return (((x & 1) * ((y >> 1) & 1)) ^ (((x >> 1) & 1) * (y & 1))
            ^ (((x >> 2) & 1) * ((y >> 3) & 1))
            ^ (((x >> 3) & 1) * ((y >> 2) & 1)))

def quotient(s):
    return min(s, s ^ 63)

def bitpair(x, y):
    return (x & y).bit_count() & 1

def permute_mask(g, s):
    return sum(1 << g[i] for i in range(6) if (s >> i) & 1)

def canon_pair(g, edge):
    return tuple(sorted(g[i] for i in edge))

def canon_matching(g, m):
    return tuple(sorted(canon_pair(g, edge) for edge in m))

def manifest_check():
    path = HERE / 'source-manifest.json'
    data = json.loads(path.read_text())
    for row in data['files']:
        f = HERE / row['path']
        require(f.is_file(), ('input missing', row['path']))
        require(hashlib.sha256(f.read_bytes()).hexdigest() == row['sha256'],
                ('input digest mismatch', row['path']))
    return hashlib.sha256(path.read_bytes()).hexdigest()

def geometry():
    labels = tuple(range(1, 7))
    duads = list(combinations(labels, 2))
    matchings = sorted(matching_generator(labels))
    independent_matchings = sorted(tuple(sorted(t)) for t in combinations(duads, 3)
                                  if len(set(sum((list(e) for e in t), []))) == 6)
    require(matchings == independent_matchings and len(matchings) == 15,
            'two complete matching generators')
    lagrangians = sorted({tuple(sorted((x, y, x ^ y)))
                         for x, y in combinations(range(1, 16), 2)
                         if symplectic(x, y) == 0})
    require(len(lagrangians) == 15, 'binary symplectic Lagrangian count')
    for line in lagrangians:
        require(len(line) == 3 and line[0] ^ line[1] ^ line[2] == 0,
                'Lagrangian additive closure')
        require(all(symplectic(x, y) == 0 for x in line for y in line),
                'Lagrangian isotropy')
    require(all(sum(v in line for line in lagrangians) == 3 for v in range(1, 16)),
            'three Lagrangians through each nonzero vector')
    even = [s for s in range(64) if s.bit_count() % 2 == 0]
    qspace = sorted({quotient(s) for s in even})
    require(len(even) == 32 and len(qspace) == 16, 'even-complement quotient dimensions')
    radical = [s for s in even if all(bitpair(s, t) == 0 for t in even)]
    require(radical == [0, 63], 'radical is the all-ones line')
    e = next(s for s in qspace if s)
    f = next(s for s in qspace if bitpair(e, s))
    orth = [s for s in qspace if not bitpair(e, s) and not bitpair(f, s)]
    g = next(s for s in orth if s)
    h = next(s for s in orth if bitpair(g, s))
    basis = (e, f, g, h)
    coord = {}
    for v in range(16):
        s = 0
        for i, b in enumerate(basis):
            if (v >> i) & 1:
                s ^= b
        require(quotient(s) not in coord, 'quotient symplectic basis independent')
        coord[quotient(s)] = v
    for s in even:
        for t in even:
            require(bitpair(s, t) == symplectic(coord[quotient(s)], coord[quotient(t)]),
                    'quotient pairing agrees with F2^4')
            require(coord[quotient(s ^ t)] == coord[quotient(s)] ^ coord[quotient(t)],
                    'quotient addition agrees with F2^4')
    to_vec = {edge: coord[quotient(sum(1 << (i - 1) for i in edge))] for edge in duads}
    require(set(to_vec.values()) == set(range(1, 16)), 'duads give every nonzero vector')
    images = [tuple(sorted(to_vec[e] for e in m)) for m in matchings]
    require(sorted(images) == lagrangians, 'matching-Lagrangian bijection')
    neighbors = [{j for j, n in enumerate(matchings)
                  if i != j and set(m) & set(n)} for i, m in enumerate(matchings)]
    require(all(len(ns) == 6 for ns in neighbors), 'syntheme graph degree')
    for i, j in combinations(range(15), 2):
        require(len(neighbors[i] & neighbors[j]) == (1 if j in neighbors[i] else 3),
                'SRG adjacent and nonadjacent common-neighbor counts')
    # Independent matrix identity: A^2 = 3I - 2A + 3J.
    adjacency = [[int(j in neighbors[i]) for j in range(15)] for i in range(15)]
    for i in range(15):
        for j in range(15):
            require(sum(adjacency[i][k] * adjacency[k][j] for k in range(15))
                    == 3 * (i == j) - 2 * adjacency[i][j] + 3,
                    'SRG matrix polynomial')
    pentads = [tuple(t) for t in combinations(range(15), 5)
               if len(set(e for i in t for e in matchings[i])) == 15]
    require(len(pentads) == 6, 'six complete one-factorizations of K6')
    pentad_pairs = [tuple(i for i, p in enumerate(pentads) if m in p) for m in range(15)]
    require(all(len(x) == 2 for x in pentad_pairs), 'each matching in two pentads')
    require(len(set(pentad_pairs)) == 15, 'synthemes biject with pentad duads')
    for i, j in combinations(range(15), 2):
        require((j in neighbors[i]) == (not (set(pentad_pairs[i]) & set(pentad_pairs[j]))),
                'explicit graph isomorphism to KG(6,2)')
        require(sum(i in p and j in p for p in pentads) == (0 if j in neighbors[i] else 1),
                'unique pentad for edge-disjoint matchings')
    permutations6 = list(permutations(range(6)))
    mid = {m: i for i, m in enumerate(matchings)}
    pid = {p: i for i, p in enumerate(pentads)}
    matching_actions, phi = {}, {}
    character = {}
    for p in permutations6:
        # Matching objects are labelled 1..6; permutations are indexed 0..5.
        p1 = (0,) + tuple(x + 1 for x in p)
        ma = tuple(mid[canon_matching(p1, m)] for m in matchings)
        matching_actions[p] = ma
        pa = tuple(pid[tuple(sorted(ma[i] for i in pentad))] for pentad in pentads)
        phi[p] = pa
        for i, m in enumerate(matchings):
            require(tuple(sorted(pa[j] for j in pentad_pairs[i])) == pentad_pairs[ma[i]],
                    'twisted equivariance of graph isomorphism')
            moved_vectors = tuple(sorted(coord[quotient(permute_mask(p,
                              quotient(sum(1 << (a - 1) for a in edge))))]
                                         for edge in m))
            require(moved_vectors == images[ma[i]], 'natural matching-Lagrangian equivariance')
        fixed_duads = sum(canon_pair(p1, e) == e for e in duads)
        fixed_matchings = sum(ma[i] == i for i in range(15))
        key = cycle_type(p)
        row = (fixed_duads, fixed_matchings, cycle_type(pa))
        require(key not in character or character[key] == row, 'cycle-type character consistency')
        character[key] = row
    require(len(set(phi.values())) == 720, 'pentad action is faithful and onto S6')
    adjacent = []
    for i in range(5):
        p = list(range(6))
        p[i], p[i+1] = p[i+1], p[i]
        adjacent.append(tuple(p))
    for p in permutations6:
        for a in adjacent:
            require(phi[compose(p, a)] == compose(phi[p], phi[a]), 'homomorphism along generators')
    transposition = adjacent[0]
    require(character[cycle_type(transposition)][:2] == (7, 3), 'fixed-count non-equivariance hostile')
    require(cycle_type(phi[transposition]) == (2, 2, 2), 'outer action witness')
    m0 = matchings[0]
    stab_matching = [p for p in permutations6 if matching_actions[p][0] == 0]
    stab_duad = [p for p in permutations6 if set(p[:2]) == {0, 1}]
    require(len(stab_matching) == len(stab_duad) == 48, 'two stabilizer orders')
    orbit_m = sorted(map(len, point_orbits(stab_matching, 6)))
    orbit_d = sorted(map(len, point_orbits(stab_duad, 6)))
    require(orbit_m == [6] and orbit_d == [2, 4], 'nonconjugate stabilizers in natural S6')
    # Directly generate the wreath stabilizer: swap one pair and permute pairs.
    wreath_generators = [(1,0,2,3,4,5),(2,3,0,1,4,5),(0,1,4,5,2,3)]
    require(group_closure(wreath_generators) == stab_matching, 'S2 wr S3 exact stabilizer')
    return {
        'counts': {'nonzero_vectors': 15, 'lagrangians': 15, 'flags': 45,
                   'synthemes': 15, 'pentads': 6, 'S6_permutations': 720},
        'symplectic_basis_even_subset_masks': basis,
        'duad_to_F2_4': [{'duad': e, 'vector': to_vec[e]} for e in duads],
        'matching_to_lagrangian_and_pentad_duad': [
            {'matching_index': i, 'matching': m, 'nonzero_vectors': images[i],
             'pentad_duad': pentad_pairs[i]} for i, m in enumerate(matchings)],
        'pentads_as_matching_indices': pentads,
        'syntheme_graph': {'srg_parameters': [15,6,1,3], 'edges':45,
                            'spectrum_eigenvalue_multiplicity': [[6,1],[1,9],[-3,5]],
                            'adjacency': adjacency},
        'action_character_table': [{'source_cycle_type': key, 'fixed_duads': row[0],
              'fixed_synthemes': row[1], 'pentad_action_cycle_type': row[2]}
             for key, row in sorted(character.items())],
        'stabilizer_hostile': {'matching': m0, 'matching_stabilizer_order':48,
              'duad_stabilizer_order':48, 'natural_point_orbits_matching':orbit_m,
              'natural_point_orbits_duad':orbit_d},
        'outer_transposition_image_zero_based': phi[transposition]
    }

def bolza():
    data = json.loads((HERE / 'inputs' / 'bolza-kernel-action.json').read_text())
    generators = []
    for name, cycles in data['generators_cycles'].items():
        p = list(range(15))
        used = set()
        for cyc in cycles:
            require(not (used & set(cyc)), ('disjoint cycle transcription', name))
            used.update(cyc)
            for a, b in zip(cyc, cyc[1:]+cyc[:1]):
                p[a-1] = b-1
        require(sorted(p) == list(range(15)), 'source generator is a permutation')
        generators.append(tuple(p))
    group = group_closure(generators)
    require(len(group) == 24, 'published Bolza action group order')
    orbits = [[i+1 for i in orbit] for orbit in point_orbits(group, 15)]
    require(sorted(orbits) == sorted(data['expected_orbits']), 'published Bolza kernel orbit partition')
    records = []
    for orbit in data['expected_orbits']:
        stab = sum(g[orbit[0]-1] == orbit[0]-1 for g in group)
        require(stab * len(orbit) == 24, 'Bolza orbit-stabilizer')
        if len(orbit) == 4:
            require(len({tuple(g[i-1]+1 for i in orbit) for g in group}) == 24,
                    'faithful four-point action identifies S4')
        records.append({'kernel_indices':orbit, 'weight':len(orbit), 'stabilizer_order':stab})
    require(sum(row['weight'] for row in records) == 15, 'weighted out-degree not edge count')
    require(len(records) == 4, 'four automorphism edge classes')
    return {'status':'FINITE-EXACT group replay of CITED geometric input; not a field or curve computation',
            'reduced_group_order':24, 'orbits':records,
            'edge_class_count':4, 'weighted_out_degree':15,
            'codomain_identification_status':'CITED only; see source manifest and source input'}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--receipt', default=str(HERE / 'receipt.json'))
    args = parser.parse_args()
    manifest_hash = manifest_check()
    geo = geometry()
    bol = bolza()
    result = {'status':'PROVED elementary geometry; FINITE-EXACT full universes; CITED Richelot/Bolza interpretation',
              'source_manifest_sha256':manifest_hash,
              'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'geometry':geo, 'bolza':bol,
              'boundaries':['No bijection to the Ogg prime set is defined.',
                            'No preferred modulus 15 is obtained.',
                            'No formula for or correction to the MRB constant is obtained.',
                            'The syntheme graph is not the Richelot isogeny graph.',
                            'No rationality, superspeciality, or field specialization is proved by this finite checker.',
                            'Geometric codomain types are cited, not recomputed.'],
              'checks':CHECKS}
    encoded=json.dumps(result, indent=2, sort_keys=True)+'\n'
    Path(args.receipt).write_text(encoded)
    print(encoded,end='')

if __name__ == '__main__':
    main()
