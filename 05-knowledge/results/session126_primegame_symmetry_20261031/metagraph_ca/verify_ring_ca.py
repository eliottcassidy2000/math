#!/usr/bin/env python3
"""Exact n=5 fixed-path carrier, perimeter CA, and observer quotients.

Standard library only. All verification gates survive python -O. No producer
imports, no inherited atlas input, no network, and no repository operations.
Outputs beside this file unless --output-dir is supplied.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations
from pathlib import Path
import argparse
import csv
import json

N = 5
FREE = [(5, 1), (4, 1), (3, 1), (5, 2), (4, 2), (5, 3)]
PERIMETER = [0, 1, 2, 4, 5, 3]
REFLECTION = [0, 3, 5, 1, 4, 2]
PAIRS = list(combinations(range(N), 2))
PERMS = list(permutations(range(N)))


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def matrix(mask):
    """Rows encode outgoing arcs. Mask zero is ascending transitive order."""
    rows = [sum(1 << j for j in range(i+1, N)) for i in range(N)]
    for k, (a, b) in enumerate(FREE):
        if mask >> k & 1:
            a, b = a-1, b-1
            rows[b] ^= 1 << a
            rows[a] ^= 1 << b
    return rows


def canonical(rows):
    return min(sum(((rows[p[i]] >> p[j]) & 1) << k
                   for k, (i, j) in enumerate(PAIRS)) for p in PERMS)


def converse(rows):
    return [((1 << N)-1) ^ row ^ (1 << i) for i, row in enumerate(rows)]


def reflect(mask):
    return sum(((mask >> i) & 1) << j for i, j in enumerate(REFLECTION))


def complement(mask):
    return mask ^ 63


def local_rule(mask, center):
    return sum((((mask >> PERIMETER[(i-1) % 6]) ^
                 (center * (mask >> PERIMETER[i])) ^
                 (mask >> PERIMETER[(i+1) % 6])) & 1) << PERIMETER[i]
               for i in range(6))


def a90(mask):
    return local_rule(mask, 0)


def b150(mask):
    return local_rule(mask, 1)


def projection(mask):
    e = ((mask >> 0) ^ (mask >> 2) ^ (mask >> 5)) & 1
    o = ((mask >> 1) ^ (mask >> 3) ^ (mask >> 4)) & 1
    return e*37 ^ o*26


def canonical_partition(signatures):
    lookup, ids, blocks = {}, [], []
    for mask, signature in enumerate(signatures):
        if signature not in lookup:
            lookup[signature] = len(blocks)
            blocks.append([])
        idx = lookup[signature]
        ids.append(idx)
        blocks[idx].append(mask)
    return ids, blocks


def continuation_quotient(observations, controls):
    """Coarsest observation-preserving deterministic quotient by refinement."""
    ids, blocks = canonical_partition(observations)
    rounds = [len(blocks)]
    while True:
        new, newblocks = canonical_partition([
            (ids[x], tuple(ids[f(x)] for f in controls.values()))
            for x in range(64)])
        if len(newblocks) == len(blocks):
            need(new == ids, 'canonical stable partition')
            break
        ids, blocks = new, newblocks
        rounds.append(len(blocks))
    transitions = {name: [ids[f(block[0])] for block in blocks]
                   for name, f in controls.items()}
    for block in blocks:
        need(len({observations[x] for x in block}) == 1, 'observer preserved')
        for name, f in controls.items():
            need(len({ids[f(x)] for x in block}) == 1,
                 'well-defined quotient transition ' + name)
    return {'states': len(blocks), 'refinement_counts': rounds,
            'block_size_histogram': dict(sorted(Counter(map(len, blocks)).items())),
            'blocks': blocks, 'class_by_mask': ids,
            'observations': [observations[b[0]] for b in blocks],
            'transitions': transitions}


def dynamics(fn):
    cycles, depths = set(), Counter()
    for start in range(64):
        path, pos, x = [], {}, start
        while x not in pos:
            pos[x] = len(path)
            path.append(x)
            x = fn(x)
        cycle = path[pos[x]:]
        rotations = [tuple(cycle[i:] + cycle[:i]) for i in range(len(cycle))]
        cycles.add(min(rotations))
        depths[pos[x]] += 1
    return {'rank': (len({fn(x) for x in range(64)})).bit_length()-1,
            'image': sorted({fn(x) for x in range(64)}),
            'cycles': sorted(cycles),
            'cycle_length_histogram': dict(sorted(Counter(map(len, cycles)).items())),
            'transient_depth_histogram': dict(sorted(depths.items()))}


def build():
    mats = [matrix(x) for x in range(64)]
    canonical_classes = [canonical(m) for m in mats]
    cids, classblocks = canonical_partition(canonical_classes)
    merged_keys = [min(canonical(m), canonical(converse(m))) for m in mats]
    pi, blocks = canonical_partition(merged_keys)
    need(len(classblocks) == 12 and len(blocks) == 10, '12 classes -> 10 nodes')
    # Independent all-labelled census, not merely fixed-path coverage.
    all_canon = set()
    for word in range(1024):
        rows = [0]*N
        for k, (i, j) in enumerate(PAIRS):
            if word >> k & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
        all_canon.add(canonical(rows))
    need(len(all_canon) == 12, '1024 labelled tournaments -> 12 classes')
    need(set(canonical_classes) == all_canon, 'fixed path covers every class')

    blue = [reflect(x) == x for x in range(64)]
    pure_blue = [i for i, b in enumerate(blocks) if all(blue[x] for x in b)]
    need(pure_blue == [0, 7, 8], 'pure-blue nodes')
    need([blocks[i] for i in pure_blue] == [[0], [11, 46, 53], [16]], 'PB fibres')
    edges = {'blue': Counter(), 'black': Counter()}
    for x in range(32):
        edges['blue' if blue[x] else 'black'][tuple(sorted((pi[x], pi[x^63])))] += 1
    need(sum(edges['blue'].values()) == 8 and sum(edges['black'].values()) == 24,
         'complement-line counts')

    for x in range(64):
        need(reflect(reflect(x)) == x, 'R involution')
        need(reflect(complement(x)) == complement(reflect(x)), 'RC=CR')
        need(pi[reflect(x)] == pi[x], 'R invisible to merged observer')
        need(a90(a90(a90(x))) == a90(x), 'A^3=A')
        need(b150(b150(b150(x))) == b150(b150(x)), 'B^3=B^2')
        need(b150(b150(x)) == projection(x), 'explicit B^2 projection')
        for f in [a90, b150]:
            need(f(reflect(x)) == reflect(f(x)), 'reflection equivariance')
        need(a90(complement(x)) == a90(x), 'A erases complement')
        need(b150(complement(x)) == complement(b150(x)), 'B complement equivariance')
    need([REFLECTION[PERIMETER[i]] for i in range(6)] ==
         [PERIMETER[(-i) % 6] for i in range(6)], 'R is perimeter reflection')

    da, db = dynamics(a90), dynamics(b150)
    need(da['rank'] == 4 and da['cycle_length_histogram'] == {1: 4, 2: 6},
         'A rank and cycles')
    need(da['transient_depth_histogram'] == {0: 16, 1: 48}, 'A transient counts')
    need(db['rank'] == 4 and db['cycle_length_histogram'] == {1: 4}, 'B dynamics')
    need(db['transient_depth_histogram'] == {0: 4, 1: 12, 2: 48}, 'B depths')
    need(sorted({projection(x) for x in range(64)}) == [0,26,37,63], 'rank2 image')
    need(Counter(projection(x) for x in range(64)) == {0:16,26:16,37:16,63:16},
         'projection fibres')
    jordan_basis = [37,26,46,10,29,5]
    combinations_of_basis = {0}
    for v in jordan_basis:
        combinations_of_basis |= {x^v for x in list(combinations_of_basis)}
    need(len(combinations_of_basis)==64, 'explicit Jordan basis spans carrier')
    need([b150(v) for v in jordan_basis] == [37,26,0,46,0,29],
         'B = I2 plus two nilpotent size2 blocks')

    observers = {'pi': pi, 'pi_blue': [(pi[x],blue[x]) for x in range(64)]}
    controls = {'CR': {'C':complement,'R':reflect}, 'A': {'A':a90}, 'B':{'B':b150},
                'AB':{'A':a90,'B':b150},
                'ACR':{'A':a90,'C':complement,'R':reflect},
                'BCR':{'B':b150,'C':complement,'R':reflect},
                'ABCR':{'A':a90,'B':b150,'C':complement,'R':reflect},
                'coordinate_flips': {str(i): lambda x,i=i:x^(1<<i) for i in range(6)}}
    quotients = {on: {cn:continuation_quotient(obs,fs) for cn,fs in controls.items()}
                 for on,obs in observers.items()}
    expected = {'pi': {'CR':38,'A':29,'B':30,'AB':37,'ACR':38,'BCR':40,'ABCR':40,
                        'coordinate_flips':64},
                'pi_blue': {'CR':38,'A':29,'B':36,'AB':37,'ACR':38,'BCR':40,
                            'ABCR':40,'coordinate_flips':64}}
    for on, cs in expected.items():
        for cn,n in cs.items():
            need(quotients[on][cn]['states'] == n, 'quotient count '+on+'/'+cn)
    signatures = [(pi[x], pi[x^63]) for x in range(64)]
    need(canonical_partition(signatures)[0] == quotients['pi']['CR']['class_by_mask'],
         'minimal CR signature')
    need(canonical_partition([(pi[x],pi[a90(x)],pi[a90(a90(x))])
                               for x in range(64)])[0] ==
         quotients['pi']['A']['class_by_mask'], 'minimal A signature')
    need(canonical_partition([(pi[x],pi[b150(x)],pi[b150(b150(x))])
                               for x in range(64)])[0] ==
         quotients['pi']['B']['class_by_mask'], 'minimal B signature')
    reflection_ids = canonical_partition([min(x,reflect(x)) for x in range(64)])[0]
    need(reflection_ids == quotients['pi']['BCR']['class_by_mask'], 'BCR=R orbit quotient')
    need(quotients['pi']['CR']['class_by_mask'] == quotients['pi_blue']['CR']['class_by_mask'],
         'CR signature already determines blue')

    # Coordinate-labelled flip words separate every distinct pair by translating
    # the first mask to the singleton transitive fibre {0}.
    for x in range(64):
        for y in range(x):
            need(pi[x^x] != pi[y^x], 'full coordinate-control separation')

    # Enumerate all eight homogeneous nearest-neighbour F2-linear rules.
    selected = []
    for left in (0,1):
        for center in (0,1):
            for right in (0,1):
                f = lambda x: sum((((left*(x>>PERIMETER[(i-1)%6])) ^
                                     (center*(x>>PERIMETER[i])) ^
                                     (right*(x>>PERIMETER[(i+1)%6])))&1)
                                   << PERIMETER[i] for i in range(6))
                is_symmetric = all(f(reflect(x)) == reflect(f(x)) for x in range(64))
                is_c_equivariant = all(f(x^63) == (f(x)^63) for x in range(64))
                if left and right and is_symmetric and is_c_equivariant:
                    selected.append([left,center,right])
    need(selected == [[1,1,1]], 'unique active-neighbour linear symmetric C-equivariant rule')

    rows = []
    for x in range(64):
        rows.append({'mask':x,'explorer_bits_lsb_first':''.join(str(x>>i&1) for i in range(6)),
                     'perimeter_bits':''.join(str(x>>i&1) for i in PERIMETER),
                     'tournament_class':cids[x],'merged_node':pi[x],
                     'blue':int(blue[x]),'C':x^63,'R':reflect(x),
                     'A':a90(x),'A2':a90(a90(x)),'B':b150(x),'B2':projection(x),
                     'pi_A':pi[a90(x)],'pi_B':pi[b150(x)],
                     'pi_C':pi[x^63]})
    report = {'status':'PASS','conventions':{'n':N,'explorer_cells':FREE,
                'mask_zero':'ascending transitive tournament; set bit reverses named nonconsecutive arc',
                'perimeter_order':PERIMETER,'reflection_permutation':REFLECTION,
                'C':'xor63','A':'rule90 on C6, T+T^-1','B':'rule150 on C6, I+T+T^-1'},
              'census':{'labelled_tournaments':1024,'classes':12,'merged_nodes':10,
                        'carrier_size':64,'class_by_mask':cids,'pi':pi,'node_fibres':blocks,
                        'pure_blue_nodes':pure_blue,'pure_blue_fibres':[blocks[i] for i in pure_blue],
                        'blue_fixed_masks':[x for x in range(64) if blue[x]],
                        'edges':{col:[[a,b,n] for (a,b),n in sorted(es.items())]
                                 for col,es in edges.items()}},
              'A':da,'B':db,
              'B2':{'rank':2,'fixed_masks':[0,26,37,63],
                     'basins':{str(v):[x for x in range(64) if projection(x)==v]
                               for v in [0,26,37,63]},
                     'even_perimeter_explorer_bits':[0,2,5],
                     'odd_perimeter_explorer_bits':[1,3,4]},
              'linear_decomposition':{'basis_explorer_masks':jordan_basis,
                  'B_images_of_basis':[b150(v) for v in jordan_basis],
                  'A_images_of_basis':[a90(v) for v in jordan_basis],
                  'B_block_form':'I_2 direct_sum J_2(0) direct_sum J_2(0)',
                  'A_block_form':'0_2 direct_sum J_2(1) direct_sum J_2(1)',
                  'complementary_projections':'P=B^2=I+A^2 onto ker A; Q=A^2 onto im A'},
              'requested_masks':[{k:rows[x][k] for k in ['mask','merged_node','blue','A','A2','B','B2']}
                                 for x in [7,17,55]],
              'PB_descent_counterexamples':[{k:rows[x][k] for k in ['mask','merged_node','A','pi_A','B','pi_B']}
                                           for x in [11,46,53]],
              'quotient_state_counts':expected,
              'CR_four_element_blocks':[b for b in quotients['pi']['CR']['blocks'] if len(b)==4],
              'unique_linear_rule_coefficients':selected}
    return report, quotients, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output-dir', type=Path, default=Path(__file__).parent)
    args = ap.parse_args()
    report, quotients, rows = build()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name,data in [('ring_ca_results.json',report),('observer_quotients.json',quotients)]:
        (args.output_dir/name).write_text(json.dumps(data,separators=(',',':'),sort_keys=True)+'\n')
    with (args.output_dir/'transitions.tsv').open('w',newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(json.dumps({'status':report['status'],'census':{k:report['census'][k]
                      for k in ['labelled_tournaments','classes','merged_nodes','carrier_size',
                                'pure_blue_fibres']},
                      'A':report['A'],'B':report['B'],
                      'requested_masks':report['requested_masks'],
                      'quotient_state_counts':report['quotient_state_counts']},
                     indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
