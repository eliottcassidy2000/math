#!/usr/bin/env python3
"""Independent finite audit: support partitions, triangular determinant,
integer combinations and local-replacement union-find. No discovery imports.
"""
import argparse
from collections import defaultdict, Counter
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys

CHECKS = 0


def need(v, msg):
    global CHECKS
    CHECKS += 1
    if not v:
        raise RuntimeError(msg)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--producer', type=Path, required=True)
    parser.add_argument('--supply', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    producer, repo, supply = args.producer.resolve(), args.repo.resolve(), args.supply.resolve()
    original = repo/'05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json'
    electrical = json.loads(original.read_text())['order8']
    X = electrical['higher_layers'][-1]['lifts'][0]
    Y = electrical['radical_lifts_mod4'][:4]
    records = json.loads(supply.read_text())['records']
    sources = [r['source'] for r in records]
    need(len(sources) == 28 and sources == sorted(set(sources)), 'Exactly 28 named distinct canonical source slots')
    labels = []
    for row in records:
        s = row['source']
        q = [v[s] % 2 for v in Y] + [X[s] % 128]
        need(q == row['source_label'], 'Actual marked label matches original electrical input')
        labels.append((sum(q[j] << j for j in range(4)),q[-1]))
    need(len(set(labels)) == 28, 'Distinct group labels forbid a support-two trade')

    def key(indices):
        parity = mass = 0
        for i in indices:
            parity ^= labels[i][0]
            mass += labels[i][1]
        return parity,mass % 128

    # Independent enumeration: first choose the whole changed support, then
    # every balanced sign partition. The parity guard is a support-only sieve.
    direct = {}
    support_counts = {}
    for k in (2,3):
        vectors = set()
        tested = parity_survivors = 0
        for support in itertools.combinations(range(28),2*k):
            tested += 1
            all_parity, all_mass = key(support)
            if all_parity:
                continue
            parity_survivors += 1
            for tail in itertools.combinations(support[1:],k-1):
                positive = (support[0],)+tail
                if (2*key(positive)[1]-all_mass) % 128:
                    continue
                negative = tuple(i for i in support if i not in positive)
                need(key(positive) == key(negative), 'Support partition satisfies full marked charge')
                vectors.add(tuple(int(i in positive)-int(i in negative) for i in range(28)))
        direct[k] = vectors
        support_counts[k] = {'supports_tested': tested, 'parity_surviving_supports': parity_survivors, 'balanced_oriented_moves':len(vectors)}
    need(len(direct[2]) == 27 and len(direct[3]) == 1917, 'Independent exact support-four/six census')

    # Recover the documented lexical order only for the integer witness.
    ordered = {}
    buckets = {}
    for k in (1,2,3):
        groups = defaultdict(list)
        for s in itertools.combinations(range(28),k):
            groups[key(s)].append(s)
        buckets[k] = {q:[sum(1 << i for i in s) for s in sets] for q,sets in groups.items()}
        if k == 1:
            continue
        rows = []
        for sets in groups.values():
            for a,b in itertools.combinations(sets,2):
                if set(a).isdisjoint(b):
                    rows.append(tuple(int(i in a)-int(i in b) for i in range(28)))
        need(set(rows) == direct[k], 'Independent support census agrees with witness ordering')
        ordered[k] = rows

    # Independent exact RREF, eliminating above AND below pivots.
    matrix = [list(map(Fraction,r)) for r in ordered[2]]
    pivots = []
    for col in range(28):
        row = next((j for j in range(len(pivots),len(matrix)) if matrix[j][col]),None)
        if row is None:
            continue
        p = len(pivots)
        matrix[p],matrix[row] = matrix[row],matrix[p]
        denominator = matrix[p][col]
        matrix[p] = [x/denominator for x in matrix[p]]
        for j in range(len(matrix)):
            if j != p and matrix[j][col]:
                factor = matrix[j][col]
                matrix[j] = [x-factor*y for x,y in zip(matrix[j],matrix[p])]
        pivots.append(col)
    need(len(pivots) == 23, 'Exact independent support-four rank is 23')
    i,j = sources.index(1638),sources.index(1974)
    need(all(v[i]+v[j] == 0 for v in direct[2]), 'All support-four moves preserve x1638+x1974 over the integers')

    witness = json.loads((producer/'integer_basis_witness.json').read_text())
    need(witness['prefix_moves'] == 48 and witness['dimension'] == 27, 'Witness declared ambient dimensions')
    moves = (ordered[2]+ordered[3])[:48]
    columns = witness['basis_columns']
    combinations = witness['integer_combinations']
    need(len(columns) == len(combinations) == 27, '27 retained basis witnesses')
    for column,coefficients in zip(columns,combinations):
        need(len(column) == 27 and len(coefficients) == 48, 'Individual witness dimensions')
        need(all(type(v) is int for v in column+coefficients), 'Witness coefficients are exact integers')
        full_column = column+[-sum(column)]
        recovered = [sum(c*v[r] for c,v in zip(coefficients,moves)) for r in range(28)]
        need(recovered == full_column, 'Integer combination verifies all 28 coordinates, including dropped margin coordinate')
        need(sum(full_column) == 0, 'Certified basis lies in exact zero-margin hyperplane')
        need(sum(full_column[r]*labels[r][1] for r in range(28)) % 128 == 0, 'Basis cyclic guard')
        need(all(sum(full_column[r]*((labels[r][0] >> b)&1) for r in range(28)) % 2 == 0 for b in range(4)), 'Basis binary guards')
    # The supplied basis is upper triangular in column form. This determinant
    # proof is independent of the producer's Bareiss elimination.
    need(all(columns[col][row] == 0 for row in range(27) for col in range(row)), 'Exact upper-triangular integer basis')
    diagonal = [columns[i][i] for i in range(27)]
    determinant = 1
    for x in diagonal:
        determinant *= x
    need(abs(determinant) == 2048, 'Independent product-of-diagonal determinant gives index 2048')

    # Explicit subgroup closure, using differences against the first source
    # rather than the producer's final-source basis.
    differences = [(q^labels[0][0],(r-labels[0][1]) % 128) for q,r in labels[1:]]
    image = {(0,0)}
    for parity,mass in differences:
        previous = list(image)
        image = {(q^(parity if t % 2 else 0),(r+t*mass) % 128) for q,r in previous for t in range(128)}
    need(len(image) == 2048, 'Independent finite image has the full authenticated Q order')

    fiber_reports = []
    hostile = None
    for cardinality in (2,3,4):
        subsets = list(itertools.combinations(range(28),cardinality))
        masks = [sum(1 << i for i in s) for s in subsets]
        index = {mask:i for i,mask in enumerate(masks)}
        parents = list(range(len(masks)))

        def find(a):
            while a != parents[a]:
                parents[a] = parents[parents[a]]
                a = parents[a]
            return a

        def union(a,b):
            a,b = find(a),find(b)
            if a != b:
                parents[b] = a

        replacements = 0
        for state_index,(selected,mask) in enumerate(zip(subsets,masks)):
            for k in range(1,min(cardinality,3)+1):
                for removed in itertools.combinations(selected,k):
                    remainder = mask ^ sum(1 << i for i in removed)
                    for inserted in buckets[k][key(removed)]:
                        if inserted & remainder:
                            continue
                        target = remainder | inserted
                        need(target in index, 'Legal local replacement remains in cardinality fiber')
                        need((target^mask).bit_count() <= 6, 'Constructed local event changes at most six slots')
                        union(state_index,index[target])
                        replacements += 1
        groups = defaultdict(list)
        for state_index,selected in enumerate(subsets):
            groups[key(selected)].append(state_index)
        distribution = Counter(len({find(i) for i in state_ids}) for state_ids in groups.values())
        fiber_reports.append({'selected':cardinality,'states':len(subsets),'fibers':len(groups),
                              'component_distribution':dict(sorted(distribution.items())),
                              'legal_local_replacements_checked':replacements})
        if cardinality == 4:
            state_ids = groups[(7,15)]
            hostile = [masks[i] for i in state_ids]
            need(len(state_ids) == 2 and len({find(i) for i in state_ids}) == 2, 'Entire 28-source hostile fiber consists of two disconnected states')
    need(fiber_reports[0]['component_distribution'] == {1:351}, 'Cardinality-two independent components')
    need(fiber_reports[1]['component_distribution'] == {1:1646}, 'Cardinality-three independent components')
    need(fiber_reports[2]['component_distribution'] == {1:1851,2:175,3:21,4:1}, 'Cardinality-four independent full component distribution')
    need(sorted(hostile) == [32901,136318978] and (hostile[0]^hostile[1]).bit_count() == 8, 'Eight-slot hostile is exact')
    hostile_sources = [[sources[i] for i in range(28) if mask >> i & 1] for mask in sorted(hostile)]

    # Replay the producer separately, without importing its implementation.
    replay_path = args.out.parent/('producer.optimized.json' if sys.flags.optimize else 'producer.normal.json')
    command = [sys.executable,'-B']+(['-O'] if sys.flags.optimize else [])
    command += [str(producer/'check.py'),'--repo',str(repo),'--supply',str(supply),'--out',str(replay_path)]
    completed = subprocess.run(command,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    replay = json.loads(completed.stdout)
    need(replay['status'] == 'PASS' and replay['checks'] == 10184, 'Producer replay passes its exact checks')
    need(replay_path.read_text() == completed.stdout, 'Producer output file agrees with stdout')
    need(replay['binary_hostile']['selected_sources'] == hostile_sources, 'Independent hostile agrees with producer')

    result = {'status':'PASS','checks':CHECKS,'independent_support_census':support_counts,
              'four_move_rank':len(pivots),'basis_diagonal':diagonal,'integer_basis_index':abs(determinant),
              'zero_margin_charge_image':len(image),'exact_lattice_equality_justification':'basis lattice contained in move lattice contained in charge kernel; outer indices both 2048',
              'independent_local_replacement_fibers':fiber_reports,
              'hostile_selected_sources':hostile_sources,'producer_checks':replay['checks'],
              'pins':{'producer_checker':digest(producer/'check.py'),'integer_witness':digest(producer/'integer_basis_witness.json'),
                      'supply':digest(supply),'original_electrical_receipt':digest(original)},
              'scope':'Independent arithmetic and full binary-fiber audit; native pin-domain/resource authentication is a separate retained prerequisite'}
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')


if __name__ == '__main__':
    main()
