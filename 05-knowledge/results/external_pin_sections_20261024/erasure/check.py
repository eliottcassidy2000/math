"""Exact grade coverage after external pins freeze up to three receiver slots.

The subtraction recurrence removes whole source slots; it never treats a
second kernel record or a charge copy as a second occupied resource.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
COUNT = 0


def need(ok, message):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(message)


def add(a, b):
    return ((((a >> 4) + (b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def sub(a, b):
    return ((((a >> 4) - (b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def census(labels):
    d = [[0] * 2048 for _ in range(len(labels) + 1)]
    d[0][0] = 1
    for i, label in enumerate(labels):
        for k in range(i + 1, 0, -1):
            previous, target = d[k - 1], d[k]
            for q, n in enumerate(previous):
                if n:
                    target[add(q, label)] += n
    return d


def remove(d, permutation):
    out = [d[0][:]]
    for row in d[1:-1]:
        previous = out[-1]
        out.append([n - previous[permutation[q]] for q, n in enumerate(row)])
    return out


def main():
    old_path = RESULTS / 'native_packet_decomposition_20261023/native/domains.json'
    old = json.loads(old_path.read_text())['sources']
    rows = []
    for r in old:
        paths = [s['path'] for s in r['exact_all_shared_clocks_target_domain']]
        if sum(a != b for a, b in zip(*paths)) == 2:
            code = r['code']
            rows.append((r['source'], (code[0] << 4) | sum(v << j for j, v in enumerate(code[1:]))))
    rows.sort()
    ids, labels = map(list, zip(*rows))
    need(len(ids) == 20, 'twenty distinct native two-clock core source slots')
    repair_row = next(r for r in old if r['source'] == 1841)
    repair_code = repair_row['code']
    repair_label = (repair_code[0] << 4) | sum(v << j for j,v in enumerate(repair_code[1:]))
    repair_paths = [r['path'] for r in repair_row['exact_all_shared_clocks_target_domain']]
    need(sum(a != b for a,b in zip(*repair_paths)) == 3 and 1841 not in ids,
         'receiver extension is an existing distinct three-clock core source')
    base = census(labels)
    need([sum(n > 0 for n in row) for row in base] ==
         [1,20,187,910,1899] + [2048]*11 + [1899,910,187,20,1],
         'graded bank rebuilt in the actual Y1,...,Y4 paired coordinates')
    permutations = [[sub(q, a) for q in range(2048)] for a in labels]
    levels = []
    witnesses = []
    repair_rows = []
    for e in range(4):
        profiles = Counter()
        universal = set(range(21-e))
        bad_near_center = []
        nonmonotone = []
        for erased in combinations(range(20), e):
            d = base
            for i in erased:
                d = remove(d, permutations[i])
            coverage = tuple(sum(n > 0 for n in row) for row in d)
            good = {k for k, n in enumerate(coverage) if n == 2048}
            universal &= good
            profiles[coverage] += 1
            need(all(min(row) >= 0 for row in d), 'deletion recurrence retains nonnegative counts')
            need(sum(map(sum, d)) == 1 << (20-e), 'each remaining binary assignment counted exactly once')
            need(coverage == coverage[::-1], 'complementation preserves each deletion coverage profile')
            if not good:
                bad_near_center.append([ids[i] for i in erased])
            if e == 3:
                repaired = [sum(d[k][q] + d[k-1][sub(q, repair_label)] > 0
                                for q in range(2048)) for k in (8,9,10)]
                need(repaired == [2048]*3, 'existing core1841 repairs every three-slot deletion at grades8..10')
                if coverage[8] < 2048 or coverage[9] < 2048:
                    holes = {str(k): [q for q,n in enumerate(d[k]) if not n] for k in (7,8,9,10)}
                    repair_rows.append(dict(erased=[ids[i] for i in erased], complete_grades=sorted(good),
                                            original_holes=holes, repaired_coverage=repaired))
                    if good:
                        nonmonotone.append([ids[i] for i in erased])
                    fresh = census([a for i,a in enumerate(labels) if i not in erased] + [repair_label])
                    need(all(all(v > 0 for v in fresh[k]) for k in (8,9,10)),
                         'independent fresh eighteen-slot census validates exceptional repair')
            if erased == tuple(range(e)):
                fresh = census([a for i, a in enumerate(labels) if i not in erased])
                need(d == fresh, 'independent fresh subset DP matches slot subtraction')
                witnesses.append(dict(erased=[ids[i] for i in erased], coverage=coverage))
        levels.append(dict(erasure_count=e, total_pin_sets=sum(profiles.values()),
                           complete_grades_for_every_pin_set=sorted(universal),
                           no_complete_grade_pin_sets=bad_near_center,
                           complete_noncentral_but_incomplete_central_pin_sets=nonmonotone,
                           profiles=[dict(coverage=p, pin_sets=n) for p,n in sorted(profiles.items())]))
        print('erased', e, 'universal grades', sorted(universal), 'no complete grade', len(bad_near_center), flush=True)
    out = dict(status='FINITE-EXACT whole-slot erasure coverage, not arbitrary relational pins',
               sources=ids, labels=labels, levels=levels, independent_deletion_controls=witnesses,
               adaptive_repair=dict(source=1841, label=repair_label, cost=3,
                                    all_triple_pin_sets=1140, complete_remaining_grades=[8,9,10],
                                    exceptional_original_centers=repair_rows),
               checks=COUNT, input_sha256=sha256(old_path.read_bytes().replace(b'\r\n',b'\n')).hexdigest())
    out = json.loads(json.dumps(out))
    target = HERE/'receipt.json'
    if '--write' in sys.argv:
        target.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    else:
        need(out == json.loads(target.read_text()), 'frozen full erasure receipt')
    print('PASS', COUNT)


if __name__ == '__main__':
    main()
