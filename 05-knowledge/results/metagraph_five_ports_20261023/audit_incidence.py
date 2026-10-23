"""Independent canonical-gauge audit of the residual rooted-five records.

Uses the separately audited little-endian flag canonicalizer, not probe.py.
The producer uses the opposite bit gauge; compare types through an explicit
bijection before comparing every rooted incidence coordinate.
"""
import importlib.util
import json
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'flag_tables', HERE.parent/'metagraph_flag_transport_20261023/check.py')
flags = importlib.util.module_from_spec(spec)
spec.loader.exec_module(flags)
checks = 0


def need(ok, why):
    global checks
    checks += 1
    if not ok:
        raise ValueError(why)


def main():
    data = json.loads((HERE/'probe.json').read_text())
    _, _, transforms = flags.tables(5)
    byroot = [[t for t in transforms if t[0][0] == r] for r in range(5)]
    canonical = [[min(tm[m] for _, _, tm in byroot[r]) for r in range(5)]
                 for m in range(1024)]
    reverse = lambda x: int(f'{x:010b}'[::-1], 2)
    oldtypes = data['rooted_type_representatives']
    newtypes = [canonical[reverse(c)][0] for c in oldtypes]
    need(len(set(newtypes)) == 48, 'opposite-gauge type bijection')
    inverse = {c: i for i, c in enumerate(newtypes)}
    lut = [[inverse[c] for c in row] for row in canonical]
    all_records = {}
    colliding_pairs = []
    for n in (8, 9):
        for rec in data['orders'][str(n)]['records']:
            code = rec['code']
            edges = list(combinations(range(n), 2))
            a = [[0]*n for _ in range(n)]
            for bit, (i, j) in enumerate(reversed(edges)):
                a[i][j] = (code >> bit) & 1
                a[j][i] = 1-a[i][j]
            counts = [[0]*48 for _ in range(n)]
            for subset in combinations(range(n), 5):
                local = sum(a[i][j] << bit
                            for bit, (i, j) in enumerate(combinations(subset, 2)))
                for r, v in enumerate(subset):
                    counts[v][lut[local][r]] += 1
            for v in range(n):
                need(counts[v] == rec['rooted5'][v], 'literal rooted incidence')
            all_records[n, code] = counts
        menu = json.loads((HERE/f'menu_n{n}.json').read_text())
        base = next(x for x in menu['menus'] if x['name'] == 'v4' and x['round'] == 1)
        need(base['all_collisions_recorded'], 'complete residual buckets')
        colliding_pairs.extend(((n, a), (n, b))
                               for g in base['collisions'] for a, b in combinations(g, 2))
    def histogram(key, cols):
        return sorted(tuple(row[c] for c in cols) for row in all_records[key])
    def succeeds(cols):
        return all(histogram(a, cols) != histogram(b, cols) for a, b in colliding_pairs)
    singles = [k for k in range(48) if succeeds((k,))]
    pairs = [list(p) for p in combinations(range(48), 2) if succeeds(p)]
    need(not singles and len(pairs) == 14, 'coordinate cardinality minimum')
    need(pairs == data['coordinate_minimum']['successful_coordinate_pairs'], 'all good pairs')
    conv = [inverse[canonical[reverse(c) ^ 1023][0]] for c in oldtypes]
    need(all(conv[conv[c]] == c for c in range(48)), 'converse involution')
    need(oldtypes[conv[6]] == 457 and conv[24] == 24, 'marked converse triple')
    need(not any({conv[c] for c in p} == set(p) for p in pairs), 'no closed two-coordinate repair')
    triple = {6, 24, conv[6]}
    need(len(triple) == 3 and {conv[c] for c in triple} == triple, 'closed triple')
    for records in all_records.values():
        for row in records:
            plus, minus = row[6]+row[conv[6]], row[6]-row[conv[6]]
            need((plus+minus)//2 == row[6] and (plus-minus)//2 == row[conv[6]]
                 and (plus-minus) % 2 == 0 and abs(minus) <= plus, 'lossless parity coordinates')
    print(json.dumps({'status': 'PASS', 'checks': checks,
                      'classes': len(all_records), 'collision_pairs': len(colliding_pairs),
                      'successful_pairs': len(pairs), 'raw_minimum': 2,
                      'converse_closed_raw_minimum': 3}))


if __name__ == '__main__':
    main()
