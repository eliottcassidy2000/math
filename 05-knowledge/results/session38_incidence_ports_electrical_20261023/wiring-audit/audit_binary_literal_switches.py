"""Verify every binary quotient's literal reflected switch against our own ledger.

Only the witness matching is taken from the producer. Class assignment, native
port counts, and merged edge counts come from the independent FIRST C++ replay.
"""
import argparse
import collections
import csv
import hashlib
import itertools
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
FULL = (1 << 21) - 1
TILES = [(u, v) for u in range(6) for v in range(u + 2, 8)]
INDEX = {pair: i for i, pair in enumerate(TILES)}
REFLECTION = [INDEX[7 - v, 7 - u] for u, v in TILES]
CHECKS = collections.Counter()


def need(test, label):
    CHECKS[label] += 1
    if not test:
        raise ValueError(label)


def clean(counter):
    return {key: value for key, value in counter.items() if value}


def reflect(mask):
    return sum(((mask >> i) & 1) << j for i, j in enumerate(REFLECTION))


def adjacency(mask):
    out = [1 << (i - 1) if i else 0 for i in range(8)]
    for i, (u, v) in enumerate(TILES):
        if mask >> i & 1:
            out[u] |= 1 << v
        else:
            out[v] |= 1 << u
    return out


def main():
    parser = argparse.ArgumentParser()
    default = HERE.parent / 'binary-counterfeits' / 'receipt.json'
    if not default.exists():
        default = HERE.parent / 'session38-binary-counterfeits' / 'receipt.json'
    parser.add_argument('--witness', type=pathlib.Path, default=default)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    witness = json.loads(args.witness.read_text())
    partitions = {'0' + format(i, '08b') for i in range(1, 256)}
    records = witness['results']
    need(len(records) == 255 and {r['partition'] for r in records} == partitions, 'all 255 partitions exactly once')
    classes = [{k: int(v) for k, v in r.items()} for r in csv.DictReader((HERE / 'independent_classes.tsv').open(), delimiter='\t')]
    need(len(classes) == 6880 and [r['cid'] for r in classes] == list(range(6880)), 'independent class ordering')
    raw = (HERE / 'independent_class_of_mask.u16le').read_bytes()
    need(len(raw) == 2 * (1 << 21), 'independent full mask map')

    def cid(mask):
        need(0 <= mask <= FULL, 'literal mask range')
        return raw[2 * mask] + 256 * raw[2 * mask + 1]

    data = []
    for cls in classes:
        A = [[0] * 8 for _ in range(8)]
        for k, (u, v) in enumerate(itertools.combinations(range(8), 2)):
            A[u][v] = (cls['canon'] >> (27 - k)) & 1
            A[v][u] = 1 - A[u][v]
        c = [0] * 8
        e = [[0] * 8 for _ in range(8)]
        for u, v, w in itertools.combinations(range(8), 3):
            if A[u][v] == A[v][w] == A[w][u]:
                for x in (u, v, w):
                    c[x] += 1
                for x, y in ((u, v), (u, w), (v, w)):
                    e[x][y] += 1
                    e[y][x] += 1
        symbols = [[3 * (c[v] % 3) + e[u][v] % 3 for v in range(8) if A[u][v]] for u in range(8)]
        data.append((sum(c) // 3, symbols))

    def sig(i, partition):
        return tuple(sorted((len(row), sum(partition[z] == '1' for z in row)) for row in data[i][1]))

    def mirror(cell):
        i, j, a, b = cell
        return classes[i]['converse_cid'], classes[j]['converse_cid'], 7 - b, 7 - a

    requested = {tuple(row[:4]) for r in records for row in r['ordinary_cell_delta']}
    baseline = {k: 0 for k in requested}
    for row in csv.DictReader((HERE / 'independent_ports.tsv').open(), delimiter='\t'):
        if int(row['blue']):
            continue
        k = tuple(int(row[x]) for x in ('source_cid', 'target_cid', 'start_score', 'end_score'))
        if k in baseline:
            baseline[k] = int(row['count'])
    merged_baseline = {}
    for row in csv.DictReader((HERE / 'independent_directed_edges.tsv').open(), delimiter='\t'):
        if not int(row['blue']):
            k = int(row['source_node']), int(row['target_node'])
            value = int(row['literal_count'])
            need(value % 2 == 0, 'all native merged black multiplicities even')
            merged_baseline[k] = value
    summaries = []
    all_switches = set()
    for r in records:
        part = r['partition']
        matching_rows = r['literal_matching']
        matching = {z['closed_source_mask']: z['counterfeit_open_target_mask'] for z in matching_rows}
        need(len(matching) == len(matching_rows) == 4, 'four distinct actual closed ports')
        need(len(set(matching.values())) == 4, 'four distinct actual open ports')
        need(set(matching.values()) == {m ^ FULL for m in matching}, 'exact open inventory preserved bijectively')
        delta = collections.Counter()
        for z in matching_rows:
            m, t = z['closed_source_mask'], z['counterfeit_open_target_mask']
            old = m ^ FULL
            need(z['native_open_target_mask'] == old, 'native complement partner')
            need(m >> INDEX[0, 7] & 1 and not (t >> INDEX[0, 7] & 1), 'closed to open ports')
            need(t != old, 'literal matching genuinely changes')
            rm, rt = reflect(m), reflect(t)
            need(m != rm and t != rt, 'actual black colors unchanged')
            need(matching.get(rm) == rt, 'literal bijection commutes with reflection')
            i, j, k = cid(m), cid(old), cid(t)
            am, at = adjacency(m), adjacency(t)
            a, b = am[7].bit_count(), am[0].bit_count()
            need((at[7].bit_count(), at[0].bit_count()) == (8 - a, 6 - b), 'literal target endpoint compatibility')
            need(sig(j, part) == sig(k, part), 'actual target quotient signature preserved')
            need(sig(classes[j]['converse_cid'], part) == sig(classes[k]['converse_cid'], part), 'target signature also preserved at mirror')
            need(cid(rm) == classes[i]['converse_cid'] and cid(rt) == classes[k]['converse_cid'], 'literal reflection class labels')
            need(data[k][0] - data[i][0] == a - b - 1, 'new pointwise triangle law')
            delta[i, j, a, b] -= 1
            delta[i, k, a, b] += 1
        delta = clean(delta)
        expected = {tuple(row[:4]): row[5] for row in r['ordinary_cell_delta']}
        need(len(expected) == len(r['ordinary_cell_delta']) and delta == expected, 'literal matching induces complete reported cell delta')
        row_change = collections.Counter()
        target_change = collections.Counter()
        directed = collections.Counter()
        undirected = collections.Counter()
        degree = collections.Counter()
        current = collections.Counter()
        for row in r['ordinary_cell_delta']:
            i, j, a, b, before, change, after = row
            cell = i, j, a, b
            need(baseline[cell] == before and before + change == after >= 0, 'native independent before counts and nonnegative after counts')
            need(delta.get(mirror(cell), 0) == change, 'full cell reflection equality')
            need(cell != mirror(cell) or after % 2 == 0, 'reflection-fixed ordinary black cells even')
            need(data[j][0] - data[i][0] == a - b - 1, 'every changed cell triangle law')
            row_change[i, a, b, sig(j, part)] += change
            target_change[j, 8 - a, 6 - b] += change
            u, v = classes[i]['node_id'], classes[j]['node_id']
            directed[u, v] += change
            undirected[min(u, v), max(u, v)] += change
            degree[u] += change
            degree[v] += change
            current[u] += change * (data[j][0] - data[i][0])
            current[v] += change * (data[i][0] - data[j][0])
        need(not clean(row_change), 'all source endpoint color signature margins unchanged')
        need(not clean(target_change), 'all target endpoint color capacities unchanged')
        need(not clean(degree), 'all merged multigraph degrees unchanged')
        need(not clean(current), 'all merged cyclic-triangle currents unchanged')
        directed, undirected = clean(directed), clean(undirected)
        need(directed and undirected, 'directed and undirected projections both change')
        need([[*k, v] for k, v in sorted(directed.items())] == r['directed_merged_delta'], 'exact directed projection')
        need([[*k, v] for k, v in sorted(undirected.items())] == r['undirected_merged_delta'], 'exact undirected projection')
        need(all(v % 2 == 0 and merged_baseline.get(k, 0) + v >= 0 for k, v in directed.items()), 'all changed black twins even and nonnegative')
        need(all(u != v for u, v in undirected), 'all loop counts untouched')
        all_switches.add(tuple(sorted(matching.items())))
        summaries.append({'partition': part, 'changed_ordinary_cells': len(delta), 'changed_directed_pairs': len(directed), 'changed_undirected_pairs': len(undirected), 'status': 'PASS'})
    need(len(all_switches) == 61, 'distinct literal switch count')
    ternary = json.loads((HERE / 'independent_ternary_quotient_receipt.json').read_text())
    need(ternary['selected']['cycle_rank'] == 0 and ternary['selected']['unpeeled_classes'] == 0, 'independent ternary upper certificate')
    files = ['independent_class_of_mask.u16le', 'independent_classes.tsv', 'independent_ports.tsv', 'independent_directed_edges.tsv', 'independent_ternary_quotient_receipt.json']
    result = {
        'status': 'PASS',
        'scope': 'Exact order-eight minimum of three colors among all set-partition quotients of the nine E33 symbols for the specified reflected endpoint/color/signature port-matching model. All 255 binary quotients have actual literal port-bijection counterfeits with altered merged wiring. No claim that these alternate matchings are canonical-complement tournament realizations.',
        'binary_partitions': len(records),
        'verified_literal_matchings': len(records),
        'distinct_literal_switches': len(all_switches),
        'matched_source_ports_counted_with_repetitions': 4 * len(records),
        'minimum_colors_in_stated_model': 3,
        'all_blue_entries_and_all_loops_unchanged': True,
        'unary_lower_bound': 'Any binary witness remains valid after the final two colors are merged.',
        'independent_baseline_sha256': {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in files},
        'producer_witness_sha256': hashlib.sha256(args.witness.read_bytes()).hexdigest(),
        'checks': dict(sorted(CHECKS.items())),
        'results': summaries,
    }
    dest = HERE / 'independent_binary_literal_receipt.json'
    if args.write:
        dest.write_text(json.dumps(result, indent=2) + '\n')
    else:
        need(json.loads(dest.read_text()) == result, 'frozen receipt mismatch')
    print(json.dumps({k: v for k, v in result.items() if k != 'results'}, indent=2))


if __name__ == '__main__':
    main()
