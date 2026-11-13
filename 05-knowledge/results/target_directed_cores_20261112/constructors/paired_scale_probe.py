"""Targeted synchronized two-module probes in inherited complete ready fibres."""
import gzip
import json
from itertools import combinations
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
BANK = '05-knowledge/results/multicore_metagraph_20261113/inflation_descent.json.gz'


def inflate(rows, sizes):
    groups = []
    offset = 0
    for size in sizes:
        groups.append(tuple(range(offset, offset + size)))
        offset += size
    out = [0] * offset
    for u, group in enumerate(groups):
        for i, vertex in enumerate(group):
            out[vertex] |= sum(1 << w for w in group[i + 1:])
            for v, other in enumerate(groups):
                if rows[u] >> v & 1:
                    out[vertex] |= sum(1 << w for w in other)
    return tuple(out), groups


def evaluate(cases):
    payload = ''.join(f'{len(rows)} {root} ' + ' '.join(map(str, rows)) + '\n'
                      for rows, root in cases)
    answer = subprocess.check_output(['/tmp/target_directed_count_profiles'],
                                     input=payload, text=True)
    return [tuple(map(int, line.split())) for line in answer.splitlines()]


def main():
    bank = json.loads(gzip.decompress(subprocess.check_output(['git', 'show', 'HEAD:' + BANK])))
    stride = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    if stride not in (2, 3):
        raise ValueError('Only the declared stride2/3 probes are implemented')
    source_values = (1691, 3311) if stride == 2 else (3311, 6551, 13031)
    sources = [r for r in bank['records'] if r['H'] in source_values]
    cases, metadata = [], []
    for source in sources:
        rows = source['rows']; n = len(rows)
        menu = [((v, stride + 1),) for v in range(1, n - 1)]
        if stride == 2:
            menu += [((u, 2), (v, 2)) for u, v in combinations(range(1, n - 1), 2)]
        else:
            menu += [((u, 3), (v, 2)) for u in range(1, n - 1)
                     for v in range(1, n - 1) if u != v]
            menu += [tuple((v, 2) for v in vertices)
                     for vertices in combinations(range(1, n - 1), 3)]
        for substitution in menu:
            sizes = [1] * n
            for v, size in substitution:
                sizes[v] = size
            child, groups = inflate(rows, sizes)
            cases.append((child, -1))
            metadata.append({'source_H': source['H'], 'source_rows': rows,
                             'parent_order': source['parent_order'],
                             'parent_rank': source['parent_rank'],
                             'substitution_at_m2': substitution,
                             'target': 2**stride * source['H'] - (2**stride - 1) * 71})
    hits = []
    near = []
    for data, (h, q, s) in zip(metadata, evaluate(cases)):
        if h == data['target']:
            hits.append({**data, 'H': h, 'readiness': 'pending literal root count'})
        elif abs(h - data['target']) <= 20:
            near.append({**data, 'H': h, 'q': q, 's': s})
    out = {'scope': 'Synchronized transitive-module increments totaling the stated stride; unchanged exterior root/predecessor',
           'stride': stride,
           'inherited_bank': BANK, 'source_count': len(sources), 'tested_graphs': len(cases),
           'hits': hits, 'near_within20': near}
    receipt = 'paired_scale_receipt.json' if stride == 2 else 'triple_scale_receipt.json'
    (HERE / receipt).write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out))


if __name__ == '__main__':
    main()
