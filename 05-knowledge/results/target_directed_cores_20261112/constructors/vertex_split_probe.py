"""Protected common induced graphs between consecutive inherited ready fibres."""
import gzip
import json
from itertools import combinations
from pathlib import Path
import subprocess
import networkx as nx

HERE = Path(__file__).resolve().parent
BANK = '05-knowledge/results/multicore_metagraph_20261113/inflation_descent.json.gz'


def graph(rows):
    g = nx.DiGraph()
    n = len(rows)
    for v in range(n):
        g.add_node(v, port='root' if v == n - 1 else 'predecessor' if v == 0 else 'other')
        for w in range(n):
            if rows[v] >> w & 1:
                g.add_edge(v, w)
    return g


def matches(left, right, removed_old=1):
    m, n = len(left), len(right)
    removed_new = n - m + removed_old
    if removed_new < 0:
        return []
    source, target = graph(left), graph(right)
    found = []
    for old in combinations(range(1, m - 1), removed_old):
        subsource = source.subgraph([v for v in source if v not in old])
        degree_type = sorted((subsource.nodes[v]['port'], subsource.in_degree(v)) for v in subsource)
        for new in combinations(range(1, n - 1), removed_new):
            subtarget = target.subgraph([v for v in target if v not in new])
            if degree_type != sorted((subtarget.nodes[v]['port'], subtarget.in_degree(v)) for v in subtarget):
                continue
            matcher = nx.algorithms.isomorphism.DiGraphMatcher(
                subsource, subtarget, node_match=lambda a, b: a['port'] == b['port'])
            if matcher.is_isomorphic():
                found.append({'delete_old': old, 'insert_new': new,
                              'retained_map': sorted(matcher.mapping.items())})
    return found


def main():
    bank = json.loads(gzip.decompress(subprocess.check_output(['git', 'show', 'HEAD:' + BANK])))
    old = [r for r in bank['records'] if r['H'] == 881]
    new = [r for r in bank['records'] if r['H'] == 1691]
    records = []
    for a in old:
        for b in new:
            hit = matches(a['rows'], b['rows'])
            records.append({'old_order': len(a['rows']), 'old_rank': a['parent_rank'],
                            'new_order': len(b['rows']), 'new_rank': b['parent_rank'],
                            'old_rows': a['rows'], 'new_rows': b['rows'], 'matches': hit})
    out = {'scope': 'Delete exactly one old nonport vertex; insert enough nonports to reach target order; preserve both physical ports and every retained arc',
           'source_H': 881, 'target_H': 1691, 'records': records}
    (HERE / 'vertex_split_receipt.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'tested_source_pairs': len(records),
                      'matched_pairs': [(r['old_rank'], r['new_rank'], len(r['matches']))
                                        for r in records if r['matches']]}))


if __name__ == '__main__':
    main()
