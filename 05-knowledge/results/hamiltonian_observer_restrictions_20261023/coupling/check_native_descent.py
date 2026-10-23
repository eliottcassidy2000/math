"""Constructive all-size native-component descent: exhaustive n=3..6 controls."""
from collections import Counter
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('restricted_coupling', HERE/'check.py')
q = importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)
CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def run(n):
    m = n-1
    pairs, idx, perms, actions = q.old.setup(m)
    npairs, nidx, _, nactions = q.old.setup(n)
    size = 1 << len(pairs)
    gc, _, _ = q.old.orbit_table(size, actions, False)
    tc, _, _ = q.old.orbit_table(1 << len(npairs), nactions, True)
    native_positions = [nidx[i, j+1] for i, j in pairs]

    def reorder_graph(g, order):
        new = sum(1 << k for k, (i, j) in enumerate(pairs)
                  if g >> idx[tuple(sorted((order[i], order[j])))] & 1)
        need(gc[new] == gc[g], 'retained source graph isomorphism')
        return new

    def path_change(g, path):
        t = q.old.relabel(g, native_positions)
        def arc(a, b):
            bit = t >> nidx[tuple(sorted((a, b)))] & 1
            return bool(bit) if a > b else not bit
        need(sorted(path) == list(range(n)), 'actual owner permutation')
        need(all(arc(a, b) for a, b in zip(path, path[1:])), 'actual Hamiltonian junctions')
        tnew = sum(1 << k for k, (i, j) in enumerate(npairs) if not arc(path[i], path[j]))
        need(tc[tnew] == tc[t], 'retained native target tournament')
        new = sum(1 << k for k, (i, j) in enumerate(pairs) if not arc(path[i], path[j+1]))
        need(q.old.relabel(new, native_positions) == tnew, 'lossless recompression')
        return new

    cases = Counter()
    max_steps = 0
    max_matching_stalls = 0
    for start in range(1, size):
        g = start
        steps = stalls = 0
        while g.bit_count() > 1:
            old_weight = g.bit_count()
            adj = [{j for j in range(m) if j != i and g >> idx[tuple(sorted((i,j)))] & 1}
                   for i in range(m)]
            mid = next((u for u in range(m) if 2 <= len(adj[u]) <= m-2), None)
            if mid is not None:
                u = mid
                a, b = sorted(adj[u])[:2]
                w = next(v for v in range(m) if v != u and v not in adj[u])
                first = [a,u,b,w]
                order = first + [v for v in range(m) if v not in first]
                g = reorder_graph(g, order)
                triangle = bool(g >> idx[0,2] & 1)
                prefix = [2,3,0,1] if triangle else [2,0,3,1]
                g = path_change(g, prefix + list(range(4,n)))
                need(g.bit_count() == old_weight-(2 if triangle else 1), 'guarded strong-four drop')
                case = 'guarded_triangle' if triangle else 'guarded_path'
            elif max(map(len, adj)) <= 1:
                edges = [(a,b) for a,b in pairs if g >> idx[a,b] & 1]
                first = list(edges[0]+edges[1])
                order = first + [v for v in range(m) if v not in first]
                g = reorder_graph(g, order)
                g = path_change(g, [1,2,0]+list(range(3,n)))
                need(g.bit_count() == old_weight, 'matching stall preserves weight')
                need(g >> idx[0,1] & 1 and g >> idx[1,3] & 1, 'matching creates adjacent edges')
                need(2 <= sum(g >> idx[tuple(sorted((1,j)))] & 1 for j in range(m) if j != 1) <= m-2,
                     'next matching state has an intermediate-degree vertex')
                stalls += 1
                case = 'matching_prepare'
            elif g == size-1:
                path = [v for a in range(n-2, -1, -2) for v in (a,a+1)]
                if n % 2:
                    path.append(0)
                g = path_change(g, path)
                need(g.bit_count() == (n-1)//2, 'complete descending-pair count')
                case = 'complete_pairs'
            else:
                centers = [u for u in range(m) if len(adj[u]) == m-1]
                need(len(centers) == 1 and all(len(adj[v]) == 1 for v in range(m) if v != centers[0]),
                     'only remaining graph is spanning star')
                u = centers[0]
                g = reorder_graph(g, [u]+[v for v in range(m) if v != u])
                g = path_change(g, list(range(1,n))+[0])
                need(g.bit_count() == 1, 'spanning star collapses to one edge')
                case = 'star_rotate'
            cases[case] += 1
            steps += 1
            need(g != 0, 'nontransitivity retained')
            need(steps <= 2*(start.bit_count()-1), 'constructive radius bound during descent')
        need(g.bit_count() == 1, 'one-edge terminal')
        max_steps = max(max_steps, steps)
        max_matching_stalls = max(max_matching_stalls, stalls)
    return {'n': n, 'nonempty_labelled_sources': size-1, 'descent_steps_by_case': dict(sorted(cases.items())),
            'maximum_steps': max_steps, 'maximum_matching_preparations': max_matching_stalls}


if __name__ == '__main__':
    reports = [run(n) for n in range(3,7)]
    result = {'status': 'FINITE-EXACT', 'checks': CHECKS, 'reports': reports,
              'scope': 'all labelled compressed graphs through native order six; every retained actual Hamiltonian junction and isomorphism independently checked'}
    receipt = HERE/'native_descent_receipt.json'
    if '--write' in sys.argv:
        receipt.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    else:
        need(result == json.loads(receipt.read_text(encoding='utf-8')), 'saved receipt')
    print(json.dumps({'status': 'PASS', 'checks': CHECKS, 'reports': reports}, indent=2))
