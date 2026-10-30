"""Independent VF2 and maximal-clique audit of the native finite certificates."""
from pathlib import Path
from itertools import combinations
from collections import Counter
import argparse
import json
import networkx as nx

HERE = Path(__file__).resolve().parent


def require(value, why):
    if not value:
        raise RuntimeError(why)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    reference = json.loads((HERE/'verification.json').read_text())
    report = {'networkx_version': nx.__version__, 'banks': {}, 'witnesses': {}}
    for name, s, t in [('r35_13.g6', 3, 5), ('r36_17.g6', 3, 6), ('r44_17.g6', 4, 4)]:
        graphs = nx.read_graph6(HERE/'inputs'/name)
        if not isinstance(graphs, list):
            graphs = [graphs]
        entries = []
        for i, G in enumerate(graphs):
            omega = max(map(len, nx.find_cliques(G)))
            alpha = max(map(len, nx.find_cliques(nx.complement(G))))
            require(omega < s and alpha < t, 'Independent clique audit')
            aut = list(nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter())
            orbits = {tuple(sorted(p[u] for p in aut)) for u in G}
            # Remove multiplicities before measuring each orbit.
            orbit_sizes = sorted(len(set(O)) for O in orbits)
            expected = reference['banks'][name]['entries'][i]
            require(len(aut) == expected['aut_order'] and orbit_sizes == expected['root_orbits'], 'Independent rooted inventory')
            pairs = {(u, v) for u in G for v in G}
            pair_count = 0
            while pairs:
                u, v = min(pairs)
                pairs -= {(p[u], p[v]) for p in aut}
                pair_count += 1
            require(pair_count == expected['ordered_pair_orbits'], 'Pair orbits by direct enumeration vs Burnside')
            entries.append({'omega': omega, 'alpha': alpha, 'aut_order': len(aut), 'root_orbits': orbit_sizes,
                            'ordered_pair_orbits': pair_count})
        report['banks'][name] = entries
        if name == 'r35_13.g6':
            H = graphs[0]
            deletions, additions = Counter(), Counter()
            for u, v in combinations(H, 2):
                edited = H.copy()
                if H.has_edge(u, v):
                    edited.remove_edge(u, v)
                    sizes = [len(S) for S in nx.find_cliques(nx.complement(edited))]
                    require(max(sizes) == 5, 'Edge-deletion independence number')
                    deletions[sizes.count(5)] += 1
                else:
                    edited.add_edge(u, v)
                    additions[sum(nx.triangles(edited).values())//3] += 1
            require(deletions == {3: 26} and additions == {1: 26, 2: 26}, 'Independent saturation counts')
            report['H13_saturation'] = {'deletion': dict(deletions), 'addition': dict(additions)}
    for name in ['balanced_ports', 'unequal_ports']:
        X = json.loads((HERE/(name+'.json')).read_text())['cross_colors']
        colors = [nx.Graph() for _ in range(3)]
        for G in colors:
            G.add_nodes_from(range(27))
        for i in range(13):
            colors[0].add_edge(0, i+1)
            colors[1].add_edge(0, i+14)
            for j in range(13):
                colors[X[i][j]].add_edge(i+1, j+14)
        for i, j in combinations(range(13), 2):
            h = (i-j) % 13 in {1, 5, 8, 12}
            colors[1 if h else 2].add_edge(i+1, j+1)
            colors[0 if h else 2].add_edge(i+14, j+14)
        maxima = [max(map(len, nx.find_cliques(G))) for G in colors]
        require(maxima == [2, 2, 4], 'Independent 27-vertex witness check')
        require(sum(G.number_of_edges() for G in colors) == 351, 'Complete edge partition')
        report['witnesses'][name] = {'max_clique_sizes': maxima, 'edges': [G.number_of_edges() for G in colors]}
    text = json.dumps(report, indent=2, sort_keys=True)+'\n'
    target = HERE/'independent_audit.json'
    if args.write:
        target.write_text(text, encoding='utf-8', newline='\n')
    else:
        require(json.loads(target.read_text()) == json.loads(text), 'Independent audit reproduces')
    print(text, end='')


if __name__ == '__main__':
    main()
