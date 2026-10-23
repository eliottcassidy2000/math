"""Exact finite incidence-flow controls and a rational mixture collision."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import json
import sys

from check import setup, orbit_table

checks = 0


def require(value, why):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(why)


def null_vector(rows):
    """A rational vector annihilating the supplied row matrix, if any."""
    width = len(rows)
    matrix = [[F(rows[j][i]) for j in range(width)] for i in range(len(rows[0]))]
    pivots = []
    rank = 0
    for c in range(width):
        pivot = next((r for r in range(rank, len(matrix)) if matrix[r][c]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = matrix[rank][c]
        matrix[rank] = [x / scale for x in matrix[rank]]
        for r in range(len(matrix)):
            if r != rank and matrix[r][c]:
                a = matrix[r][c]
                matrix[r] = [x - a*y for x, y in zip(matrix[r], matrix[rank])]
        pivots.append(c)
        rank += 1
    free = next((c for c in range(width) if c not in pivots), None)
    if free is None:
        return None
    vector = [F(0)] * width
    vector[free] = F(1)
    for r, c in reversed(list(enumerate(pivots))):
        vector[c] = -sum(matrix[r][j]*vector[j] for j in range(c + 1, width))
    multiple = lcm(*(x.denominator for x in vector))
    integer = [int(x * multiple) for x in vector]
    common = gcd(*integer)
    return [x // common for x in integer]


def transport(weights):
    r = list(map(sum, weights))
    c = [sum(row[j] for row in weights) for j in range(len(weights[0]))]
    mass = sum(r)
    forward = [[F(w, r[i]) for w in row] for i, row in enumerate(weights)]
    P = [[sum(F(weights[i][t]*weights[j][t], r[i]*c[t])
              for t in range(len(c))) for j in range(len(r))] for i in range(len(r))]
    stationary = [F(x, mass) for x in r]
    for i in range(len(r)):
        require(sum(P[i]) == 1 and P[i][i] > 0, 'stochastic positive diagonal')
        for j in range(len(r)):
            require(stationary[i]*P[i][j] == stationary[j]*P[j][i], 'detailed balance')
    reached = {0}
    while True:
        added = {j for i in reached for j in range(len(r)) if P[i][j]}
        if added <= reached:
            break
        reached |= added
    require(len(reached) == len(r), 'two-step support connected')
    # Independently compare the edge and conditional-variance quadratic forms.
    vectors = [[int(i == a) for i in range(len(r))] for a in range(len(r))]
    vectors += [[int(i == a or i == b) for i in range(len(r))]
                for a, b in combinations(range(len(r)), 2)]
    for f in vectors:
        edge = sum(stationary[i]*P[i][j]*(f[i]-f[j])**2
                   for i in range(len(r)) for j in range(len(r))) / 2
        variance = sum(sum(weights[i][t]*f[i]**2 for i in range(len(r)))
                       - F(sum(weights[i][t]*f[i] for i in range(len(r)))**2, c[t])
                       for t in range(len(c))) / mass
        require(edge == variance, 'electrical/conditional-variance identity')
    return forward, P, stationary


def main():
    inherited = json.loads(Path(__file__).with_name('receipt.json').read_text(encoding='utf-8'))
    census = []
    witness = None
    for n in range(1, 5):
        pairs, index, perms, actions = setup(n)
        size = 1 << len(pairs)
        gc, greps, ga = orbit_table(size, actions, False)
        tc, treps, ta = orbit_table(size, actions, True)
        gi, ti = {g: i for i, g in enumerate(greps)}, {t: i for i, t in enumerate(treps)}
        B = [[0] * len(treps) for _ in greps]
        for word in range(size):
            B[gi[gc[word]]][ti[tc[word]]] += 1
        require(B == inherited['incidence_small'][str(n)]['matrix'], 'literal words vs ordered incidence')
        for i, g in enumerate(greps):
            require(sum(B[i])*ga[g] == len(perms), 'graph row orbit size')
        for j, t in enumerate(treps):
            require(sum(row[j] for row in B)*ta[t] == len(perms), 'tournament column orbit size')
        transport(B)
        merged_g = sorted({min(g, gc[size-1-g]) for g in greps})
        merged_t = sorted({min(t, tc[size-1-t]) for t in treps})
        mg, mt = {g: i for i, g in enumerate(merged_g)}, {t: i for i, t in enumerate(merged_t)}
        W = [[0] * len(merged_t) for _ in merged_g]
        for word in range(size):
            g = min(gc[word], gc[size-1-word])
            t = min(tc[word], tc[size-1-word])
            W[mg[g]][mt[t]] += 1
        forward, P, stationary = transport(W)
        require(len({tuple(row) for row in forward}) == len(forward), 'distinct merged rows')
        census.append({'n': n, 'word_mass': size, 'graph_types': len(greps),
                       'tournament_types': len(treps), 'merged_graph_types': len(merged_g),
                       'merged_tournament_types': len(merged_t),
                       'incidence_edges_merged': sum(w > 0 for row in W for w in row)})
        if n == 4:
            for support_size in range(2, len(merged_g) + 1):
                for support in combinations(range(len(merged_g)), support_size):
                    relation = null_vector([forward[i] for i in support])
                    if relation is None:
                        continue
                    require(sum(relation) == 0 and any(x > 0 for x in relation), 'balanced signed mixture')
                    total = sum(max(0, x) for x in relation)
                    positive = [F(max(0, x), total) for x in relation]
                    negative = [F(max(0, -x), total) for x in relation]
                    left = [sum(positive[k]*forward[i][j] for k, i in enumerate(support))
                            for j in range(len(merged_t))]
                    right = [sum(negative[k]*forward[i][j] for k, i in enumerate(support))
                             for j in range(len(merged_t))]
                    require(left == right and positive != negative, 'distinct probability mixtures same output')
                    witness = {'graph_pair_representatives': [merged_g[i] for i in support],
                               'edge_lists': [[[u, v] for k, (u, v) in enumerate(pairs)
                                               if merged_g[i] >> k & 1] for i in support],
                               'integer_relation': relation,
                               'positive_mixture': list(map(str, positive)),
                               'negative_mixture': list(map(str, negative)),
                               'tournament_pair_representatives': merged_t,
                               'common_distribution': list(map(str, left)),
                               'support_minimal_in_this_census': support_size,
                               'merged_incidence_matrix': W,
                               'stationary_distribution': list(map(str, stationary))}
                    break
                if witness:
                    break
    return {'status': 'FINITE-EXACT', 'checks': checks, 'census': census,
            'mixture_collision_order4': witness}


if __name__ == '__main__':
    result = main()
    receipt = Path(__file__).with_name('flow_receipt.json')
    if '--write' in sys.argv:
        receipt.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    else:
        require(result == json.loads(receipt.read_text(encoding='utf-8')), 'saved flow receipt')
    print(json.dumps({'status': 'PASS', 'checks': checks,
                      'mixture_collision_order4': result['mixture_collision_order4']}))
