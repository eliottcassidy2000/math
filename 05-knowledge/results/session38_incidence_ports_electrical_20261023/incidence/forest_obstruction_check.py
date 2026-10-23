"""Direct dual-incidence seed checks for the all-order forest obstruction.

The induction is proved in FOREST_OBSTRUCTION.md. Finite padding controls
are regression tests, not the source of the all-order quantifier.
"""
from itertools import combinations
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
checks = 0


def need(condition, message):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def graph(S):
    return [[int((v - u) % 15 in S) for v in range(15)]
            for u in range(15)]


def opposite(A):
    return [list(row) for row in zip(*A)]


def incidence_signature(A):
    n = len(A)
    c = [0] * n
    edge = [[0] * n for _ in range(n)]
    for u, v, w in combinations(range(n), 3):
        cyclic = int(A[u][v] == A[v][w] == A[w][u])
        for x in (u, v, w):
            c[x] += cyclic
        for x, y in ((u, v), (u, w), (v, w)):
            edge[x][y] += cyclic
            edge[y][x] += cyclic
    rows = []
    for u in range(n):
        row = tuple(sorted((c[v], edge[u][v]) for v in range(n)
                           if A[u][v]))
        need(sum(e for _, e in row) == c[u], 'outgoing triangle balance')
        rows.append(row)
    return tuple(sorted(rows))


def add_source(A):
    n = len(A)
    return [row + [0] for row in A] + [[1] * n + [0]]


def add_sink(A):
    n = len(A)
    return [row + [1] for row in A] + [[0] * (n + 1)]


def source_signature(E):
    new_row = tuple(sorted((sum(e for _, e in row), 0) for row in E))
    return tuple(sorted(E + (new_row,)))


def sink_signature(E):
    return tuple(sorted([tuple(sorted(row + ((0, 0),))) for row in E] + [()]))


def hamiltonian_count(A):
    n = len(A)
    outgoing = [sum(A[u][v] << v for v in range(n)) for u in range(n)]
    dp = [[0] * n for _ in range(1 << n)]
    for u in range(n):
        dp[1 << u][u] = 1
    for mask in range(1, 1 << n):
        for u, value in enumerate(dp[mask]):
            if not value:
                continue
            choices = outgoing[u] & ~mask
            while choices:
                bit = choices & -choices
                choices -= bit
                v = bit.bit_length() - 1
                dp[mask | bit][v] += value
    return sum(dp[-1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    sets = [[14, 13, 12, 11, 5, 6, 8], [14, 13, 12, 4, 10, 9, 8]]
    seeds = [graph(S) for S in sets]
    signatures = []
    reverse_signatures = []
    counts = []
    for A in seeds:
        need(all(A[u][u] == 0 for u in range(15)), 'seed has a loop')
        need(all(A[u][v] + A[v][u] == 1 for u, v in combinations(range(15), 2)),
             'seed is not a tournament')
        signatures.append(incidence_signature(A))
        reverse_signatures.append(incidence_signature(opposite(A)))
        counts.append(hamiltonian_count(A))
    need(signatures[0] == signatures[1], 'forward exact E differs')
    need(reverse_signatures[0] == reverse_signatures[1], 'converse exact E differs')
    need(counts == [196264905, 196277235], 'Hamiltonian witnesses differ')
    need(counts[0] != counts[1], 'seed classes not separated by H')
    current = seeds
    orders_checked = []
    for n in range(15, 26):
        forward = [incidence_signature(A) for A in current]
        backward = [incidence_signature(opposite(A)) for A in current]
        need(forward[0] == forward[1], 'source-padded forward E differs')
        need(backward[0] == backward[1], 'source-padded converse E differs')
        orders_checked.append(n)
        if n == 25:
            break
        for A, E, E_star in zip(current, forward, backward):
            padded = add_source(A)
            need(incidence_signature(padded) == source_signature(E),
                 'source signature formula')
            need(opposite(padded) == add_sink(opposite(A)),
                 'converse of source padding is sink padding')
            need(incidence_signature(opposite(padded)) == sink_signature(E_star),
                 'sink signature formula')
        current = [add_source(A) for A in current]
    output = {
        'scope': 'Exact dual-E seed collision and finite controls of the proved padding induction; no wiring-nonuniqueness claim.',
        'seed_order': 15,
        'connection_sets': sets,
        'hamiltonian_counts': counts,
        'forward_exact_E_equal': True,
        'converse_exact_E_equal': True,
        'common_forward_E': signatures[0],
        'common_converse_E': reverse_signatures[0],
        'source_padding_orders_checked': orders_checked,
        'all_order_statement': 'Every n >= 15 has two distinct tournament isotypes with identical ordered exact-E/converse-E pair; therefore every E33 quotient signature incidence multigraph has parallel class edges.',
        'checks_before_receipt_comparison': checks,
    }
    output = json.loads(json.dumps(output))
    receipt = HERE / 'forest_obstruction_receipt.json'
    if args.write:
        receipt.write_text(json.dumps(output, indent=2) + '\n')
    else:
        need(json.loads(receipt.read_text()) == output, 'frozen receipt differs')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
