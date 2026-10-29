"""Joint two-run endpoint modules, integral run-swap glue and hostiles."""
import argparse
from collections import Counter
from itertools import combinations, permutations
import json
from math import prod
from pathlib import Path

from atom_precision import automorphisms, determinant, hostile
from paley_pencil import diagonal, local_layers, valuation
from sympy import Matrix, Poly, factorint, rem, symbols

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise ValueError(why)


def subset_paths(rows):
    n = len(rows)
    table = [[[0]*n for _ in range(n)] for _ in range(1 << n)]
    for v in range(n):
        table[1 << v][v][v] = 1
    for mask in range(1, 1 << n):
        for s in range(n):
            for e, value in enumerate(table[mask][s]):
                if not value:
                    continue
                todo = rows[e] & ~mask
                while todo:
                    bit = todo & -todo
                    todo -= bit
                    table[mask | bit][s][bit.bit_length()-1] += value
    return table


def two_runs(rows):
    n = len(rows)
    ports = list(permutations(range(n), 2))
    index = {p: i for i, p in enumerate(ports)}
    paths = subset_paths(rows)
    full = (1 << n)-1
    response = [[0]*len(ports) for _ in ports]
    for mask in range(1, full):
        left = [(i, j, x) for i, row in enumerate(paths[mask]) for j, x in enumerate(row) if x]
        right = [(i, j, x) for i, row in enumerate(paths[full ^ mask]) for j, x in enumerate(row) if x]
        for a, b, x in left:
            for c, d, y in right:
                response[index[a, c]][index[b, d]] += x*y
    return ports, response, paths[full]


def walk_pair_inclusion_exclusion(rows, starts, ends):
    """Independent walks with n total vertex occurrences, then full support."""
    n, result = len(rows), 0
    forced = sum(1 << v for v in set(starts+ends))
    other = [v for v in range(n) if not(forced >> v & 1)]
    for subset in range(1 << len(other)):
        mask = forced | sum(1 << v for i, v in enumerate(other) if subset >> i & 1)
        counts = []
        for s, t in zip(starts, ends):
            vector = [int(v == s) for v in range(n)]
            values = [vector[t]]
            for _ in range(n-2):
                new = [0]*n
                for u, value in enumerate(vector):
                    todo = rows[u] & mask if value else 0
                    while todo:
                        bit = todo & -todo
                        todo -= bit
                        new[bit.bit_length()-1] += value
                vector = new
                values.append(vector[t])
            counts.append(values)
        term = sum(counts[0][k]*counts[1][n-2-k] for k in range(n-1))
        result += (-1)**(n-mask.bit_count())*term
    return result


def certify_smith(matrix):
    ds = diagonal(matrix)
    need(all(ds) and prod(ds) == abs(determinant(matrix)), 'independent Bareiss determinant')
    primes = sorted({int(p) for d in set(ds) for p in factorint(d)})
    certificates = {}
    for p in primes:
        exponents = [valuation(d, p) for d in ds]
        precision = max(exponents)+1
        seen = local_layers(matrix, p, precision)
        counts = Counter(exponents)
        need(seen == [counts[k] for k in range(precision)], 'independent local Smith layers')
        certificates[str(p)] = seen
    return dict(smith=ds, local_layers=certificates)


def certify(label, rows, full_smith=True):
    n = len(rows)
    ports, M, H = two_runs(rows)
    index = {p: i for i, p in enumerate(ports)}
    if n <= 7:
        brute = [[0]*len(ports) for _ in ports]
        for path in permutations(range(n)):
            mandatory = [i+1 for i in range(n-1) if not(rows[path[i]] >> path[i+1] & 1)]
            if len(mandatory) > 1:
                continue
            for cut in mandatory if mandatory else range(1, n):
                brute[index[path[0], path[cut]]][index[path[cut-1], path[-1]]] += 1
        need(M == brute, 'complete permutation-and-cut census')
    nonzero = [(i, j) for i, row in enumerate(M) for j, x in enumerate(row) if x]
    probes = {(0, 0), (0, len(ports)-1)}
    probes.update(nonzero[k] for k in (0, len(nonzero)//2, len(nonzero)-1))
    for i, j in sorted(probes):
        need(M[i][j] == walk_pair_inclusion_exclusion(rows, ports[i], ports[j]),
             'independent two-walk inclusion exclusion')
    group = automorphisms(rows)
    for g in group:
        moved = [index[g[a], g[b]] for a, b in ports]
        need(all(M[moved[i]][moved[j]] == M[i][j]
                 for i in range(len(ports)) for j in range(len(ports))), 'full transporter covariance')
    swapped = [index[b, a] for a, b in ports]
    need(all(M[swapped[i]][swapped[j]] == M[i][j]
             for i in range(len(ports)) for j in range(len(ports))), 'simultaneous run swap')
    # This also checks the baseline hostile endpoint result against its saved receipt.
    inherited = json.loads((HERE/'atom_precision.json').read_text(encoding='utf-8'))
    baseline = next(x for x in inherited['objects'] if x['label'] == label)
    need(H == baseline['endpoints'], 'inherited single-run matrix unchanged')
    answer = dict(label=label, rows=rows, ports=ports, response=M,
                  total_two_run_covers=sum(map(sum, M)), total_hamilton_paths=sum(map(sum, H)),
                  full_automorphism_order=len(group), independent_probe_count=len(probes))
    if full_smith:
        answer.update(certify_smith(M))
        unordered = list(combinations(range(n), 2))
        Bplus = [[M[index[x]][index[y]]+M[index[x]][index[y[::-1]]] for y in unordered] for x in unordered]
        Bminus = [[M[index[x]][index[y]]-M[index[x]][index[y[::-1]]] for y in unordered] for x in unordered]
        need(all((a-b) % 2 == 0 for u, v in zip(Bplus, Bminus) for a, b in zip(u, v)),
             'the two sector maps preserve parity glue')
        answer['symmetric'] = dict(matrix=Bplus, **certify_smith(Bplus))
        answer['antisymmetric'] = dict(matrix=Bminus, **certify_smith(Bminus))
        need(abs(determinant(M)) == abs(determinant(Bplus)*determinant(Bminus)),
             'rational splitting preserves determinant')
    else:
        layers = local_layers(M, 3, 12)
        # This finite digit census is declared complete only if every pivot is seen.
        need(sum(layers) == len(ports), 'complete 3-primary valuation census')
        answer['three_layers'] = layers
    return answer


def six_vertex_control():
    pairs = list(combinations(range(4), 2))
    child = [0]*4
    for k, (i, j) in enumerate(pairs):
        u, v = (i, j) if 7 >> k & 1 else (j, i)
        child[u] |= 1 << v
    parent = child+[0, 0]
    for u in range(4):
        parent[u] |= 1 << 4
        parent[5] |= 1 << u
    parent[4] |= 1 << 5
    def complemented(rows, paths):
        size = len(rows)
        out = [((1 << size)-1) ^ row ^ (1 << i) for i, row in enumerate(rows)]
        for path in paths:
            for u, v in zip(path, path[1:]):
                need(rows[u] >> v & 1, 'legal directed component')
                out[u] ^= 1 << v
                out[v] ^= 1 << u
        return out
    def canonical(rows, marks=()):
        size, results = len(rows), []
        for order in permutations(range(size)):
            inverse = {v: i for i, v in enumerate(order)}
            code = sum((rows[order[i]] >> order[j] & 1) << k
                       for k, (i, j) in enumerate(combinations(range(size), 2)))
            results.append((tuple(inverse[v] for v in marks), code))
        return min(results)
    C, D = ((0, 1), (3, 2)), ((0, 2, 1), (3,))
    outer = [canonical(complemented(child, cover), (cover[0][0], cover[-1][-1])) for cover in (C, D)]
    need(outer[0] == outer[1], 'the two outer endpoints agree as rooted target types')
    all_ports = [canonical(complemented(child, cover), tuple(x for path in cover for x in (path[0], path[-1]))) for cover in (C, D)]
    need(all_ports[0] != all_ports[1], 'joint ordered ports retain the distinction')
    targets = [canonical(complemented(parent, [path]))[1]
               for path in ((0, 1, 4, 5, 3, 2), (0, 2, 1, 4, 5, 3))]
    need(targets == [81, 83], 'literal complement revisit hostile unchanged')
    return dict(outer_types=outer, joint_port_types=all_ports, complement_target_codes=targets)


def affine_convolution(record):
    M = record['response']
    index = {tuple(p): i for i, p in enumerate(record['ports'])}
    group = [(b, a) for a in (1, 2, 4) for b in range(7)]
    gi = {g: i for i, g in enumerate(group)}
    z = symbols('z')
    phi = sum(z**i for i in range(7))
    exponents = (1, 2, 4)
    output = []
    for sign, name in ((1, 'symmetric'), (-1, 'antisymmetric')):
        kernel = [M[index[0, 1]][index[b, (b+a) % 7]]
                  +sign*M[index[0, 1]][index[(b+a) % 7, b]] for b, a in group]
        sector = []
        for b, a in group:
            inv = pow(a, -1, 7)
            row = []
            for c, d in group:
                value = M[index[b, (b+a) % 7]][index[c, (c+d) % 7]]
                value += sign*M[index[b, (b+a) % 7]][index[(c+d) % 7, c]]
                need(value == kernel[gi[((c-b)*inv % 7, d*inv % 7)]], 'affine convolution law')
                row.append(value)
            sector.append(row)
        totals = [sum(kernel[k*7:(k+1)*7]) for k in range(3)]
        trivial = sum(totals)
        cubic_norm = sum(x*x for x in totals)-sum(totals[i]*totals[j] for i, j in combinations(range(3), 2))
        point_matrix = [[0]*7 for _ in range(7)]
        block = Matrix.zeros(3)
        for weight, (b, a) in zip(kernel, group):
            inv = pow(a, -1, 7)
            for x in range(7):
                point_matrix[(a*x+b) % 7][x] += weight
            for j, e in enumerate(exponents):
                block[exponents.index(e*inv % 7), j] += weight*z**((-e*b*inv) % 7)
        determinant_polynomial = Poly(rem(block.det(), phi, z), z)
        c = int(determinant_polynomial.nth(0))
        d = int(determinant_polynomial.nth(1))
        need(rem(determinant_polynomial.as_expr()-c-d*(z+z**2+z**4), phi, z) == 0,
             'quadratic fixed field determinant')
        need(rem(determinant_polynomial.as_expr().subs(z, z*z)-determinant_polynomial.as_expr(), phi, z) == 0,
             'order-three Galois phase preserves the block determinant')
        norm = c*c-c*d+2*d*d
        need(determinant(point_matrix) == trivial*norm, 'independent seven-point norm determinant')
        need(determinant(sector) == trivial*cubic_norm*norm**3, 'regular representation determinant')
        output.append(dict(sector=name, group=group, kernel=kernel, quotient_totals=totals,
                           trivial=trivial, cubic_norm=cubic_norm,
                           quadratic_determinant=[c, d], quadratic_norm=norm))
    need([x['quadratic_determinant'] for x in output] == [[24, -24], [248, 608]],
         'explicit quadratic determinant values')
    return output


def run():
    records = []
    for n, label in ((3, 'Paley3'), (5, 'T5'), (7, 'Paley7')):
        steps = {x*x % n for x in range(1, n)} if n != 5 else {1, 2}
        rows = [sum(1 << ((u+d) % n) for d in steps) for u in range(n)]
        records.append(certify(label, rows))
    records.extend(certify(label, hostile(flip), False)
                   for label, flip in (('81-control', False), ('27-control', True)))
    need([x['three_layers'][0] for x in records[-2:]] == [24, 22], 'joint-port hostile separator')
    need(records[-2]['total_two_run_covers'] == records[-1]['total_two_run_covers'] == 9011898,
         'total-cover hostile stays matched')
    return dict(objects=records, six_vertex_revisit=six_vertex_control(),
                paley7_affine_convolution=affine_convolution(records[2]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    target = HERE/'two_run.json'
    # JSON normalizes the tuple port keys and lists used in the pure checker.
    encoded = json.dumps(result, indent=2)+'\n'
    if args.write:
        target.write_text(encoded, encoding='utf-8')
    else:
        need(json.loads(encoded) == json.loads(target.read_text(encoding='utf-8')), 'receipt replay')
    print(json.dumps(dict(status='PASS', objects=[dict(label=x['label'], mass=x['total_two_run_covers'],
                      three_layers=x.get('three_layers')) for x in result['objects']])))


if __name__ == '__main__':
    main()
