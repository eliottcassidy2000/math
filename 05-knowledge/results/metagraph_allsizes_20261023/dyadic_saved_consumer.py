"""Independent consumer of the order-eight packet's literal dyadic lifts.

This imports no code from the producing packet. It constructs a sparse
integer reduced Laplacian, uses highest-bit binary elimination, recomputes
the entire half-form, and consumes the saved higher lifts without finding
any corrections itself.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import csv
import gzip
import io
import json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def lf(path):
    raw = path.read_bytes()
    if path.suffix == '.gz':
        raw = gzip.decompress(raw)
    return raw.replace(b'\r\n', b'\n')


def bits(x):
    while x:
        p = x.bit_length() - 1
        yield p
        x ^= 1 << p


def echelon(rows):
    pivots = {}
    for row in rows:
        while row:
            p = row.bit_length() - 1
            if p in pivots:
                row ^= pivots[p]
            else:
                pivots[p] = row
                break
    return pivots


def nullspace(rows, width):
    need(all(type(r) is int and 0 <= r < (1 << width) for r in rows),
         'binary matrix type and width')
    pivots = echelon(rows)
    basis = []
    order = sorted(pivots)
    for free in range(width):
        if free in pivots:
            continue
        x = 1 << free
        for p in order:
            if (pivots[p] & x).bit_count() & 1:
                x ^= 1 << p
        need(all(not ((r & x).bit_count() & 1) for r in rows),
             'independent binary nullspace vector')
        basis.append(x)
    need(len(basis) + len(pivots) == width, 'rank-nullity')
    return basis, len(pivots)


def combine(coefficients, basis):
    x = 0
    for i in bits(coefficients):
        x ^= basis[i]
    return x


def same_span(left, right, name):
    rl = len(echelon(left))
    rr = len(echelon(right))
    need(rl == len(left) and rr == len(right), name + ': both bases independent')
    need(rl == rr == len(echelon(left + right)), name + ': exactly equal spans')


def read_graph(receipt):
    graph = RESULTS / 'metagraph_order8_20261023'
    raw_edges = lf(graph / 'merged_edges.tsv.gz')
    raw_nodes = lf(graph / 'merged_nodes.tsv')
    hashes = {k: sha256(v).hexdigest() for k, v in
              [('edges', raw_edges), ('nodes', raw_nodes)]}
    need(all(hashes[k] == receipt['input_sha256_lf'][k] for k in hashes),
         'complete graph input hashes')
    nodes = list(csv.DictReader(io.StringIO(raw_nodes.decode()), delimiter='\t'))
    n = len(nodes)
    need([int(x['merged_id']) for x in nodes] == list(range(n)), 'literal vertex IDs')
    root = receipt['order8']['root']
    need(n == receipt['order8']['vertices'] and root == n - 1,
         'order-eight root and dimension')
    # Sparse rows also act as columns because every insertion is symmetric.
    A = [{} for _ in range(n - 1)]
    adjacency = [set() for _ in range(n)]
    edges = list(csv.DictReader(io.StringIO(raw_edges.decode()), delimiter='\t'))
    for e in edges:
        u, v, blue, black = (int(e[k]) for k in ('source', 'target', 'blue', 'black'))
        need(0 <= u <= v < n and blue >= 0 and black >= 0 and black % 2 == 0,
             'typed colored edge')
        w = blue + black
        if u == v or not w:
            continue
        adjacency[u].add(v)
        adjacency[v].add(u)
        if u != root:
            A[u][u] = A[u].get(u, 0) + w
        if v != root:
            A[v][v] = A[v].get(v, 0) + w
        if u != root and v != root:
            A[u][v] = A[u].get(v, 0) - w
            A[v][u] = A[v].get(u, 0) - w
    reached = {root}
    stack = [root]
    while stack:
        for v in adjacency[stack.pop()]:
            if v not in reached:
                reached.add(v)
                stack.append(v)
    need(len(reached) == n, 'connected graph gives nonsingular reduced Laplacian')
    need(all(A[j].get(i, 0) == w for i, row in enumerate(A) for j, w in row.items()),
         'literal sparse matrix symmetric')
    return A, n, root, hashes


def sparse_image(A, x):
    # Ordinary row multiplication: independent of the producer's edge flow.
    return [sum(w * x[j] for j, w in row.items()) for row in A]


def binary_potential_image(A, potential):
    # Fast column addition for sparse binary potentials, using symmetry.
    out = [0] * len(A)
    for j in bits(potential):
        for i, w in A[j].items():
            out[i] += w
    return out


def independent_half_form(A):
    rows = [sum((w % 2) << j for j, w in row.items()) for row in A]
    V1, rankA = nullspace(rows, len(A))
    half_images = []
    for b in V1:
        image = binary_potential_image(A, b)
        need(all(v % 2 == 0 for v in image), 'actual binary lift has even image')
        half_images.append(sum(((v // 2) % 2) << j for j, v in enumerate(image)))
    # q1 evaluated directly in the independently chosen full binary kernel.
    Q = [sum(((u & image).bit_count() % 2) << j for j, u in enumerate(V1))
         for image in half_images]
    need(all(((Q[i] >> j) & 1) == ((Q[j] >> i) & 1)
             for i in range(len(Q)) for j in range(i)), 'independent full half-form symmetric')
    radical, rankQ = nullspace(Q, len(V1))
    actual_radical = [combine(c, V1) for c in radical]
    return V1, rankA, Q, rankQ, actual_radical


def consume_layer(A, n, root, level, lifts, saved_form, expected_span):
    need(type(lifts) is list and type(saved_form) is list, 'literal saved lists')
    need(all(type(v) is int for v in saved_form), 'saved form uses literal integer rows')
    modulus = 1 << level
    parities = []
    images = []
    for x in lifts:
        need(type(x) is list and len(x) == n and
             all(type(v) is int for v in x) and x[root] == 0,
             'literal integer lift with root zero at level ' + str(level))
        parities.append(sum((v % 2) << j for j, v in enumerate(x[:-1])))
        image = sparse_image(A, x[:-1])
        need(all(v % modulus == 0 for v in image), 'actual image divisibility at level ' + str(level))
        # The omitted root image is minus the sum of all reduced images.
        need((-sum(image)) % modulus == 0, 'root image divisibility')
        images.append(image)
    same_span(parities, expected_span, 'nested radical at level ' + str(level))
    form = []
    for x in lifts:
        row = 0
        for j, image in enumerate(images):
            energy = sum(v * a for v, a in zip(x[:-1], image))
            need(energy % modulus == 0, 'integer energy quotient')
            row |= ((energy // modulus) % 2) << j
        form.append(row)
    need(form == saved_form, 'saved full form agrees at level ' + str(level))
    need(all(((form[i] >> j) & 1) == ((form[j] >> i) & 1)
             for i in range(len(form)) for j in range(i)), 'saved form symmetric')
    radical, rank = nullspace(form, len(lifts))
    return [combine(c, parities) for c in radical], {
        'level': level, 'dimension': len(lifts), 'rank': rank,
        'next_dimension': len(radical), 'form_rows': form}


def hostile_controls(A, n, root, lifts, form, expected):
    def rejected(xs, rows, message):
        try:
            consume_layer(A, n, root, 2, xs, rows, expected)
        except ValueError:
            need(True, message)
        else:
            need(False, message)

    # A full modulus change is a valid alternative lift, not a bad certificate.
    gauge = [x[:] for x in lifts]
    gauge[0][0] += 4
    _, entry = consume_layer(A, n, root, 2, gauge, form, expected)
    need(entry['rank'] == 4, 'noncanonical lift accepted with unchanged form')
    bad = [x[:] for x in lifts]
    bad[0][0] = float(bad[0][0])
    rejected(bad, form, 'float alias rejected')
    bad = [x[:] for x in lifts]
    j = next(j for j, row in enumerate(A) if any(v % 4 for v in row.values()))
    bad[0][j] += 1
    rejected(bad, form, 'actual image divisibility mutation rejected')
    rejected(lifts[:-1], form, 'omitted independent mode rejected')
    badform = form[:]
    badform[0] ^= 1
    rejected(lifts, badform, 'energy bit mutation rejected')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    producer = RESULTS / 'metagraph_order8_electrical_20261023'
    raw_receipt = lf(producer / 'receipt.json')
    receipt = json.loads(raw_receipt)
    claim = receipt['order8']
    A, n, root, inputs = read_graph(receipt)
    V1, rankA, Q, rankQ, expected = independent_half_form(A)
    need(rankA == claim['blue_reduced_rank_mod2'], 'independent Laplacian mod-two rank')
    need(len(V1) == claim['critical_two_rank'] and rankQ == claim['half_form_rank'],
         'independent entire half-form dimension and rank')
    need(len(expected) == claim['four_divisible_factors'], 'independent first radical dimension')
    hostile_controls(A, n, root, claim['radical_lifts_mod4'], claim['quarter_form_rows'], expected)
    layers = [{'level': 1, 'dimension': len(V1), 'rank': rankQ, 'next_dimension': len(expected)}]
    saved = [(2, claim['radical_lifts_mod4'], claim['quarter_form_rows']),
             (3, claim['eighth_radical_lifts_mod8'], claim['eighth_form_rows'])]
    saved += [(s['level'], s['lifts'], s['form_rows']) for s in claim['higher_layers']]
    need([s[0] for s in saved] == list(range(2, 9)), 'consecutive saved layers through eight')
    for level, lifts, form in saved:
        expected, entry = consume_layer(A, n, root, level, lifts, form, expected)
        if level >= 4:
            declared = next(s for s in claim['higher_layers'] if s['level'] == level)
            need(all(entry[k] == declared[k] for k in ('dimension', 'rank', 'next_dimension')),
                 'higher saved dimensions and rank')
        layers.append(entry)
    need(not expected and claim['complete_two_primary'] is True, 'filtration genuinely terminates')
    multiplicities = {str(1 << s['level']): s['rank'] for s in layers if s['rank']}
    valuation = sum(s['level'] * s['rank'] for s in layers)
    need(multiplicities == {'2': 3358, '4': 4, '256': 1} and valuation == 3374,
         'exact finite 2-primary group and tree valuation')
    result = {
        'status': 'FINITE-EXACT independent saved-lift consumer',
        'checks': CHECKS,
        'input_sha256_lf': inputs,
        'producer_receipt_sha256_lf': sha256(raw_receipt).hexdigest(),
        'producer_code_sha256_lf': sha256(lf(producer / 'check.py')).hexdigest(),
        'consumer_code_sha256_lf': sha256(lf(Path(__file__))).hexdigest(),
        'independent_half_form_sha256': sha256(json.dumps(Q).encode()).hexdigest(),
        'vertices': n, 'root': root, 'reduced_mod_two_rank': rankA,
        'layers': layers, 'two_primary_cyclic_multiplicities': multiplicities,
        'tree_count_two_valuation': valuation,
        'scope': 'Consumes the authenticated finite graph export; does not re-enumerate tournaments.'}
    target = HERE / 'dyadic_saved_receipt.json'
    if args.write:
        target.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    else:
        need(json.loads(target.read_text(encoding='utf-8')) == result, 'independent saved receipt')
    print(json.dumps({'checks': CHECKS, 'layers': layers,
                      'two_primary_cyclic_multiplicities': multiplicities,
                      'tree_count_two_valuation': valuation}, indent=2))


if __name__ == '__main__':
    main()
