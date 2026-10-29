"""Native three-run responses and the exact S3 integral Fourier guard."""
import argparse
from itertools import combinations, permutations
import json
from pathlib import Path

from atom_precision import automorphisms, determinant
from paley_pencil import diagonal, local_layers
from two_run import subset_paths
from sympy import Matrix

HERE = Path(__file__).resolve().parent
GROUP = list(permutations(range(3)))


def need(ok, why):
    if not ok:
        raise ValueError(why)


def multiply(A, B):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*B)] for row in A]


def fourier():
    columns, matrices = [], {}
    for g in GROUP:
        sign = (-1)**sum(g[i] > g[j] for i in range(3) for j in range(i+1, 3))
        images = []
        for source in (0, 1):
            v = [0, 0, 0]
            v[g[source]] += 1
            v[g[2]] -= 1
            images.append(v[:2])
        rho = [[images[j][i] for j in range(2)] for i in range(2)]
        matrices[g] = rho
        columns.append([1, sign]+sum(rho, []))
    F = [list(row) for row in zip(*columns)]
    inverse = Matrix(F).inv()*6
    need(all(x.q == 1 for x in inverse), 'inverse denominator divides six')
    inverse = [[int(x) for x in row] for row in inverse.tolist()]
    sixI = [[6*int(i == j) for j in range(6)] for i in range(6)]
    need(multiply(F, inverse) == multiply(inverse, F) == sixI, 'exact inverse on both sides')
    need(abs(determinant(F)) == 54 and diagonal(F) == [1, 1, 1, 3, 3, 6], 'index54 integral character map')
    need(local_layers(inverse[:3], 3, 1) == local_layers(inverse, 3, 1) == [3],
         'first three congruences span every ternary guard')
    need(all([x % 2 for x in row] == [1, 1, 0, 0, 0, 0] for row in inverse),
         'the single binary guard')
    for g in GROUP:
        for h in GROUP:
            composed = tuple(g[h[i]] for i in range(3))
            need(multiply(matrices[g], matrices[h]) == matrices[composed], 'standard representation multiplication')
    modular = {}
    for modulus in (2, 3, 4, 6):
        reached = set()
        # This exhaustive universe is deliberately only 6^6 at its largest.
        from itertools import product
        for source in product(range(modulus), repeat=6):
            reached.add(tuple(sum(a*b for a, b in zip(row, source)) % modulus for row in F))
        from math import gcd
        expected = modulus**6//(gcd(modulus, 3)**2*gcd(modulus, 6))
        need(len(reached) == expected, 'independent finite image count')
        if modulus == 6:
            guarded = {y for y in product(range(6), repeat=6)
                       if (y[0]-y[1]) % 2 == 0
                       and all(sum(a*b for a, b in zip(row, y)) % 3 == 0 for row in inverse[:3])}
            need(guarded == reached, 'one binary and three ternary guards exactly reconstruct the source lattice')
        modular[str(modulus)] = expected
    return F, inverse, modular


def native(rows):
    n = len(rows)
    triples = list(combinations(range(n), 3))
    ports = [tuple(triple[g[i]] for i in range(3)) for triple in triples for g in GROUP]
    index = {port: i for i, port in enumerate(ports)}
    M = [[0]*len(ports) for _ in ports]
    for path in permutations(range(n)):
        required = {i+1 for i in range(n-1) if not(rows[path[i]] >> path[i+1] & 1)}
        if len(required) > 2:
            continue
        optional = [i for i in range(1, n) if i not in required]
        for extra in combinations(optional, 2-len(required)):
            bounds = [0]+sorted(required | set(extra))+[n]
            starts = tuple(path[bounds[i]] for i in range(3))
            ends = tuple(path[bounds[i+1]-1] for i in range(3))
            M[index[starts]][index[ends]] += 1
    # Independent disjoint-subset convolution, with no permutation cuts.
    table = subset_paths(rows)
    endpoints = [[(i, j, count) for i, row in enumerate(value) for j, count in enumerate(row) if count]
                 for value in table]
    full = (1 << n)-1
    other = [[0]*len(ports) for _ in ports]
    for first in range(1, full):
        remaining = full ^ first
        second = (remaining-1) & remaining
        while second:
            third = remaining ^ second
            for a, b, u in endpoints[first]:
                for c, d, v in endpoints[second]:
                    for e, f, w in endpoints[third]:
                        other[index[a, c, e]][index[b, d, f]] += u*v*w
            second = (second-1) & remaining
    need(M == other, 'complete independent three-subset path convolution')
    for g in automorphisms(rows):
        moved = [index[tuple(g[v] for v in port)] for port in ports]
        need(all(M[moved[i]][moved[j]] == M[i][j]
                 for i in range(len(ports)) for j in range(len(ports))), 'native transporter action')
    for g in GROUP:
        moved = [index[tuple(port[g[i]] for i in range(3))] for port in ports]
        need(all(M[moved[i]][moved[j]] == M[i][j]
                 for i in range(len(ports)) for j in range(len(ports))), 'joint run permutation action')
    return triples, ports, M


def split_matrix(M, F, inverse):
    size = len(M)
    N = [[0]*size for _ in range(size)]
    for start in range(0, size, 6):
        for end in range(0, size, 6):
            block = [M[start+i][end:end+6] for i in range(6)]
            image = multiply(multiply(F, block), inverse)
            need(all(x % 6 == 0 for row in image for x in row), 'integer action on the ambient split sectors')
            for i in range(6):
                for j in range(6):
                    N[start+i][end+j] = image[i][j]//6
    # Direct intertwiner check, avoiding any conclusion from spectra alone.
    for start in range(0, size, 6):
        for end in range(0, size, 6):
            lhs = multiply(F, [M[start+i][end:end+6] for i in range(6)])
            rhs = multiply([N[start+i][end:end+6] for i in range(6)], F)
            need(lhs == rhs, 'source-preserving Fourier intertwiner')
    return N


def boundary_kernel(triples, M):
    # Oriented tetrahedron boundaries in the ten sorted triangles on five vertices.
    boundary = []
    for triple in triples:
        boundary.append([sum((-1)**j for j in range(4) if tetra[:j]+tetra[j+1:] == triple)
                         for tetra in combinations(range(5), 4)])
    need(local_layers(boundary, 5, 1) == [4], 'four independent boundary relations')
    lifted = [row[:] for row in boundary for _ in GROUP]
    need(all(x == 0 for row in multiply(M, lifted) for x in row), 'all five tetrahedron boundaries lie in the kernel')
    need(local_layers(M, 5, 1) == [56], 'rank56 lower bound plus four independent null vectors')
    return boundary


def complement_probe(record):
    rows, M = record['rows'], record['native']
    B = [row[:] for row in record['tetrahedron_boundary'] for _ in GROUP]
    converse = [31 ^ row ^ (1 << i) for i, row in enumerate(rows)]
    _, _, reverse_matrix = native(converse)
    need(reverse_matrix == [list(row) for row in zip(*M)], 'global converse transposes joint response')
    need(all(x == 0 for row in multiply(reverse_matrix, B) for x in row), 'boundary kernel survives global converse')
    cases = []
    for path in permutations(range(5)):
        if not all(rows[u] >> v & 1 for u, v in zip(path, path[1:])):
            continue
        target = converse[:]
        for u, v in zip(path, path[1:]):
            target[u] ^= 1 << v
            target[v] ^= 1 << u
        for start in range(5):
            reached = 1 << start
            while True:
                updated = reached
                for u in range(5):
                    if reached >> u & 1:
                        updated |= target[u]
                if updated == reached:
                    break
                reached = updated
            need(reached == 31, 'literal complement target is strong')
        _, _, N = native(target)
        need(local_layers(M+N, 5, 1) == [60], 'one legal complement kills all four hidden directions')
        cases.append(dict(path=path, target_rows=target, stacked_rank=60))
        if path == (0, 1, 2, 3, 4):
            # Boundary of tetrahedron0124: 124 + 014 - 024 - 012.
            column, source = 1, 0
            plus = [i for i, row in enumerate(B) if row[column] == 1]
            minus = [i for i, row in enumerate(B) if row[column] == -1]
            before = [sum(M[source][i] for i in side) for side in (plus, minus)]
            after = [sum(N[source][i] for i in side) for side in (plus, minus)]
            need(before == [1, 1] and after == [2, 0] and len(plus) == len(minus) == 12,
                 'positive equal-mass port-law witness')
            witness = dict(source=record['ports'][source], plus_columns=plus, minus_columns=minus,
                           before=before, after=after, target_rows=target)
    need(len(cases) == 15, 'every legal Hamilton complement was tested')
    # The nullspace is not present after imposing full vertex-automorphism invariance.
    ports = [tuple(p) for p in record['ports']]
    index = {p: i for i, p in enumerate(ports)}
    group, unseen, orbits = automorphisms(rows), set(range(len(ports))), []
    while unseen:
        seed = min(unseen)
        orbit = {index[tuple(g[v] for v in ports[seed])] for g in group}
        orbits.append(sorted(orbit))
        unseen -= orbit
    invariant = [[sum(row[i] for i in orbit) for orbit in orbits] for row in M]
    need(len(orbits) == 12 and local_layers(invariant, 7, 1) == [12], 'vertex-unnamed invariant laws have no hidden directions')
    return dict(all_legal_complements=cases, witness=witness,
                global_converse_stacked_rank=56, vertex_invariant_orbits=orbits, vertex_invariant_rank=12)


def run():
    F, inverse, modular = fourier()
    objects = []
    for n, label in ((3, 'C3'), (5, 'T5'), (7, 'Paley7')):
        steps = {x*x % n for x in range(1, n)} if n != 5 else {1, 2}
        rows = [sum(1 << ((u+d) % n) for d in steps) for u in range(n)]
        triples, ports, M = native(rows)
        N = split_matrix(M, F, inverse)
        rank = len(ports) if n != 5 else 56
        if n != 5:
            need(local_layers(M, 7, 1) == [rank], 'full rational rank certified at a good prime')
        item = dict(label=label, rows=rows, triples=triples, ports=ports,
                    total_covers=sum(map(sum, M)), rational_rank=rank,
                    native=M, ambient_split=N, primary={})
        for p in (2, 3):
            native_layers, split_layers = local_layers(M, p, 8), local_layers(N, p, 8)
            need(sum(native_layers) == sum(split_layers) == rank, 'complete nonzero primary layers')
            if rank == len(ports):
                need(sum(e*v for e, v in enumerate(native_layers)) == sum(e*v for e, v in enumerate(split_layers)),
                     'equal full-rank determinant valuations')
            item['primary'][str(p)] = dict(native=native_layers, ambient_split=split_layers)
        if n == 5:
            item['tetrahedron_boundary'] = boundary_kernel(triples, M)
            need(sum(1 for x in diagonal(M) if x) == rank, 'independent exact integer Smith rank')
        objects.append(item)
    need(objects[2]['primary']['2']['native'][0] == 165 and objects[2]['primary']['2']['ambient_split'][0] == 160,
         'actual dyadic three-run glue witness')
    need(objects[2]['primary']['3']['native'][0] == 201 and objects[2]['primary']['3']['ambient_split'][0] == 197,
         'actual ternary three-run glue witness')
    return dict(group=GROUP, fourier=F, six_times_inverse=inverse, modular_image_counts=modular, objects=objects,
                t5_actual_complement=complement_probe(objects[1]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    # Keep the large matrices compact; the receipt remains plain JSON.
    encoded = json.dumps(result, separators=(',', ':'))+'\n'
    target = HERE/'three_run.json'
    if args.write:
        target.write_text(encoded, encoding='utf-8')
    else:
        need(json.loads(encoded) == json.loads(target.read_text(encoding='utf-8')), 'receipt replay')
    print(json.dumps(dict(status='PASS', objects=[dict(label=x['label'], mass=x['total_covers'],
                                                     rank=x['rational_rank']) for x in result['objects']])) )


if __name__ == '__main__':
    main()
