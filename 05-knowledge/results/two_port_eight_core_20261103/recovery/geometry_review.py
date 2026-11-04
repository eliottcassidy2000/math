"""Independent field distances, obstruction witnesses and cube inclusion-exclusion."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
ROOT = HERE.parents[3]
GEOMETRY = PACKET / 'geometry'
CHECKS = 0
ZERO = (F(0),) * 4
ONE = (F(1), F(0), F(0), F(0))


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def plus(a, b):
    return tuple(x + y for x, y in zip(a, b))


def minus(a, b):
    return tuple(x - y for x, y in zip(a, b))


def times(x, y):
    # Explicit basis multiplication, independent of the producer's XOR loop.
    a, b, c, d = x
    e, f, g, h = y
    return (a*e + 3*b*f + 11*c*g + 33*d*h,
            a*f + b*e + 11*c*h + 11*d*g,
            a*g + c*e + 3*b*h + 3*d*f,
            a*h + d*e + b*g + c*f)


def scalar(k):
    return (F(k), F(0), F(0), F(0))


def difference(a, b):
    return minus(a[0], b[0]), minus(a[1], b[1])


def norm(v):
    return plus(times(v[0], v[0]), times(v[1], v[1]))


def distance(a, b):
    return norm(difference(a, b))


def edge_set(n, mask):
    return {edge for i, edge in enumerate(combinations(range(n), 2)) if mask >> i & 1}


def chromatic_by_edge_subsets(n, edges, palettes):
    totals = {q: 0 for q in palettes}
    for mask in range(1 << len(edges)):
        parent = list(range(n))
        def find(v):
            while parent[v] != v:
                v = parent[v]
            return v
        for k, (u, v) in enumerate(edges):
            if mask >> k & 1:
                parent[find(u)] = find(v)
        components = len({find(v) for v in range(n)})
        sign = (-1) ** mask.bit_count()
        for q in palettes:
            totals[q] += sign * q ** components
    return totals


def main():
    paths = [HERE / 'geometry_review.py', HERE / 'GEOMETRY_REVIEW.md',
             GEOMETRY / 'PROOF.md', GEOMETRY / 'verify.py', GEOMETRY / 'receipt.json',
             GEOMETRY / 'SOURCE_PINS.json', PACKET / 'cores/order8.json',
             PACKET / 'cores/independent.py', PACKET / 'cores/independent.cpp',
             PACKET / 'cores/independent.json', PACKET / 'cores/independent.optimized.json',
             HERE / 'normal.json', HERE / 'optimized.json', HERE / 'geometry_replay.json']
    pins = {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}
    for pin in json.loads((GEOMETRY / 'SOURCE_PINS.json').read_text()):
        need(sha256((ROOT / pin['path']).read_bytes()).hexdigest() == pin['sha256'],
             'inherited source pin')
    receipt = json.loads((GEOMETRY / 'receipt.json').read_text())
    need(receipt == json.loads((HERE / 'geometry_replay.json').read_text()), 'producer replay matches frozen receipt')
    need(sha256((GEOMETRY / 'verify.py').read_bytes()).hexdigest() == receipt['verifier_sha256'], 'producer code pin')
    for name, digest in receipt['artifacts_sha256'].items():
        need(sha256((GEOMETRY / name).read_bytes()).hexdigest() == digest, 'artifact pin')
    core_audit = json.loads((PACKET / 'cores/independent.json').read_text())
    need(core_audit == json.loads((PACKET / 'cores/independent.optimized.json').read_text()), 'universe audit normal/O')
    for name, digest in core_audit['source_sha256'].items():
        need(sha256((PACKET / 'cores' / name).read_bytes()).hexdigest() == digest, 'universe audit source pin')
    need(core_audit['eligible_graphs'] == 724080 and core_audit['orbit_union_equals_complete_survivor_set'], 'complete alternate universe')
    need((HERE / 'normal.json').read_bytes() == (HERE / 'optimized.json').read_bytes(), 'coefficient audit normal/O')

    raw = json.loads((GEOMETRY / 'coordinates.json').read_text())
    points = [tuple(tuple(F(c) for c in coordinate) for coordinate in point)
              for point in raw['old_spindle'] + raw['private_cube']]
    need(len(points) == 15 and len(set(points)) == 15, 'fifteen distinct exact points')
    edges = set()
    for i, j in combinations(range(15), 2):
        d = distance(points[i], points[j])
        need(d != ZERO, 'distinct pair distance nonzero')
        if d == ONE:
            edges.add((i, j))
    need(sorted(map(list, edges)) == sorted(raw['full_unit_edges']), 'all105 exact pair decisions')
    need(len(edges) == 25 and 105 - len(edges) == 80, 'unit/nonunit cardinalities')
    old_edges = {(u, v) for u, v in edges if v < 7}
    need(len(old_edges) == 11, 'old Moser edges')
    private = sorted((u - 7, v - 7) for u, v in edges if u >= 7)
    cube = [(u, v) for u, v in combinations(range(8), 2) if (u ^ v).bit_count() == 1]
    need(private == cube, 'private graph is exactly cube')
    cross = [(u, v - 7) for u, v in sorted(edges) if u < 7 <= v]
    need(cross == [(0, 0), (2, 0)] and (0, 2) in old_edges, 'two different boundary colors forbidden')

    unit_triples = []
    for triple in combinations(range(7), 3):
        a, b, c = [points[v] for v in triple]
        v, w = difference(b, a), difference(c, a)
        determinant = minus(times(v[0], w[1]), times(v[1], w[0]))
        radius_one = determinant != ZERO and times(times(norm(v), norm(w)), distance(b, c)) == times(scalar(4), times(determinant, determinant))
        centers = [i for i in range(7) if all(distance(points[i], points[j]) == ONE for j in triple)]
        need(radius_one == bool(centers), 'unit circumradius iff old center exists')
        need(len(centers) <= 1, 'unique noncollinear circumcenter')
        if radius_one:
            unit_triples.append({'triple': list(triple), 'old_center': centers[0]})
    need(len(unit_triples) == 10, 'ten unit triples and no exterior triple center')
    need(unit_triples == json.loads((GEOMETRY / 'spindle_contact_audit.json').read_text())['unit_triples'], 'independent circle-center audit')

    exclusions = json.loads((GEOMETRY / 'metric_exclusions.json').read_text())
    for row in exclusions:
        graph = edge_set(8, row['mask'])
        if row['reason'] == 'rigid_Moser_exterior_capacity':
            embedding = row['moser_embedding_from_standard_coordinates']
            outside = row['exterior_vertex']
            need(len(set(embedding)) == 7 and outside not in embedding, 'distinct Moser and exterior labels')
            need(all(tuple(sorted((embedding[u], embedding[v]))) in graph for u, v in old_edges), 'eleven mandatory Moser edges')
            neighbors = [v for v in range(8) if tuple(sorted((outside, v))) in graph]
            need(neighbors == row['exterior_neighbours'] and len(neighbors) >= 3, 'exterior contact contradiction')
        else:
            hub, path, z = row['hub'], row['ordered_neighbour_path'], row['forced_duplicate_vertex']
            need(len({hub, z, *path}) == 6, 'six distinct fan vertices')
            required = {tuple(sorted(edge)) for edge in
                        [(hub, v) for v in path] + list(zip(path, path[1:])) + [(z, path[0]), (z, path[-1])]}
            need(len(required) == 9 and required <= graph, 'nine mandatory fan edges')
    valid_walks = []
    for mask in range(8):
        walk = [0]
        for k in range(3):
            walk.append((walk[-1] + (1 if mask >> k & 1 else -1)) % 6)
        if len(set(walk)) == 4:
            valid_walks.append(walk)
    need(valid_walks == [[0, 5, 4, 3], [0, 1, 2, 3]], 'only injective four-point hexagon paths have antipodal endpoints')

    chromatic = chromatic_by_edge_subsets(8, cube, (3, 4, 5, 6))
    need([chromatic[q] for q in (4, 5, 6)] == [2652, 29660, 198030], 'independent cube chromatic evaluations')
    tables = json.loads((GEOMETRY / 'boundary_tables.json').read_text())
    for q in (4, 5, 6):
        need(chromatic[q] % q == 0, 'fixed color symmetry integral')
        extensions = chromatic[q] // q * (q - 2)
        rows = tables[str(q)]
        need(len(rows) == {4: 16, 5: 42, 6: 52}[q], 'retained boundary table size')
        total = 0
        for row in rows:
            pi, tail = row['partition'], row['private_witness']
            word = pi + tail
            need(all(word[u] != word[v] for u, v in edges), 'stored witness obeys every edge')
            need(all(0 <= color < q for color in word), 'witness palette valid')
            need(row['extension_count'] == extensions, 'every boundary weight from independent inclusion-exclusion')
            named = factorial(q) // factorial(q - max(pi) - 1)
            need(row['named_boundary_assignments'] == named, 'named boundary multiplicity')
            total += named * extensions
        need(total == receipt['boundary_controls'][str(q)]['full_graph_colourings'], 'full named color count')
    need(pins == {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}, 'review inputs unchanged')
    print(json.dumps({'status': 'PASS INDEPENDENT GEOMETRY AUDIT', 'checks': CHECKS,
                      'source_sha256': pins, 'pair_count': 105, 'unit_edges': 25,
                      'triple_count': 35, 'old_unit_circle_centers': 10,
                      'exclusions': len(exclusions), 'cube_edge_subsets': 4096,
                      'cube_chromatic_evaluations': chromatic}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
