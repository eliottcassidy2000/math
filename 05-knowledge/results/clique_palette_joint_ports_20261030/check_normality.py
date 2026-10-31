"""Exact independent validation of the Paley palette normality certificate."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0
def need(ok, text):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(text)

def determinant(matrix):
    """Fraction Gaussian elimination, independent of the producer's Bareiss routine."""
    a = [[F(x) for x in row] for row in matrix]
    out = F(1)
    for col in range(len(a)):
        pivot = next((i for i in range(col, len(a)) if a[i][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[pivot], a[col] = a[col], a[pivot]
            out = -out
        scale = a[col][col]
        out *= scale
        for i in range(col+1, len(a)):
            ratio = a[i][col]/scale
            for j in range(col+1, len(a)):
                a[i][j] -= ratio*a[col][j]
    need(out.denominator == 1, 'integer determinant')
    return int(out)

def geometry():
    cert = json.loads((HERE/'normality_certificate.json').read_text())
    rows = [sum(1 << j for j in range(7) if (j-i) % 7 in (1,2,4)) for i in range(7)]
    cycles = []
    for triple in combinations(range(7), 3):
        mask = sum(1 << v for v in triple)
        if all((rows[v] & mask).bit_count() == 1 for v in triple):
            cycles.append(mask)
    expected = [s for s in range(128) if all(s & c != c for c in cycles)]
    verts = cert['vertex_masks']
    need(verts == expected and len(verts) == 50 and len(cycles) == 14, 'complete literal transitive subsets')
    need(max(s.bit_count() for s in verts) == 3, 'no transitive four-set')
    gens = [[1]+[(v >> i) & 1 for i in range(7)] for v in verts]
    inequalities = []
    for i in range(7):
        for sign, rhs in ((-1, 0), (1, 1)):
            a = [0]*7
            a[i] = sign
            inequalities.append((a, rhs))
    inequalities += [([(c >> i) & 1 for i in range(7)], 2) for c in cycles]
    inequalities.append(([1]*7, 3))
    for v in gens:
        need(all(sum(x*y for x,y in zip(a, v[1:])) <= rhs for a,rhs in inequalities), 'all vertices feasible')
    simplices = [tuple(s) for s in cert['simplices']]
    need(len(simplices) == len(set(simplices)), 'no duplicate simplex')
    faces, determinants = defaultdict(list), []
    det_hist = Counter()
    for j, s in enumerate(simplices):
        need(len(s) == 8 and s == tuple(sorted(set(s))) and min(s) >= 0 and max(s) < len(verts), 'simplex vertex indices')
        matrix = [[gens[v][i] for v in s] for i in range(8)]
        determinant_value = determinant(matrix)
        need(determinant_value != 0, 'full dimensional simplex')
        determinants.append(determinant_value)
        det_hist[abs(determinant_value)] += 1
        for pos, vertex in enumerate(s):
            face = s[:pos]+s[pos+1:]
            # Orientation of [face in sorted order, opposite vertex last].
            side = determinant_value*(-1)**(7-pos)
            faces[face].append((j, side))
    boundary = 0
    facet_witnesses = set()
    for face, incidents in faces.items():
        if len(incidents) == 2:
            need(incidents[0][1]*incidents[1][1] < 0, 'paired facets have opposite sides')
        else:
            need(len(incidents) == 1, 'facet incidence count')
            need(any(all(sum(x*y for x,y in zip(a, gens[v][1:])) == rhs for v in face)
                     for a,rhs in inequalities), 'unpaired facet lies on a declared boundary')
            facet_witnesses.update(i for i,(a,rhs) in enumerate(inequalities)
                                  if all(sum(x*y for x,y in zip(a,gens[v][1:]))==rhs for v in face))
            boundary += 1
    need(len(facet_witnesses)==len(inequalities),'every declared inequality has a full facet witness')
    records = {(r['simplex'], tuple(r['residue'])): r for r in cert['fundamental_points']}
    need(len(records) == len(cert['fundamental_points']), 'no repeated fundamental point key')
    used, heights = set(), Counter()
    for j, s in enumerate(simplices):
        p = abs(determinants[j])
        need(p in (1,2,3), 'declared determinant range')
        if p == 1:
            continue
        matrix = [[gens[v][i] for v in s] for i in range(8)]
        residues = []
        for r in product(range(p), repeat=8):
            numerator = [sum(x*y for x,y in zip(row, r)) for row in matrix]
            if any(x % p for x in numerator):
                continue
            residues.append(r)
            if not any(r):
                continue
            key = (j, r)
            need(key in records, 'every nonzero fundamental class supplied')
            record = records[key]
            point = [x//p for x in numerator]
            need(record['point'] == point, 'exact fundamental point coordinates')
            colors = record['colors']
            need(len(colors) == point[0] and all(v in verts for v in colors), 'actual legal colors at exact height')
            need([sum(v >> i & 1 for v in colors) for i in range(7)] == point[1:], 'exact demand decomposition')
            used.add(key)
            heights[point[0]] += 1
        need(len(residues) == p, 'all lattice classes in fundamental parallelepiped')
    need(used == set(records), 'no unverified extra records')
    return {'vertices': len(verts), 'inequalities': len(inequalities), 'simplices': len(simplices),
            'determinants': dict(sorted(det_hist.items())), 'boundary_facets': boundary,
            'internal_facets': len(faces)-boundary, 'fundamental_points': len(records),
            'fundamental_heights': dict(sorted(heights.items())),
            'all_size_consequence': 'The 29-inequality polytope equals the transitive-set polytope and has integer decomposition.'}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'geometry': geometry(), 'checks': CHECKS,
              'scope': 'Exact certificate plus displayed finite-to-all-size proof; no floating point or producer imports.'}
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text, encoding='utf-8', newline='\n')
    else:
        print(text, end='')

if __name__ == '__main__':
    main()
