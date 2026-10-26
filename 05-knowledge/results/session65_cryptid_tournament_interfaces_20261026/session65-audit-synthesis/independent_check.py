"""Independent finite audit; no producer imports or whole-tournament census."""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
BASE = HERE.parent

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
for row in pins['files']:
    path = BASE / row['path']
    need(path.stat().st_size == row['bytes'], 'pinned byte count: ' + row['path'])
    need(digest(path) == row['sha256'], 'pinned source: ' + row['path'])

# Revalidate existing audit seals, not their long computational runs.
seal_rows = 0
def verify_one(root, rel, expected):
    global seal_rows
    need(digest(root / rel) == expected, 'inherited seal: ' + str(root / rel))
    seal_rows += 1

def verify_sums(folder):
    root = BASE / folder
    for line in (root / 'SHA256SUMS').read_text().splitlines():
        expected, rel = line.split(maxsplit=1)
        verify_one(root, rel.lstrip('*'), expected)

for folder in ('session65-audit-cryptids', 'session65-audit-goldbach',
               'session65-cryptids', 'session65-twentyfive',
               'session65-root/multiperfect', 'session65-root/fifteen-kernel'):
    verify_sums(folder)
for folder, name, fields in (
    ('session65-audit-fifteens', 'audit-manifest.json',
      (('audit_sha256', 'session65-audit-fifteens'), ('frozen_lane_sha256', 'session65-fifteens'))),
    ('session65-audit-weakcollatz', 'manifest.json',
      (('audit_sha256', 'session65-audit-weakcollatz'), ('target_sha256', 'session65-weakcollatz'))),
    ('session65-fifteens', 'source-manifest.json',
      (('inputs_sha256', 'session65-fifteens'), ('frozen_artifacts_sha256', 'session65-fifteens')))):
    manifest = json.loads((BASE / folder / name).read_text())
    for field, root in fields:
        for rel, expected in manifest[field].items():
            verify_one(BASE / root, rel, expected)
for folder, filename, field, root in (
    ('session65-audit-reductions', 'MANIFEST.json', 'files', 'session65-audit-reductions'),
    ('session65-reductions', 'MANIFEST.json', 'files', 'session65-reductions'),
    ('session65-audit-goldbach', 'MANIFEST.json', 'payload_files', 'session65-audit-goldbach'),
    ('session65-audit-goldbach', 'INPUTS.json', 'files', '')):
    manifest = json.loads((BASE / folder / filename).read_text())
    for row in manifest[field]:
        verify_one(BASE / root, row['path'], row['sha256'])
need(digest(BASE / 'session65-cryptids/REPORT.md') ==
     json.loads((BASE / 'session65-audit-cryptids/AMENDMENT.json').read_text())['new_sha256'],
     'cryptids audited report amendment')
need(digest(BASE / 'session65-reductions/REPORT.md') ==
     json.loads((BASE / 'session65-audit-reductions/MANIFEST.json').read_text())['producer_report_sha256'],
     'reductions audited report identity')

# Generate all matchings by filtering every three-element subset of the full
# duad carrier, not by the producer's recursive partner construction.
D = list(combinations(range(6), 2))
M = [es for es in combinations(D, 3)
     if sorted(x for e in es for x in e) == list(range(6))]
need(len(D) == len(set(D)) == len(M) == len(set(M)) == 15, 'complete carriers')
need(len(list(combinations(D, 3))) == 455, 'complete matching candidate universe')

# Independent finite symplectic dictionary: canonical representatives for
# even-subset classes modulo the full six-set; duads give all nonzero classes.
canon = lambda mask: min(mask, mask ^ 63)
classes = {canon(mask) for mask in range(64) if mask.bit_count() % 2 == 0}
duad_masks = [sum(1 << v for v in d) for d in D]
need({canon(x) for x in duad_masks} == classes - {0}, 'all nonzero quotient vectors')
planes = []
for triple in combinations(duad_masks, 3):
    if canon(triple[0] ^ triple[1] ^ triple[2]) == 0 and all(
            (a & b).bit_count() % 2 == 0 for a, b in combinations(triple, 2)):
        planes.append(frozenset(canon(x) for x in triple))
matching_planes = [frozenset(canon(sum(1 << v for v in e)) for e in m) for m in M]
need(len(planes) == 15 and set(planes) == set(matching_planes), 'syntheme Lagrangian carrier')
need(sum(canon(d) in p for d in duad_masks for p in planes) == 45, 'point plane incidence')

E = list(combinations(range(15), 2))
Eindex = {e: i for i, e in enumerate(E)}
receipt = json.loads((BASE / 'session65-root/fifteen-kernel/normal.json').read_text())
need((BASE / 'session65-root/fifteen-kernel/normal.json').read_bytes() ==
     (BASE / 'session65-root/fifteen-kernel/optimized.json').read_bytes(),
     'root receipt modes identical')

def convolution(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result

def defects(rows, perm):
    need(len(rows) == 15, 'literal 15 rows')
    need(all(type(row) is int and 0 <= row < (1 << 15) for row in rows), 'literal row range')
    need(all((rows[i] >> i) & 1 == 0 for i in range(15)), 'no loops')
    need(all(((rows[a] >> b) & 1) + ((rows[b] >> a) & 1) == 1 for a, b in E),
         'complete antisymmetric tournament')
    return [e for e in E if ((rows[e[0]] >> e[1]) & 1) !=
            ((rows[perm[e[0]]] >> perm[e[1]]) & 1)]

def audit(name, carrier, label_perm, expected):
    moved = []
    for v in carrier:
        if carrier is D:
            image = tuple(sorted(label_perm[x] for x in v))
        else:
            image = tuple(sorted(tuple(sorted(label_perm[x] for x in d)) for d in v))
        moved.append(carrier.index(image))
    need(sorted(moved) == list(range(15)), 'actual permutation')
    need(all(moved[moved[i]] == i for i in range(15)), 'actual involution')
    f = sum(moved[i] == i for i in range(15))
    b = (15-f)//2
    need((f, b) == expected, 'fixed-point character')
    # Compare each edge bit with its forward image bit, reversing the image
    # bit when image endpoints occur in descending canonical order.
    target, reverse = [], []
    for a, c in E:
        u, v = moved[a], moved[c]
        target.append(Eindex[tuple(sorted((u, v)))])
        reverse.append(int(u > v))
    unseen, orbit_data, polynomial, global_bits = set(range(105)), [], [1], {}
    local_assignments = 0
    while unseen:
        start = min(unseen)
        orbit, e = [], start
        while e not in orbit:
            orbit.append(e)
            e = target[e]
        need(e == start, 'closed permutation orbit')
        unseen.difference_update(orbit)
        L = len(orbit)
        sign = -1 if sum(reverse[e] for e in orbit) % 2 else 1
        local = [0] * (L+1)
        best = None
        for bits in product((0, 1), repeat=L):
            local_assignments += 1
            assignment = dict(zip(orbit, bits))
            d = sum(assignment[e] != (assignment[target[e]] ^ reverse[e]) for e in orbit)
            local[d] += 1
            if best is None or d < best[0]:
                best = (d, assignment)
        global_bits.update(best[1])
        polynomial = convolution(polynomial, local)
        orbit_data.append({'edges': [list(E[e]) for e in orbit], 'length': L,
                           'sign': sign, 'local_coefficients': local})
    m = b*(f+b-1)
    exponent = comb(f,2)+b+m
    expected_poly = [0]*106
    for j in range(m+1):
        expected_poly[b+2*j] = (1 << exponent)*comb(m,j)
    polynomial += [0]*(106-len(polynomial))
    need(polynomial == expected_poly, 'complete exact polynomial coefficients')
    need(sum(polynomial) == 1 << 105, 'all orientations counted algebraically')
    counts = Counter((o['length'], o['sign']) for o in orbit_data)
    need(counts == {(1,1): comb(f,2), (1,-1): b, (2,1): m}, 'all pair orbit types')
    rows = [0]*15
    for e, bit in global_bits.items():
        a, c = E[e]
        rows[a if bit else c] |= 1 << (c if bit else a)
    independent_defects = defects(rows, moved)
    need(len(independent_defects) == b, 'independent literal minimizer')
    r = receipt[name]
    need(r['permutation'] == moved, 'receipt actual carrier action')
    need((r['fixed_vertices'], r['swapped_pairs'], r['minimum_defect']) == (f,b,b),
         'receipt fixed points and minimum')
    root_defects = defects(r['one_minimizer_rows'], moved)
    need(len(root_defects) == b, 'root literal orientation attains minimum')
    need([list(e) for e in root_defects] == r['forced_defect_pairs'], 'root forced defect list')
    need(all(moved[a] == c and moved[c] == a for a, c in root_defects), 'forced swapped pairs')
    expected_text = f'2^{exponent} z^{b} (1+z^2)^{m}'
    need(r['polynomial'] == expected_text, 'root polynomial text')
    need({(x['length'], x['sign']): x['count'] for x in r['pair_orbits']} == counts,
         'root orbit inventory')
    # One edge flip in a positive doubleton must add two failures; neither
    # ordinary orientation complementation nor swapping a forced edge fixes it.
    doubleton = next(o for o in orbit_data if o['length'] == 2)
    a, c = doubleton['edges'][0]
    mutant = rows.copy(); mutant[a] ^= 1 << c; mutant[c] ^= 1 << a
    need(len(defects(mutant, moved)) == b+2, 'doubleton one-arc hostile')
    return {'fixed_vertices': f, 'swapped_vertex_pairs': b, 'unordered_pairs': 105,
            'unordered_pair_orbits': len(orbit_data), 'local_assignments_checked': local_assignments,
            'polynomial': expected_text, 'polynomial_coefficients_0_to_105': polynomial,
            'actual_permutation': moved, 'independent_minimizer_rows': rows,
            'root_literal_minimizer_valid': True, 'minimum_defect': b,
            'all_pair_orbits': orbit_data, 'doubleton_flip_defect': b+2}

results = {
    'duads_under_transposition': audit('duads_under_transposition', D, [1,0,2,3,4,5], (7,4)),
    'synthemes_under_transposition': audit('synthemes_under_transposition', M, [1,0,2,3,4,5], (3,6)),
    'duads_under_triple_transposition': audit('duads_under_triple_transposition', D, [1,0,3,2,5,4], (3,6))
}
print(json.dumps({'status': 'PASS', 'source_pins_checked': len(pins['files']), 'inherited_seal_rows_checked': seal_rows,
                  'carriers': {'duads': D, 'synthemes': M, 'matching_candidates': 455,
                               'symplectic_classes': 16, 'lagrangians': 15, 'flags': 45},
                  'results': results,
                  'scope': 'Complete edge-orbit factorization and literal minimizers, not enumeration of 2^105 tournaments; no full outer automorphism or arithmetic Richelot proof.'},
                 sort_keys=True, indent=2))
