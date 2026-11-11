"""Complete fixed transitive-module return classification, with exact controls.

Run from any directory. The inherited order<=8 corpus is certified here by
degree-permutation canonicalization, automorphism mass, and fresh path DP.
The module-response identity uses the inherited positive placeholder DP;
literal inflated graphs are independently counted. No assert statements.
"""
from collections import Counter, defaultdict
from itertools import permutations, product
from pathlib import Path
import gzip, hashlib, json, math, sys, time

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent
OLD = RESULTS / 'chanpak_five_core_bridges_20261110/alternatives'
sys.path.insert(0, str(OLD))
import check as inherited
import two_birth


def need(value, message):
    if not value:
        raise RuntimeError(message)


def H(rows):
    """Fresh last-vertex recurrence, opposite traversal to inherited.dp."""
    n = len(rows)
    if not n:
        return 1
    end = [[0] * n for _ in range(1 << n)]
    back = [sum((rows[u] >> v & 1) << u for u in range(n)) for v in range(n)]
    for mask in range(1, 1 << n):
        todo = mask
        while todo:
            bit = todo & -todo
            todo -= bit
            v = bit.bit_length() - 1
            if mask == bit:
                end[mask][v] = 1
                continue
            prev = mask ^ bit
            avail = back[v] & prev
            total = 0
            while avail:
                b = avail & -avail
                avail -= b
                total += end[prev][b.bit_length() - 1]
            end[mask][v] = total
    return sum(end[-1])


def strong(rows):
    if not rows:
        return False
    for reverse in (False, True):
        reached = 1
        while True:
            nxt = reached
            for v in range(len(rows)):
                if reached >> v & 1:
                    nxt |= rows[v] if not reverse else sum((rows[u] >> v & 1) << u for u in range(len(rows)))
            if nxt == reached:
                break
            reached = nxt
        if reached != (1 << len(rows)) - 1:
            return False
    return True


def canonical_and_aut(rows):
    """Each lexicographically minimal adjacency code occurs |Aut| times."""
    buckets = defaultdict(list)
    for v, r in enumerate(rows):
        buckets[r.bit_count()].append(v)
    best, multiplicity = None, 0
    for pieces in product(*(permutations(vs) for _, vs in sorted(buckets.items()))):
        order = sum(pieces, ())
        code = tuple(sum((rows[u] >> v & 1) << j for j, v in enumerate(order)) for u in order)
        if best is None or code < best:
            best, multiplicity = code, 1
        elif code == best:
            multiplicity += 1
    return best, multiplicity


def converse(rows):
    return tuple(((1 << len(rows)) - 1) ^ r ^ (1 << v) for v, r in enumerate(rows))


def pointed_iso(rows, vertex, other, point):
    n = len(rows)
    if n != len(other):
        return None
    domains = [[w for w in range(n) if other[w].bit_count() == rows[v].bit_count()
                and (w == point) == (v == vertex)] for v in range(n)]
    for p in product(*domains):
        if len(set(p)) != n:
            continue
        if all((rows[u] >> v & 1) == (other[p[u]] >> p[v] & 1) for u in range(n) for v in range(n)):
            return p
    return None


def frontier(source, bound):
    current = [{'rows': source, 'path': []}]
    layers = []
    while current:
        children, cuts = [], 0
        for item in current:
            values = two_birth.all_birth_counts(item['rows'])
            cuts += len(values)
            for out, h in enumerate(values):
                if h > bound:
                    continue
                child = inherited.attach(item['rows'], out)
                if not strong(child):
                    continue
                need(H(child) == h, 'accepted frontier child fresh DP')
                children.append({'rows': child, 'h': h, 'path': item['path'] + [out]})
        layers.append({'cuts': cuts, 'count': len(children),
                       'histogram': dict(sorted(Counter(c['h'] for c in children).items())),
                       'children': children})
        current = children
    return layers


def main():
    started = time.perf_counter()
    source = RESULTS / 'metagraph_aspect_sequences_20261111/records.jsonl.gz'
    records = [json.loads(line) for line in gzip.open(source, 'rt')]
    need(len(records) == 7412, 'declared complete order1..8 universe')
    audit = []
    candidates = []
    for n in range(1, 9):
        group = [r for r in records if r['n'] == n]
        seen, mass = set(), 0
        for rec in group:
            rows = tuple(rec['rows'])
            need(all(not (rows[v] >> v & 1) for v in range(n)), 'no loops')
            need(all(((rows[v] >> w & 1) + (rows[w] >> v & 1)) == 1
                     for v in range(n) for w in range(v + 1, n)), 'tournament type')
            code, aut = canonical_and_aut(rows)
            need(code not in seen, 'ordinary class distinctness')
            seen.add(code)
            need(aut == rec['aut'], 'fresh automorphism order')
            mass += math.factorial(n) // aut
            h, sc = H(rows), strong(rows)
            need(h == rec['H'] and sc == bool(rec['strong']), 'fresh H and strongness')
            if h == 73 and sc:
                candidates.append(rec)
        need(mass == 2 ** (n * (n - 1) // 2), 'full labelled mass certifies completeness')
        audit.append({'n': n, 'classes': len(group), 'labelled_mass': mass})
        print('Complete universe order', n, 'PASS', flush=True)

    refs = {'A': inherited.ready_source((2, 3, 1), 1),
            'B': inherited.ready_source((2, 2, 1), 1),
            'C': (126, 236, 120, 16, 226, 72, 8, 109)}
    hits = []
    for rec in candidates:
        rows = tuple(rec['rows'])
        for v in range(len(rows)):
            runs = two_birth.module_runs(rows, v)
            if runs[1] != 73 or runs[2] != 34 or any(runs[3:]):
                continue
            identification = []
            for name, ref in refs.items():
                for side, other in [('forward', ref), ('converse', converse(ref))]:
                    p = pointed_iso(rows, v, other, 0)
                    if p is not None:
                        identification.append({'class': name, 'side': side, 'map': p})
            need(len(identification) == 1, 'unique old pointed family identification')
            counts = [H(two_birth.inflate(rows, v, k)) for k in range(1, 5)]
            need(counts == [34 * 2 ** k + 5 for k in range(1, 5)], 'literal uniform controls')
            hits.append({'n': rec['n'], 'rep': rec['rep'], 'rows': rows,
                         'module': v, 'run_coefficients': runs,
                         'identification': identification[0], 'literal_k1_to4': counts})
    need(len(candidates) == 19 and len(hits) == 6, 'complete template counts')
    # A deliberately close wrong template passes at k2 but not the identity.
    hostile = (254, 508, 504, 240, 480, 448, 384, 256, 9)
    need(H(hostile) == 141, 'positive141 hostile count')
    hostile_runs = [two_birth.module_runs(hostile, v) for v in range(len(hostile))]
    need(not any(c[2] == 68 and not any(c[3:]) for c in hostile_runs),
         'finite successful reset is not an affine dyadic module family starting k2')

    frontiers = {}
    for name, ref in refs.items():
        T = two_birth.inflate(ref, 0, 2)
        need(H(T) == 141 and strong(T), 'frontier source')
        layers = frontier(T, 425)
        need(all(c['h'] != 425 for layer in layers for c in layer['children']), 'no425 at any frontier layer')
        frontiers[name] = {'source': T, 'layers': layers}
        print('All-order425 source frontier', name, [x['count'] for x in layers], 'PASS', flush=True)
    report = {'status': 'PASS', 'universe': audit,
              'input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'strong_H73_by_order': dict(sorted(Counter(r['n'] for r in candidates).items())),
              'strong_H73_bases': [{'n': r['n'], 'rep': r['rep'], 'rows': r['rows']} for r in candidates],
              'templates': hits, 'finite_reset_hostile_run_coefficients': hostile_runs,
              'frontiers425': frontiers,
              'scope': 'Complete single transitive-module families at H(k)=34*2^k+5; no arbitrary variable-incidence reset exclusion.'}
    (HERE / 'module_classification.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print('PASS; seconds', round(time.perf_counter() - started, 3), flush=True)


if __name__ == '__main__':
    main()
