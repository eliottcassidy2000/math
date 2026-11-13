"""Independent literal-path audit of the root's all-order induced-count unions."""
from collections import Counter
from functools import lru_cache
from pathlib import Path
import csv
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PACKET = HERE.parent
BANK = PACKET.parent / 'actual_hamilton_multiplicities_20261021/session24-multiplicities'

def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)

def walk_inventory(rows):
    # Enumerate individual simple paths, with no subset recurrence or endpoint sums.
    counts = [0] * (1 << len(rows))
    counts[0] = 1
    walks = 0
    def visit(last, used):
        nonlocal walks
        counts[used] += 1
        walks += 1
        choices = rows[last] & ~used
        while choices:
            bit = choices & -choices
            choices -= bit
            visit(bit.bit_length() - 1, used | bit)
    for start in range(len(rows)):
        visit(start, 1 << start)
    return counts, walks

def strong(rows):
    n = len(rows)
    for start in range(n):
        reached = 1 << start
        while True:
            new = reached
            for v in range(n):
                if reached >> v & 1:
                    new |= rows[v]
            if new == reached:
                break
            reached = new
        if reached != (1 << n) - 1:
            return False
    return True

need(5 ** 8 > 73 ** 3, 'Busch all-order cutoff')
records = json.loads((BANK / 'census.json').read_text(encoding='utf-8'))['records']
cores = [(record['H'], tuple(record['adjacency_outmasks'])) for record in records
         if record['strong'] and record['H'] <= 73]
with (BANK / 'n8_classes.csv').open(newline='', encoding='utf-8') as stream:
    for record in csv.DictReader(stream):
        if record['strong'] == '1' and int(record['H']) <= 73:
            code = int(record['code'])
            cores.append((int(record['H']), tuple((code // 256 ** i) % 256 for i in range(8))))
atom_values = {}
subset_count = 0
literal_paths = 0
by_order = Counter()
for h, rows in cores:
    need(strong(rows), 'strong flag')
    values, paths = walk_inventory(rows)
    need(values[-1] == h, 'literal full count')
    atom_values.setdefault(h, set()).update(values)
    subset_count += len(values)
    literal_paths += paths
    by_order[len(rows)] += 1

@lru_cache(None)
def factors(target, least=3):
    # All nondecreasing multisets of nontrivial SCC counts, independently of
    # the producer's ascending pair-convolution recurrence.
    if target == 1:
        return frozenset((1,))
    result = set()
    for d in sorted(atom_values):
        if d < least or target % d:
            continue
        for a in atom_values[d]:
            for b in factors(target // d, d):
                result.add(a * b)
    return frozenset(result)

actual = {str(h): sorted(factors(h)) for h in range(1, 74, 2)}
expected = json.loads((PACKET / 'containment_edges.json').read_text(encoding='utf-8'))
need(actual == expected['all_order_induced_count_unions'], 'every induced-count union')

def raw_parent(h):
    # Derive via the ordinary binary-heap index m=(h-1)/2.
    return 2 * (((h - 1) // 2) // 2) + 1

failures = []
for h in range(3, 74, 2):
    if h in (7, 21):
        need(actual[str(h)] == [], 'forbidden target')
        continue
    parent = raw_parent(h)
    while parent in (7, 21):
        parent = raw_parent(parent)
    need(parent < h and parent % 2 == 1, 'strict odd descent')
    if parent not in actual[str(h)]:
        failures.append([parent, h])
need(failures == expected['missing_arithmetic_tree_edges'], 'all missing edges')
need(len(cores) == 148 and subset_count == 23512, 'universe sizes')
sources = [BANK / 'census.json', BANK / 'n8_classes.csv', BANK / 'PROOF.md',
           PACKET / 'containment_edges.py', PACKET / 'containment_edges.json',
           PACKET / 'tree_blueprint.py',
           PACKET.parent / 'metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md']
out = dict(status='PASS independent literal-path and factor-multiset audit',
           scope='All-order conclusion inherits complete strong census and Busch bound; class census not regenerated',
           cores=len(cores), core_orders=dict(sorted(by_order.items())),
           subset_counts=subset_count, individual_simple_paths_enumerated=literal_paths,
           all37_count_unions_match=True, missing_edges=failures,
           source_pins=[dict(path=str(p.relative_to(ROOT)).replace('\\', '/'),
                             sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources])
suffix = '.optimized' if sys.flags.optimize else ''
(HERE / ('containment_audit' + suffix + '.json')).write_text(json.dumps(out, indent=2) + '\n', encoding='utf-8')
print(json.dumps({k: v for k, v in out.items() if k != 'source_pins'}))
