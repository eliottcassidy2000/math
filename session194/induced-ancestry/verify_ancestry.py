#!/usr/bin/env python3
"""Independent literal positive controls for the induced ancestry theorem."""
import hashlib
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def need(b, message):
    if not b:
        raise ValueError(message)

def literal(rows):
    h = 0
    for word in itertools.permutations(range(len(rows))):
        if all(rows[u] >> v & 1 for u, v in zip(word, word[1:])):
            h += 1
    return h

def birth(rows, mask):
    n = len(rows)
    return [row | (0 if mask >> i & 1 else 1 << n)
            for i, row in enumerate(rows)] + [mask]

def strong(rows):
    n = len(rows)
    for i in range(n):
        reached = {i}
        while True:
            new = reached | {j for u in reached for j in range(n) if rows[u] >> j & 1}
            if new == reached:
                break
            reached = new
        if len(reached) != n:
            return False
    return True

def main():
    receipt = json.loads((ROOT / 'ancestry_full.json').read_text())
    need(receipt['mode'] == 'full' and receipt['status'] == 'PASS', 'full receipt required')
    expected = dict(source_first_hits=0, source_second_hits=10, source_direct_checks=1696,
                    first_tested=1280, first_retained=272, second_tested=69632,
                    second_retained=1220, third_tested=624640,
                    final_direct_checks=624640, target_hits=[0, 0, 0])
    for key, value in expected.items():
        need(receipt[key] == value, 'exact finite universe: ' + key)
    need(len(receipt['seeds']) == 10, 'ten source records')
    D = [2, 12, 9, 17, 7]
    need(literal(D) == 13 and strong(D), 'strong H13 core')
    seeds = []
    for record in receipt['seeds']:
        r = birth(birth(D, record['first_mask']), record['second_mask'])
        need(r == record['rows'], 'literal source birth masks')
        need(literal(r) == 41 and strong(r), 'literal strong H41 source')
        seeds.append(r)
    pre = birth(D, 31)
    A = birth(pre, 47)
    B = birth(pre, 55)
    a69 = birth(A, 31)
    a125 = birth(a69, 31)
    b131 = birth(B, 10)
    b133 = birth(B, 29)
    skip29 = birth(D, 29)
    skip75 = birth(skip29, 7)
    skip127 = birth(skip75, 1)
    controls = [('D', D, 13, True), ('source', pre, 13, False),
                ('A', A, 41, True), ('B', B, 43, True),
                ('A69', a69, 69, True), ('A125', a125, 125, True),
                ('B131', b131, 131, True), ('B133', b133, 133, True),
                ('bypass29', skip29, 29, True), ('bypass75', skip75, 75, True),
                ('bypass127', skip127, 127, True)]
    for name, rows, h, is_strong in controls:
        need(literal(rows) == h, 'literal control ' + name)
        need(strong(rows) == is_strong, 'strongness ' + name)
    for exponent, h in [(5, 13), (7, 41), (10, 127)]:
        need(5 ** exponent > h ** 3, 'strict order cutoff')
    out = dict(status='PASS', finite_receipt=expected,
               literal_H41_source_records=len(seeds),
               positive_controls=[dict(name=name, rows=rows, H=h, strong=s)
                                  for name, rows, h, s in controls],
               cited_classification='Unique strong H13 core D; inherited, not re-proved here',
               cited_order_bound='Busch H(U)>=5^((|U|-1)/3) for strong U')
    raw = (json.dumps(out, indent=2, sort_keys=True) + '\n').encode()
    (ROOT / 'ancestry_controls.json').write_bytes(raw)
    print(json.dumps(dict(status='PASS', literal_source_records=10,
                         literal_positive_controls=len(controls),
                         receipt_sha256=hashlib.sha256(raw).hexdigest()), sort_keys=True))

if __name__ == '__main__':
    main()
