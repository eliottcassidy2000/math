def _publication_require(value, message):
    if not value:
        raise RuntimeError(message)
if '--negative-check' in __import__('sys').argv:
    _publication_require(False, 'injected publication guard failure')
import json, math
from pathlib import Path
from collections import deque
P = Path(__file__).resolve().parent.parent / 'observer-profile-algebra'
O = Path(__file__).resolve().parent.parent / 'observer-profile-audit'
source = json.loads((P / 'walsh_certificate.json').read_text())
answers = []
for n in range(1, 7):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    indices = {p: k for k, p in enumerate(pairs)}
    N = 1 << len(pairs)
    gens = []
    for i in range(n - 1):
        p = list(range(n))
        p[i], p[i + 1] = (p[i + 1], p[i])
        dest = [indices[tuple(sorted((p[a], p[b])))] for a, b in pairs]
        table = [0] * N
        for x in range(1, N):
            bit = x & -x
            table[x] = table[x ^ bit] ^ 1 << dest[bit.bit_length() - 1]
        gens.append((table, 1 << indices[i, i + 1]))
    seen = set()
    rows = []
    allgraphs = 0
    for x in range(N):
        if x in seen:
            continue
        allgraphs += 1
        signs = {x: 1}
        q = deque([x])
        coherent = True
        while q:
            a = q.popleft()
            for table, bit in gens:
                b = table[a]
                s = signs[a] * (-1 if a & bit else 1)
                if b in signs:
                    coherent &= signs[b] == s
                else:
                    signs[b] = s
                    q.append(b)
        seen.update(signs)
        if coherent and x.bit_count() % 2 == 0:
            aut = math.factorial(n) // len(signs)
            rows.append({'mask': x, 'edges': [list(pairs[k]) for k in range(len(pairs)) if x >> k & 1], 'aut': aut, 'signed_order_sum': aut * sum(signs.values()), 'orbit_size': len(signs)})
    actual = source[n - 1]
    _publication_require(rows == actual['rows'], n)
    rank = sum((r['signed_order_sum'] != 0 for r in rows))
    _publication_require(rank == actual['rank'] and len(rows) == actual['dimension_converse_invariants'], 'independent audit assertion failed')
    answers.append({'n': n, 'graphs': allgraphs, 'dimension': len(rows), 'rank': rank, 'zero_masks': [r['mask'] for r in rows if not r['signed_order_sum']], 'all_walsh_rows_exact_match': True})
(O / 'independent_walsh_certificate.json').write_text(json.dumps(answers, indent=2) + '\n')
print(json.dumps(answers, indent=2))
