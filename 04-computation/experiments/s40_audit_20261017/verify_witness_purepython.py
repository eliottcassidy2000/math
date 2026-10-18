"""Pure-Python (no numba) recomputation of witness-order excesses: Legendre symbols by pow(), pair loop.
Usage: python verify_witness_purepython.py JSON [JSON ...]  (rows with keys p and order or best_order)"""
import sys, json
def excess(order, p):
    leg = [0] + [1 if pow(a, (p - 1) // 2, p) == 1 else -1 for a in range(1, p)]
    s = 0
    for i in range(p):
        u = order[i]
        for j in range(i + 1, p):
            s += leg[(order[j] - u) % p]
    return s
for fn in sys.argv[1:]:
    for r in json.load(open(fn)):
        o = r.get('order') or r.get('best_order')
        p = r['p']
        assert sorted(o) == list(range(p))
        e = excess(o, p); nat = excess(list(range(p)), p)
        print(fn.split('/')[-1], 'p=%d witness excess %d natural %d beaten %s' % (p, e, nat, e > nat))
