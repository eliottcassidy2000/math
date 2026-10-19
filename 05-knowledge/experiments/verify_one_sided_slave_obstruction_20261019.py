#!/usr/bin/env python3
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import json

def count(p):
    n = len(p)

    @lru_cache(None)
    def f(m):
        if m == (1 << n) - 1:
            return 1
        return sum((f(m | 1 << i) for i in range(n) if not m >> i & 1 and (not p[i] & ~m)))
    return f(0)

def main():
    rows = []
    for r in range(1, 21):
        p = [0, 0, 0] + [1 | sum((1 << 3 + j for j in range(i))) for i in range(r)]
        p[-1] |= 4
        E = count(p)
        q = p[:]
        q[2] |= 1
        Nxz = count(q)
        q = p[:]
        q[1] |= 4
        Nzy = count(q)
        if not E == (r + 1) * (r + 3):
            raise AssertionError('exact certificate check failed')
        if not Nxz == r * (r + 3):
            raise AssertionError('exact certificate check failed')
        if not Nzy == (r + 1) * (r + 4) // 2:
            raise AssertionError('exact certificate check failed')
        if not 4 * E - 3 * (Nxz + Nzy) == -(r + 4) * (r - 3) // 2:
            raise AssertionError('exact certificate check failed')
        rows.append({'upper_slave_count': r, 'extensions': E, 'x_before_z': Nxz, 'z_before_y': Nzy, 'aggregate': str(Fraction(Nxz + Nzy, E)), 'deficit': 4 * E - 3 * (Nxz + Nzy)})
    result = {'status': 'PASS', 'actual_posets': rows}
    out = Path(__file__).with_name('one_sided_slave_aggregate_obstruction_audit_20261019.json')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
