#!/usr/bin/env python3
"""Exact range certificates for the cost-carry univalence question."""
import hashlib
import json

checks = 0
rows = []
for L in range(1, 17):
    b = (3 ** L).bit_length() - 1
    Bmin = 3 ** L - 2 ** L
    max_width = 0
    for A in range(L, b + 1):
        prefix = [min((3 ** i).bit_length() - 1, A - L + i)
                  for i in range(L + 1)]
        checks += 1
        if not (prefix[0] == 0 and prefix[-1] == A and
                all(prefix[i+1] > prefix[i] for i in range(L)) and
                all(2 ** prefix[i] < 3 ** i for i in range(1, L + 1))):
            raise ArithmeticError('upper-prefix envelope is not an admissible word')
        Bmax = sum(3 ** (L - 1 - i) * 2 ** prefix[i] for i in range(L))
        max_width = max(max_width, Bmax - Bmin)
    sufficient = max_width < 2 * 3 ** L
    checks += 1
    if sufficient != (L <= 12):
        raise ArithmeticError('range threshold mismatch in the declared universe')
    rows.append({'length': L, 'total_cost_ceiling': b, 'B_min': Bmin,
                 'B_max_at_max_cost': Bmin + max_width,
                 'maximum_interval_width': max_width,
                 'minimum_nonzero_same_cost_carry_difference': 2 * 3 ** L,
                 'range_forces_univalence': sufficient})
result = {'status': 'FINITE-EXACT range certificates; no claim about univalence beyond the range bound',
          'checks': checks, 'lengths': [1, 16], 'rows': rows,
          'scope': 'The general extremal formula is proved in the companion roadmap.'}
result['payload_sha256'] = hashlib.sha256(json.dumps(result, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
print(json.dumps(result, indent=2, sort_keys=True))
