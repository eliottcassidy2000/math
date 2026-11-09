def require(condition, message):
    if not condition:
        raise RuntimeError(message)
from itertools import product
import json
out = {}
for m in range(2, 7):
    vals = set()
    for q in product((0, 1), repeat=m - 2):
        w = (1,) + q + (0,)
        c = sum((1 if j - i == 1 else 2 ** (j - i - 2) for i in range(m) for j in range(i + 1, m) if w[i] and (not w[j])))
        vals.add(c)
    require(3 not in vals and 10 not in vals, 'check_transitive.py:10: validation failed')
    out[m] = sorted(vals)
require(2 ** (7 - 3) > 10, 'check_transitive.py:12: validation failed')
print(json.dumps(out, indent=2))
