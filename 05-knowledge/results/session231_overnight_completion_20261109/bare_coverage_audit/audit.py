import sys

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

if "--false-gate" in sys.argv:
    check(False, "deliberate-false gate sentinel")

import json, hashlib, time
from pathlib import Path
ROOT = Path(__file__).resolve().parent
N = 100000000
flags = bytearray((N + 15) // 16)
reps = 0
targets = set()
visits = 0
maxP = (N - 577 + 25) // 818
for p in range(maxP + 1):
    a = t = e = 0
    for bit in bin(p)[2:] if p else '':
        c = int(bit)
        a, t, e = (2 * a + c, t + c * (a + e), c)
    A = 1440 * a + 878
    T = 577 * a + 45 * t + 332 + 25 * e
    E = 14
    f = A - T + 45 - E
    if f <= N:
        check(f > 0 and f % 2, 'correctness gate failed')
        targets.add(f)
        reps += 1
    stack = [(A, T, E)]
    while stack:
        A, T, E = stack.pop()
        visits += 1
        if A > N - 45:
            continue
        f = 2 * A - T + 45
        if f <= N:
            check(f > 0 and f % 2, 'correctness gate failed')
            targets.add(f)
            reps += 1
        stack.append((2 * A, T, 0))
        stack.append((2 * A + 45, T + A + E, 45))
for f in targets:
    idx = f // 2
    flags[idx // 8] |= 1 << idx % 8
ordered = sorted(targets)
gaps = sorted(((b - a - 2, a, b) for a, b in zip(ordered, ordered[1:])), reverse=True)[:10]
result = {'method': 'integer prefix DP and pruned suffix transfer tree; no producer imports', 'cap': N, 'prefix_bound': maxP, 'tree_visits': visits, 'representations_k0': reps, 'count_le_1e7': sum((f <= 10000000 for f in targets)), 'count_le_1e8': len(targets), 'largest_internal_gaps': gaps, 'bitset_encoding': 'LSB first; bit floor(n/2) for positive odd n; byte length 6250000', 'bitset_sha256': hashlib.sha256(flags).hexdigest()}
(ROOT / 'union.bitset').write_bytes(flags)
(ROOT / 'RESULT.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
