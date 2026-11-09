import sys

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

if "--false-gate" in sys.argv:
    check(False, "deliberate-false gate sentinel")

import json, hashlib
from pathlib import Path
import runpy
transfer = runpy.run_path(str(Path(__file__).resolve().parent / 'controls.py'))['transfer']
ROOT = Path(__file__).resolve().parent
producer = ROOT.parent / 'bare_coverage'
if not producer.is_dir():
    producer = ROOT.parent / 'session231-bare-coverage'
r = json.loads((producer / 'RESULT.json').read_text())
flags = (ROOT / 'union.bitset').read_bytes()
X = tuple(map(int, '101000111000'))
Y = tuple(map(int, '111011100110'))

def bits(p, b=None):
    return tuple(map(int, f'{p:0{b}b}')) if b else tuple(map(int, bin(p)[2:])) if p else ()
checked = 0
for n, w in r['witnesses'].items():
    k, j, P, R, b, order = w
    n = int(n)
    p = bits(P) + X * j
    s = Y * (k - j) + bits(R, b)
    check(transfer(p, s) == n, (n, w))
    check(order == len(p) + len(s) + 9, 'correctness gate failed')
    check(flags[n // 2 // 8] >> n // 2 % 8 & 1, 'correctness gate failed')
    if k:
        checked += 1
for k in range(4):
    for j in range(k + 1):
        z = (13617152 * 4096 ** k + 1602510 * k + 717259) // 24843
        check(transfer(X * j, Y * (k - j)) == z - 2 * j, 'correctness gate failed')
for p, s, expected in [((), (), [2245227, 2245225]), ((1,), (), [5759665, 5759643]), ((), (1, 0), [9781153, 9781151]), ((1, 0), (0,), [24838551, 24838563])]:
    check([transfer(p + X * j, Y * (1 - j) + s) for j in range(2)] == expected, 'correctness gate failed')
from functools import lru_cache
C = [4, 1, 18, 7, 75, 31, 175, 63]
mask = 124

@lru_cache(None)
def endcounts(s, v):
    if s == 1 << v:
        return 1
    rest = s ^ 1 << v
    return sum((endcounts(rest, u) for u in range(8) if rest >> u & 1 and C[u] >> v & 1))

@lru_cache(None)
def startcounts(s, v):
    if s == 1 << v:
        return 1
    rest = s ^ 1 << v
    return sum((startcounts(rest, u) for u in range(8) if rest >> u & 1 and C[v] >> u & 1))
h = [1] + [sum((endcounts(s, v) for v in range(8) if s >> v & 1)) for s in range(1, 256)]
e = [sum((endcounts(s, v) for v in range(8) if s >> v & mask >> v & 1)) for s in range(256)]
a = [sum((startcounts(s, v) for v in range(8) if s >> v & mask >> v & 1)) for s in range(256)]
state = (h[255], sum((h[s] * h[255 ^ s] for s in range(256))), sum((h[s] * a[255 ^ s] for s in range(256))), sum((e[s] * h[255 ^ s] for s in range(256))), sum((e[s] * a[255 ^ s] for s in range(256))), a[255], e[255])
check(state == (45, 1440, 878, 577, 332, 25, 14), 'correctness gate failed')
result = {'all_receipt_witnesses_verified': len(r['witnesses']), 'k1_witnesses_verified': checked, 'core_state': state, 'bare_formula_through_k': 3, 'endpoint_controls': 'PASS', 'absorption': 'PASS', 'producer_normal_optimized_identical': (producer / 'RESULT.json').read_bytes() == (producer / 'RESULT.normal.json').read_bytes(), 'sha256': hashlib.sha256(flags).hexdigest()}
(ROOT / 'CONTROLS.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
