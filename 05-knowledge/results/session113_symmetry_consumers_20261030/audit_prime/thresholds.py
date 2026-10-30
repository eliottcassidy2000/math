import sys, itertools, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'prime_ranking'))
from prime_rank import PrimeBank
B = object.__new__(PrimeBank)
B.p = 3
checks = 0
for q in (0, 1, 2):
    stride = 3 ** q
    reps = []
    for w in itertools.product((0, 1), repeat=9):
        if w != min((w[s:] + w[:s] for s in range(0, 9, stride))):
            continue
        v = next((v for v in range(3) if all((w[i] == w[i % 3 ** v] for i in range(9)))))
        reps.append((int(''.join(map(str, w)), 2), v))
    for t in range(513):
        want = [sum((x >= t and v == e for x, v in reps)) for e in range(3)]
        got = B.ge_periods(2, q, t)
        if got != want:
            raise RuntimeError((q, t, got, want))
        checks += 1
print(json.dumps(dict(status='PASS', thresholds=checks, length=9, strides=[1, 3, 9], scope='abstract period arithmetic only; no p=3 ring theorem'), indent=2))
