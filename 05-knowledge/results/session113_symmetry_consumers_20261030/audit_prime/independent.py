import sys, json, itertools, random
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'prime_ranking'))
from prime_rank import PrimeBank
from literal_ring import Ring
B = PrimeBank(5, 2)
R = Ring(5, 2)
mp = []
seen = {}
for d, (_, j) in zip(B.depths, B.blocks):
    u = seen.get(d, 0) + 1
    seen[d] = u
    z = (u, 0)
    for _ in range(d):
        z = R.mul(z, (0, 1))
    for i in range(5 ** j):
        s = pow(6, i, 25)
        mp.append(R.lookup[z[0] * s % 25, z[1] * s % 25])
inv = {old: (i, e) for i, (old, e) in enumerate(mp)}
cuts = [[(inv[a][0], inv[b][0], e ^ inv[a][1] ^ inv[b][1]) for a, b, e in cut] for cut in R.cuts]
count = 0
all_count = 0
samples = []
targets = {0, 1, 2, 3, 4, 15, 16, 17, 255, 256, 257, 332927, 332928, 665853, 665854, 665855}
rng = random.Random(2026113)
targets.update((rng.randrange(665856) for _ in range(30)))

def visit(i, group, flat):
    global count, all_count
    if i == len(B.blocks):
        all_count += 1
        prime = not any((all((flat[a] ^ flat[b] == e for a, b, e in cut)) for cut in cuts))
        if prime:
            if count in targets:
                samples.append((count, flat))
            count += 1
        return
    kind, j = B.blocks[i]
    L = 5 ** j
    shifts = {(a + (b if kind == 'B' else 0)) % L for a, b in group}
    for w in itertools.product((0, 1), repeat=L):
        if any((w > w[s:] + w[:s] for s in shifts)):
            continue
        subgroup = [(a, b) for a, b in group if w == w[(a + (b if kind == 'B' else 0)) % L:] + w[:(a + (b if kind == 'B' else 0)) % L]]
        visit(i + 1, subgroup, flat + w)
visit(0, list(itertools.product(range(5), repeat=2)), ())
if not (all_count == 692224 and count == 665856 == B.prime_total()):
    raise RuntimeError('independent prime audit failed')
for r, x in samples:
    if not B.flat(B.unrank(r)) == x:
        raise RuntimeError('independent prime audit failed')
    if not B.rank(B.split(x)) == r:
        raise RuntimeError('independent prime audit failed')
for r in (-1, B.prime_total(), B.prime_total() + 1, 1.5, '0', None):
    try:
        B.unrank(r)
    except (ValueError, TypeError):
        pass
    else:
        raise AssertionError(('invalid rank accepted', r))
for t in [0, 1, 2, 2 ** 25 - 2, 2 ** 25 - 1, 2 ** 25] + [rng.randrange(2 ** 25) for _ in range(30)]:
    vals = []
    for v in (0, 1):
        n = 5 ** v
        words = [w for w in itertools.product((0, 1), repeat=n) if v == 0 or len(set(w)) > 1]
        vals.append(sum((int(''.join(map(str, w * (25 // n))), 2) >= t for w in words)))
    if not B.ge_periods(2, 2, t)[:2] == vals:
        raise RuntimeError('independent prime audit failed')
out = dict(status='PASS', all_orbits_explicit=all_count, prime_orbits_literal=count, independent_rank_controls=len(samples), invalid_rank_controls=6, hidden_original_period_thresholds=36)
print(json.dumps(out, indent=2))
