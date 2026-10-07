"""Audit A: which Polymath-ladder fields Frac Z[w1, w3, (w4), wt] are excluded by Lemma R (c-stable prime with known kappa <= 5)."""
import itertools
from auditA_multiquad import analyse, sqf
def field(ts):
    gens = []
    for t in [1] + list(ts):
        m = sqf(-(4 * t - 1))
        cand = gens + [m]
        ok = all(sqf(eval('*'.join(map(str, sub))) if len(sub) > 1 else sub[0]) != 1
                 for r in range(1, len(cand) + 1) for sub in itertools.combinations(cand, r))
        if ok: gens.append(m)
    return gens
def excluded(gens):
    res = analyse(gens, pmax=400)
    good = [r for r in res if isinstance(r[2], int) and r[2] <= 5]
    return (min(good, key=lambda r: r[2]) if good else None), [r[0] for r in res[:3]]
for base, label in (([3], "Z[w1,w3]"), ([3, 4], "Z[w1,w3,w4]")):
    print(f"one extra rotation over {label}:")
    for t in range(2, 16):
        if t in base: continue
        g = field(base + [t]); ex, first = excluded(g)
        print(f"  + w{t}: gens {g}: " + (f"excluded chi<={ex[2]} via p={ex[0]}" if ex else f"NOT excluded; first c-stable primes {first}"))
