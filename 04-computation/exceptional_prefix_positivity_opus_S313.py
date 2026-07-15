# opus-2026-07-15-S313 -- HYP-6920: positivity sweep for the 4 exceptional
# prefixes (B1/B2 mechanisms silent). For every single-cluster pattern dvec
# (residues from R, diameter <= 42, sorted consecutive gaps <= 7), find an
# exact point t0 in E with a_dvec(t0) > 0 (proves integral_E a > 0 with an
# explicit triangle floor). Report any pattern with NO positive sample
# (candidate uniform-coverage counterexample => genuine gap).
from fractions import Fraction
import itertools, time

W = Fraction(1, 13)

def safe_set(P):
    ivs = [(Fraction(0), Fraction(1))]
    for q in P:
        bands = [(Fraction(13*k+1, 13*q), Fraction(13*(k+1)-1, 13*q)) for k in range(q)]
        new = []
        for (a, b) in ivs:
            for (c, d) in bands:
                lo, hi = max(a, c), min(b, d)
                if lo < hi: new.append((lo, hi))
        ivs = sorted(new)
    return ivs

def a_at(ds, t):
    arcs = []
    for d in ds:
        c = (-d * t) % 1
        lo, hi = c - W, c + W
        if lo < 0: arcs.extend([(lo % 1, Fraction(1)), (Fraction(0), hi)])
        elif hi > 1: arcs.extend([(lo, Fraction(1)), (Fraction(0), hi % 1)])
        else: arcs.append((lo, hi))
    arcs.sort()
    tot, cur = Fraction(0), Fraction(0)
    for (lo, hi) in arcs:
        if hi <= cur: continue
        tot += hi - max(lo, cur)
        cur = hi
    return 1 - tot

CASES = [
    ([2, 4, 6, 8, 10], [1, 3, 5, 7, 9, 11, 12]),
    ([2, 5, 7, 10, 12], [1, 3, 4, 6, 8, 9, 11]),
    ([4, 5, 6, 11, 12], [1, 2, 3, 7, 8, 9, 10]),
    ([8, 9, 10, 11, 12], [1, 2, 3, 4, 5, 6, 7]),
]

t0 = time.time()
for P, R in CASES:
    E = safe_set(P)
    # sample points: cell midpoints from a moderately fine exact grid per component
    samples = []
    for (a, b) in E:
        n = 40
        samples.extend(a + (b - a) * Fraction(2*i+1, 2*n) for i in range(n))
    npat = nfail = 0
    fails = []
    for base in R:
        cs = sorted(((r - base) % 13) for r in R if r != base)
        # lift choices per difference
        choices = [[c + 13*k for k in range(4) if c + 13*k <= 42] for c in cs]
        for combo in itertools.product(*choices):
            ds = sorted(set([0] + list(combo)))
            if len(ds) != 7: continue
            # single-cluster filter: consecutive sorted gaps <= 7
            if any(y - x > 7 for x, y in zip(ds, ds[1:])): continue
            npat += 1
            ok = False
            for t in samples:
                if a_at(ds, t) > 0: ok = True; break
            if not ok:
                nfail += 1
                fails.append(ds)
    print(f"P={P}: patterns={npat}, positivity failures={nfail} "
          f"[{time.time()-t0:.0f}s]", flush=True)
    for f in fails[:5]: print("    FAIL dvec:", f, flush=True)
print(f"done {time.time()-t0:.0f}s")
