"""Exhaustive, dependency-free controls for two-sided dyadic transport lifts."""
from itertools import product
from collections import Counter
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def mul(d, u):
    return tuple(sum(x*y for x, y in zip(row, u)) for row in d)


def trans(d):
    return tuple(zip(*d))


def matrix(flat, r, c):
    return tuple(tuple(flat[i*c:(i+1)*c]) for i in range(r))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b)) % 2


def predicted(u, v, a, b):
    if dot(v, a) != dot(u, b) or (not any(u) and any(a)) or (not any(v) and any(b)):
        return 0
    r, c = len(v), len(u)
    rank = r + c - 1 if any(u) and any(v) else r if any(u) else c if any(v) else 0
    return 1 << (r*c-rank)


def main():
    cases = lifts = 0
    for r, c in product(range(1, 4), repeat=2):
        ds = [matrix(z, r, c) for z in product(range(2), repeat=r*c)]
        for u in product(range(2), repeat=c):
            for v in product(range(2), repeat=r):
                table = Counter((tuple(t % 2 for t in mul(d, u)),
                                 tuple(t % 2 for t in mul(trans(d), v))) for d in ds)
                for a in product(range(2), repeat=r):
                    for b in product(range(2), repeat=c):
                        need(table[a, b] == predicted(u, v, a, b), 'binary image/count')
                        cases += 1
    # Exhaust all integral vectors modulo four and binary base matrices in 1x1,
    # 1x2, 2x1 and 2x2. Test lifts at levels one and two separately.
    for r, c in product(range(1, 3), repeat=2):
        ds = [matrix(z, r, c) for z in product(range(2), repeat=r*c)]
        for u in product(range(4), repeat=c):
            for v in product(range(4), repeat=r):
                for base in ds:
                    cu, cv = mul(base, u), mul(trans(base), v)
                    for h in (1, 2):
                        q = 1 << h
                        if any(z % q for z in cu+cv):
                            continue
                        a = tuple((-z//q) % 2 for z in cu)
                        b = tuple((-z//q) % 2 for z in cv)
                        ub, vb = tuple(x % 2 for x in u), tuple(x % 2 for x in v)
                        need(dot(vb, a) == dot(ub, b), 'shared residual corner')
                        actual = 0
                        for d in ds:
                            lifted = tuple(tuple(x+q*y for x, y in zip(br, dr)) for br, dr in zip(base, d))
                            actual += not any(z % (2*q) for z in mul(lifted, u)+mul(trans(lifted), v))
                        need(actual == predicted(ub, vb, a, b), 'direct dyadic lift count')
                        lifts += 1
    need(predicted((1,), (1,), (0,), (1,)) == 0, 'incompatible corner hostile')
    need(all((2*(1+2*d)) % 4 != 0 for d in (0, 1)), 'even differences hostile')
    candidates = [tuple(1+2*d for d in ds) for ds in product(range(2), repeat=2)]
    admitted = [c for c in candidates if sum(c) % 4 == 0]
    need(admitted == [(1, 3), (3, 1)], 'exact resource hostile lifts')
    need(not any(all(0 <= z <= 1 for z in c) for c in admitted), 'resource box has no lift')
    out = {'status': 'FINITE-EXACT controls; general proof in TRANSPORT_LIFTS.md',
           'binary_image_cases': cases, 'dyadic_lift_cases': lifts, 'checks': CHECKS,
           'resource_hostile': {'base': [1, 1], 'unrestricted_lifts_mod4': [list(c) for c in admitted], 'box': [0, 1]},
           'scope': 'No path existence, no native electrical input, no all-size wiring claim.'}
    target = HERE / 'receipt.json'
    if target.exists():
        need(json.loads(target.read_text()) == out, 'saved receipt')
    else:
        target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(out, sort_keys=True))


if __name__ == '__main__':
    main()
