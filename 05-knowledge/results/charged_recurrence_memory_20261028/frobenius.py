"""F_27 Paley affine/semilinear actions and their positive recurrence supports."""
from collections import Counter
from pathlib import Path
import json
from memory import sum_power

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def digits(x):
    return (x % 3, (x // 3) % 3, x // 9)


def encode(a):
    return sum((v % 3) * 3 ** i for i, v in enumerate(a))


def add(x, y):
    return encode([a + b for a, b in zip(digits(x), digits(y))])


def negative(x):
    return encode([-a for a in digits(x)])


def multiply(x, y):
    c = [0] * 5
    for i, a in enumerate(digits(x)):
        for j, b in enumerate(digits(y)):
            c[i + j] += a * b
    # t^3=t+1 in F_3[t]/(t^3-t-1).
    for i in (4, 3):
        c[i - 3] += c[i]
        c[i - 2] += c[i]
    return encode(c[:3])


def power(x, n):
    a = 1
    for _ in range(n):
        a = multiply(a, x)
    return a


def main(write=True):
    global CHECKS
    CHECKS = 0
    for x in range(1, 27):
        need(power(x, 26) == 1, 'nonzero field element unit')
    H = {multiply(x, x) for x in range(1, 27)}
    need(len(H) == 13 and negative(1) not in H, 'square and negative-square partition')
    T = tuple(sum(1 << y for y in range(27) if add(y, negative(x)) in H) for x in range(27))
    affine, semi, central = Counter(), Counter(), Counter()
    elements = set()
    for j in range(3):
        for a in H:
            for b in range(27):
                g = tuple(add(multiply(a, power(x, 3 ** j)), b) for x in range(27))
                need(len(set(g)) == 27, 'semilinear bijection')
                need(all(((T[x] >> y) & 1) == ((T[g[x]] >> g[y]) & 1)
                         for x in range(27) for y in range(27)), 'actual Paley tournament automorphism')
                fixed = sum(i == x for i, x in enumerate(g))
                semi[fixed] += 1
                if j == 0:
                    affine[fixed] += 1
                if b == 0:
                    central[fixed] += 1
                elements.add(g)
    need(len(elements) == 1053, 'specified semilinear subgroup has 1053 distinct elements')
    need(affine == {0: 26, 1: 324, 27: 1}, 'affine fixed-count histogram')
    need(semi == {0: 494, 1: 324, 3: 234, 27: 1}, 'Frobenius changes the offspring support')
    need(central == {1: 12, 3: 26, 27: 1}, 'negation centralizer fixed-count histogram')
    frob = tuple(power(x, 3) for x in range(27))
    need(all(frob[frob[frob[x]]] == x for x in range(27)) and sum(frob[x] == x for x in range(27)) == 3,
         'Frobenius order three with literal prime-subfield fixed points')
    supports = {0, 1, 3, 27}
    rows = []
    for k in range(1, 4):
        if k > 1:
            # Every possible root fixed count is padded to 27 by zero offspring.
            supports = sum_power(supports, 27)
        bound = 27 ** k
        coin = [False] * (bound + 1)
        coin[0] = True
        for i in range(1, bound + 1):
            coin[i] = any(i >= g and coin[i - g] for g in (24, 26, 27))
        if k >= 2:
            need({bound - x for x in supports} == {i for i, ok in enumerate(coin) if ok},
                 'bounded support equals reflected three-generator semigroup')
            need(coin.count(False) == 109 and len(supports) == bound + 1 - 109,
                 'exact semilinear arity-rank deficit')
        rows.append({'depth': k, 'semilinear_arity_order': len(supports),
                     'affine_arity_order': 3 if k == 1 else bound + 1 - 325})
    # Elementary interval proof's finite initial gap certificate.
    coin = [False] * 300
    coin[0] = True
    for i in range(1, 300):
        coin[i] = any(i >= g and coin[i - g] for g in (24, 26, 27))
    gaps = [i for i, ok in enumerate(coin) if not ok]
    need(len(gaps) == 109 and max(gaps) == 217, 'three-coin genus and conductor')
    out = {'status': 'FINITE-EXACT specified automorphism subgroups; no full-Aut claim',
           'checks': CHECKS, 'field_polynomial': 't^3-t-1 over F3',
           'affine_histogram': sorted(affine.items()), 'semilinear_histogram': sorted(semi.items()),
           'negation_centralizer_histogram': sorted(central.items()), 'rows': rows,
           'deficit_semigroup_generators': [24, 26, 27], 'deficit_gaps': gaps}
    if write:
        (HERE / 'frobenius.json').write_text(json.dumps(out, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'checks': CHECKS, 'semilinear_histogram': sorted(semi.items()), 'rows': rows}))
    return out


if __name__ == '__main__':
    main()
