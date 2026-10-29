"""Integer, multiplier-resolved cycle inventories; no assertion-only guards."""
from collections import Counter
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'complement_port_crt_20261028'))
import paley

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def generator(p):
    h = (p - 1) // 2
    return next(a for a in paley.squares(p)
                if len({pow(a, i, p) for i in range(h)}) == h)


def order(a, p):
    return next(d for d in range(1, p) if pow(a, d, p) == 1)


def multiply(a, b, h):
    out = Counter()
    for (u, x), m in a.items():
        for (v, y), n in b.items():
            out[(tuple((i + j) % h for i, j in zip(u, v)),
                 tuple(sorted(x + y)))] += m * n
    return out


def power(a, n, h, k):
    out = Counter({((0,) * k, ()): 1})
    for _ in range(n):
        out = multiply(out, a, h)
    return out


def transform(a, h, length=1, copies=1, charge=1):
    out = Counter()
    for (u, x), m in a.items():
        key = (tuple(charge * v % h for v in u),
               tuple(sorted(tuple(length * v for v in x) * copies)))
        out[key] += m
    return out


def prepend(out, a, root_charge, multiplicity=1):
    for (u, x), m in a.items():
        out[((root_charge,) + u, x)] += multiplicity * m


@lru_cache(None)
def inventory(p, k, central=False):
    if k == 0:
        return Counter({((), (1,)): 1})
    h = (p - 1) // 2
    g = generator(p)
    a = inventory(p, k - 1)
    size = sum(a.values())
    out = Counter()
    if central:
        c = inventory(p, k - 1, True)
        for e in range(h):
            d = order(pow(g, e, p), p)
            paired = transform(a, h, length=d, copies=2, charge=2)
            term = multiply(c, power(paired, h // d, h, k - 1), h)
            prepend(out, term, e, size ** (h - h // d))
        expected = h * sum(c.values()) * size ** h
    else:
        prepend(out, power(a, p, h, k - 1), 0)
        prepend(out, transform(a, h, length=p), 0,
                (p - 1) * size ** (p - 1))
        for e in range(1, h):
            d = order(pow(g, e, p), p)
            term = multiply(a, power(transform(a, h, length=d),
                                     (p - 1) // d, h, k - 1), h)
            prepend(out, term, e, p * size ** (p - 1 - (p - 1) // d))
        expected = p * h * size ** p
    need(sum(out.values()) == expected, 'cycle inventory group order')
    need(all(sum(cycles) == p ** k and len(u) == k
             for u, cycles in out), 'cycle inventory types')
    return out


def cycle_type(g):
    seen = set()
    out = []
    for i in range(len(g)):
        if i in seen:
            continue
        j = i
        length = 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = g[j]
        out.append(length)
    return tuple(sorted(out))


def literal_central_2(p):
    """Independent element enumeration, without cycle-index recurrence."""
    h = (p - 1) // 2
    H = paley.squares(p)
    A = paley.affine(p)
    t = paley.neg(p, 1)
    logs = {pow(generator(p), e, p): e for e in range(h)}
    counts = Counter()
    for root in H:
        for zero in H:
            for gs in product(A, repeat=h):
                children = [None] * p
                children[0] = tuple(zero * x % p for x in range(p))
                for i, g in zip(H, gs):
                    children[i] = g
                    children[-i % p] = paley.compose(t, paley.compose(g, t))
                c = paley.lift(p, root, 0, children)
                charges = tuple(logs[x] for x in paley.characters(c, p, 2))
                counts[(charges, cycle_type(c))] += 1
    return counts


def summary(p, k, central, kernel=False):
    data = inventory(p, k, central)
    if kernel:
        data = Counter({key: m for key, m in data.items() if not any(key[0])})
    total = sum(data.values())
    fixed = Counter()
    for (_, cycles), m in data.items():
        fixed[cycles.count(1)] += m
    moments = [sum(m * x ** r for x, m in fixed.items()) // total for r in range(7)]
    palettes = []
    for q in (2, 3):
        numerator = sum(m * q ** len(cycles) for (_, cycles), m in data.items())
        need(numerator % total == 0, 'integral palette orbit count')
        palettes.append(numerator // total)
    for r, value in enumerate(moments):
        need(value * total == sum(m * x ** r for x, m in fixed.items()),
             'integral ordered port orbit count')
    return {'prime': p, 'depth': k, 'group': 'Hol' if kernel else ('C' if central else 'A'),
            'order': total, 'inventory_cells': len(data), 'fixed_support': sorted(fixed),
            'fixed_histogram': sorted(fixed.items()), 'ordered_ports_0_to_6': moments,
            'palette_orbits_q2_q3': palettes}


def main(write=True):
    global CHECKS
    CHECKS = 0
    inventory.cache_clear()
    rows = []
    for p in (3, 7, 11):
        for k in (1, 2):
            for central, kernel in ((True, False), (True, True)):
                rows.append(summary(p, k, central, kernel))
            h = (p - 1) // 2
            c = inventory(p, k, True)
            sectors = Counter()
            for (charge, _), m in c.items():
                sectors[charge] += m
            need(len(sectors) == h ** k and len(set(sectors.values())) == 1,
                 'equal multiplier sectors')
        literal = Counter()
        logs = {pow(generator(p), e, p): e for e in range((p - 1) // 2)}
        for a in paley.affine(p):
            literal[((logs[paley.characters(a, p, 1)[0]],), cycle_type(a))] += 1
        need(literal == inventory(p, 1), 'base affine inventory versus literal permutations')
    for p in (3, 7):
        literal = literal_central_2(p)
        need(literal == inventory(p, 2, True), 'every depth-two charged cycle cell agrees')
        rows.append({'literal_centralizer_prime': p, 'elements': sum(literal.values()),
                     'charged_cycle_cells': len(literal)})
    # Same full cycle index does not classify joint transport (inherited hostile).
    parallel = [tuple(x if x < 3 else 3 * (x // 3) + (x % 3 + a) % 3
                      for x in range(9)) for a in range(3)]
    mirrored = [tuple(x if x < 3 else 3 * (x // 3) + (x % 3 + (a if x < 6 else -a)) % 3
                      for x in range(9)) for a in range(3)]
    need(Counter(map(cycle_type, parallel)) == Counter(map(cycle_type, mirrored)),
         'full cycle index hostile agrees')
    counts = [sum((g[3], g[6]) == (4, 7) for g in group) for group in (parallel, mirrored)]
    need(counts == [1, 0], 'joint transporter hostile separates')
    out = {'status': 'FINITE-EXACT checks of proved recurrences; not an independent researcher audit',
           'checks': CHECKS, 'rows': rows, 'transporter_hostile': counts}
    if write:
        (HERE / 'charged.json').write_text(json.dumps(out, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'checks': CHECKS, 'summaries': [(r.get('prime'), r.get('depth'),
                      r.get('group'), r.get('ordered_ports_0_to_6')) for r in rows]}))
    return out


if __name__ == '__main__':
    main()
