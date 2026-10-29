"""Exact fixed-port Fourier jets, grouped by deepest nontrivial character."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json
import charged

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def add(*terms):
    return [sum(row[i] for row in terms) for i in range(len(terms[0]))]


def scale(a, c):
    return [c * x for x in a]


def mul(a, b):
    out = [Q(0)] * len(a)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:len(a) - i]):
                if y:
                    out[i + j] += x * y
    return out


def power(a, n):
    out = [Q(1)] + [Q(0)] * (len(a) - 1)
    for _ in range(n):
        out = mul(out, a)
    return out


def double(a):
    return [x * 2 ** i for i, x in enumerate(a)]


def laws(p, depth, R):
    h = (p - 1) // 2
    unit = [Q(1)] + [Q(0)] * R
    ordinary = [Q(1, factorial(i)) for i in range(R + 1)]
    central = ordinary[:]
    af, cf = [], []
    rows = []
    for k in range(depth + 1):
        hol = add(central, *[scale(c, h ** j * (h - 1)) for j, c in enumerate(cf)])
        rows.append((ordinary, central, af, cf, hol))
        a0 = scale(add(power(ordinary, p), scale(ordinary, -p), scale(unit, p - 1)), Q(1, p * h))
        c0 = scale(mul(central, add(power(double(ordinary), h), scale(unit, -1))), Q(1, h))
        new_af = [a0] + [scale(power(a, p), Q(1, p * h)) for a in af]
        new_cf = [c0] + [scale(mul(c, power(double(a), h)), Q(1, h))
                         for a, c in zip(af, cf)]
        central = scale(mul(central, add(power(double(ordinary), h), scale(unit, h - 1))), Q(1, h))
        ordinary = add(scale(power(ordinary, p), Q(1, p * h)),
                       scale(ordinary, Q(h - 1, h)), scale(unit, Q(p - 1, p * h)))
        af, cf = new_af, new_cf
    return rows


def moment(a, r):
    return a[r] * factorial(r)


def main(write=True):
    global CHECKS
    CHECKS = 0
    output = []
    for p in (7, 11, 19):
        h = (p - 1) // 2
        rows = laws(p, 6, p + 1)
        for k, (a, c, af, cf, hol) in enumerate(rows):
            for r in range(p + 2):
                need(moment(hol, r).denominator == 1, 'integral holonomy moment jet')
            if k:
                need(moment(hol, 1) == p + 2 * k - 2, 'all-depth one-port holonomy count')
                need(moment(hol, 2) == 8 * k * k + (10 * p - 30) * k + p * p - 10 * p + 22,
                     'all-depth two-port holonomy count')
            if k >= 2:
                d = add(c, scale(cf[0], h - 1))
                need(hol[:p] == d[:p], 'all arities below p miss deeper charge')
                delta = moment(hol, p) - moment(d, p)
                need(delta == (h - 1) * 2 ** p * factorial(p), 'sharp depth-independent p-port gap')
            if k == 1 or (k == 2 and p in (7, 11)):
                hist = charged.inventory(p, k, True)
                ker = {key: m for key, m in hist.items() if not any(key[0])}
                for r in range(p + 2):
                    literal = Q(sum(m * cyc.count(1) ** r for (_, cyc), m in ker.items()), sum(ker.values()))
                    need(moment(hol, r) == literal, 'jet versus complete charged cycle inventory')
            if k:
                output.append({'prime': p, 'depth': k,
                               'holonomy_ports_0_to_4': [int(moment(hol, r)) for r in range(5)],
                               'first_deeper_charge_gap_at_p': (h - 1) * 2 ** p * factorial(p) if k >= 2 else None})
    # The next threshold 49 is checked as a complete rational Taylor jet.
    rows = laws(7, 4, 49)
    thresholds = []
    for k in (3, 4):
        for j, (a, c) in enumerate(zip(rows[k][2], rows[k][3])):
            r = 7 ** j
            if r > 49:
                need(not any(c), 'character deeper than truncation is exactly zero in jet')
                continue
            need(not any(c[:r]) and c[r] > 0, 'first nonzero character moment at p-power')
            exponent = ((7 ** j - 1) // 6 - j) // 2
            expected = Q(2 ** r, 3 ** j * 21 ** exponent)
            need(c[r] == expected, 'explicit leading character coefficient')
            if 2 * r <= 49:
                need(not any(a[:2 * r]) and a[2 * r] == Q(1, 21 ** ((r - 1) // 6)),
                     'automorphism character valuation and leading coefficient')
            delta = factorial(r) * 3 ** j * 2 * c[r]
            need(delta.denominator == 1, 'integer first sector split')
            thresholds.append({'prime': 7, 'depth': k, 'charge_depth': j,
                               'first_visible_arity': r, 'first_orbit_gap': int(delta)})
    out = {'status': 'FINITE-EXACT jets of PROVED recurrences; same-agent alternate representation',
           'checks': CHECKS, 'rows': output, 'thresholds': thresholds}
    if write:
        (HERE / 'jets.json').write_text(json.dumps(out, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'checks': CHECKS, 'thresholds': thresholds}))
    return out


if __name__ == '__main__':
    main()
