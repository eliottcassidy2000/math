"""Exact modular selectors for the proved physical module-repair family."""
from pathlib import Path
from math import gcd
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def count_mod(n, m, modulus):
    return (15 * pow(100, n, modulus) * pow(2, m, modulus)
            - pow(11, n, modulus)) % modulus


def select(n, k, target):
    if n < 1 or k < 1 or target % 11 == 0:
        raise ValueError('n,k positive and target a unit modulo 11')
    m = next(j for j in range(10) if count_mod(n, j, 11) == target % 11)
    for j in range(1, k):
        power = 11 ** j
        error = (target - count_mod(n, m, power * 11)) // power % 11
        slope = 93 * 15 * pow(100, n, 11) * pow(2, m, 11) % 11
        digit = error * pow(slope, -1, 11) % 11
        m += digit * 10 * 11 ** (j - 1)
    return m


def select_exact_seven(n, k, target):
    period = 10 * 11 ** (k - 1)
    m = select(n, k, target)
    m += next(t for t in range(3) if (m + t * period - n) % 3 == 0) * period
    if (2 + 5 * n + (m - n) // 3) % 7 == 0:
        m += 3 * period
    return m


def valuation(n, p):
    n = abs(n)
    if not n:
        raise ValueError('zero has no finite valuation')
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def main():
    tables = []
    for n in (1, 2, 3):
        for k in (1, 2, 3):
            modulus = 11 ** k
            period = 10 * 11 ** (k - 1)
            image = [count_mod(n, m, modulus) for m in range(period)]
            units = {a for a in range(modulus) if gcd(a, modulus) == 1}
            need(len(set(image)) == period and set(image) == units,
                 'independent complete exponent image')
            for a in sorted(units):
                m = select(n, k, a)
                need(0 <= m < period and image[m] == a, 'digit inverse')
                seven = select_exact_seven(n, k, a)
                need(0 <= seven < 6 * period, 'joint exponent bound')
                need(count_mod(n, seven, modulus) == a, 'joint target digit')
                residue = count_mod(n, seven, 49)
                need(residue % 7 == 0 and residue != 0, 'exact seven valuation')
            tables.append(dict(n=n, k=k, unit_targets=len(units)))
    need(set(count_mod(0, m, 11) for m in range(10)) == set(range(11)) - {10},
         'n=0 hostile boundary')
    pairs = 0
    for n in (1, 2, 3):
        for delta in range(1, 243):
            difference = 15 * 100 ** n * (2 ** delta - 1)
            predicted = 0 if delta % 10 else 1 + valuation(delta // 10, 11)
            need(valuation(difference, 11) == predicted, 'count-difference ancestry')
            pairs += 1
    leading = []
    for j in range(1, 7):
        power = 11 ** j
        residue = pow(2, 10 * 11 ** (j - 1), power * 11)
        need((residue - 1) // power == 93 % 11, 'leading unit')
        leading.append(dict(precision=j, leading_unit=(residue - 1) // power))
    fixtures = []
    for n, k, a in ((1, 2, 1), (1, 3, 7), (2, 2, 63), (3, 3, 881)):
        fixtures.append(dict(n=n, k=k, target=a, exponent=select(n,k,a),
                             exact_seven_exponent=select_exact_seven(n,k,a)))
    result = dict(status='PASS', scope='Exact residues in a proved physical family; not absolute count coverage',
                  tables=tables, unit_targets=sum(t['unit_targets'] for t in tables),
                  valuation_controls=pairs, leading_units=leading, fixtures=fixtures)
    Path(__file__).with_name('residue_receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print('PASS digit-lift residue selector, simultaneous exact factor 7, actual-count valuations')


if __name__ == '__main__':
    main()
