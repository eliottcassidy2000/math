"""Native exact periods for all short bounded programs; independent trace partitions."""
from fractions import Fraction
from functools import reduce
from itertools import product
from math import gcd, lcm
from pathlib import Path
import json

GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def observe(program, n, m):
    for i, (a, b) in enumerate(program):
        if n % b == 0:
            return i, (n//b*a) % m
    return ('HALT',)


def least_period(word):
    size = len(word)
    for d in range(1, size+1):
        if size % d == 0 and all(word[i] == word[i % d] for i in range(size)):
            return d
    raise RuntimeError('no period')


def period_formula(program, m):
    active = [(a, b) for i, (a, b) in enumerate(program)
              if not any(b % bj == 0 for _, bj in program[:i])]
    d = lcm(*(b for _, b in active))
    g = reduce(gcd, (a for a, _ in active))
    direct = lcm(*(m*b//gcd(m, a) for a, b in active))
    need(gcd(d, g) == 1, 'active denominator/numerator coprime aggregates')
    need(direct == d*m//gcd(m, g), 'primewise simplified period')
    return direct


def trace(n, h):
    out = []
    halted = False
    for _ in range(h):
        if halted:
            out.append(0)
        elif n % 2 == 0:
            n //= 2
            out.append(1)
        elif n % 3 == 0:
            n = n//3*2
            out.append(2)
        else:
            halted = True
            out.append(0)
    return tuple(out)


def valuation(n, p):
    result = 0
    while n % p == 0:
        n //= p
        result += 1
    return result


def main():
    fractions = [(a, b) for a in range(1, 5) for b in range(1, 5) if gcd(a, b) == 1]
    programs = [(f,) for f in fractions] + list(product(fractions, repeat=2))
    cases, max_period = 0, 0
    for program in programs:
        # Known sufficient period uses every denominator, including inactive rules.
        for m in range(1, 13):
            sup = m*lcm(*(b for _, b in program))
            word = [observe(program, n, m) for n in range(1, sup+1)]
            actual, expected = least_period(word), period_formula(program, m)
            need(actual == expected, 'minimal native one-step observation period')
            need(all(observe(program, n, m) == observe(program, n+expected, m)
                     for n in range(1, sup+1)), 'native periodicity beyond first window')
            cases += 1
            max_period = max(max_period, actual)
    rows = []
    for h in range(10):
        modulus = 2**h*3**((h+1)//2)
        word = [trace(n, h) for n in range(1, modulus+1)]
        period = least_period(word)
        classes = len(set(word))
        expected_classes = 1+h+(h+1)**2//4
        need(period == modulus, 'least trace modulus')
        need(classes == expected_classes, 'least response class count')
        # A different enumeration uses abstract valuations, not residues.
        abstract = set()
        for x in range(h+1):
            for y in range(h+1):
                abstract.add(tuple(([1]*x+[2, 1]*y+[0]*h)[:h]))
        need(set(word) == abstract, 'residue and valuation trace universes agree')
        rows.append({'horizon': h, 'minimal_modulus': period, 'response_classes': classes,
                     'naive_modulus': 6**h})
    for n in range(1, 10001):
        rank = valuation(n, 2)+2*valuation(n, 3)
        state = n
        for i in range(rank):
            if state % 2 == 0:
                state //= 2
            elif state % 3 == 0:
                state = state//3*2
            else:
                raise RuntimeError('early halt contradicts rank')
            need(valuation(state, 2)+2*valuation(state, 3) == rank-i-1, 'rank drops one')
        need(state % 2 != 0 and state % 3 != 0, 'native halt at exact rank')
    need(period_formula(((1, 2), (2, 3)), 1) == 6, 'one-step index clock')
    need(rows[2]['minimal_modulus'] == 12 and rows[2]['naive_modulus'] == 36,
         'priority makes naive precision nonminimal')
    need(period_formula(((2, 1), (1, 3)), 6) == 3, 'inactive division excluded')
    crt_program = ((1, 2), (2, 3))
    chart5, chart7 = period_formula(crt_program, 5), period_formula(crt_program, 7)
    compatible = {(a, b) for a in range(chart5) for b in range(chart7) if (a-b) % 6 == 0}
    native_pairs = {(n % chart5, n % chart7) for n in range(210)}
    need(compatible == native_pairs and len(compatible) == 210, 'compatible CRT input charts')
    for m in range(1, 13):
        for n in range(1, 13):
            lm, ln = period_formula(crt_program, m), period_formula(crt_program, n)
            need(lcm(lm, ln) == period_formula(crt_program, lcm(m, n)), 'clock lcm compatibility')
            need(gcd(lm, ln) == period_formula(crt_program, gcd(m, n)), 'clock gcd compatibility')
    # Same current parity cannot support one more exact division/parity query.
    need(observe(((1, 2),), 2, 2) != observe(((1, 2),), 4, 2), 'closed parity quotient rejected')
    report = {'status': 'PASS', 'gates': GATES,
              'universe': 'Every reduced a/b with 1<=a,b<=4, lists of length 1 or 2; m=1..12',
              'programs': len(programs), 'one_step_cases': cases,
              'largest_observed_minimal_one_step_period': max_period,
              'trace_table': rows, 'rank_inputs_checked': [1, 10000],
              'crt_example': {'chart_moduli': [chart5, chart7], 'shared_modulus': 6,
                              'compatible_pairs': 210, 'cartesian_pairs': chart5*chart7}}
    Path(__file__).with_name('clock_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'gates': GATES, 'one_step_cases': cases,
                      'last_trace_row': rows[-1]}))


if __name__ == '__main__':
    main()
