"""Exact nonlinear bootstrap and flexible-dual endpoint certificate."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json

GATES = 0


def need(value, why):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(why)


def main():
    r = 28
    need(sum(F(4**j, factorial(j)) for j in range(8)) == F(16319, 315) > 50,
         'initial exponential bound')
    need(F(49, 200) > F(1, 7), 'exclude small variance')
    levels = [F(49, 200), F(1, 2), F(91, 100), F(99, 100)]
    expected_mass = [F(122, 175), F(139, 56), F(7493, 1400)]
    bootstrap = []
    for (c, nxt), expected in zip(zip(levels, levels[1:]), expected_mass):
        mass = c*r/4-1-F(1, 2*r)
        need(mass == expected > 0, 'bootstrap exact mass bound')
        taylor = sum(mass**j/factorial(j) for j in range(13))
        slack = taylor*(1-nxt-F(1, 4)**r)-1
        need(slack > 0, 'strict rational nonlinear bootstrap')
        bootstrap.append({'from': str(c), 'to': str(nxt), 'mass': str(mass),
                          'taylor_degree': 12, 'rational_slack': str(slack)})
    c, t, d = levels[-1], F(1, 7), F(381, 500)
    lam, q = t*t/d, r+F(33, 2)+F(8, r)
    need(c-4*t > 0, 'retained-edge bound increases')
    need((c-4*t)*r-4-F(2, r) == F(2677, 350) > 7, 'eight heavy edges')
    tail_slack = c-d-q/(r-8)*F(3, 4)**8
    need(tail_slack == F(1751781, 458752000) > 0, 'tight geometric endpoint')
    need(q*F(3, 4)**(r-1) < c-d, 'other geometric endpoint')
    need(F(3, 4)*(1+F(1, r)) < 1, 'endpoint decreases')
    need(d-4*t > 0, 'positive light deficit')
    need(lam == F(500, 18669), 'dual constant')
    need(3*t/4-4*lam == F(1, 74676) > 0, 'local remainder positive')
    need(t-4*lam == F(667, 18669) > 0, 'total dual positive')
    need(t*r == 4, 'strict contradiction valid at endpoint')
    square = {(1, 0): F(1, 2), (0, 0): -t, (0, 1): F(1, 2)}
    rhs = {}
    for e, x in square.items():
        for f, y in square.items():
            key = tuple(a+b for a, b in zip(e, f))
            rhs[key] = rhs.get(key, 0)+x*y
    rhs[(0, 1)] += t-4*lam
    rhs[(0, 2)] -= F(1, 4)
    lhs = {(0, 0): lam*d, (1, 0): -t, (2, 0): F(1, 4),
           (1, 1): F(1, 2), (0, 1): -4*lam}
    need({e: x for e, x in rhs.items() if x} == lhs, 'exact local dual polynomial')
    for n in range(56, 1500):
        rr = n//2
        qq = rr+F(33, 2)+F(8, rr)
        need((c-4*t)*rr-4-F(2, rr) > 7, 'finite heavy-edge replay')
        need(max(qq*F(3, 4)**(rr-1), qq/(rr-8)*F(3, 4)**8) < c-d, 'finite tail replay')
        need(F(4, rr) <= t, 'finite closing coefficient replay')
    report = {'status': 'PASS', 'gates': GATES, 'threshold': 56,
              'remaining_open': [10, 55], 'closed_since_original_1500': 1444,
              'bootstrap': bootstrap, 'tight_tail_slack': str(tail_slack),
              'constants': {'c': str(c), 't': str(t), 'd': str(d), 'lambda': str(lam)}}
    Path(__file__).with_name('stability_56_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'gates': GATES, 'threshold': 56}))


if __name__ == '__main__':
    main()
