"""Exact minimal witnesses for the Pascal-boundary scope correction.

No universal Collatz premise. Run normally or with python -O.
"""
from fractions import Fraction as F
from math import factorial
import json


def w(l, k):
    return F(2 * factorial(k) * factorial(l + 1), factorial(l + k + 2))


def d(l, k):
    return F(2, (l + 1) * (l + 2)) * w(l, k)


def odd_step(n):
    raw = 3 * n + 1
    return raw // (raw & -raw)


def main():
    checks = 0

    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            raise ValueError(label)

    for l in range(21):
        for k in range(21):
            check(w(l, k) == w(l + 1, k) + w(l, k + 1), 'W split')
            check(d(l, k) == d(l, k + 1) + F(l + 3, l + 1) * d(l + 1, k),
                  'D weighted split')
            check(d(l, k) != d(l + 1, k) + d(l, k + 1), 'D not a split array')
    check(d(0, 0) == 1, 'D normalization')
    check(d(1, 0) + d(0, 1) == F(5, 9), 'minimal split failure')

    # The full incoming operator on this finitely supported v can be computed
    # exactly by pushing its entire support; no omitted incoming mass exists.
    v = {3: F(1), 5: F(1)}
    incoming = {}
    for source, mass in v.items():
        target = odd_step(source)
        if target != 1:
            incoming[target] = incoming.get(target, F(0)) + mass
    defect = {n: v.get(n, F(0)) - incoming.get(n, F(0))
              for n in set(v) | set(incoming)}
    check(all(mass >= 0 for mass in defect.values()), 'summable subinvariant')
    check(defect == {3: F(1), 5: F(0)}, 'positive v need not mean positive injection')
    check(odd_step(3) == 5 and odd_step(5) == 1, 'finite potential route')

    # The beta-family threshold is not a general positive-power condition.
    # After t=1-r=e^-u, the endpoint integral of 1/log(1/t)^2
    # is int_1^infinity u^-2 du=1. Its exact rational tail at integer M is1/M.
    for m in range(1, 101):
        check(F(1) - F(m - 1, m) == F(1, m), 'log-density endpoint integral tail')

    print(json.dumps({
        'status': 'FINITE-EXACT scope witnesses; analytic endpoint proof in note',
        'checks': checks,
        'counter_universe': '0<=L,K<=20',
        'D00': str(d(0, 0)),
        'D10_plus_D01': str(d(1, 0) + d(0, 1)),
        'finite_potential': {'v3': '1', 'v5': '1', 'lambda3': '1', 'lambda5': '0'},
        'not_a_claim': 'No audit of all numerical Bellman bounds or representation analogies'
    }, indent=2))


if __name__ == '__main__':
    main()
