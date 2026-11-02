"""Exact certificate for the combined Session139/baseline proof; no floats."""
from fractions import Fraction as F
from pathlib import Path
import json

GATES = 0


def need(value, why):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(why)


def kappa(n):
    r = n//2
    return F(4, r)*(1-F(3, 4)**r-(F(1, 4)**r if n % 2 == 0 else 0))


def main():
    c, t, d = F(999, 4000), F(25001, 1000000), F(667, 5000)
    lam, r = t*t/d, 160
    q = F(r)+F(33, 2)+F(8, r)
    need(F(1, 40) < c/2, 'V<=2 excluded at new cutoff')
    need(c-4*t > 0, 'retained matching bound increases')
    need((c-4*t)*r-16-F(8, r) == F(98867, 12500) > 7, 'eight heavy edges endpoint')
    slack = (c-d)-q/F(r-8)*F(3, 4)**8
    need(slack == F(1681793, 24903680000) > 0, 'tight geometric endpoint')
    need(q*F(3, 4)**(r-1) < c-d, 'other geometric endpoint')
    need(F(3, 4)*(1+F(1, r)) < 1, 'first endpoint monotonic bound')
    need(d-4*t > 0, 'light deficit positive')
    need(3*t/4-4*lam == F(25001, 2900000000) > 0, 'local dual nonnegative remainder')
    need(t-4*lam == F(9075363, 1450000000) > 0, 'positive total dual')
    need(lam*d == t*t and 0 < lam < 1, 'dual normalization')
    need(t-F(1, 40) == F(1, 1000000), 'strict cutoff slack')
    need(800000000*3**160 < 4**160, 'rho integer endpoint')
    need(t-F(1, 100000000) > F(1, 40), 'final contradiction endpoint')
    # Exact coefficients in S^i w^j, expanding both sides of (7).
    left = {(0, 0): lam*d, (1, 0): -t, (2, 0): F(1, 4),
            (1, 1): F(1, 2), (0, 1): -4*lam}
    square_terms = {(1, 0): F(1, 2), (0, 0): -t, (0, 1): F(1, 2)}
    right = {}
    for (i, j), x in square_terms.items():
        for (k, ell), y in square_terms.items():
            exponent = (i+k, j+ell)
            right[exponent] = right.get(exponent, 0)+x*y
    right[(0, 1)] += t-4*lam
    right[(0, 2)] -= F(1, 4)
    right = {e: x for e, x in right.items() if x}
    need(left == right, 'local dual identity by independent polynomial expansion')
    # Native boundaries: S=0, w=0, selected fair weight, and both sides of cutoff.
    for w in (F(0), t/2, t):
        for s in (F(0), t, 2*t, F(1), F(25)):
            actual = w*s/2+s*s/4-t*s+lam*(d-4*w)
            need(actual >= 0, 'finite boundary replay of dual')
    for n in range(320, 1500):
        r, k = n//2, kappa(n)
        q = F(r)+F(33, 2)+F(8, r)
        need((c-4*t)*r-16-F(8, r) > 7, 'finite heavy-matching replay')
        need(max(q*F(3, 4)**(r-1), q/F(r-8)*F(3, 4)**8) < c-d,
             'finite uniform tail replay')
        need(t-lam*(4-k*r) > k, 'finite contradiction replay')
    # All-order follows from the proof; these just check the displayed recurrences.
    for r in range(162, 200):
        odd = lambda j: F(j, 4)*kappa(2*j+1)
        even = lambda j: F(j, 4)*kappa(2*j)
        need(odd(r) == F(3, 4)*odd(r-1)+F(1, 4), 'odd coefficient recurrence')
        need(even(r) == even(r-1)-F(3, 16)*even(r-2)+F(3, 16), 'even recurrence')
        need(kappa(2*r+1)-kappa(2*r) == F(4, r)*F(1, 4)**r, 'parity correction')
    report = {'status': 'PASS', 'gates': GATES, 'threshold': 320,
              'scope': 'Exact endpoint certificate; infinite proof in STABILITY_320.md',
              'constants': {'c': str(c), 't': str(t), 'd': str(d), 'lambda': str(lam)},
              'tight_tail_slack': str(slack), 'remaining_open': [10, 319],
              'closed_since_original_1500_threshold': 1180,
              'closed_beyond_concurrent_526_threshold': 206}
    Path(__file__).with_name('stability_320_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'gates': GATES, 'threshold': 320}))


if __name__ == '__main__':
    main()
