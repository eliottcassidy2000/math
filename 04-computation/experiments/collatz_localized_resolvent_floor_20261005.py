"""Source-local signed kernel floors and a finite ROOT-extraction interface.

Exact arithmetic only. Moment intervals are conditional premises; controls
using a completed ROOT head do not create new rooted sources.
"""
from fractions import Fraction as F
from math import comb
import json
from collatz_atom_polynomial_dual_20261005 import injection_atom
from collatz_floor_transport_deadlines_20261005 import threshold_receipt, replay


def natural(n):
    if type(n) is not int or n < 0:
        raise ValueError('exact natural required')


def kernel(m, j):
    natural(m)
    natural(j)
    t = 1 << abs(m-j)
    return F(4*t, (1+t)**2)


def selector(m, degree, j):
    natural(degree)
    h = kernel(m, j)
    return (9*h-8)*h**degree


def interval(pair):
    if type(pair) is not tuple or len(pair) != 2:
        raise ValueError('exact interval pair required')
    if any(type(x) not in (int, F) for x in pair):
        raise ValueError('exact rational endpoints required')
    lo, hi = map(F, pair)
    if not 0 <= lo <= hi <= 1:
        raise ValueError('bounded kernel moment interval required')
    return lo, hi


def dual_bounds(degree, hd, hd1):
    """Conditional bounds, assuming BOTH packets describe the same actual law."""
    natural(degree)
    lo, hi = interval(hd)
    lo1, hi1 = interval(hd1)
    if degree == 0 and (lo, hi) != (F(1), F(1)):
        raise ValueError('zeroth moment is exactly one')
    if lo1 > hi:
        raise ValueError('consecutive bounded moments cannot increase')
    return (9*lo1-8*hi, 9*hi1-8*lo+8*F(16, 25)**degree)


def head_interval(m, degree, head):
    natural(m)
    natural(degree)
    if type(head) is not tuple or len(head) <= m:
        raise ValueError('declared complete head must pass target')
    if any(type(p) is not F or p < 0 for p in head) or sum(head) > 1:
        raise ValueError('exact subprobability head required')
    if degree == 0:
        return F(1), F(1)
    lo = sum((p*kernel(m, j)**degree for j, p in enumerate(head)), F(0))
    hi = lo + (1-sum(head))*kernel(m, len(head))**degree
    return lo, hi


def dyadic_floor(value):
    if type(value) is not F or value <= 0:
        raise ValueError('proved positive rational required')
    exponent = max(0, value.denominator.bit_length()-value.numerator.bit_length())
    result = F(1, 1 << exponent)
    while result > value:
        exponent += 1
        result /= 2
    return result


def extract_receipt(m, degree, hd, hd1):
    """Attempt literal extraction; actual replay prevents false positive output.

    A positive readout is conditional, not validation of the moment premise.
    If a claimed floor fails actual bounded replay, reject it. Very small
    floors can require a large finite computation; no complexity cap is hidden.
    """
    natural(m)
    lower, upper = dual_bounds(degree, hd, hd1)
    if lower <= 0:
        return {'status': 'no_positive_floor'}
    if lower > F(1, 3):
        raise ValueError('claimed nonroot weight exceeds its universal ceiling')
    epsilon = dyadic_floor(lower)
    result = threshold_receipt(6*m+3, epsilon)
    if result['status'] != 'met' or replay(6*m+3, result['word']) != 1:
        raise ValueError('claimed moment floor failed actual finite ROOT extraction')
    return {'status': 'literal_ROOT', 'epsilon': epsilon, 'degree': degree,
            'source': 6*m+3, 'word': result['word'],
            'counter_ceiling': result['counter_ceiling']}


def main():
    checks = 0

    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            raise ValueError(label)

    for m in range(13):
        for degree in range(9):
            for j in range(41):
                h = kernel(m, j)
                q = selector(m, degree, j)
                check(0 < h <= 1, 'bounded kernel')
                check(q <= int(m == j), 'global-sign discrete minorant')
                check(int(m == j)-q <= 8*F(16, 25)**degree, 'uniform error')
                check(q <= selector(m, degree+1, j), 'monotone refinement')
                check(kernel(m+3, j+3) == h, 'source shift preserves exact kernel')
                check(kernel(2*m, 2*j) == h*h/(2-h)**2, 'quadratic scale map')
    for distance in range(13):
        t = F(1 << distance)
        trace = t+1/t
        check(kernel(0, distance) == 4/(trace+2), 'trace coordinate')
        check(t*t+1/(t*t) == trace*trace-2, 'Chebyshev trace doubling')
        for degree in range(9):
            r = 1/(1+t)
            expanded = 4**degree*sum(((-1)**j*comb(degree, j)*r**(degree+j)
                                      for j in range(degree+1)), F(0))
            check(expanded == kernel(0, distance)**degree, 'resolvent expansion')
            check(4**degree*sum(comb(degree, j) for j in range(degree+1)) == 8**degree,
                  'primitive-conversion conditioning bill')

    # Independent finite probability controls, including an absent target.
    law = {0:F(1, 7), 1:F(2, 7), 3:F(4, 7)}
    for m in range(6):
        old_lower = F(-100)
        for degree in range(21):
            a = sum((p*kernel(m, j)**degree for j, p in law.items()), F(0))
            b = sum((p*kernel(m, j)**(degree+1) for j, p in law.items()), F(0))
            lo, hi = dual_bounds(degree, (a,a), (b,b))
            check(lo <= law.get(m, F(0)) <= hi, 'independent law enclosure')
            check(old_lower <= lo, 'exact readout refinement')
            old_lower = lo
            eta = F(1, 10**6)
            aw = (max(F(0),a-eta), min(F(1),a+eta)) if degree else (F(1),F(1))
            bw = (max(F(0),b-eta), min(F(1),b+eta))
            noisy_lo, _ = dual_bounds(degree, aw, bw)
            check(lo-17*eta <= noisy_lo <= lo, 'signed absolute-error budget')

    # Completed head is only a pipeline control, explicitly not new coverage.
    head = tuple(injection_atom(m)[0] for m in range(16))
    controls = []
    for m in (0, 1, 4):
        for degree in range(1, 101):
            a, b = head_interval(m, degree, head), head_interval(m, degree+1, head)
            low, high = dual_bounds(degree, a, b)
            check(low <= head[m] <= high, 'actual known-head enclosure')
            if low > 0:
                result = extract_receipt(m, degree, a, b)
                check(result['epsilon'] <= low <= head[m], 'retained dyadic floor')
                check(replay(result['source'], result['word']) == 1, 'literal pipeline output')
                controls.append({'source': result['source'], 'degree': degree,
                                 'floor': str(result['epsilon']),
                                 'actual_weight': str(head[m]),
                                 'deadline': result['counter_ceiling'],
                                 'actual_odd_steps': len(result['word'])})
                break
        else:
            raise ValueError('declared known-head control did not give positive readout')

    # An unsupported but numerically positive packet never becomes a false receipt.
    try:
        extract_receipt(4, 1, (F(1, 4), F(1, 4)), (F(1, 4), F(1, 4)))
    except ValueError:
        check(True, 'false source27 floor rejected by actual ROOT check')
    else:
        raise ValueError('unsupported high floor accepted')
    for bad in (True, 1.0, -1):
        try:
            kernel(bad, 0)
        except ValueError:
            check(True, 'malformed source index')
        else:
            raise ValueError('malformed index accepted')
    print(json.dumps({'status': 'PROVED conditional floor; universal positivity OPEN',
                      'checks': checks, 'direct_kernel_coefficient_norm': 17,
                      'error_bound': '8*(16/25)^d',
                      'known_head_controls': controls,
                      'coverage_scope': 'all example ROOT sources were already in the head'}, indent=2))


if __name__ == '__main__':
    main()
