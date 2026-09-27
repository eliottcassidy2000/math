"""Independent audit of incoming coefficient-descent cylinder claims."""
from fractions import Fraction
from itertools import product
import json
import math


def need(ok, message):
    if not ok:
        raise ValueError(message)


def v2(n):
    need(n > 0, 'positive valuation')
    return (n & -n).bit_length() - 1


def word(n, j):
    vals = []
    for _ in range(j):
        v = v2(3 * n + 1)
        n = (3 * n + 1) >> v
        vals.append(v)
    return tuple(vals), n


def gap_audit():
    p3 = 1
    worst = Fraction()
    witness = None
    for j in range(1, 5001):
        p3 *= 3
        A = p3.bit_length()
        need((1 << (A - 1)) < p3 < (1 << A), 'minimal exact exponent')
        ratio = Fraction(p3 - (1 << j), ((1 << A) - p3) * (1 << j))
        need(ratio < 1, 'one-member gap inequality')
        if ratio > worst:
            worst, witness = ratio, (j, A)
    for j in range(1, 156):
        exact = (3 ** j).bit_length()
        for A in range(exact + 3):
            need((A > j * math.log2(3)) == (A >= exact), 'reported finite float cutoff')
    return dict(jmax=5000, largest_ratio=str(worst), witness=witness,
                float_coefficient_cutoffs_independently_checked_through_step=155)


def cylinder_controls():
    controls = 0
    for j in range(1, 5):
        for vals in product(range(1, 6), repeat=j):
            A = S = 0
            for v in vals:
                S = 3 * S + (1 << A)
                A += v
            p3 = 3 ** j
            coarse = (-S * pow(p3, -1, 1 << A)) % (1 << A)
            exact = (((1 << A) - S) * pow(p3, -1, 1 << (A + 1))) % (1 << (A + 1))
            for lift in range(4):
                n = coarse + (lift << A)
                actual_vals, actual_end = word(n, j)
                Q = (p3 * n + S) >> A
                need(actual_vals[:-1] == vals[:-1] and actual_vals[-1] >= vals[-1],
                     'coarse cylinder meaning')
                need(actual_end == Q >> v2(Q), 'coarse endpoint oddpart')
                n = exact + (lift << (A + 1))
                actual_vals, actual_end = word(n, j)
                need(actual_vals == vals and actual_end == (p3 * n + S) >> A,
                     'exact cylinder meaning')
                controls += 2
    return dict(controls=controls, minimal_hostile=dict(word=[1], A=1, source=1,
                claimed_affine_endpoint=2, actual_endpoint=1))


def census(jmax=14):
    p3 = [3 ** j for j in range(jmax + 1)]
    checked = 0
    coarse_exceptions = []
    exact_exceptions = []

    def visit(j, A, S, vals):
        nonlocal checked
        next_j = j + 1
        next_S = 3 * S + (1 << A)
        critical_A = p3[next_j].bit_length()
        max_A = 40 if next_j == 1 else 41
        for next_A in range(critical_A, max_A + 1):
            final_v = next_A - A
            need(final_v >= 1, 'positive final valuation')
            denominator = (1 << next_A) - p3[next_j]
            coarse = (-next_S * pow(p3[next_j], -1, 1 << next_A)) % (1 << next_A)
            exact = (((1 << next_A) - next_S) * pow(p3[next_j], -1, 1 << (next_A + 1))) % (1 << (next_A + 1))
            need(next_S < (1 << next_A) * denominator, 'at most one coarse candidate')
            if coarse * denominator <= next_S:
                coarse_exceptions.append((coarse, vals + (final_v,), str(Fraction(next_S, denominator))))
            if exact * denominator <= next_S:
                exact_exceptions.append((exact, vals + (final_v,), str(Fraction(next_S, denominator))))
            checked += 1
        if next_j < jmax:
            for final_v in range(1, critical_A - A):
                visit(next_j, A + final_v, next_S, vals + (final_v,))

    visit(0, 0, 0, ())
    for j in range(1, jmax + 1):
        need(j * 3 ** (j - 1) < (1 << 41) - 3 ** j,
             'omitted total-precision tail has threshold below1')
    need(checked == 606746, 'incoming census reproduction')
    need(coarse_exceptions == [(1, (2,), '1')], 'coarse exception census')
    need(exact_exceptions == [(1, (2,), '1')], 'exact exception census')
    return dict(jmax=jmax, labeled_coarse_cylinders=checked,
                exact_word_cylinders=checked, coarse_exceptions=coarse_exceptions,
                exact_exceptions=exact_exceptions,
                actual_cutoff='j=1 through totalA40; j>=2 through totalA41',
                uniform_omitted_tail='A>=41 implies N<1 using S<=j*3^(j-1)')


def density_audit():
    counts = {0: 1}
    densities = []
    for j in range(1, 61):
        largest_A = (3 ** j).bit_length() - 1
        following = {}
        for A, count in counts.items():
            for new_A in range(A + 1, largest_A + 1):
                following[new_A] = following.get(new_A, 0) + count
        counts = following
        densities.append(sum((Fraction(count, 1 << A) for A, count in counts.items()), Fraction()))
    return {str(j): dict(exact=str(densities[j - 1]), decimal=float(densities[j - 1]))
            for j in (1, 2, 3, 4, 41, 60)}


def bank_hostile():
    rows = []
    for q in range(1, 342, 2):
        c = 3 * q
        A = J = 0
        while c >= 2 * q:
            need(J < 1000, 'finite old core bank')
            previous_A = A
            a = v2(3 * c + 1)
            c = (3 * c + 1) >> a
            A += a
            J += 1
        for R in range(previous_A + 1, A + 1):
            if (1 << (R + 1)) - 3 ** (J + 1) > max(0, 3 * (1 << (A - R)) * c - 6 * q + 1):
                break
        else:
            raise ValueError('missing old bank row')
        K = R + 1
        residue = ((6 * q - 1) * pow(3, -1, 1 << K)) % (1 << K)
        rows.append((residue, K))
    source = 4091
    need(not any(source % (1 << K) == residue for residue, K in rows), 'outside old bank')
    n = source
    orbit = []
    A = 0
    coefficient_time = None
    while n >= source:
        need(len(orbit) < 100, 'small hostile cap')
        a = v2(3 * n + 1)
        n = (3 * n + 1) >> a
        A += a
        orbit.append(n)
        if coefficient_time is None and 3 ** len(orbit) < (1 << A):
            coefficient_time = len(orbit)
    need(len(orbit) == 8 and n == 1639 and coefficient_time == 8, 'bank residual witness')
    return dict(source=source, old_bank_member=False, orbit=orbit,
                first_descent=8, first_coefficient_descent=8)


if __name__ == '__main__':
    print(json.dumps(dict(status='Independent exact finite audit; global claims separately scoped',
                         gap=gap_audit(), cylinders=cylinder_controls(),
                         census=census(), densities=density_audit(), bank_residual=bank_hostile()), indent=2))
