"""Exact carry-defect bounds and all-height finite-offset certificates.

No global Collatz conclusion. Checks use explicit exceptions under python -O.
"""
from fractions import Fraction
from hashlib import sha256
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def valuation(n):
    need(n > 0, "positive valuation argument")
    return (n & -n).bit_length() - 1


def odd_step(n):
    m = 3 * n + 1
    a = valuation(m)
    return m >> a, a


def jet(n):
    return sum(s * (1 << s) for s in range(n.bit_length()) if (n >> s) & 1)


def carry(n):
    state, total, slack, s = 1, 0, 0, 0
    transitions = []
    while (n >> s) or state:
        bit = (n >> s) & 1
        after = (3 * bit + state) // 2
        output = 3 * bit + state - 2 * after
        charge = (6 if state == 0 and bit else 2 if state == 2 and bit else 0)
        need(charge == 8 * bit + 2 * (state == 2) - 4 * (after == 2) - 2 * after,
             "local state potential")
        transitions.append((s, state, bit, output, after, charge))
        total += after << (s + 1)
        slack += charge << s
        state, s = after, s + 1
    need(slack == 8 * n - total, "telescoped slack")
    return total, slack, transitions


def lower_equality(n):
    return n % 8 == 1 and not ((n >> 3) & (n >> 4))


def first_descent(n, cap=10000):
    y = n
    for count in range(1, cap + 1):
        y, _ = odd_step(y)
        if y < n:
            return count, y
    raise RuntimeError(f"descent not found within finite cap for {n}")


def trace_one(n, cap=10000):
    count, total = 0, 0
    while n != 1 and count < cap:
        n, a = odd_step(n)
        total += a
        count += 1
    need(n == 1, "core verification did not terminate within cap")
    return count, total


def certificate(d):
    """Sufficient universal certificate for every integer n=(2^t+d-1)/3."""
    a = valuation(d)
    core = d >> a
    j0, a0 = trace_one(core)
    j, total, appended = j0, a0, 0
    while 3 ** (j + 1) >= 1 << (a + total):
        j += 1
        total += 2
        appended += 1
    denominator = 1 << (a + total)
    numerator = 3 ** (j + 1)
    cutoff = max(a + total + 1, 1)
    while (1 << cutoff) * (denominator - numerator) <= (4 - d) * denominator:
        cutoff += 1
    exceptions = []
    for t in range(1, cutoff):
        if d >= (1 << t) or ((1 << t) + d - 1) % 3:
            continue
        n = ((1 << t) + d - 1) // 3
        if n > 1:
            steps, y = first_descent(n)
            exceptions.append((t, steps, y))
    # Replay two genuinely high representatives against the symbolic lift.
    lifted = []
    for t in (cutoff, cutoff + 1):
        if ((1 << t) + d - 1) % 3:
            continue
        n = ((1 << t) + d - 1) // 3
        y = n
        for _ in range(j + 1):
            y, _ = odd_step(y)
        expected = 3 ** j * (1 << (t - a - total)) + 1
        need(y == expected and y < n, "high-lift affine certificate")
        lifted.append((t, y))
    need(lifted, "one parity of t is admissible")
    return dict(d=d, core=core, initial_valuation=a, core_steps=j0,
                core_divisions=a0, appended_unit_steps=appended,
                large_steps=j + 1, cutoff=cutoff,
                lambda_numerator=numerator, lambda_denominator=denominator,
                exceptions=exceptions, lifted=lifted)


def main():
    maximum = 1 << 18
    eq_lower = eq_upper = eq_mantissa = 0
    for n in range(1, maximum, 2):
        c, e, _ = carry(n)
        m = 3 * n + 1
        d = m - (1 << (m.bit_length() - 1))
        need(c == jet(m) - 3 * jet(n), "independent z=2 derivative identity")
        need(2 * n + 6 <= c <= 8 * n, "sharp carry bounds")
        need((c == 2 * n + 6) == lower_equality(n), "lower equality language")
        need((e == 0) == (odd_step(n)[0] == 1), "upper equality basin")
        need(e >= 2 * d, "ordinary mantissa bound")
        no_two_zero = "00" not in bin(n)[2:]
        need((e == 2 * d) == no_two_zero, "mantissa equality language")
        eq_lower += c == 2 * n + 6
        eq_upper += e == 0
        eq_mantissa += e == 2 * d
        h = n.bit_length() - 1
        mu = Fraction(jet(n), n)
        need(h - 1 < mu <= h, "normalized first-jet shell bound")
        if n < 4096:
            c2, e2, _ = carry(4 * n + 1)
            need(c2 == 4 * c + 8 and e2 == 4 * e, "inverse prolongation")
            y, a = odd_step(n)
            need(e % (1 << (a + 1)) == 0, "valuation-scaled integral defect")

    family_checks = 0
    for k in range(3, 401):
        base = ((1 << (2 * k)) - 1) // 3
        n = base + 2
        h = base.bit_length() - 1
        need(n.bit_length() == base.bit_length(), "common dyadic shell")
        need(n ^ base == 2, "single low-bit perturbation")
        need(carry(base)[1] == 0 and carry(n)[1] == 12, "near-boundary defects")
        y1, a1 = odd_step(n)
        need(a1 == 1 and y1 == (1 << (2 * k - 1)) + 3, "first bad reset")
        need(carry(y1)[1] == 3 * (1 << (2 * k)) + 4, "macroscopic defect creation")
        y2, a2 = odd_step(y1)
        y3, a3 = odd_step(y2)
        need(a2 == 1 and y2 == 3 * (1 << (2 * k - 2)) + 5, "second rise")
        if k >= 4:
            need(a3 == 4 and y3 == 9 * (1 << (2 * k - 6)) + 1,
                 "three-step return")
            need(y3 < n, "three-step strict descent")
        else:
            need(y3 == 5 and a3 == 5, "small carry collision")
        if k >= 5:
            need(carry(y3)[1] == 54 * (1 << (2 * k - 6)),
                 "defect need not decrease at descending return")
        for r in range(5):
            old = sum((h - s) ** r * (1 << s)
                      for s in range(h + 1) if (base >> s) & 1)
            new = sum((h - s) ** r * (1 << s)
                      for s in range(h + 1) if (n >> s) & 1)
            need(new == old + 2 * (h - 1) ** r, "all fixed normalized jet collisions")
        family_checks += 1

    family_27 = []
    for k in range(3, 81):
        n = ((1 << (2 * k)) + 17) // 3
        need(carry(n)[1] == 36, "27 family constant defect")
        steps, endpoint = first_descent(n)
        need(steps == (37 if k == 3 else 28 if k == 4 else 6),
             "27 family exact first-descent times")
        if k >= 6:
            y = n
            for _ in range(6):
                y, _ = odd_step(y)
            need(y == 243 * (1 << (2 * k - 10)) + 5,
                 "27 family explicit all-height descent formula")
        if k <= 8:
            family_27.append((k, n, steps, endpoint))

    growth_controls = 0
    for length in range(1, 65):
        for k in (1, 2, 5, 20, 80):
            n = ((1 << (length + 2 * k)) + (1 << (length + 1)) - 3) // 3
            need(bin(n)[2:] == '10' * (k - 1) + '1' * (length + 1),
                 "growth family binary language")
            need(carry(n)[1] == (1 << (length + 2)) - 4,
                 "growth family fixed absolute defect")
            y = n
            for j in range(1, length + 1):
                y, a = odd_step(y)
                need(a == 1 and y > n and y == 3 ** j * ((n + 1) >> j) - 1,
                     "growth family exact rising prefix")
            growth_controls += 1

    # PROVED from these finite certificates and the written high-lift argument:
    # every odd n>1 with E(n)<=4096 has a finite smaller odd iterate.
    bound = 4096
    certs = [certificate(d) for d in range(2, bound // 2 + 1, 2) if d % 3 != 1]
    exception_rows = [(c['d'], *row) for c in certs for row in c['exceptions']]
    max_horizon = max([c['large_steps'] for c in certs] + [r[2] for r in exception_rows])
    # Exact E<=12 classification, including its sole four-step exception n=7.
    for n in range(3, maximum, 2):
        if carry(n)[1] <= 12:
            steps, _ = first_descent(n)
            need(steps <= 4 and (steps <= 3 or n == 7), "E<=12 first-descent bound")
    summary = dict(
        status="FINITE-EXACT controls; all-height theorem uses written lift proof",
        odd_control_count=maximum // 2,
        lower_equality_count=eq_lower,
        upper_equality_count=eq_upper,
        mantissa_equality_count=eq_mantissa,
        large_family_count=family_checks,
        defect_bound=bound,
        positive_offset_count=len(certs),
        exception_count=len(exception_rows),
        maximum_core_steps=max(c['core_steps'] for c in certs),
        maximum_cutoff=max(c['cutoff'] for c in certs),
        maximum_selected_horizon=max_horizon,
        maximum_large_horizon=max(c['large_steps'] for c in certs),
        maximum_exception_source_bits=max((((1 << r[1]) + r[0] - 1) // 3).bit_length()
                                          for r in exception_rows),
        maximal_exception_witness=next(r for r in exception_rows if r[2] == max_horizon),
        family_27_first_descent=family_27,
        normalized_near_boundary_growth_controls=growth_controls,
        certificate_sha256=sha256(json.dumps(certs, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
        example_certificates=[{k: v for k, v in c.items() if k not in ('exceptions', 'lifted')}
                              for c in certs if c['d'] in (2, 6, 18, 54, 162, 2048)],
    )
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
