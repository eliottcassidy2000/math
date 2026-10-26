#!/usr/bin/env python3
"""Exact odd-valuation coding, Antihydra barrier and hostile controls.

No floating-point arithmetic, convergence oracle, or random-orbit assertion.
The native Turing-machine clock is deliberately outside this checker.
"""
from fractions import Fraction as Q
from itertools import product
from hashlib import sha256
import argparse
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def valuation(x, p=2):
    x = Q(x)
    if x == 0:
        raise ValueError('valuation at zero is infinite, not a finite letter')
    def vi(n):
        n = abs(n)
        count = 0
        while n % p == 0:
            count += 1
            n //= p
        return count
    return vi(x.numerator) - vi(x.denominator)


def inverse(x, k, kind):
    if not isinstance(k, int) or k < 1:
        raise ValueError('valuation letters must be positive integers')
    if valuation(x) != 0:
        raise ValueError('endpoint must be a 2-adic odd unit')
    if kind == 'collatz':
        return (2 ** k * Q(x) - 1) / 3
    if kind == 'twisted':
        return Q(1, 3) + Q(2, 3) ** k * x
    raise ValueError(kind)


def forward(x, kind):
    x = Q(x)
    if valuation(x) != 0:
        raise ValueError('source must be a 2-adic odd unit')
    if kind == 'collatz':
        t = 3 * x + 1
    elif kind == 'twisted':
        t = 3 * x - 1
    else:
        raise ValueError(kind)
    k = valuation(t)
    need(k >= 1, 'odd source did not yield positive valuation')
    y = t / (2 ** k)
    if kind == 'twisted':
        y *= 3 ** (k - 1)
    need(valuation(y) == 0, 'odd return lost oddness')
    return y, k


def decode_prefix(word, tail, kind):
    x = Q(tail)
    for k in reversed(word):
        x = inverse(x, k, kind)
    return x


def prefix_affine(word, kind):
    """Return exact (M,B) with decoded source = B + M * tail."""
    m, b = Q(1), Q(0)
    for k in reversed(word):
        if not isinstance(k, int) or k < 1:
            raise ValueError('valuation letters must be positive integers')
        if kind == 'collatz':
            factor, offset = Q(2 ** k, 3), Q(-1, 3)
        elif kind == 'twisted':
            factor, offset = Q(2, 3) ** k, Q(1, 3)
        else:
            raise ValueError(kind)
        m, b = factor * m, offset + factor * b
    return m, b


def periodic(word, kind):
    if not word:
        raise ValueError('period must be nonempty')
    m, b = prefix_affine(word, kind)
    need(m != 1, 'powers of two and three cannot coincide')
    return b / (1 - m)


def residue(x, bits):
    x = Q(x)
    mod = 1 << bits
    if x.denominator % 2 == 0:
        raise ValueError('not a 2-adic integer')
    return x.numerator * pow(x.denominator, -1, mod) % mod


def debt_summary(word):
    delta, minimum = 0, 0
    for k in word:
        delta += 2 * k - 3
        minimum = min(minimum, delta)
    return delta, minimum


def compose_debt(left, right):
    delta, minimum = left
    change, trough = right
    return delta + change, min(minimum, delta + trough)


def common_future(x, y, kind, max_steps):
    """Finite witnessed common future; None means unestablished, not absent."""
    left = {}
    state = Q(x)
    for i in range(max_steps + 1):
        if state in left:
            break
        left[state] = i
        state, _ = forward(state, kind)
    state, seen = Q(y), set()
    for j in range(max_steps + 1):
        if state in left:
            return state, left[state], j
        if state in seen:
            break
        seen.add(state)
        state, _ = forward(state, kind)
    return None


def reverse_barrier(word):
    """Least a ensuring a_j >= 1 at EVERY vertex, including terminal."""
    requirement = 1
    for k in reversed(word):
        requirement = max(1, requirement + 3 - 2 * k)
    return requirement


def prefix_barrier(word):
    delta = 0
    worst = 0
    for k in word:
        delta += 2 * k - 3
        worst = max(worst, -delta)
    return 1 + worst


def admissible(word, a):
    if a < 1:
        return False
    for k in word:
        a += 2 * k - 3
        if a < 1:
            return False
    return True


def twisted_integer_inverse(y):
    if not isinstance(y, int) or y % 2 == 0:
        raise ValueError('expected odd integer')
    r = valuation(y, 3)
    m = y // 3 ** r
    k = r + 1
    numerator = 2 ** k * m + 1
    return None if numerator % 3 else (numerator // 3, k)


def macro_step(a, b):
    """Parameter transition only; no claimed native machine-step clock."""
    if a < 0 or b < 0:
        raise ValueError('invalid tape-parameter state')
    if b % 2 == 0:
        return a + 2, 3 * (b // 2) + 2
    if a == 0:
        return None
    return a - 1, 3 * ((b + 1) // 2)


def simulate_blank(steps):
    a, b = 0, 4
    initial = [(a, b)]
    while (b + 4) % 2 == 0:
        a, b = macro_step(a, b)
        initial.append((a, b))
    need((a, b + 4) == (6, 27), 'blank-to-odd section')
    x, total, min_a, min_at = b + 4, 0, a, 0
    letters = []
    digest = sha256()
    for j in range(steps):
        if a == 0:
            return {'status': 'HIT_BARRIER', 'odd_returns': j,
                    'a': a, 'x_hex': hex(x), 'initial_macro_states': initial}
        t = 3 * x - 1
        k = (t & -t).bit_length() - 1
        y = 3 ** (k - 1) * (t >> k)
        need(y > x, 'positive integer twisted orbit must grow')
        # Independently replay every native-parameter macro in this return.
        aa, bb = a, x - 4
        for _ in range(k):
            state = macro_step(aa, bb)
            need(state is not None, 'legal odd guard failed within return')
            aa, bb = state
        a += 2 * k - 3
        need((aa, bb + 4) == (a, y), 'odd-return compression mismatch')
        total += k
        need(a == 6 + 2 * total - 3 * (j + 1), 'barrier ledger')
        if a < min_a:
            min_a, min_at = a, j + 1
        x = y
        letters.append(k)
        digest.update(str(k).encode('ascii') + b',')
    return {'status': 'FINITE_PREFIX_ONLY', 'odd_returns': steps,
            'valuation_sum': total, 'a': a, 'min_a': min_a,
            'min_a_first_index': min_at, 'x_bit_length': x.bit_length(),
            'x_sha256_hex_encoding': sha256(hex(x).encode('ascii')).hexdigest(),
            'valuation_csv_sha256': digest.hexdigest(),
            'first_35_valuations': letters[:35],
            'reverse_requirement_for_prefix': reverse_barrier(letters),
            'initial_macro_states': initial}


def finite_ruin_mass(a, steps):
    """Exact mass of ruin by horizon under iid P(k)=2^-k.

    At time t, resource > steps-t cannot hit zero by the horizon because
    every increment is >= -1. Aggregate all those infinitely many letters
    into a safe sink instead of dropping a tail probability.
    """
    alive = {a: Q(1)}
    dead, safe = Q(0), Q(0)
    for t in range(steps):
        nxt = {}
        remaining = steps - t - 1
        for height, mass in alive.items():
            if height == 0:
                dead += mass
                continue
            # New height = height + 2*k - 3. Relevant iff <= remaining.
            kmax = max(0, (remaining - height + 3) // 2)
            for k in range(1, kmax + 1):
                h = height + 2 * k - 3
                weight = mass / 2 ** k
                if h == 0:
                    dead += weight
                else:
                    nxt[h] = nxt.get(h, Q(0)) + weight
            safe += mass / 2 ** kmax
        alive = nxt
        need(dead + safe + sum(alive.values(), Q(0)) == 1,
             'finite probability accounting lost an infinite tail')
    need(not alive, 'horizon accounting should have no active states')
    return dead


def check():
    counts = {'prefix_words': 0, 'branch_round_trips': 0,
              'cylinder_checks': 0, 'barrier_checks': 0,
              'isometry_pairs': 0, 'integer_inverse_checks': 0,
              'periodic_words': 0, 'exceptional_prefixes': 0}
    tails = [Q(1), Q(3), Q(-7, 5)]
    for length in range(6):
        for word in product(range(1, 5), repeat=length):
            counts['prefix_words'] += 1
            A = sum(word)
            R = reverse_barrier(word)
            need(R == 1 - debt_summary(word)[1], 'inherited debt-monoid requirement')
            need(R == prefix_barrier(word), 'reverse/prefix barrier mismatch')
            for a in range(0, 9):
                need(admissible(word, a) == (a >= R), 'minimal barrier requirement')
                counts['barrier_checks'] += 1
            for split in range(length + 1):
                u, v = word[:split], word[split:]
                delta = 2 * sum(u) - 3 * len(u)
                need(compose_debt(debt_summary(u), debt_summary(v)) == debt_summary(word),
                     'chronological debt composition')
                need(R == max(reverse_barrier(u), reverse_barrier(v) - delta),
                     'barrier concatenation identity')
            for kind in ['collatz', 'twisted']:
                m, b = prefix_affine(word, kind)
                need(valuation(m) == A, 'prefix linear multiplier valuation')
                decoded = [decode_prefix(word, tail, kind) for tail in tails]
                for tail, x in zip(tails, decoded):
                    need(x == b + m * tail, 'affine decoding identity')
                    y = x
                    for expected in word:
                        y, actual = forward(y, kind)
                        need(actual == expected, 'exact prefix letter')
                        counts['branch_round_trips'] += 1
                    need(y == tail, 'exact prefix tail')
                need(residue(decoded[0], A + 1) == residue(decoded[1], A + 1),
                     'same-cylinder congruence')
                need(valuation(decoded[1] - decoded[0]) == A + 1,
                     'cylinder radius sharpness')
                counts['cylinder_checks'] += 1
                exceptional = Q(-1 if kind == 'collatz' else 1, 3)
                y = decode_prefix(word, exceptional, kind)
                for expected in word:
                    y, actual = forward(y, kind)
                    need(actual == expected, 'finite exceptional prefix')
                try:
                    forward(y, kind)
                except ValueError as err:
                    need('infinite' in str(err), 'wrong singular error')
                else:
                    raise RuntimeError('singular valuation must not be silently coded')
                counts['exceptional_prefixes'] += 1
                if word:
                    cyc = periodic(word, kind)
                    y = cyc
                    for expected in word:
                        y, actual = forward(y, kind)
                        need(actual == expected, 'periodic word coding')
                    need(y == cyc, 'periodic closure')
                    if kind == 'twisted':
                        need(Q(1, 3) < cyc <= 1, 'twisted periodic real range')
                    counts['periodic_words'] += 1
            if length <= 3:
                for p in range(1, 5):
                    for q in range(p + 1, 5):
                        values = []
                        for kind in ['collatz', 'twisted']:
                            fixed_one = Q(-1 if kind == 'collatz' else 1)
                            x = decode_prefix(word + (p,), fixed_one, kind)
                            y = decode_prefix(word + (q,), fixed_one, kind)
                            values.append(valuation(x - y))
                        need(values == [A + p, A + p], 'cross-map valuation isometry')
                        counts['isometry_pairs'] += 1
    images = {}
    for x in range(-999, 1000, 2):
        y, k = forward(x, 'twisted')
        need(y.denominator == 1, 'integer twisted closure')
        y = int(y)
        need(twisted_integer_inverse(y) == (x, k), 'forced integer inverse')
        need(y not in images, 'twisted integer merger detected')
        images[y] = x
        need(y > x if x > 1 else (y == x if x == 1 else y < x),
             'integer drift direction')
        counts['integer_inverse_checks'] += 1
    for y in range(-999, 1000, 2):
        answer = twisted_integer_inverse(y)
        candidates = [(inverse(y, k, 'twisted'), k) for k in range(1, 13)]
        candidates = [(int(x), k) for x, k in candidates if x.denominator == 1]
        need(candidates == ([] if answer is None else [answer]),
             'integer inverse uniqueness enumeration')
        counts['integer_inverse_checks'] += 1
    need(common_future(3, 13, 'collatz', 8) == (Q(5), 1, 1),
         'witnessed odd-step shared ancestry')
    need(common_future(1, Q(-5, 3), 'collatz', 8) is None,
         'common-prefix hostile cycles are disjoint')
    need(common_future(27, 135, 'twisted', 8) == (Q(135), 1, 0),
         'twisted same-ray common future')
    need(forward(3, 'collatz')[0] == forward(13, 'collatz')[0] == 5,
         'Collatz merger hostile')
    need(forward(1, 'collatz') == (1, 2), 'positive Collatz fixed point')
    need(forward(Q(-5, 3), 'collatz') == (-1, 2), 'same-prefix disconnected tails')
    need(periodic((2,), 'twisted') == Q(3, 5), 'integer preservation hostile')
    need(periodic((1,), 'collatz') == -1, 'negative integer coding hostile')
    for a in range(1, 16):
        x = 1 + 2 ** (a + 1)
        y = x
        for _ in range(a):
            y, k = forward(y, 'twisted')
            need(k == 1, 'positive integer barrier hostile word')
        need(y == 1 + 2 * 3 ** a, 'positive integer barrier hostile endpoint')
        need(not admissible((1,) * a, a), 'counter must hit zero')
    # Real-normalized affine identities, every stage, without logarithms.
    for kind, x0 in [('collatz', Q(27)), ('twisted', Q(27))]:
        x, word, previous = x0, [], x0
        for _ in range(60):
            x, k = forward(x, kind)
            word.append(k)
            m, b = prefix_affine(word, kind)
            current = m * x
            need(current == x0 - b, 'normalized affine ledger')
            if kind == 'twisted':
                need(x0 - 1 < current < previous, 'twisted real potential bound')
            else:
                need(current > previous, 'Collatz real normalized ledger')
            previous = current
    ruin = {str(h): str(finite_ruin_mass(6, h)) for h in [1, 6, 12, 24, 48]}
    for value in ruin.values():
        p = Q(value)
        need(9 - p > 0 and (9 - p) ** 2 > 80,
             'finite ruin mass must be less than exact limit 9-4*sqrt(5)')
    need(ruin['6'] == '1/64', 'first possible ruin requires six ones')
    return {'status': 'PASS', 'controls': counts, 'finite_horizon_ruin_mass_a6': ruin,
            'hostiles': ['zero/infinite valuation rejected',
                         'finite prefix collision is not common future',
                         'Collatz integer siblings do not transport to twisted integers',
                         'twisted positive integer barrier-hitting words 1^a',
                         'periodic words produce negative or noninteger 2-adic points',
                         'coefficient drift is not positive integer convergence']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--steps', type=int, default=10000)
    parser.add_argument('--output')
    args = parser.parse_args()
    if args.steps < 0:
        raise ValueError('steps must be nonnegative')
    report = check()
    report['blank_antihydra_macro_prefix'] = simulate_blank(args.steps)
    text = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        with open(args.output, 'w') as stream:
            stream.write(text)
    print(text, end='')
