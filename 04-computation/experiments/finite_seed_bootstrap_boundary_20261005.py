"""Exact limits and constructive survivors of finite-seed Collatz bootstraps.

No universal convergence assertion. Run normally and with -O; checks do not
depend on assertions. The literal replay and backward reader are independent
of the affine-carrier compiler.
"""
from dataclasses import dataclass
from itertools import product
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def positive_odd(n):
    need(type(n) is int and n > 0 and n % 2 == 1, 'exact positive odd integer')


def valuation(n, p):
    need(type(n) is int and n > 0 and p in (2, 3), 'positive valuation input')
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def step(n):
    positive_odd(n)
    a = valuation(3*n+1, 2)
    return (3*n+1) >> a, a


def replay(n, word):
    states = [n]
    for a in word:
        n, actual = step(n)
        need(actual == a, 'actual valuation, including the last oddness bit')
        states.append(n)
    return tuple(states)


def backward(hub, word):
    out = hub
    for a in reversed(word):
        numerator = (1 << a)*out-1
        need(numerator > 0 and numerator % 3 == 0, 'literal inverse guard')
        out = numerator//3
        positive_odd(out)
    return out


def carrier(word):
    need(type(word) is tuple and all(type(a) is int and a >= 1 for a in word),
         'positive integer valuation tuple')
    p, q, b = 1, 1, 0
    for a in word:
        p, q, b = 3*p, (1 << a)*q, 3*b+q
    return p, q, b


def log2_mod3(target, depth):
    """Unique exponent modulo 2*3**(depth-1), by one ternary digit at a time."""
    need(type(depth) is int and depth >= 1 and type(target) is int and target % 3,
         'a unit and positive ternary depth')
    exponent = 0 if target % 3 == 1 else 1
    modulus, period = 3, 2
    for _ in range(1, depth):
        modulus *= 3
        candidates = [exponent+j*period for j in range(3)
                      if pow(2, exponent+j*period, modulus) == target % modulus]
        need(len(candidates) == 1, 'unique lifted ternary address')
        exponent = candidates[0]
        period *= 3
    return exponent, period


@dataclass(frozen=True)
class TailReceipt:
    depth: int
    residue: int
    cutoff: int
    hub: int
    template: tuple
    extra: int
    exponent_period: int
    height_threshold: int
    source: int

    @property
    def word(self):
        return self.template[:-1]+(self.template[-1]+self.extra,)


def dyadic_tail_member(depth, residue, cutoff, hub, lift=0):
    """Construct a member of the tail that reaches the supplied hub strictly.

    A root proof for the hub is a separate input obligation, not manufactured
    by this compiler. Increasing lift chooses a later exponent in the same
    exact ternary phase and preserves the dyadic tail.
    """
    need(type(depth) is int and depth >= 1, 'positive binary depth')
    need(type(residue) is int and 0 < residue < 1 << depth and residue % 2,
         'canonical odd residue')
    need(type(cutoff) is int and cutoff >= 0, 'nonnegative tail cutoff')
    need(type(lift) is int and lift >= 0, 'nonnegative exponent lift')
    positive_odd(hub)
    need(hub % 3, 'supplied hub is a ternary unit')
    x, letters = residue, []
    for _ in range(depth):
        x, a = step(x)  # U(1)=1 is allowed only in this finite template.
        letters.append(a)
    word = tuple(letters)
    p, q, b = carrier(word)
    extra, period = log2_mod3(b*pow(q*hub, -1, p) % p, depth)
    if extra == 0:
        extra = period
    threshold = max(cutoff, hub)
    for i in range(depth):
        pi, qi, bi = carrier(word[:i])
        threshold = max(threshold, (qi*hub-bi)//pi)
    while (q*(1 << extra)*hub-b)//p <= threshold:
        extra += period
    extra += lift*period
    numerator = q*(1 << extra)*hub-b
    need(numerator % p == 0, 'full ternary guard retains Q')
    source = numerator//p
    return TailReceipt(depth, residue, cutoff, hub, word, extra, period, threshold, source)


def verify_tail(receipt):
    need(type(receipt) is TailReceipt, 'typed receipt')
    need(all(type(v) is int for v in (receipt.depth, receipt.residue, receipt.cutoff,
                                     receipt.hub, receipt.extra, receipt.exponent_period,
                                     receipt.height_threshold, receipt.source)),
         'exact integer receipt fields')
    need(receipt.depth >= 1 and 0 < receipt.residue < 1 << receipt.depth
         and receipt.residue % 2 and receipt.cutoff >= 0 and receipt.extra >= 1,
         'receipt domain')
    positive_odd(receipt.hub)
    need(receipt.hub % 3 and type(receipt.template) is tuple
         and len(receipt.template) == receipt.depth, 'hub and exact word depth')
    carrier(receipt.template)
    need(receipt.exponent_period == 2*3**(receipt.depth-1), 'declared exponent period')
    n = receipt.source
    positive_odd(n)
    need(n > receipt.height_threshold and n > receipt.cutoff,
         'strict tail and preterminal height bound')
    need(n % (1 << receipt.depth) == receipt.residue, 'specified dyadic source')
    states = replay(n, receipt.word)
    need(states[-1] == receipt.hub and all(x > receipt.hub for x in states[:-1]),
         'first hit of the supplied hub; no root padding')
    need(backward(receipt.hub, receipt.word) == n, 'independent backward identity')
    return states


def l_step(n):
    positive_odd(n)
    need(n % 256 == 219, 'native L guard')
    return (9*n-3)//16


def l_fuel(n):
    positive_odd(n)
    return max(0, (valuation(7*n+3, 2)-4)//4)


def l_ancestor_depth(u):
    positive_odd(u)
    return valuation(7*u+3, 3)//2 if valuation(7*u+3, 2) >= 4 else 0


def l_ancestors(u):
    positive_odd(u)
    out = [u]
    for _ in range(l_ancestor_depth(u)):
        numerator = 16*out[-1]+3
        need(numerator % 9 == 0, 'inverse L integrality')
        parent = numerator//9
        need(l_step(parent) == out[-1], 'inverse retains native source guard')
        out.append(parent)
    return tuple(out)


def root_route(n):
    """Finite control only; callers supply an explicit finite universe."""
    positive_odd(n)
    out = [n]
    for _ in range(10000):
        if n == 1:
            return tuple(out)
        n, _ = step(n)
        out.append(n)
    raise ValueError('finite root control exceeded its explicit cap')


def compositions(total):
    if total == 0:
        yield ()
    else:
        for first in range(1, total+1):
            for tail in compositions(total-first):
                yield (first,)+tail


def contracting_boundary(word):
    p, q, b = carrier(word)
    need(q > p, 'contracting coefficient')
    c = (-b*pow(p, -1, q)) % q
    e = (p*c+b)//q
    first = max(0, (e-c)//(q-p)+1)
    return p, q, b, c, e, first


def checks():
    counts = {}
    count = 0
    for depth in range(1, 7):
        modulus = 3**depth
        period = 2*3**(depth-1)
        independent = {pow(2, e, modulus): e for e in range(period)}
        need(len(independent) == period, 'exhaustive primitive-power control')
        for target, expected in independent.items():
            need(log2_mod3(target, depth) == (expected, period), 'digit log vs power table')
            count += 1
    counts['ternary_log_units_depth_1_to_6'] = count
    count, max_bits = 0, 0
    for depth in range(1, 6):
        for residue in range(1, 1 << depth, 2):
            for hub, cutoff in product((1, 5, 7, 11, 13), (0, 10000)):
                previous = 0
                for lift in (0, 1):
                    row = dyadic_tail_member(depth, residue, cutoff, hub, lift)
                    verify_tail(row)
                    need(row.source > previous, 'later exponent gives a new tail member')
                    previous = row.source
                    max_bits = max(max_bits, row.source.bit_length())
                    count += 1
    counts['tail_receipts'] = count
    counts['largest_constructed_source_bits'] = max_bits
    count = 0
    for n in range(1, 20000, 2):
        x, repeated = n, 0
        while x % 256 == 219:
            y = l_step(x)
            need(0 < y < x and 7*y+3 == 9*(7*x+3)//16, 'native paid conjugacy')
            k = (81*x+53)//128
            target, a = step(k)
            left = replay(x, (1, 2, 1, 1, a+2))
            right = replay(y, (1, 2, a))
            need(left[-1] == target == right[-1], 'literal common-future receipt')
            need(all(z > 1 for z in left[:-1]+right[:-1]), 'no padded first root')
            x, repeated = y, repeated+1
        need(repeated == l_fuel(n), 'complete literal L fuel check')
        # Read inverse parents directly without consulting the valuation formula.
        x, independent = n, [n]
        while (16*x+3) % 9 == 0:
            y = (16*x+3)//9
            if y % 256 != 219:
                break
            independent.append(y)
            x = y
        need(tuple(independent) == l_ancestors(n), 'exact finite ancestry depth')
        count += 1
    counts['literal_L_fuel_and_ancestry_sources'] = count
    seeds = (1, 123, 555)
    closure = sorted({n for seed in seeds for n in l_ancestors(seed)})
    need(closure == [1, 123, 219, 555, 987, 1755], 'finite certified seed closure')
    for seed in seeds:
        need(root_route(seed)[-1] == 1, 'literal finite seed discharge')
    for n in closure:
        need(root_route(n)[-1] == 1, 'independent finite closure control')
    counts['L_only_seeds'] = list(seeds)
    counts['L_only_closure'] = closure
    count, threshold_controls, exceptions = 0, 0, []
    for total in range(1, 11):
        for word in compositions(total):
            p, q, b = carrier(word)
            if q <= p:
                continue
            p, q, b, c, e, first = contracting_boundary(word)
            need(first <= 1, 'finite control of cited one-member theorem')
            if first:
                exceptions.append((word, c))
            for k in range(4):
                n = c+q*k
                x = n
                for i, a in enumerate(word):
                    x, actual = step(x)
                    need(actual == a if i+1 < len(word) else actual >= a,
                         'coarse guard keeps earlier exact valuations')
                nominal = e+p*k
                need(x == nominal >> valuation(nominal, 2), 'coarse oddpart endpoint')
                need((nominal < n) == (k >= first), 'sharp nominal height threshold')
                if k >= first:
                    need(x < n, 'sufficient actual descent')
                threshold_controls += 1
            count += 1
    counts['contracting_words_total_cost_at_most_10'] = count
    counts['threshold_parameter_controls_0_to_3'] = threshold_controls
    counts['words_with_exceptional_least_representative'] = len(exceptions)
    # Decisive scope controls: paying a smaller child is not closure of the domain.
    need(9 % 4 == 1 and step(9)[0] == 7 and 7 % 4 != 1, 'domain closure hostile')
    p, q, b, c, e, first = contracting_boundary((2,))
    need((c, e, first) == (1, 1, 1), 'root fixed-point boundary')
    need(step(1)[0] == 1 and carrier((1,)) == (3, 2, 1), 'coarse equality hostile')
    for bad in (True, 1.0, 0, 2, -1):
        try:
            l_fuel(bad)
        except ValueError:
            pass
        else:
            raise ValueError('invalid positive odd input accepted')
    try:
        dyadic_tail_member(3, 3, 0, 3)
    except ValueError:
        pass
    else:
        raise ValueError('3-divisible hub accepted')
    example = dyadic_tail_member(3, 3, 10000, 5)
    counts['example_tail_member'] = dict(source=example.source, hub=example.hub,
                                       word=example.word, extra=example.extra)
    print(json.dumps(counts, indent=2, sort_keys=True))
    print('PASS: exact finite controls; all-height claims are proved in the companion note.')


if __name__ == '__main__':
    checks()
