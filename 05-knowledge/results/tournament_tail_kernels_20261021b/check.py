"""Event-driven future expiration and bounded-divisor multiplicity ladders."""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import factorial, gcd, comb
from pathlib import Path
import argparse
import json
import sys
import time
import statistics
import tracemalloc

sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / 'tournament_enumeration_20261021a'
CHECKS = Counter()


def need(condition, label):
    CHECKS[label] += 1
    if not condition:
        raise ValueError(label)


def digest(path):
    return sha256(Path(path).read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def load(path, name):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def natural(x):
    if type(x) is not int or x < 0:
        raise ValueError('literal nonnegative integer required')
    return x


def predecessor(ceiling, start):
    """Largest allowed odd length no greater than ceiling; zero if absent."""
    if ceiling < start:
        return 0
    return ceiling - 1 + ceiling % 2


def count(bound, *, core=False, no_fixed=False, calendar=True):
    natural(bound)
    if any(type(x) is not bool for x in (core, no_fixed, calendar)):
        raise ValueError('literal mode flags required')
    start = 3 if no_fixed else 1
    lengths = list(range(start, bound + 1, 2))[::-1]
    width = max(1, bound.bit_length())
    digit_mask = (1 << width) - 1
    ncoords = max(0, (bound - start) // 2 + 1)
    masks = [(1 << (width * n)) - 1 for n in range(ncoords + 1)]
    fac = factorial(bound)
    common = fac // (fac & -fac)
    outputs = [0] * (bound + 1)
    if not lengths:
        outputs[0] = common
        return [Fraction(x, common) for x in outputs], {'transitions': 0, 'rows': [],
                    'ladder_divisions': 0, 'forced_zero_layers_omitted': 0,
                    'active_states': 0, 'peak_queued_states': 1}
    queues = {ell: defaultdict(int) for ell in lengths}
    queues[lengths[0]][0, 0] = common
    queued = 1
    peak_queued = 1
    rows = []
    divisions = 0
    skipped = 0
    for ell in lengths:
        states = queues.pop(ell)
        queued -= len(states)
        idx = (ell - start) // 2
        grow = sum(gcd(ell, j) << (width * i)
                   for i, j in enumerate(range(start, ell, 2)))
        row_transitions = 0
        for (degree, response), weight in states.items():
            charge = (response >> (width * idx)) & digit_mask
            limit = (bound - degree) // ell
            factor = weight
            for k in range(limit + 1):
                if k:
                    two = (k & -k).bit_length() - 1
                    shift = (ell - 1) // 2 + ell * (k - 1) + charge - two
                    if shift < 0:
                        raise ValueError('negative ladder shift')
                    factor, remainder = divmod(factor << shift, ell * (k >> two))
                    if remainder:
                        raise ValueError('nonintegral multiplicity ladder')
                    divisions += 1
                updated = response + k * grow
                excesses = range(limit - k + 1) if core and k else (0,)
                for excess in excesses:
                    d = degree + ell * (k + excess)
                    newell = predecessor(min(ell - 2, bound - d), start) if calendar else ell - 2
                    coefficient = (-1) ** excess * comb(k + excess - 1, excess) if k else 1
                    value = factor * coefficient
                    if newell < start:
                        outputs[d] += value
                    else:
                        digits = max(0, (min(newell, bound - d) - start) // 2 + 1)
                        key = d, updated & masks[digits]
                        target = queues[newell]
                        old = target.get(key, 0)
                        new = old + value
                        if new:
                            target[key] = new
                        elif key in target:
                            del target[key]
                        queued += bool(new) - bool(old)
                    if calendar:
                        # An ordinary stagewise history would traverse each of
                        # these intervening layers solely by the k=0 branch.
                        stopped = newell if newell >= start else start - 2
                        skipped += (ell - 2 - stopped) // 2
                    row_transitions += 1
                peak_queued = max(peak_queued, queued)
        rows.append({'cycle_length': ell, 'states': len(states),
                     'transitions': row_transitions, 'queued_states': queued})
    if queued or queues:
        raise ValueError('unconsumed scheduled state')
    return [Fraction(x, common) for x in outputs], {
        'transitions': sum(r['transitions'] for r in rows),
        'active_states': sum(r['states'] for r in rows),
        'peak_active_states': max([1] + [r['states'] for r in rows]),
        'peak_queued_states': peak_queued,
        'ladder_divisions': divisions,
        'forced_zero_layers_omitted': skipped,
        'rows': rows}


def inherited():
    spec = spec_from_file_location('old_counter', PRIOR / 'check.py')
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def poly_count(bound, *, core=False, no_fixed=False):
    """Gaussian-normalized cycle-count polynomials, encoded with signed digits."""
    natural(bound)
    if type(core) is not bool or type(no_fixed) is not bool:
        raise ValueError('literal flags')
    start = 3 if no_fixed else 1
    lengths = list(range(start, bound + 1, 2))[::-1]
    width = max(1, bound.bit_length())
    digit_mask = (1 << width) - 1
    ncoords = max(0, (bound - 3) // 2 + 1)
    masks = [(1 << (width * n)) - 1 for n in range(ncoords + 1)]
    common = factorial(bound)
    coefficient_bits = comb(bound, 2) + 2 * bound + common.bit_length() + 2
    coefficient_mask = (1 << coefficient_bits) - 1
    half_radix = 1 << (coefficient_bits - 1)
    # Payload (lowest occupied cycle count, packed signed polynomial).
    states = {(0, 0): (0, common)}
    rows = []
    divisions = 0
    for ell in lengths:
        idx = (ell - 3) // 2
        grow = sum((gcd(ell, j) - 1) << (width * i)
                   for i, j in enumerate(range(3, ell, 2)))
        new = {}
        transitions = 0
        for (degree, response), (offset, weight) in states.items():
            charge = (response >> (width * idx)) & digit_mask if ell >= 3 else 0
            limit = (bound - degree) // ell
            factor = weight
            for k in range(limit + 1):
                if k:
                    shift = ((ell - 1) // 2) * (2 * k - 1) + charge
                    factor, rem = divmod(factor << shift, ell * k)
                    if rem:
                        raise ValueError('nonintegral normalized cycle polynomial')
                    divisions += 1
                updated = response + k * grow
                excesses = range(limit - k + 1) if core and k else (0,)
                for excess in excesses:
                    d = degree + ell * (k + excess)
                    digits = max(0, (min(ell - 2, bound - d) - 3) // 2 + 1)
                    key = d, updated & masks[digits]
                    coefficient = (-1) ** excess * comb(k + excess - 1, excess) if k else 1
                    value = factor * coefficient
                    off = offset + k
                    if key in new:
                        prior_offset, prior_value = new[key]
                        base = min(prior_offset, off)
                        value = (prior_value << ((prior_offset - base) * coefficient_bits)) + (
                            value << ((off - base) * coefficient_bits))
                        off = base
                    if value:
                        leading_zeros = ((value & -value).bit_length() - 1) // coefficient_bits
                        new[key] = off + leading_zeros, value >> (leading_zeros * coefficient_bits)
                    elif key in new:
                        del new[key]
                    transitions += 1
        states = new
        rows.append({'cycle_length': ell, 'states': len(states), 'transitions': transitions,
                     'payload_bits': sum(abs(v).bit_length() for _, v in states.values())})
    answer = [Fraction(0)] * (bound + 1)
    for (degree, response), (offset, value) in states.items():
        if response:
            raise ValueError('unfinished normalized response')
        total = 0
        while value:
            digit = value & coefficient_mask
            if digit >= half_radix:
                digit -= 1 << coefficient_bits
            total += digit << comb(offset, 2)
            value = (value - digit) >> coefficient_bits
            offset += 1
        answer[degree] += Fraction(total, common)
    return answer, {'transitions': sum(r['transitions'] for r in rows),
                    'peak_states': max([1] + [r['states'] for r in rows]),
                    'peak_payload_bits': max([common.bit_length()] + [r['payload_bits'] for r in rows]),
                    'coefficient_bits': coefficient_bits, 'ladder_divisions': divisions, 'rows': rows}


def list_poly_count(bound, *, core=False, no_fixed=False, terminal=True):
    """Same normalized polynomial module with ordinary sparse coefficient maps."""
    natural(bound)
    if any(type(x) is not bool for x in (core, no_fixed, terminal)):
        raise ValueError('literal flags')
    start = 3 if no_fixed else 1
    lengths = list(range(start, bound + 1, 2))[::-1]
    width = max(1, bound.bit_length())
    digit_mask = (1 << width) - 1
    masks = [(1 << (width * n)) - 1 for n in range(max(0, (bound - 3) // 2 + 1) + 1)]
    common = factorial(bound)
    states = {(0, 0): {0: common}}
    outputs = [0] * (bound + 1)
    rows = []
    divisions = 0
    for ell in lengths:
        idx = (ell - 3) // 2
        grow = sum((gcd(ell, j) - 1) << (width * i)
                   for i, j in enumerate(range(3, ell, 2)))
        new = {}
        transitions = 0
        for (degree, response), weight in states.items():
            charge = (response >> (width * idx)) & digit_mask if ell >= 3 else 0
            limit = (bound - degree) // ell
            factor = weight
            for k in range(limit + 1):
                if k:
                    shift = ((ell - 1) // 2) * (2 * k - 1) + charge
                    f = {}
                    for t, c in factor.items():
                        v, rem = divmod(c << shift, ell * k)
                        if rem:
                            raise ValueError('coefficient division')
                        f[t] = v
                        divisions += 1
                    factor = f
                updated = response + k * grow
                if terminal and ell == start:
                    scalar = sum(c << comb(t + k, 2) for t, c in factor.items())
                    for excess in (range(limit - k + 1) if core and k else (0,)):
                        d = degree + ell * (k + excess)
                        coefficient = (-1) ** excess * comb(k + excess - 1, excess) if k else 1
                        outputs[d] += coefficient * scalar
                        transitions += 1
                    continue
                for excess in (range(limit - k + 1) if core and k else (0,)):
                    d = degree + ell * (k + excess)
                    digits = max(0, (min(ell - 2, bound - d) - 3) // 2 + 1)
                    key = d, updated & masks[digits]
                    coefficient = (-1) ** excess * comb(k + excess - 1, excess) if k else 1
                    target = new.setdefault(key, {})
                    for t, c in factor.items():
                        target[t + k] = target.get(t + k, 0) + coefficient * c
                    transitions += 1
        states = {key: {t: c for t, c in p.items() if c} for key, p in new.items()}
        states = {key: p for key, p in states.items() if p}
        rows.append({'cycle_length': ell, 'states': len(states), 'transitions': transitions,
                     'coefficient_cells': sum(len(p) for p in states.values()),
                     'payload_bits': sum(abs(c).bit_length() for p in states.values() for c in p.values())})
    answer = [Fraction(v, common) for v in outputs]
    for (degree, response), weight in states.items():
        if response:
            raise ValueError('unfinished response')
        answer[degree] = Fraction(sum(c << comb(t, 2) for t, c in weight.items()), common)
    return answer, {'transitions': sum(r['transitions'] for r in rows),
                    'peak_states': max([1] + [r['states'] for r in rows]),
                    'peak_coefficient_cells': max([1] + [r['coefficient_cells'] for r in rows]),
                    'peak_payload_bits': max([1] + [r['payload_bits'] for r in rows]),
                    'ladder_divisions': divisions, 'rows': rows}


def module_dp(bound, *, core=False, no_fixed=False):
    """Degree-local factorial polynomial module, with final functional fold."""
    natural(bound)
    if type(core) is not bool or type(no_fixed) is not bool:
        raise ValueError('literal flags')
    return _module_engine(bound, core, no_fixed)


def _module_engine(bound, core, no_fixed, initial=None, highest=None, settled=None):
    start = 3 if no_fixed else 1
    width = max(1, bound.bit_length())
    digit_mask = (1 << width) - 1
    masks = [(1 << (width * n)) - 1 for n in range(bound // 2 + 2)]
    facts = [factorial(d) for d in range(bound + 1)]
    states = {(0, 0): {0: 1}} if initial is None else initial
    outputs = [0] * (bound + 1) if settled is None else settled
    rows = []
    precomputed = 0
    coefficient_updates = 0
    for ell in range((bound if bound % 2 else bound - 1) if highest is None else highest, start - 1, -2):
        idx = (ell - 3) // 2
        grow = sum((gcd(ell, j) - 1) << (width * i)
                   for i, j in enumerate(range(3, ell, 2)))
        tables = {}
        for degree, _ in states:
            if degree in tables:
                continue
            entries = []
            limit = (bound - degree) // ell
            for k in range(limit + 1):
                base_degree = degree + ell * k
                multiplier, rem = divmod(facts[base_degree], facts[degree] * ell ** k * facts[k])
                if rem:
                    raise ValueError('local normalized transfer nonintegral')
                precomputed += 1
                extra = []
                for excess in (range(limit - k + 1) if core and k else (0,)):
                    d = base_degree + ell * excess
                    coefficient = facts[d] // facts[base_degree]
                    if excess:
                        coefficient *= (-1) ** excess * comb(k + excess - 1, excess)
                    digits = max(0, (min(ell - 2, bound - d) - 3) // 2 + 1)
                    extra.append((d, coefficient, masks[digits]))
                entries.append((k, ((ell - 1) // 2) * k * k, multiplier, k * grow, extra))
            tables[degree] = entries
        new = {}
        transitions = 0
        for (degree, response), weight in states.items():
            charge = (response >> (width * idx)) & digit_mask if ell >= 3 else 0
            for k, base_shift, multiplier, offset, extra in tables[degree]:
                shift = base_shift + k * charge
                if ell == start:
                    scalar = multiplier * sum(c << comb(t + k, 2) for t, c in weight.items())
                    scalar <<= shift
                    for d, coefficient, mask in extra:
                        outputs[d] += scalar * coefficient
                        transitions += 1
                    coefficient_updates += len(weight)
                    continue
                factor = [(t + k, (c * multiplier) << shift) for t, c in weight.items()]
                updated = response + offset
                for d, coefficient, mask in extra:
                    target = new.setdefault((d, updated & mask), {})
                    for t, c in factor:
                        target[t] = target.get(t, 0) + coefficient * c
                    coefficient_updates += len(factor)
                    transitions += 1
        states = {key: {t: c for t, c in p.items() if c} for key, p in new.items()}
        states = {key: p for key, p in states.items() if p}
        rows.append({'cycle_length': ell, 'states': len(states), 'transitions': transitions,
                     'coefficient_cells': sum(len(p) for p in states.values()),
                     'payload_bits': sum(abs(c).bit_length() for p in states.values() for c in p.values())})
    answer = [Fraction(v, facts[d]) for d, v in enumerate(outputs)]
    for (degree, response), weight in states.items():
        if response:
            raise ValueError('unfinished response')
        answer[degree] += Fraction(sum(c << comb(t, 2) for t, c in weight.items()), facts[degree])
    return answer, {'transitions': sum(r['transitions'] for r in rows),
                    'peak_states': max([1] + [r['states'] for r in rows]),
                    'peak_coefficient_cells': max([1] + [r['coefficient_cells'] for r in rows]),
                    'peak_payload_bits': max([1] + [r['payload_bits'] for r in rows]),
                    'coefficient_updates': coefficient_updates,
                    'precomputed_transfers': precomputed, 'rows': rows}


def hybrid_dp(bound, *, core=False, switch=5):
    """Incoming scalar degree-local states, followed by a polynomial tail.

    Fixed points must be present: their live R_1 coordinate supplies t.
    """
    natural(bound)
    if type(core) is not bool or type(switch) is not int or switch < 1 or switch % 2 != 1:
        raise ValueError('literal flag and positive odd switch required')
    width = max(1, bound.bit_length())
    digit_mask = (1 << width) - 1
    masks = [(1 << (width * n)) - 1 for n in range(bound // 2 + 2)]
    facts = [factorial(d) for d in range(bound + 1)]
    twos = [(f & -f).bit_length() - 1 for f in facts]
    odds = [f >> v for f, v in zip(facts, twos)]
    states = {(0, 0): 1}
    rows = []
    precomputed = 0
    for ell in range(bound if bound % 2 else bound - 1, switch, -2):
        idx = (ell - 1) // 2
        grow = sum(gcd(ell, j) << (width * i) for i, j in enumerate(range(1, ell, 2)))
        tables = {}
        for degree, _ in states:
            if degree in tables:
                continue
            entries = []
            limit = (bound - degree) // ell
            for k in range(limit + 1):
                base_shift = k * (ell - 1) // 2 + ell * comb(k, 2) - twos[k]
                for excess in (range(limit - k + 1) if core and k else (0,)):
                    d = degree + ell * (k + excess)
                    multiplier, rem = divmod(odds[d], odds[degree] * ell ** k * odds[k])
                    if rem:
                        raise ValueError('incoming local transfer division')
                    if excess:
                        multiplier *= (-1) ** excess * comb(k + excess - 1, excess)
                    digits = max(0, (min(ell - 2, bound - d) - 1) // 2 + 1)
                    entries.append((k, d, base_shift, multiplier, masks[digits], k * grow))
                    precomputed += 1
            tables[degree] = entries
        new = defaultdict(int)
        transitions = 0
        for (degree, response), weight in states.items():
            charge = (response >> (width * idx)) & digit_mask
            for k, d, shift, multiplier, mask, offset in tables[degree]:
                new[d, (response + offset) & mask] += (weight * multiplier) << (shift + k * charge)
                transitions += 1
        states = {key: c for key, c in new.items() if c}
        rows.append({'cycle_length': ell, 'states': len(states), 'transitions': transitions})
    converted = {}
    outputs = [0] * (bound + 1)
    top = min(switch, bound if bound % 2 else bound - 1)
    for (d, response), weight in states.items():
        if d == bound:
            outputs[d] += weight << twos[d]
            continue
        t = response & digit_mask
        charge_count = max(0, (min(top, bound - d) - 3) // 2 + 1)
        residual = 0
        for i in range(charge_count):
            charge = ((response >> (width * (i + 1))) & digit_mask) - t
            if charge < 0:
                raise ValueError('negative residual gcd charge')
            residual += charge << (width * i)
        exponent = twos[d] - comb(t, 2)
        if exponent >= 0:
            coefficient = weight << exponent
        else:
            coefficient, rem = divmod(weight, 1 << -exponent)
            if rem:
                raise ValueError('Gaussian normalization not integral')
        target = converted.setdefault((d, residual), {})
        target[t] = target.get(t, 0) + coefficient
    answer, tail = _module_engine(bound, core, False, converted, top, outputs)
    return answer, {'transitions': sum(r['transitions'] for r in rows) + tail['transitions'],
                    'precomputed_transfers': precomputed + tail['precomputed_transfers'],
                    'peak_states': max([tail['peak_states']] + [r['states'] for r in rows]),
                    'tail_statistics': {k: v for k, v in tail.items() if k != 'rows'},
                    'switch': switch, 'rows': rows + tail['rows']}


def probe(bound, core):
    old = inherited()
    funcs = [('old_packed_factorials', lambda: old.packed_dp(bound, core=core)),
             ('stagewise_ladder', lambda: count(bound, core=core, calendar=False)),
             ('calendar_ladder', lambda: count(bound, core=core)),
             ('polynomial_payload', lambda: poly_count(bound, core=core)),
             ('list_polynomial', lambda: list_poly_count(bound, core=core)),
             ('local_polynomial', lambda: module_dp(bound, core=core))]
    answer = None
    records = []
    for name, fn in funcs:
        stamp = time.perf_counter()
        result, stats = fn()
        elapsed = time.perf_counter() - stamp
        if answer is None:
            answer = result
        elif result != answer:
            raise ValueError('coefficient disagreement')
        records.append({'algorithm': name, 'seconds': elapsed, 'statistics': stats})
        print(name, elapsed, stats['transitions'], flush=True)
    return {'bound': bound, 'core': core, 'runs': records}


def controls():
    old = inherited()
    unmerged = old.inherited()
    local = load(HERE / 'incoming_local.py', 'incoming_local_control')
    for bound in range(25):
        for core in (False, True):
            for no_fixed in (False, True):
                reference, _ = old.packed_dp(bound, core=core, no_fixed=no_fixed)
                for kernel in (list_poly_count, module_dp):
                    answer, stats = kernel(bound, core=core, no_fixed=no_fixed)
                    need(answer == reference, 'polynomial kernels all four modes')
                plain, _ = list_poly_count(bound, core=core, no_fixed=no_fixed, terminal=False)
                need(plain == reference, 'terminal functional before excess equals full polynomial')
                incoming, _ = local.local_dp(bound, core=core, no_fixed=no_fixed)
                need(incoming == reference, 'frozen incoming local-scale control')
                if bound <= 18:
                    need(reference == unmerged.direct_series(bound, core=core, no_fixed=no_fixed),
                         'independent unmerged rational partition expansion')
                if bound <= 12:
                    radix, _ = poly_count(bound, core=core, no_fixed=no_fixed)
                    need(radix == reference, 'signed radix prototype exact small control')
                if bound <= 16:
                    for calendar in (False, True):
                        scheduled, stats = count(bound, core=core, no_fixed=no_fixed, calendar=calendar)
                        need(scheduled == reference, 'bounded-divisor ladder and scheduled zero layers')
                        need(stats['forced_zero_layers_omitted'] >= 0, 'zero-layer shortcut counter boundary')
                if not no_fixed:
                    for switch in (1, 3, 5, 9):
                        mixed, _ = hybrid_dp(bound, core=core, switch=switch)
                        need(mixed == reference, 'exact local scalar-to-polynomial change of basis')
    # Independent raw-history arithmetic and every possible future insertion.
    for mass in range(31):
        for parts in unmerged.odd_partitions(mass):
            t = len(parts)
            e = sum((a - 1) // 2 for a in parts) + sum(
                gcd(a, b) for i, a in enumerate(parts) for b in parts[i + 1:])
            g = e - comb(t, 2)
            need(g >= 0, 'Gaussian-normalized exponent is nonnegative')
            weight = unmerged.cycle_weight(parts)
            need((factorial(mass) * weight / (1 << comb(t, 2))).denominator == 1,
                 'local full factorial clears normalized history')
            if all(a >= 3 for a in parts):
                need(2 * g <= t * (mass - t), 'nontrivial normalized exponent bound')
                need(9 * g <= mass * mass, 'sharp quadratic bound before unit stage')
            for ell in (1, 3, 5, 7):
                response = sum(gcd(a, ell) for a in parts)
                residual = response - t
                need(0 <= residual <= mass, 'residual response digit bound')
                for k in range(4):
                    inc = k * (ell - 1) // 2 + ell * comb(k, 2) + k * response
                    normalized = (ell - 1) * k * k // 2 + k * residual
                    need(inc - comb(t + k, 2) + comb(t, 2) == normalized,
                         'exact normalized insertion identity')
                    new_mass = mass + ell * k
                    need(factorial(new_mass) % (factorial(mass) * ell ** k * factorial(k)) == 0,
                         'integer local normalized transfer')
    for k in range(16):
        need((3 * k) ** 2 // 9 == k * k, 'equal three-cycle histories attain normalized bound')
    # The final Gaussian functional is linear, not evaluation or a ring map.
    need((1 << comb(2, 2)) != (1 << comb(1, 2)) ** 2, 'Gaussian functional is not multiplicative')
    a, b = (9,), (3, 3, 3)
    normal = [unmerged.cycle_weight(x) / (1 << comb(len(x), 2)) for x in (a, b)]
    need(sum(normal) != sum(unmerged.cycle_weight(x) for x in (a, b)),
         'dropping retained cycle polynomial corrupts equal-degree unit-tail state')
    h1, h2 = (5, 13), (9, 9)
    need(sum(h1) == sum(h2) and len(h1) == len(h2), 'same mass and cycle count hostile')
    need(sum(gcd(a, 3) - 1 for a in h1) == 0 and
         sum(gcd(a, 3) - 1 for a in h2) == 4,
         'cycle polynomial does not replace residual gcd response')
    # At k=1,s=1 a core excess shifts degree twice, but only one actual cycle.
    need(-comb(1, 0) * (1 << comb(1, 2)) == -1 and
         -(1 << comb(2, 2)) == -2, 'core excess is not an actual cycle')
    need(list_poly_count(0)[0] == [1] and hybrid_dp(0)[0] == [1], 'empty universe and switch boundary')
    observer, _ = list_poly_count(3, no_fixed=True)
    need(observer[3] == Fraction(2, 3), 'no-fixed observer remains rational')
    for kernel in (module_dp, list_poly_count, hybrid_dp, count, poly_count):
        for bad in (-1, True, 2.0, '2'):
            try:
                kernel(bad)
            except (ValueError, TypeError):
                need(True, 'literal input rejection')
            else:
                need(False, 'literal input rejection')
        for bad in (0, 1, 'false'):
            try:
                kernel(4, core=bad)
            except (ValueError, TypeError):
                need(True, 'literal mode rejection')
            else:
                need(False, 'literal mode rejection')
    return unmerged


def run():
    CHECKS.clear()
    unmerged = controls()
    previous = json.loads((PRIOR / 'receipt.json').read_text())
    series = {}
    cases = [('tournaments', 200, False, False), ('pair_free_cores', 100, True, False),
             ('no_fixed_core_mass', 100, True, True), ('no_fixed_mass', 100, False, True)]
    local = load(HERE / 'incoming_local.py', 'incoming_large_control')
    for name, bound, core, no_fixed in cases:
        answer, stats = list_poly_count(bound, core=core, no_fixed=no_fixed)
        second, second_stats = local.local_dp(bound, core=core, no_fixed=no_fixed)
        for x, y in zip(answer, second):
            need(x == y, 'large independent incoming scalar kernel coefficient')
        if name in previous['series']:
            need([str(x) for x in answer] == previous['series'][name]['coefficients'],
                 'previous frozen exact coefficient vector')
        series[name] = {'bound': bound, 'core': core, 'no_fixed': no_fixed,
                        'coefficients': [str(x) for x in answer], 'statistics': stats,
                        'incoming_statistics': second_stats}
        if name == 'tournaments':
            for line in (Path(unmerged.__file__).parent / 'a000568-reference.txt').read_text().splitlines():
                if line and not line.startswith('#'):
                    n, value = map(int, line.split())
                    need(answer[n] == value, 'independent stored OEIS coefficients through 76')
            need(answer[:61] == unmerged.direct_series(60), 'independent unmerged partitions through 60')
    return {'schema': 1, 'method': 'Gaussian-normalized cycle polynomial and terminal functional',
            'checks': sum(CHECKS.values()), 'checks_by_kind': dict(sorted(CHECKS.items())),
            'dependencies': {'previous_script_lf_sha256': digest(PRIOR / 'check.py'),
                             'previous_receipt_lf_sha256': digest(PRIOR / 'receipt.json'),
                             'incoming_local_lf_sha256': digest(HERE / 'incoming_local.py'),
                             'incoming_commit': '622e445138',
                             'incoming_source_path': '05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py'},
            'series': series}


def benchmarks():
    old = inherited()
    local = load(HERE / 'incoming_local.py', 'incoming_benchmark')
    functions = [('previous_common_odd', old.packed_dp), ('incoming_degree_local', local.local_dp),
                 ('global_polynomial', list_poly_count), ('local_polynomial', module_dp),
                 ('scalar_prefix_unit_tail', lambda n, **kw: hybrid_dp(n, switch=1, **kw))]
    records = []
    for bound, core in ((200, False), (100, True)):
        baseline = None
        for repeat in range(3):
            for name, fn in (functions if repeat % 2 == 0 else functions[::-1]):
                stamp = time.perf_counter()
                values, stats = fn(bound, core=core)
                elapsed = time.perf_counter() - stamp
                if baseline is None:
                    baseline = values
                if values != baseline:
                    raise ValueError('timed coefficient mismatch')
                records.append({'bound': bound, 'core': core, 'repeat': repeat,
                                'algorithm': name, 'seconds': elapsed,
                                'transitions': stats['transitions'],
                                'peak_profile_keys': stats.get('peak_states'),
                                'coefficient_sha256': sha256('\n'.join(map(str, values)).encode()).hexdigest()})
                print(json.dumps(records[-1]), flush=True)
    memories = []
    for bound, core in ((100, False), (60, True)):
        baseline = None
        for name, fn in functions[:4]:
            tracemalloc.start()
            values, stats = fn(bound, core=core)
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            if baseline is None:
                baseline = values
            if values != baseline:
                raise ValueError('memory coefficient mismatch')
            memories.append({'bound': bound, 'core': core, 'algorithm': name, 'traced_peak_bytes': peak})
    medians = []
    for bound, core in ((200, False), (100, True)):
        for name, _ in functions:
            samples = [r['seconds'] for r in records if r['bound'] == bound and r['core'] == core and r['algorithm'] == name]
            medians.append({'bound': bound, 'core': core, 'algorithm': name, 'median_seconds': statistics.median(samples)})
    return {'python': sys.version, 'shared_host': True, 'alternating_order': True,
            'setup_and_output_conversion_included': True, 'runs': records, 'medians': medians,
            'memory_runs_separate_from_timings': memories,
            'source_hashes': {'check.py': digest(__file__), 'incoming_local.py': digest(HERE / 'incoming_local.py'),
                              'previous_check.py': digest(PRIOR / 'check.py')}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--benchmark', action='store_true')
    parser.add_argument('--probe', type=int)
    parser.add_argument('--core', action='store_true')
    args = parser.parse_args()
    if args.probe is not None:
        result = probe(args.probe, args.core)
        for row in result['runs']:
            row['statistics'].pop('rows', None)
        print(json.dumps(result, indent=2))
        return
    result = benchmarks() if args.benchmark else run()
    target = HERE / ('benchmarks.json' if args.benchmark else 'receipt.json')
    content = json.dumps(result, sort_keys=True, indent=2) + '\n'
    if args.verify:
        if target.read_text(encoding='utf-8') != content:
            raise ValueError('saved receipt differs')
    else:
        target.write_text(content, encoding='utf-8', newline='\n')
    print(json.dumps({'checks': result.get('checks'), 'path': target.name,
                      'lf_sha256': digest(target), 'verified': args.verify}, sort_keys=True))


if __name__ == '__main__':
    main()
