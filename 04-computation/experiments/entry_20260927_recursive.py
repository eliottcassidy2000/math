"""Source-relative guarded Collatz certificate calculus and exact controls."""
from dataclasses import dataclass
from hashlib import sha256
from fractions import Fraction
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def v2(n):
    need(n > 0, 'positive valuation argument')
    return (n & -n).bit_length() - 1


def U(n):
    return (3 * n + 1) >> v2(3 * n + 1)


@dataclass(frozen=True)
class Transfer:
    multiplier: int = 1
    denominator: int = 1
    carry: int = 0
    steps: int = 0

    def then(self, following):
        return Transfer(following.multiplier * self.multiplier,
                        following.denominator * self.denominator,
                        following.multiplier * self.carry + self.denominator * following.carry,
                        self.steps + following.steps)

    def apply(self, n):
        numerator = self.multiplier * n + self.carry
        need(numerator % self.denominator == 0, 'nonintegral transfer')
        return numerator // self.denominator


def exact_word(n, count):
    result = Transfer()
    for _ in range(count):
        a = v2(3 * n + 1)
        step = Transfer(3, 1 << a, 1, 1)
        result = result.then(step)
        n = step.apply(n)
    return n, result


def bank():
    result = {}
    for q in range(1, 342, 2):
        c = 3 * q
        J = A = 0
        while c >= 2 * q:
            need(J < 1000, 'finite core cap')
            P = A
            A += v2(3 * c + 1)
            J += 1
            c = U(c)
        for R in range(P + 1, A + 1):
            if (1 << (R + 1)) - 3 ** (J + 1) > max(0, 3 * (1 << (A - R)) * c - 6 * q + 1):
                break
        else:
            raise ValueError('missing finite-bank threshold')
        K = R + 1
        residue = ((6 * q - 1) * pow(3, -1, 1 << K)) % (1 << K)
        for small in range(residue, 2 * q + 1, 1 << K):
            endpoint, _ = exact_word(small, J + 1)
            need(1 < small and endpoint < small, 'small-source bank exception')
        result[q] = dict(J=J, A=A, P=P, b=c, R=R, K=K, residue=residue)
    return result


def primitive(n, node, entries):
    kind, argument = node
    need(n > 0 and n & 1, 'positive odd source')
    if kind == 'step':
        return exact_word(n, 1)
    if kind == 'ones':
        length = argument
        need(length >= 1 and v2(n + 1) >= length + 1, 'ones guard')
        p, d = 3 ** length, 1 << length
        t = Transfer(p, d, p - d, length)
    elif kind == 'pairs':
        length = argument
        need(length >= 1 and v2(n + 5) >= 3 * length + 1, 'pair guard')
        p, d = 9 ** length, 8 ** length
        t = Transfer(p, d, 5 * (p - d), 2 * length)
    elif kind == 'to_one':
        numerator = 3 * n + 1
        need(numerator & (numerator - 1) == 0, 'terminal power-of-two guard')
        t = Transfer(3, numerator, 1, 1)
    elif kind == 'to_fixed':
        target = argument
        need(target > 0 and target & 1, 'positive odd terminal target')
        numerator = 3 * n + 1
        need(numerator % target == 0, 'terminal target divisibility')
        quotient = numerator // target
        need(quotient > 1 and quotient & (quotient - 1) == 0,
             'terminal quotient power-of-two guard')
        t = Transfer(3, quotient, 1, 1)
    elif kind == 'core':
        row = entries[argument]
        need(n % (1 << row['K']) == row['residue'], 'core-call source class')
        endpoint, t = exact_word(n, row['J'] + 1)
        need(endpoint < n, 'core-call own-source conclusion')
        return endpoint, t
    else:
        raise ValueError('unknown certificate rule')
    return t.apply(n), t


def verify(source, nodes, entries, require_descent=True):
    """A child proof may discharge its own target; caller target stays source."""
    n = source
    total = Transfer()
    for node in nodes:
        if node[0] == 'call':
            child_source, body = node[1]
            need(child_source == n, 'call substituted a different source')
            endpoint, transfer = verify(child_source, body, entries, True)
        else:
            endpoint, transfer = primitive(n, node, entries)
        total = total.then(transfer)
        need(total.apply(source) == endpoint, 'affine composition lost original source')
        n = endpoint
    if require_descent:
        need(n < source, 'enclosing target not discharged')
        need((total.denominator - total.multiplier) * source > total.carry,
             'terminal affine comparison')
    return n, total


def policy_certificate(source, entries, cap=20000):
    n = source
    nodes = []
    for _ in range(cap):
        if n < source:
            endpoint, transfer = verify(source, nodes, entries)
            return endpoint, transfer.steps, len(nodes), nodes
        if (3 * n + 1) & (3 * n) == 0:
            node = ('to_one', None)
        else:
            matches = [(row['J'], q) for q, row in entries.items()
                       if n % (1 << row['K']) == row['residue']]
            if matches:
                node = ('core', min(matches)[1])
            elif n % 4 == 1:
                node = ('step', None)
            elif v2(n + 5) >= 4:
                node = ('pairs', (v2(n + 5) - 1) // 3)
            else:
                node = ('ones', v2(n + 1) - 1)
        n, _ = primitive(n, node, entries)
        nodes.append(node)
    raise ValueError('finite policy cap exhausted; no completeness conclusion')


def completed_exponent(k, target=1):
    need(target in (1, 47), 'supported independently checked completion targets')
    h = 0
    for level in range(1, 2 * k + 1):
        candidates = [h + digit * 3 ** (level - 1) for digit in range(3)]
        modulus = 3 ** (level + 1)
        good = [t for t in candidates if
                ((pow(4, t, modulus) + 14) if target == 1 else
                 (47 * pow(4, t, modulus) + 7)) % modulus == 0]
        need(len(good) == 1, 'ternary exponent lift')
        h = good[0]
    return h


def return_cylinder(k):
    """All-height sufficient final valuation, preserving the original source."""
    need(k >= 0, 'nonnegative pair count')
    m = k + 1
    t = 1
    while (1 << t) * (8 ** m - 5) <= 9 ** m - 5:
        t += 1
    odd_residue = (5 * pow(9 ** m, -1, 1 << t)) % (1 << t)
    power = 3 * m + t
    residue = (odd_residue << (3 * m)) - 5
    return dict(k=k, t=t, b=odd_residue, K=power, residue=residue)


def disjoint_cylinders(rows):
    kept = []
    for residue, power in sorted(set(rows), key=lambda row: (row[1], row[0])):
        if not any(power >= old_power and residue % (1 << old_power) == old_residue
                   for old_residue, old_power in kept):
            kept.append((residue, power))
    return kept


def density_extension(entries, count):
    old = disjoint_cylinders((row['residue'], row['K']) for row in entries.values())
    additions = [return_cylinder(k) for k in range(count)]
    old_mass = sum((Fraction(1, 1 << power) for _, power in old), Fraction())
    # Every new k>=3 cylinder lies inside n=-5 mod4096. The old bank
    # misses this entire containing cylinder, so the comparison is all-k.
    for old_residue, old_power in old:
        intersects = ((4091 % (1 << old_power) == old_residue) if old_power <= 12
                      else (old_residue % 4096 == 4091))
        need(not intersects, 'old bank must miss entire -5 mod4096 cylinder')
    new_mass = Fraction()
    overlap = Fraction()
    wholly_old = []
    for row in additions:
        residue, power = row['residue'], row['K']
        mass = Fraction(1, 1 << power)
        new_mass += mass
        intersect = Fraction()
        for old_residue, old_power in old:
            if power >= old_power and residue % (1 << old_power) == old_residue:
                intersect += mass
            elif old_power >= power and old_residue % (1 << power) == residue:
                intersect += Fraction(1, 1 << old_power)
        need(intersect <= mass, 'disjoint old-bank overlap')
        overlap += intersect
        if intersect == mass:
            wholly_old.append(row['k'])
    augmented = disjoint_cylinders(old + [(r['residue'], r['K']) for r in additions])
    augmented_mass = sum((Fraction(1, 1 << power) for _, power in augmented), Fraction())
    need(augmented_mass == old_mass + new_mass - overlap, 'independent union-density calculation')
    return dict(count=count, old_disjoint_rows=len(old), old_density=str(old_mass),
                new_family_density=str(new_mass), overlap_density=str(overlap),
                added_density=str(new_mass - overlap), augmented_density=str(augmented_mass),
                augmented_density_decimal=float(augmented_mass), wholly_old_k=wholly_old,
                unlisted_family_tail_density_less_than=str(Fraction(1, 8 * 9 ** count)))


def main():
    entries = bank()
    accepted = rejected = 0
    for n in range(3, 2048, 2):
        for kind in ('ones', 'pairs'):
            for length in range(1, 7):
                legal = v2(n + (1 if kind == 'ones' else 5)) >= (length + 1 if kind == 'ones' else 3 * length + 1)
                try:
                    endpoint, transfer = primitive(n, (kind, length), entries)
                except ValueError:
                    need(not legal, 'rejected a valid loop')
                    rejected += 1
                else:
                    need(legal, 'accepted an invalid loop')
                    literal, literal_transfer = exact_word(n, transfer.steps)
                    need(endpoint == literal and transfer == literal_transfer, 'loop versus literal exact word')
                    accepted += 1

    # A valid child descent does not discharge its enclosing original target.
    try:
        verify(27, [('step', None), ('call', (41, [('step', None)]))], entries)
    except ValueError as error:
        need(str(error) == 'enclosing target not discharged', 'expected root-target rejection')
    else:
        raise ValueError('unsound nested descent accepted')
    try:
        verify(27, [('call', (9, [('step', None)]))], entries)
    except ValueError as error:
        need(str(error) == 'call substituted a different source', 'expected source rejection')
    else:
        raise ValueError('unjustified smaller-core substitution accepted')
    # Genuine nested proof: 7->11, child11->17->13->5, with5 below both roots.
    end, _ = verify(7, [('step', None), ('call', (11, [('step', None)] * 3))], entries)
    need(end == 5, 'valid nested source-relative proof')

    nested_rows = []
    for k in range(1, 129):
        start = 4 * 8 ** k - 5
        length = 4 + v2(k)
        endpoint, transfer = verify(start, [('pairs', k), ('ones', length)], entries, False)
        expected = 3 ** length * (9 ** k - 1) // (1 << (length - 2)) - 1
        need(endpoint == expected and endpoint > start, 'nested shadow endpoint')
        literal = start
        for _ in range(2 * k + length):
            literal = U(literal)
            need(literal > start, 'nested first-descent lower bound')
        need(literal == endpoint and transfer.steps == 2 * k + length, 'nested clock')
        if k & 1:
            need(endpoint == (81 * 9 ** k - 85) // 4, 'odd-k closed endpoint')
            a = v2(3 * endpoint + 1)
            need(a == v2(243 * 9 ** k - 251) - 2, 'next precision gate')
        if k <= 8:
            nested_rows.append((k, start, 2 * k + length, endpoint, v2(3 * endpoint + 1)))

    completed = []
    for k in range(1, 5):
        h = completed_exponent(k)
        need((4 ** h + 14) % 3 ** (2 * k + 1) == 0, 'completed family integrality')
        coefficient = (4 ** h + 14) // 3 ** (2 * k + 1)
        need(coefficient >= 2 and coefficient % 2 == 0, 'completed family positive even coefficient')
        source = coefficient * 8 ** k - 5
        endpoint, transfer = verify(source, [('pairs', k), ('to_one', None)], entries)
        need(endpoint == 1 and transfer.steps == 2 * k + 1, 'two-node terminal certificate')
        literal = source
        for _ in range(2 * k):
            literal = U(literal)
            need(literal > source, 'completed family exact first descent')
        need(literal == (4 ** h - 1) // 3, 'completed family terminal inverse')
        completed.append((k, h, source.bit_length(), transfer.steps, 2))

    # The family now contains27. Its two-node certificate reaches47, which
    # discharges every other member but emphatically does not discharge27.
    completed47 = []
    for k, shift in [(1, 0), (1, 1), (2, 0), (3, 0), (4, 0)]:
        h = completed_exponent(k, 47) + shift * 9 ** k
        numerator = 2 * (47 * 4 ** h + 7)
        divisor = 3 ** (2 * k + 1)
        need(numerator % divisor == 0, 'target47 integrality')
        coefficient = numerator // divisor
        source = coefficient * 8 ** k - 5
        nodes = [('pairs', k), ('to_fixed', 47)]
        if source == 27:
            try:
                verify(source, nodes, entries)
            except ValueError as error:
                need(str(error) == 'enclosing target not discharged', 'target47 source27 rejection')
            else:
                raise ValueError('arrival at47 incorrectly paid original27')
            nodes += [('call', (47, [('step', None)] * 34))]
            endpoint, transfer = verify(source, nodes, entries)
            need(endpoint == 23 and transfer.steps == 37, 'fixed27 completed certificate')
        else:
            endpoint, transfer = verify(source, nodes, entries)
            need(endpoint == 47 and transfer.steps == 2 * k + 1, 'target47 two-node certificate')
            need(v2(source + 5) == 3 * k + 1, 'target47 direct counter extraction')
        literal = source
        for _ in range(transfer.steps - 1):
            literal = U(literal)
            need(literal > source, 'target47 exact first-descent lower bound')
        need(U(literal) == endpoint, 'target47 literal endpoint')
        completed47.append(dict(k=k, h=h, source_bits=source.bit_length(),
                                steps=transfer.steps, top_level_nodes=len(nodes),
                                total_stored_nodes=37 if source == 27 else 2))

    return_controls = 0
    for k in range(100):
        row = return_cylinder(k)
        need(row['b'] & 1, 'odd return coefficient residue')
        for shift in range(8):
            b = row['b'] + (shift << row['t'])
            source = b * 8 ** (k + 1) - 5
            need(source % (1 << row['K']) == row['residue'], 'return source cylinder')
            nodes = ([('pairs', k)] if k else []) + [('step', None)] * 2
            endpoint, transfer = verify(source, nodes, entries)
            need(transfer.steps == 2 * k + 2, 'return-cylinder clock')
            literal = source
            for _ in range(2 * k + 1):
                literal = U(literal)
                need(literal > source, 'return-cylinder exact first descent')
            need(U(literal) == endpoint, 'return-cylinder literal endpoint')
            numerator = b * 9 ** (k + 1) - 5
            need(endpoint == numerator >> v2(numerator), 'return-cylinder terminal closed form')
            return_controls += 1
    return_hostiles = []
    for k, b, previous_t in [(5, 3, 1), (11, 1, 2)]:
        source = b * 8 ** (k + 1) - 5
        numerator = b * 9 ** (k + 1) - 5
        need(v2(numerator) == previous_t, 'exact lower-valuation hostile')
        literal = source
        for _ in range(2 * k + 2):
            literal = U(literal)
            need(literal > source, 'return budget boundary hostile')
        return_hostiles.append(dict(k=k, b=b, t=previous_t, source=source,
                                    endpoint=literal, steps=2 * k + 2))

    # This tests a particular next-step payment gate, not all possible returns.
    successful_gates = []
    valuation_rows = []
    for k in range(1, 2000, 2):
        numerator = 243 * 9 ** k - 251
        a = v2(numerator) - 2
        original = 4 * 8 ** k - 5
        after_gate = numerator >> (a + 2)
        valuation_rows.append((k, a))
        if k % 4 == 1:
            need(a == 2 and after_gate > original, 'quarter-class failed payment gate')
        if after_gate < original:
            need((1 << a) * 16 * 8 ** k > 243 * 9 ** k, 'necessary exponential payment budget')
            successful_gates.append((k, a, after_gate))
    for i, (k, a) in enumerate(valuation_rows[:128]):
        for ell, other_a in valuation_rows[i + 1:128]:
            need(v2(ell - k) >= min(a, other_a) - 1, 'valuation collision spacing')
    for (k, a, _), (ell, other_a, _) in zip(successful_gates, successful_gates[1:]):
        need(32 * (ell - k) * 8 ** k > 243 * 9 ** k, 'successful gate exponential spacing')

    ordinary = []
    for n in range(3, 8192, 2):
        endpoint, steps, count, _ = policy_certificate(n, entries)
        ordinary.append((n, endpoint, steps, count))
    family = []
    mersenne = []
    for k in range(1, 129):
        n = 4 * 8 ** k - 5
        endpoint, steps, count, _ = policy_certificate(n, entries)
        family.append((k, endpoint.bit_length(), steps, count))
        n = (1 << (2 * k)) - 1
        endpoint, steps, count, _ = policy_certificate(n, entries)
        mersenne.append((2 * k, endpoint.bit_length(), steps, count))
    summary = dict(
        status='Exact source-relative certificate controls; policy completeness OPEN',
        finite_core_rows=len(entries),
        accepted_loop_controls=accepted,
        rejected_illegal_loop_controls=rejected,
        negative_controls=['nested target27 remains unpaid at31', 'core9 cannot replace actual27'],
        nested_shadow_controls=128,
        nested_shadow_examples=nested_rows,
        two_node_completion_examples=completed,
        target47_completion_examples=completed47,
        all_height_return_cylinder_controls=return_controls,
        return_cylinders_first_twenty=[return_cylinder(k) for k in range(20)],
        return_budget_hostiles=return_hostiles,
        return_cylinder_density_extensions=[density_extension(entries, count) for count in (20, 100)],
        next_gate_odd_exponents_tested=1000,
        next_gate_successes=successful_gates,
        next_gate_pairwise_congruence_controls=128 * 127 // 2,
        ordinary_source_count=len(ordinary),
        maximum_ordinary_macro_nodes=max(r[3] for r in ordinary),
        maximum_ordinary_atomic_steps=max(r[2] for r in ordinary),
        fixed27_family_controls=len(family),
        maximum_fixed27_macro_nodes=max(r[3] for r in family),
        maximum_fixed27_atomic_steps=max(r[2] for r in family),
        mersenne_controls=len(mersenne),
        maximum_mersenne_macro_nodes=max(r[3] for r in mersenne),
        maximum_mersenne_atomic_steps=max(r[2] for r in mersenne),
        finite_policy_digest=sha256(json.dumps([ordinary, family, mersenne], separators=(',', ':')).encode()).hexdigest(),
        family_first_eight=family[:8],
    )
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
