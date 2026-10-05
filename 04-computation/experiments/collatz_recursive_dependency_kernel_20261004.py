"""Guarded repetition of a Collatz common-future dependency, not a new orbit map.

Run with python -B or python -B -O.  No import-time experiment or file mutation.
The companion note proves the all-height statements and scopes compression.
"""
from dataclasses import dataclass, replace
import inverse_ray_ternary_addresses_20261004 as codec

V = (1, 2, 1, 1, 1, 2)
W = (3, 2, 1, 1, 1, 2)
H_DATA = (729, 1024, 669)
W_DATA = (729, 1024, 2971)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def integer(value, lower, message):
    need(type(value) is int and value >= lower, message)


def odd(value):
    integer(value, 1, 'positive integer source')
    need(value % 2 == 1, 'odd source')


def v2(value):
    need(type(value) is int and value != 0, 'nonzero integer valuation argument')
    value = abs(value)
    return (value & -value).bit_length()-1


def step(source):
    odd(source)
    value = 3*source+1
    exponent = v2(value)
    return value >> exponent, exponent


def word_guard(word):
    need(type(word) is tuple and all(type(a) is int and a >= 1 for a in word),
         'tuple of positive integer valuation exponents')


def compose(first, second):
    """Chronological composition: second after first; tuples are P,Q,B."""
    p, q, b = first
    r, s, c = second
    return r*p, s*q, r*b+c*q


def stats(word):
    word_guard(word)
    result = (1, 1, 0)
    for a in word:
        result = compose(result, (3, 1 << a, 1))
    return result


def power(data, count):
    integer(count, 0, 'nonnegative repeat count')
    result, operations = (1, 1, 0), 0
    while count:
        if count & 1:
            result = compose(result, data)
            operations += 1
        count >>= 1
        if count:
            data = compose(data, data)
            operations += 1
    return result, operations


def literal(source, word, require_root=False):
    """Independent actual-step check, rejecting every root-padded suffix."""
    odd(source)
    word_guard(word)
    current = source
    for expected in word:
        need(current != 1, 'no root self-return or padded first-hit suffix')
        current, actual = step(current)
        need(actual == expected, 'exact valuation word')
    if require_root:
        need(current == 1, 'terminal ROOT witness')
    return current


def observed_word(source, cap=10000):
    """Finite test helper only: explicitly observes a terminal proof to supply."""
    odd(source)
    integer(cap, 0, 'nonnegative observation cap')
    current, word = source, []
    while current != 1:
        need(len(word) < cap, 'terminal observation cap')
        current, a = step(current)
        word.append(a)
    return tuple(word)


def fuel(source):
    odd(source)
    return v2(295*source-669)


def repeat_endpoint(source, repetitions):
    """Prefix-only theorem: this is NOT by itself a home certificate."""
    odd(source)
    integer(repetitions, 1, 'positive integer repeat count')
    need(fuel(source) >= 10*repetitions+1, 'retained exact binary fuel')
    delta = (729**repetitions*(295*source-669)) >> (10*repetitions)
    need((delta+669) % 295 == 0, 'closed-form endpoint integrality')
    target = (delta+669)//295
    need(1 < target < source and target % 2 == 1, 'strict positive odd dependency')
    return target


@dataclass(frozen=True)
class RepeatCertificate:
    source: int
    repetitions: int
    terminal: int
    terminal_word: tuple


def verify(record):
    need(type(record) is RepeatCertificate, 'typed repeat certificate')
    odd(record.terminal)
    need(record.terminal > 1, 'nonroot terminal dependency')
    target = repeat_endpoint(record.source, record.repetitions)
    need(target == record.terminal, 'supplied terminal is this source dependency')
    need(record.terminal_word != (), 'a nonroot terminal needs a first edge')
    literal(record.terminal, record.terminal_word, require_root=True)
    return True


def expand_word(record):
    verify(record)
    first, *tail = record.terminal_word
    return V+W*(record.repetitions-1)+(first+2,)+tuple(tail)


def codec_from_word(word):
    """Construct the old ROOT grammar by modular inverse checks, no orbit search."""
    word_guard(word)
    cert = codec.ROOT
    for exponent in reversed(word):
        residue = (pow(2, exponent, 9)*codec.mod3(cert, 2)-1) % 9
        need(residue % 3 == 0, 'inverse word has an integer odd source')
        row = residue//3
        least = codec.kappa(codec.mod3(cert, 2), row)
        need(exponent >= least and (exponent-least) % 6 == 0,
             'inverse word is an exact valuation edge')
        cert = codec.extend(cert, row, (exponent-least)//6)
    return cert


def expand_certificate(record):
    word = expand_word(record)
    literal(record.source, word, require_root=True)
    cert = codec_from_word(word)
    need(codec.expand(cert, bit_cap=codec.bit_bounds(cert)[1]) == record.source,
         'inherited codec retains supplied source identity')
    need(codec.ranks(cert) == (len(word), sum(word)+len(word)), 'exact first-hit ranks')
    return cert


def principal_log4(target, depth):
    integer(depth, 1, 'positive ternary precision')
    integer(target, 0, 'nonnegative principal-unit representative')
    need(target % 3 == 1, 'principal unit modulo three')
    result, period = 0, 1
    for precision in range(2, depth+1):
        modulus = 3**precision
        choices = [result+d*period for d in range(3)
                   if pow(4, result+d*period, modulus) == target % modulus]
        need(len(choices) == 1, 'one lifted ternary exponent digit')
        result = choices[0]
        period *= 3
    return result, period


@dataclass(frozen=True)
class CompletedPlan:
    repetitions: int
    root_exponent: int


def audit_plan(plan):
    need(type(plan) is CompletedPlan, 'typed symbolic completed family')
    integer(plan.repetitions, 1, 'positive integer repeat count')
    integer(plan.root_exponent, 4, 'positive-height terminal root exponent')
    need(plan.root_exponent % 2 == 0, 'even terminal root exponent')
    modulus = 3**(6*plan.repetitions+1)
    need((295*pow(2, plan.root_exponent, modulus)-2302) % modulus == 0,
         'full terminal ternary divisibility')
    return True


def completed_plan(repetitions, lift=0):
    integer(repetitions, 1, 'positive integer repeat count')
    integer(lift, 0, 'nonnegative exponent-family parameter')
    depth = 6*repetitions+1
    modulus = 3**depth
    first, period = principal_log4(2302*pow(295, -1, modulus) % modulus, depth)
    if first < 2:
        first += period
    plan = CompletedPlan(repetitions, 2*(first+lift*period))
    audit_plan(plan)
    return plan


def plan_word(plan):
    audit_plan(plan)
    return V+W*(plan.repetitions-1)+(plan.root_exponent+2,)


def plan_source_mod(plan, modulus):
    """Closed formula, exact even for noncoprime moduli, without source expansion."""
    audit_plan(plan)
    integer(modulus, 1, 'positive residue modulus')
    p = 729**plan.repetitions
    q = 1024**plan.repetitions
    denominator = 885*p
    enlarged = denominator*modulus
    numerator = (2007*p+q*(295*pow(2, plan.root_exponent, enlarged)-2302)) % enlarged
    need(numerator % denominator == 0, 'modular closed-form integer division')
    return numerator//denominator


def inverse_word_mod(word, modulus):
    """Independent modular reader: reverse literal inverse letters one by one."""
    word_guard(word)
    integer(modulus, 1, 'positive residue modulus')
    current_modulus = modulus*3**len(word)
    current = 1
    for a in reversed(word):
        value = (pow(2, a, current_modulus)*current-1) % current_modulus
        need(value % 3 == 0, 'independent modular inverse integrality')
        current = value//3
        current_modulus //= 3
    return current


def materialize(plan, bit_cap=1000000):
    audit_plan(plan)
    integer(bit_cap, 1, 'positive explicit expansion bit cap')
    need(plan.root_exponent+10*plan.repetitions+2 <= bit_cap,
         'declared source-expression expansion cap')
    terminal = ((1 << plan.root_exponent)-1)//3
    p, q = 729**plan.repetitions, 1024**plan.repetitions
    delta = 295*terminal-669
    need(delta % p == 0, 'terminal ternary fuel for inverse dependency')
    numerator = 669+q*(delta//p)
    need(numerator % 295 == 0, 'integer source from fixed-point expression')
    record = RepeatCertificate(numerator//295, plan.repetitions, terminal,
                               (plan.root_exponent,))
    verify(record)
    need(fuel(terminal) == 1 and fuel(record.source) == 10*plan.repetitions+1,
         'exactly the advertised repeats before the explicit root exit')
    return record


def power3_exponent(repetitions):
    """Unique power-of-three exponent class paying the repeat guard, not its exit."""
    integer(repetitions, 1, 'positive integer repeat count')
    bits = 10*repetitions+1
    modulus = 1 << bits
    target = 669*pow(295, -1, modulus) % modulus
    need(target % 8 == 3, 'target belongs to the powers-of-three subgroup')
    exponent, period = 1, 2
    for precision in range(4, bits+1):
        small_modulus = 1 << precision
        choices = [exponent+d*period for d in range(2)
                   if pow(3, exponent+d*period, small_modulus) == target % small_modulus]
        need(len(choices) == 1, 'unique binary exponent lift')
        exponent = choices[0]
        period *= 2
    need(period == 1 << (10*repetitions-1) and 0 < exponent < period,
         'canonical exponent and exact period')
    return exponent, period


def decode_word(data):
    """Positive-word carrier test; diagonal powers alone are insufficient."""
    need(type(data) is tuple and len(data) == 3
         and all(type(value) is int and value > 0 for value in data),
         'exact positive integer carrier tuple')
    p, q, carry = data
    length = 0
    while p > 1 and p % 3 == 0:
        p //= 3
        length += 1
    need(p == 1 and length >= 1 and q > 1 and q & (q-1) == 0, 'typed diagonal')
    cost, output = q.bit_length()-1, []
    for remaining in range(length, 1, -1):
        residual = carry-3**(remaining-1)
        need(residual > 0, 'positive tail carry')
        a = v2(residual)
        need(1 <= a < cost, 'positive first exponent and remaining cost')
        output.append(a)
        carry, cost = residual >> a, cost-a
    need(carry == 1 and cost >= 1, 'one-letter terminal carry')
    output.append(cost)
    need(stats(tuple(output)) == data, 'full ordered carrier equality')
    return tuple(output)


def rejected(function):
    try:
        function()
    except ValueError:
        return
    raise ValueError('hostile input was accepted')


def main():
    need(stats(V) == (729, 256, 925) and stats(W) == W_DATA, 'inherited words')
    need(compose(stats(V), W_DATA) == compose(H_DATA, stats(V)), 'conjugacy square')
    need(decode_word(W_DATA) == W, 'actual repeated word decoder')
    rejected(lambda: decode_word(H_DATA))
    need(repeat_endpoint(155, 1) == 111 and literal(155, V) == 445,
         'same diagonal does not make H an actual Collatz word')
    prefix_cases, summary_products = 0, 0
    for m in range(1, 33):
        modulus = 1 << (10*m+1)
        residue = 669*pow(295, -1, modulus) % modulus
        summary, operations = power(W_DATA, m-1)
        need(summary == stats(W*(m-1)), 'binary/linear composition agreement')
        summary_products += operations
        for lift in (0, 1, 7, 31):
            source = residue+lift*modulus
            target = repeat_endpoint(source, m)
            need(target*1024**m > source*729**m, 'strict fixed-point height bound')
            need(target.bit_length() <= source.bit_length()
                 <= target.bit_length()+(m+1)//2, 'retained integer bit-size bound')
            current = source
            for _ in range(m):
                need(current % 2048 == 155, 'every original-source guard')
                child = (729*current+669)//1024
                need(111 <= child < current, 'every decreasing positive dependency')
                need(fuel(child) == fuel(current)-10, 'exact binary fuel consumption')
                current = child
            need(current == target, 'linear/closed dependency agreement')
            endpoint, b = step(target)
            need(literal(source, V+W*(m-1)+(b+2,)) == endpoint,
                 'actual repeated prefix and terminal common future')
            bad = residue+(1 << (10*m))
            need(fuel(bad) == 10*m, 'sharp missing final oddness bit')
            rejected(lambda bad=bad, m=m: repeat_endpoint(bad, m))
            prefix_cases += 1
    certificates = []
    observations = 0
    for m in range(1, 13):
        modulus = 1 << (10*m+1)
        source = 669*pow(295, -1, modulus) % modulus
        terminal = repeat_endpoint(source, m)
        terminal_word = observed_word(terminal)
        observations += len(terminal_word)
        record = RepeatCertificate(source, m, terminal, terminal_word)
        cert = expand_certificate(record)
        need(codec.ranks(cert)[0] == len(terminal_word)+6*m, 'odd-rank substitution')
        need(codec.ranks(cert)[1] == sum(terminal_word)+len(terminal_word)+16*m,
             'ordinary-rank substitution')
        certificates.append(record)
    for m in (64, 128, 256, 512):
        modulus = 1 << (10*m+1)
        source = 669*pow(295, -1, modulus) % modulus
        target = repeat_endpoint(source, m)
        need(fuel(source)-fuel(target) == 10*m, 'large prefix-only compressed guard')
        _, operations = power(W_DATA, m-1)
        need(operations <= 2*(m-1).bit_length(), 'logarithmic composition-node count')
    print('PROVED-IMPLEMENTATION: H(n)=(729n+669)/1024 on n=155 mod2048.')
    print('Finite prefix universe: m=1..32; four lifts each; checked', prefix_cases,
          'prefixes, plus m=64,128,256,512 guard-only controls.')
    print('Closed/linear block products agree; binary composition products across m=1..32:',
          summary_products)
    print('Completed supplied-terminal controls:', len(certificates),
          '; separately observed terminal odd edges:', observations)
    print('Their added 468 odd letters use 12 repeat counts and shared V/W templates;')
    print('the 1622 supplied terminal letters and all integer/source bits remain explicit.')
    print('Each repeat adds odd rank6, halving cost10, ordinary rank16.')
    for m in (1, 2):
        plan = completed_plan(m)
        record = materialize(plan)
        cert = expand_certificate(record)
        print('All-height completed family literal control:',
              'm', m, 'b', plan.root_exponent, 'source_bits', record.source.bit_length(),
              'terminal_bits', record.terminal.bit_length(), 'ranks', codec.ranks(cert))
    modular_checks = 0
    for m in range(1, 17):
        for lift in (0, 1, 7):
            plan = completed_plan(m, lift)
            word = plan_word(plan)
            cert = codec_from_word(word)
            for modulus in (2**20, 3**8, 19**3, 295, 1024*729):
                need(plan_source_mod(plan, modulus) == inverse_word_mod(word, modulus),
                     'independent symbolic modular readers')
                modular_checks += 1
            need(plan_source_mod(plan, 2048) == 155, 'symbolic initial guard')
            source_fuel_residue = plan_source_mod(plan, 1 << (10*m+2))
            need(v2(295*source_fuel_residue-669) == 10*m+1, 'symbolic exact repeat count')
            need(plan_source_mod(plan, 2**20) == codec.mod2(cert, 20), 'old binary codec')
            need(plan_source_mod(plan, 3**8) == codec.mod3(cert, 8), 'old ternary codec')
            need(codec.ranks(cert) == (6*m+1, 16*m+plan.root_exponent+1),
                 'symbolic first-hit ranks')
    print('Symbolic completed family: m=1..16, lifts0/1/7;', modular_checks,
          'closed-form/inverse-letter residue comparisons; no large-source expansion.')
    phase_checks, previous = 0, None
    for m in range(1, 9):
        exponent, period = power3_exponent(m)
        modulus = 1 << (10*m+1)
        need(pow(3, period, modulus) == 1 and pow(3, period//2, modulus) != 1,
             'independent exact power-of-three clock')
        if previous is not None:
            need(exponent % previous[1] == previous[0], 'retained earlier exponent digits')
        for lift in (0, 1, 2, 7, 31):
            residue = pow(3, exponent+lift*period, modulus)
            need((295*residue-669) % modulus == 0, 'power-of-three repeat guard')
            phase_checks += 1
        previous = exponent, period
        print('Power-of-three prefix family:', 'm', m, 'exponent', exponent, 'mod', period)
    need(power3_exponent(1) == (483, 512), 'inherited first power-of-three phase')
    print('Power-of-three controls:', phase_checks,
          '; modular guards only, with no terminal home claim or huge source expansion.')
    record = certificates[0]
    hostiles = [lambda: verify(replace(record, source=record.source+2)),
                lambda: verify(replace(record, repetitions=record.repetitions+1)),
                lambda: verify(replace(record, terminal=record.terminal+2)),
                lambda: verify(replace(record, terminal_word=record.terminal_word+(2,))),
                lambda: verify(replace(record, terminal=1, terminal_word=())),
                lambda: verify(replace(record, repetitions=True)),
                lambda: repeat_endpoint(True, 1), lambda: repeat_endpoint(155.0, 1),
                lambda: repeat_endpoint(154, 1), lambda: repeat_endpoint(155, 0),
                lambda: audit_plan(CompletedPlan(1, 890)),
                lambda: audit_plan(CompletedPlan(True, 888)),
                lambda: decode_word((729.0, 1024, 2971)),
                lambda: materialize(completed_plan(2), bit_cap=1000)]
    for hostile in hostiles:
        rejected(hostile)
    print('Explicit API/type/root/fuel/phase/expansion-cap hostiles rejected:', len(hostiles),
          '; plus128 sharp fuel failures and H-not-word hostile.')
    print('Scope: repeat description compresses proof structure; source and arithmetic bits remain.')
    print('No finite guard bank, unconditional source coverage, or self-grounding cycle is claimed.')


if __name__ == '__main__':
    main()
