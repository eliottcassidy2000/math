"""Grounded terminal lifts are closed under guarded contracting inverse words.

Exact symbolic families, source-preserving point extension, and finite controls.
No trajectory discovery or universal Collatz assertion. Run from repository root.
"""
from dataclasses import dataclass, replace
from itertools import product

import collatz_uncovered_join_routes_20261007 as routes
import collatz_terminal_lifts_20261007 as lifts
import collatz_reset2_rules_20261007 as rules

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def natural(n, lower=0):
    need(type(n) is int and n >= lower, 'exact integer with stated lower bound')


@dataclass(frozen=True)
class Seed:
    source: int
    word: tuple


@dataclass(frozen=True)
class Family:
    seed: Seed
    inverse_word: tuple
    residue: int
    period: int
    least: int


@dataclass(frozen=True)
class Point:
    family: Family
    parameter: int


def audit_seed(seed):
    need(type(seed) is Seed, 'typed seed')
    routes.odd(seed.source)
    need(seed.source > 1, 'supplied nonroot seed')
    routes.letters(seed.word)
    need(bool(seed.word), 'nonempty supplied ROOT certificate')
    need(routes.replay(seed.source, seed.word)[-1] == 1,
         'first-hit ROOT certificate, not a conditional obligation')
    p, q, b = routes.carrier(seed.word)
    need(p*seed.source+b == q, 'seed carrier endpoint')
    return seed.source, len(seed.word), sum(seed.word), p, q


def quotient_residue(length, parameter, digits):
    """(4^(3^length*parameter)-1)/3^(length+1), modulo 3^digits.

    Stabilization prevents a long seed word from forcing a huge exponent.
    This standalone arithmetic helper validates even its non-seed inputs.
    """
    natural(length)
    natural(parameter)
    natural(digits)
    if digits == 0:
        return 0
    modulus = 3**digits
    if length >= digits-1:
        e = digits-1
        denominator = 3**(e+1)
        w = (pow(4, 3**e, denominator*modulus)-1)//denominator
        return parameter*w % modulus
    denominator = 3**(length+1)
    return ((pow(4, 3**length*parameter, denominator*modulus)-1)
            // denominator) % modulus


def _base_residue(ctx, parameter, radix, digits):
    a, length, cost, p, _ = ctx
    if digits == 0:
        return 0
    modulus = radix**digits
    if radix == 3:
        return (a+4*pow(2, cost, modulus)
                *quotient_residue(length, parameter, digits)) % modulus
    exponent = 2*p*parameter
    power = 0 if exponent >= digits else 1 << exponent
    return (a+4*pow(2, cost, modulus)*(power-1)
            *pow(3*p, -1, modulus)) % modulus


def base_residue(seed, parameter, radix, digits):
    ctx = audit_seed(seed)
    natural(parameter)
    natural(digits)
    need(type(radix) is int and radix in (2, 3), 'binary or ternary reader')
    return _base_residue(ctx, parameter, radix, digits)


def _compile(seed, ctx, inverse_word):
    routes.letters(inverse_word)
    p, q, b = routes.carrier(inverse_word)
    need(not inverse_word or q < p, 'contracting inverse word or identity')
    target = b*pow(q, -1, p) % p
    residue, period = 0, 1
    for j in range(1, len(inverse_word)+1):
        next_period = 3*period
        candidates = [residue+t*period for t in range(3)
                      if _base_residue(ctx, residue+t*period, 3, j)
                      == target % next_period]
        need(len(candidates) == 1, 'ternary isometry gives exactly one lifted digit')
        residue, period = candidates[0], next_period
    a, _, _, seed_p, seed_q = ctx
    # q*N_s>b iff C*4^(seed_p*s)>H; keep the strict boundary.
    c = 4*seed_q*q
    h = 3*seed_p*(b-q*a)+c
    floor_ratio = h//c
    cut = 0 if floor_ratio < 0 else ((floor_ratio.bit_length()+2*seed_p-1)
                                    // (2*seed_p))
    least = residue+max(0, (cut-residue+period-1)//period)*period
    return Family(seed, inverse_word, residue, period, least)


def compile_family(seed, inverse_word=()):
    return _compile(seed, audit_seed(seed), inverse_word)


def audit_family(family):
    need(type(family) is Family, 'typed family')
    for field in ('residue', 'period', 'least'):
        natural(getattr(family, field), 1 if field == 'period' else 0)
    ctx = audit_seed(family.seed)
    need(_compile(family.seed, ctx, family.inverse_word) == family,
         'recomputed canonical phase and exact finite-head cutoff')
    return ctx


def audit_point(point):
    need(type(point) is Point, 'typed point')
    ctx = audit_family(point.family)
    natural(point.parameter)
    f = point.family
    need(point.parameter >= f.least and point.parameter % f.period == f.residue,
         'supplied parameter obeys guard and positivity')
    return ctx


def refine_family(family, new_inverse_word):
    ctx = audit_family(family)
    routes.letters(new_inverse_word)
    pu, qu, _ = routes.carrier(new_inverse_word)
    need(not new_inverse_word or qu < pu, 'new inverse word contracts')
    result = _compile(family.seed, ctx, new_inverse_word+family.inverse_word)
    need(result.residue % family.period == family.residue,
         'new phase refines old phase without losing the original coordinate')
    need(result.least >= family.least, 'positive new child has positive old parent')
    return result


def extend_point(point, new_inverse_word):
    """Never replace an unguarded supplied parameter with a fitting one."""
    audit_point(point)
    result = refine_family(point.family, new_inverse_word)
    s = point.parameter
    if s < result.least or s % result.period != result.residue:
        return None
    child = Point(result, s)
    audit_point(child)
    return child


def root_word(point):
    _, _, _, p, _ = audit_point(point)
    w = point.family.seed.word
    tail = () if point.parameter == 0 else (2*(1+p*point.parameter),)
    return point.family.inverse_word+w+tail


def reseed_zero_phase(point):
    """Absorb the inverse word into a new supplied seed, preserving the source.

    The original point retains the transport sidecar s=p*t. Nonzero phases
    cannot use this canonical reseeding identity; they retain their old form.
    """
    audit_point(point)
    family = point.family
    need(family.residue == 0, 'canonical lossless reseeding requires zero phase')
    p, q, b = routes.carrier(family.inverse_word)
    numerator = q*family.seed.source-b
    need(numerator % p == 0 and point.parameter % p == 0, 'exact seed/parameter transport')
    seed = Seed(numerator//p, family.inverse_word+family.seed.word)
    result = Point(compile_family(seed), point.parameter//p)
    need(root_word(result) == root_word(point), 'identical actual ROOT word, including final valuation')
    return result


def source_residue(point, radix, digits):
    ctx = audit_point(point)
    natural(digits)
    need(type(radix) is int and radix in (2, 3), 'binary or ternary reader')
    p, q, b = routes.carrier(point.family.inverse_word)
    modulus = radix**digits
    if radix == 2:
        return (q*_base_residue(ctx, point.parameter, 2, digits)-b)*pow(p, -1, modulus) % modulus
    # Division by p requires ell additional ternary digits before reduction.
    high = _base_residue(ctx, point.parameter, 3, digits+len(point.family.inverse_word))
    numerator = (q*high-b) % (p*modulus)
    need(numerator % p == 0, 'retained ternary denominator precision')
    return numerator//p


def expand(point, bit_cap=10000):
    a, _, _, p, q = audit_point(point)
    natural(bit_cap, 1)
    s = point.parameter
    if s == 0:
        n = a
        need(n.bit_length() <= bit_cap, 'base fits supplied expansion cap')
    else:
        need(max(a.bit_length(), q.bit_length()+2+2*p*s)+1 <= bit_cap,
             'conservative explicit expansion cap')
        n = a+4*q*((4**(p*s)-1)//(3*p))
    pv, qv, bv = routes.carrier(point.family.inverse_word)
    need((qv*n-bv) % pv == 0, 'literal guarded division')
    child = (qv*n-bv)//pv
    routes.odd(child)
    return child


def rejects(callback):
    try:
        callback()
    except (ValueError, TypeError):
        return True
    return False


def valuation3(n):
    need(type(n) is int and n != 0, 'nonzero exact integer for valuation')
    n, e = abs(n), 0
    while n % 3 == 0:
        n //= 3
        e += 1
    return e


def main():
    # Independent direct quotient path, including both sides of stabilization.
    quotient_cases = 0
    for length in range(9):
        for digits in range(1, 9):
            for s in range(21):
                d = 3**(length+1)
                direct = (pow(4, 3**length*s, d*3**digits)-1)//d
                need(quotient_residue(length, s, digits) == direct,
                     'stable precision reader agrees with direct modular exponentiation')
                quotient_cases += 1
    for digits in range(1, 15):
        need(quotient_residue(10**6, 10**100+7, digits)
             == quotient_residue(digits-1, 10**100+7, digits),
             'word-length-independent stabilized residue')

    seeds = (Seed(3, (1, 4)), Seed(5, (4,)), Seed(13, (3, 4)))
    words = [()]+[w for length in range(1, 4)
                  for w in product((1, 2, 3), repeat=length)
                  if 2**sum(w) < 3**len(w)]
    realized = 0
    for seed in seeds:
        base = compile_family(seed)
        literal = [expand(Point(base, s)) for s in range(81)]
        for s, n in enumerate(literal):
            need(routes.replay(n, root_word(Point(base, s)))[-1] == 1,
                 'independent literal orbit of the lifted terminal')
            for radix in (2, 3):
                for digits in range(7):
                    need(base_residue(seed, s, radix, digits) == n % radix**digits,
                         'small literal base residue')
            if s:
                for t in range(s):
                    need(valuation3(n-literal[t]) == valuation3(s-t),
                         'exact difference valuation, including long shared ternary prefixes')
        for digits in range(1, 6):
            need(len({base_residue(seed, s, 3, digits) for s in range(3**digits)}) == 3**digits,
                 'complete residue universe is permuted')
        for word in words:
            family = compile_family(seed, word)
            p, q, b = routes.carrier(word)
            for s, n in enumerate(literal):
                expected = (q*n-b) > 0 and (q*n-b) % p == 0
                accepted = s >= family.least and s % family.period == family.residue
                need(accepted == expected, 'compiled guard iff actual positive inverse integer')
                if not accepted:
                    continue
                point = Point(family, s)
                child = expand(point)
                need(routes.replay(child, word)[-1] == n,
                     'final integrality suffices for every actual valuation prefix')
                need(routes.replay(child, root_word(point))[-1] == 1,
                     'complete inherited first-hit certificate')
                for radix in (2, 3):
                    for digits in range(5):
                        need(source_residue(point, radix, digits) == child % radix**digits,
                             'inverse source residue with denominator precision')
                need(lifts.recognize_completed(seed.source, seed.word, n).index == s,
                     'independent old recognizer recovers unchanged parameter')
                realized += 1

    generators = ((1,), (1, 2), (1, 1, 1, 2, 1, 1, 4))
    refinements = 0
    for seed in seeds:
        for depth in range(1, 4):
            for choices in product(generators, repeat=depth):
                family = compile_family(seed)
                for word in choices:
                    family = refine_family(family, word)
                    # Do not materialize these potentially enormous sources.
                    point = Point(family, family.least+family.period*10**25)
                    need(len(root_word(point)) == len(family.inverse_word)+len(seed.word)+1,
                         'unbounded symbolic source retains an explicit finite ROOT word')
                    need(0 <= source_residue(point, 3, 12) < 3**12, 'symbolic source reader')
                    refinements += 1

    # A fixed source cannot silently migrate to the phase supplied by refinement.
    point3 = Point(compile_family(seeds[0]), 0)
    need(extend_point(point3, (1,)) is None, 'seed 3 fails inverse G1 guard')
    point5 = Point(compile_family(seeds[1]), 0)
    need(refine_family(point5.family, ()) == point5.family,
         'empty family extension is identity, not a strict decrease')
    need(extend_point(point5, ()) == point5,
         'empty point extension preserves the same integer and parameter')
    child3 = extend_point(point5, (1,))
    need(child3 is not None and child3.parameter == 0 and expand(child3) == 3,
         'seed 5 extends to 3 at exactly the supplied parameter')
    need(extend_point(child3, (1,)) is None, 'next forced extension is rejected, not recycled')

    # Exact recursive reseeding when the native phase contains parameter zero.
    for t in range(11):
        old = Point(child3.family, 3*t)
        new = reseed_zero_phase(old)
        need(new.family.seed == seeds[0] and new.parameter == t,
             'G1 over seed 5 becomes precisely the canonical seed-3 lift')
        need(expand(old) == expand(new), 'lossless reseeding preserves each literal source')
    # Nonzero phases have a different terminal boundary, not the same family.
    nonzero = compile_family(seeds[0], (1,))
    need(nonzero.residue == 1 and nonzero.period == 3, 'nonzero phase hostile')
    first = Point(nonzero, 1)
    c, cw = expand(first), root_word(first)
    need((c, cw) == (828503, (1, 1, 4, 20)), 'small exact rebase boundary')
    canonical = compile_family(Seed(c, cw))
    need(rejects(lambda: reseed_zero_phase(first)), 'nonzero phase cannot silently change the family')
    for t in range(5):
        for u in range(5):
            old_n = expand(Point(nonzero, 1+3*t))
            new_n = expand(Point(canonical, u))
            need((old_n == new_n) == (t == u == 0),
                 'nonzero-phase naive canonical reseeding shares only the base point')

    # Finite-head positivity hostile: inverse ones at a small seed may be negative.
    family = compile_family(seeds[1], (1,)*9)
    need(rejects(lambda: audit_point(Point(family, 0))), 'nonpositive or illegal finite head excluded')
    small = Point(compile_family(seeds[1]), 0)
    for bad in (True, 0.0, -1):
        need(rejects(lambda bad=bad: audit_point(replace(small, parameter=bad))),
             'numeric aliases cannot bypass public validation')
    need(rejects(lambda: audit_family(replace(small.family, residue=False))), 'forged bool residue')
    need(rejects(lambda: audit_family(replace(small.family, period=2))), 'forged period')
    need(rejects(lambda: audit_seed(Seed(5, (4, 2)))), 'first-hit ROOT padding excluded')
    need(rejects(lambda: compile_family(seeds[1], (4,))), 'expanding inverse outside claimed scope')
    need(rejects(lambda: expand(Point(compile_family(seeds[0]), 10**100))), 'explicit allocation cap')

    exponent, supplied = rules.frozen_seed_word()
    seed = Seed((1 << exponent)-1, supplied)
    need(exponent == 1457, 'authenticated previously grounded Mersenne seed')
    basepoint = Point(compile_family(seed), 0)
    h = extend_point(basepoint, (1, 2))
    need(h is not None, 'supplied h child is legal without changing the seed parameter')
    z = extend_point(h, (1, 2)+(1,)*5)
    need(z is not None, 'supplied stripped terminal is in the same grounded closure')
    rows = []
    for name, point in (('h1', h), ('z', z)):
        n, word = expand(point), root_word(point)
        need(routes.replay(n, word)[-1] == 1, 'selected child directly replays to first ROOT')
        rows.append((name, n.bit_length(), len(word), sum(word), point.family.residue))
        for radix in (2, 3):
            for digits in (1, 2, 3, 7, 15):
                need(source_residue(point, radix, digits) == n % radix**digits,
                     'long supplied seed uses fast exact readers')
    need(rows == [('h1', 1457, 7349, 13105, 0), ('z', 1454, 7356, 13113, 0)],
         'independent concrete child certificate summary')
    far = Point(z.family, z.family.period*10**100)
    need(source_residue(far, 3, 30) < 3**30, 'unexpanded large-seed descendant family')
    need(len(root_word(far)) == 7357, 'one symbolic final valuation grounds the distant member')
    rebased = reseed_zero_phase(far)
    need(rebased.parameter == 10**100 and rebased.family.seed.source == expand(z),
         'large child family recanonicalizes at its actual grounded child')
    for radix in (2, 3):
        for digits in (1, 5, 15, 30):
            need(source_residue(rebased, radix, digits) == source_residue(far, radix, digits),
                 'lossless large-family reseeding retains exact binary and ternary residues')
    print('STATUS: PROVED symbolic guarded-family closure; FINITE-EXACT controls')
    print('QUOTIENT_DIRECT_CASES', quotient_cases)
    print('LITERAL_GROUNDED_CHILD_POINTS', realized)
    print('SYMBOLIC_REFINEMENTS', refinements)
    print('SUPPLIED_SEED', exponent, len(supplied), sum(supplied))
    print('GROUNDED_CHILD_ROWS(name,bits,steps,cost,parameter_phase)', rows)
    print('SOURCE_PRESERVING_GUARD_HOSTILE: seed 3 cannot take G1 at parameter 0')
    print('RESEEDING: zero phase preserves every member; nonzero phase shares only the base point')
    print('OPEN: covering an arbitrary supplied source by these grounded families')
    print('CHECKS', CHECKS)


if __name__ == '__main__':
    main()
