"""Nineteen explicit completed-route rays, with a sharp odd-depth-four bound.

This constructs new integers in requested residue classes, not certificates
for arbitrary supplied integers. Exact checks run unchanged under python -O.
"""
from collections import Counter
from dataclasses import dataclass
from math import gcd

from inverse_ray_ternary_addresses_20261004 import (
    ROOT, audit_certificate, bit_bounds, chain, encode_source, expand,
    exponent, extend, kappa, mod3, odd_step, ranks,
)


def need(condition, message):
    if not condition:
        raise ValueError(message)


def natural(value, name, positive=False):
    need(type(value) is int and value >= int(positive), name)


def valuation(value, prime):
    need(value != 0, 'nonzero valuation argument')
    result = 0
    while value % prime == 0:
        result += 1
        value //= prime
    return result


@dataclass(frozen=True)
class Ray:
    seed: int
    target: int
    first_exponent: int

    def word(self, parameter):
        natural(parameter, 'nonnegative integer ray parameter')
        if self.seed == 6:
            return (1, 7+18*parameter, 5, 4)
        head = self.first_exponent+18*parameter
        return (head,) if self.target == 1 else (head, 4)

    def affine_data(self):
        if self.seed == 6:
            return (256*53, 5, 9)
        return (2**self.first_exponent*self.target, 1, 3)

    def residue(self, parameter, modulus):
        natural(parameter, 'nonnegative integer ray parameter')
        natural(modulus, 'positive integer modulus', positive=True)
        need(gcd(modulus, 3) == 1, 'ray residue modulus coprime to three')
        if modulus == 1:
            return 0
        coefficient, offset, denominator = self.affine_data()
        return ((coefficient*pow(2, 18*parameter, modulus)-offset)
                * pow(denominator, -1, modulus)) % modulus

    def digit_coefficient(self):
        return 7 if self.seed == 6 else (3*self.seed+1) % 19


def make_bank():
    rays = {}
    for target, exponents in ((1, range(4, 21, 2)), (5, range(1, 18, 2))):
        for k in exponents:
            source = (2**k*target-1)//3
            seed = source % 19
            need(seed not in rays, 'disjoint noncritical residue rays')
            rays[seed] = Ray(seed, target, k)
    need(set(rays) == set(range(19))-{6}, 'eighteen noncritical seeds')
    rays[6] = Ray(6, 53, 7)
    return tuple(rays[r] for r in range(19))


BANK = make_bank()


def inverse_exponent(parent, k):
    natural(k, 'positive integer inverse exponent', positive=True)
    audit_certificate(parent)
    parent_mod9 = mod3(parent, 2)
    value = (pow(2, k, 9)*parent_mod9-1) % 9
    need(value % 3 == 0, 'inverse exponent has an integral source')
    row = value//3
    least = kappa(parent_mod9, row)
    need(k >= least and (k-least) % 6 == 0, 'inverse exponent matches row')
    return extend(parent, row, (k-least)//6)


def certificate_for(ray, parameter):
    need(isinstance(ray, Ray) and ray in BANK, 'one of the nineteen fixed rays')
    cert = ROOT
    for k in reversed(ray.word(parameter)):
        cert = inverse_exponent(cert, k)
    audit_certificate(cert)
    return cert


def certificate_residue(cert, modulus):
    natural(modulus, 'positive integer modulus', positive=True)
    need(gcd(modulus, 3) == 1, 'certificate modulus coprime to three')
    audit_certificate(cert)
    if modulus == 1:
        return 0
    current, inverse3 = 1, pow(3, -1, modulus)
    for node in reversed(chain(cert)):
        current = (pow(2, exponent(node), modulus)*current-1)*inverse3 % modulus
    return current


@dataclass(frozen=True)
class CompiledAddress:
    address: int
    precision: int
    ray: Ray
    parameter: int
    certificate: object


def compile_address(address, precision, lift=0):
    """Compile a NEW home-reaching integer in address mod19**precision.

    lift enumerates infinitely many distinct sources in the requested class.
    The canonical parameter (lift=0) has precision-1 base-19 digits.
    No source integer is expanded by this operation.
    """
    natural(precision, 'positive integer 19-adic precision', positive=True)
    natural(address, 'canonical nonnegative address')
    natural(lift, 'nonnegative integer lift quotient')
    need(address < 19**precision, 'canonical address below the modulus')
    ray, parameter = BANK[address % 19], 0
    inverse_digit = pow(ray.digit_coefficient(), -1, 19)
    for level in range(1, precision):
        power = 19**level
        current = ray.residue(parameter, 19*power)
        need((address-current) % power == 0, 'previous address digits retained')
        digit = ((address-current)//power*inverse_digit) % 19
        parameter += digit*19**(level-1)
    parameter += lift*19**(precision-1)
    cert = certificate_for(ray, parameter)
    need(certificate_residue(cert, 19**precision) == address,
         'independent certificate reader matches requested address')
    return CompiledAddress(address, precision, ray, parameter, cert)


def literal_check(ray, parameter):
    cert = certificate_for(ray, parameter)
    source = expand(cert, bit_cap=5000)
    coefficient, offset, denominator = ray.affine_data()
    need(source == (coefficient*2**(18*parameter)-offset)//denominator,
         'closed formula and backward certificate agree')
    current, word = source, []
    while current != 1:
        need(len(word) < 4, 'first-hit route has at most four odd steps')
        current, k = odd_step(current)
        word.append(k)
    need(tuple(word) == ray.word(parameter), 'independent forward exact valuations')
    need(encode_source(source) == cert, 'first-hit codec round trip')
    need(ranks(cert) == (len(word), sum(word)+len(word)), 'ordinary and odd clocks')
    return source


def main():
    print('MOD19 COMPLETED ROUTE LIFTS: PROVED families; FINITE-EXACT controls')
    need(pow(2, 9, 19) == 18 and pow(2, 6, 19) != 1,
         'two has exact order eighteen modulo nineteen')
    need(valuation(2**18-1, 19) == 1 and ((2**18-1)//19) % 19 == 3,
         'nondegenerate nineteen-adic digit coefficient')
    print('order19(2)=18; v19(2^18-1)=1; quotient mod19=3')

    print('seed : first source ; forward valuation word ; next-digit coefficient')
    for ray in BANK:
        source = literal_check(ray, 0)
        print(f'{ray.seed:2} : {source} ; {ray.word(0)} ; {ray.digit_coefficient()}')

    # Complete depth-two obstruction table: modulo27 and19 both use exponent
    # period18. A nonroot first predecessor can have exponent20, so residue2
    # is retained even though its least representative is the excluded root.
    table = []
    for a in (2, 4, 8, 10, 14, 16):
        target_mod9 = ((pow(2, a, 27)-1) % 27)//3
        target_mod19 = (pow(2, a, 19)-1)*pow(3, -1, 19) % 19
        b = next(k for k in range(1, 19)
                 if pow(2, k, 19)*target_mod19 % 19 == 1)
        numerator_mod9 = (pow(2, b, 9)*target_mod9-1) % 9
        legal = numerator_mod9 % 3 == 0
        row = numerator_mod9//3 if legal else None
        need(not legal or row == 0, 'every depth-two multiple19 is a ternary leaf')
        table.append((a, b, legal, row))
    print('depth-two multiple19 table (root exponent, next exponent, legal, row):', table)
    need(literal_check(BANK[6], 0) == 1507, 'sharp depth-four critical witness')

    # Independent forward replay, exact valuations of differences, and all
    # nineteen one-digit lifts at each of eight levels.
    literal_count = difference_count = digit_count = 0
    for ray in BANK:
        values = [literal_check(ray, t) for t in range(25)]
        literal_count += len(values)
        for s in range(25):
            for t in range(s+1, 25):
                need(valuation(values[t]-values[s], 19) == 1+valuation(t-s, 19),
                     'scaled nineteen-adic isometry')
                difference_count += 1
        for level in range(1, 9):
            modulus, t0 = 19**(level+1), 37
            base = ray.residue(t0, modulus)
            children = set()
            for digit in range(19):
                t = t0+digit*19**(level-1)
                value = ray.residue(t, modulus)
                expected = (base+digit*ray.digit_coefficient()*19**level) % modulus
                need(value == expected, 'exact next address digit')
                children.add(value)
                digit_count += 1
            need(len(children) == 19, 'nineteen distinct children')
    print('literal forward replays:', literal_count,
          '; exact difference valuations:', difference_count,
          '; digit-lift checks:', digit_count)

    compiled_count = 0
    for precision in range(1, 4):
        histogram = Counter()
        for address in range(19**precision):
            result = compile_address(address, precision)
            need(result.parameter < 19**(precision-1), 'canonical parameter length')
            need(ranks(result.certificate)[1] <= 18*19**(precision-1)+5,
                 'uniform canonical ordinary-rank bound')
            if precision > 1:
                parent = compile_address(address % 19**(precision-1), precision-1)
                need(result.parameter % 19**(precision-2) == parent.parameter,
                     'compatible recursive address prefixes')
            histogram[ranks(result.certificate)[0]] += 1
            compiled_count += 1
        print('all addresses at precision', precision, ':', 19**precision,
              '; odd-rank histogram', sorted(histogram.items()))
    need(compiled_count == 7239, 'declared exhaustive address universe')

    lift_checks = 0
    for address in range(19):
        for lift in range(3):
            result = compile_address(address, 1, lift)
            source = literal_check(result.ray, result.parameter)
            need(source % 19 == address, 'independent lifted-source residue')
            lift_checks += 1
    print('least-level sources with three lift quotients:', lift_checks)
    example = compile_address(19**80-13, 80)
    lower, upper = bit_bounds(example.certificate)
    print('symbolic precision80: parameter bits', example.parameter.bit_length(),
          '; odd rank', ranks(example.certificate)[0],
          '; source bit-length bounds', (lower, upper))
    need(lower > 10**100, 'source too large to expand; certificate remains small')

    hostile = [
        lambda: compile_address(True, 1), lambda: compile_address(0, True),
        lambda: compile_address(0, 0), lambda: compile_address(-1, 1),
        lambda: compile_address(19, 1), lambda: compile_address(0, 1, -1),
        lambda: compile_address(0, 1, 1.0),
        lambda: inverse_exponent(ROOT, 2),
        lambda: inverse_exponent(ROOT, 1),
        lambda: inverse_exponent(ROOT, True),
        lambda: BANK[0].residue(0, 3),
    ]
    for action in hostile:
        try:
            action()
        except ValueError:
            pass
        else:
            raise ValueError('hostile input was accepted')
    print('rejected type/domain/root/integrality hostiles:', len(hostile))
    supplied = 19
    constructed = expand(compile_address(supplied % 19, 1).certificate)
    need(constructed == 87381 and constructed != supplied,
         'residue coverage does not certify the supplied source')
    print('source-identity hostile: input19 replaced by certified87381, same residue0')
    print('UNIVERSE: all addresses mod19^a for a1..3; all19 rays, t0..24;')
    print('all19 digit choices at a1..8; one symbolic precision80; exact proof is all-height.')
    print('SCOPE: completed odd rank at most4, sharp; density-zero bank; Collatz coverage OPEN.')


if __name__ == '__main__':
    main()
