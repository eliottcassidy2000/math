"""Joint ternary/19-adic addresses of completed routes of odd rank at most four.

The third inverse exponent retains the division precision that a one-row
observer loses. This constructs new sources, never certifies a residue alias.
All checks survive -O; no orbit search supplies the address compiler's roots.
"""
from dataclasses import dataclass
from itertools import product
from math import gcd
import json

import inverse_ray_ternary_addresses_20261004 as codec
from mod19_route_lifts_20261004 import inverse_exponent, certificate_residue


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def natural(value, name, positive=False):
    need(type(value) is int and value >= int(positive), name)


def valuation(value, prime):
    need(value != 0, 'finite valuation')
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


def carry(word):
    result, total = 0, 0
    for exponent in word:
        result = 3*result + 2**total
        total += exponent
    return result


def inverse_literal(word, hub, allow_root=False):
    value = hub
    for exponent in reversed(word):
        numerator = 2**exponent*value-1
        if numerator % 3:
            return None
        value = numerator//3
        if value <= int(not allow_root):
            return None
    return value


@dataclass(frozen=True)
class JointRay:
    source: int
    prefix: tuple
    hub: int

    @property
    def coefficient(self):
        return self.hub*2**sum(self.prefix)

    @property
    def offset(self):
        return carry(self.prefix)

    def word(self, parameter):
        natural(parameter, 'nonnegative exact parameter')
        a,b,c = self.prefix
        return (a,b,c+18*parameter)+(() if self.hub == 1 else (4,))

    def residue(self, parameter, modulus):
        natural(parameter, 'nonnegative exact parameter')
        natural(modulus, 'positive modulus', True)
        # Keep the denominator before reducing, including at powers of three.
        numerator = (self.coefficient*pow(2,18*parameter,27*modulus)-self.offset) % (27*modulus)
        need(numerator % 27 == 0, 'three integral inverse divisions')
        return numerator//27

    def certificate(self, parameter):
        cert = codec.ROOT
        for exponent in reversed(self.word(parameter)):
            cert = inverse_exponent(cert,exponent)
        codec.audit_certificate(cert)
        return cert

    def literal(self, parameter, bit_cap=10000):
        natural(parameter, 'nonnegative exact parameter')
        natural(bit_cap, 'positive expansion cap', True)
        need(18*parameter+sum(self.prefix)+self.hub.bit_length() <= bit_cap,
             'declared expansion cap')
        return (self.coefficient*2**(18*parameter)-self.offset)//27


def make_bank():
    # Finite witness search, not an all-height minimum claim.
    choices = {}
    for hub in (1,5):
        for word in product(range(1,19),repeat=3):
            source = inverse_literal(word,hub)
            if source is None or (hub == 5 and source % 19 != 6):
                continue
            residue = source % 19
            key = (source,sum(word),word,hub)
            if residue not in choices or key < choices[residue]:
                choices[residue] = key
    need(set(choices) == set(range(19)), 'one explicit three-inverse-prefix witness per19 residue')
    rays = tuple(JointRay(choices[r][0],choices[r][2],choices[r][3]) for r in range(19))
    need(sum(ray.hub == 1 for ray in rays) == 18 and rays[6].hub == 5,
         'eighteen rank-three rays and one rank-four ray')
    for ray in rays:
        need(ray.coefficient-ray.offset == 27*ray.source, 'base affine source identity')
        need(gcd(ray.coefficient,57) == 1, 'two unit derivative coordinates')
    return rays


BANK = make_bank()


@dataclass(frozen=True)
class JointAddress:
    ternary_address: int
    ternary_depth: int
    prime_address: int
    prime_depth: int
    ray: JointRay
    parameter: int
    period: int
    certificate: object


def compile_address(ternary_address, ternary_depth, prime_address, prime_depth, lift=0):
    natural(ternary_depth, 'nonnegative ternary precision')
    natural(prime_depth, 'positive19 precision', True)
    natural(ternary_address, 'canonical ternary address')
    natural(prime_address, 'canonical19 address')
    natural(lift, 'nonnegative progression lift')
    m3, mp = 3**ternary_depth, 19**prime_depth
    need(ternary_address < m3 and prime_address < mp, 'canonical address ranges')
    ray = BANK[prime_address % 19]
    t3 = 0
    derivative3 = ray.coefficient % 3
    for level in range(ternary_depth):
        unit = 3**level
        delta = ternary_address-ray.residue(t3,3*unit)
        need(delta % unit == 0, 'retained ternary digits')
        digit = delta//unit*pow(derivative3,-1,3) % 3
        t3 += digit*unit
    tp = 0
    derivative19 = ray.coefficient*pow(9,-1,19) % 19
    for level in range(1,prime_depth):
        unit = 19**level
        delta = prime_address-ray.residue(tp,19*unit)
        need(delta % unit == 0, 'retained19 digits')
        digit = delta//unit*pow(derivative19,-1,19) % 19
        tp += digit*19**(level-1)
    pperiod = 19**(prime_depth-1)
    parameter = t3 if pperiod == 1 else t3+m3*((tp-t3)*pow(m3,-1,pperiod) % pperiod)
    period = m3*pperiod
    need(0 <= parameter < period, 'canonical joint block parameter')
    parameter += lift*period
    cert = ray.certificate(parameter)
    need(codec.mod3(cert,ternary_depth) == ternary_address,
         'independent ternary certificate evaluator')
    need(certificate_residue(cert,mp) == prime_address,
         'independent prime certificate evaluator')
    need(codec.ranks(cert)[0] == (3 if ray.hub == 1 else 4), 'uniform completed odd rank')
    return JointAddress(ternary_address,ternary_depth,prime_address,prime_depth,
                        ray,parameter,period,cert)


def literal_check(ray, parameter):
    source = ray.literal(parameter)
    current, observed = source, []
    while current != 1:
        need(len(observed) < 4, 'finite exact first-hit check')
        value, exponent = 3*current+1, 0
        while value % 2 == 0:
            value //= 2
            exponent += 1
        current = value
        observed.append(exponent)
    need(tuple(observed) == ray.word(parameter), 'independent literal first-hit valuations')
    cert = ray.certificate(parameter)
    need(codec.expand(cert) == source and codec.encode_source(source) == cert,
         'certificate source identity and independent first-hit encoder')
    return source


def main():
    need(valuation(2**18-1,3) == 3 and valuation(2**18-1,19) == 1,
         'exact two-prime generator valuations')
    bank_report = [dict(residue=r, source=ray.source, prefix=ray.prefix,hub=ray.hub,
                        coefficient=ray.coefficient,offset=ray.offset)
                   for r,ray in enumerate(BANK)]
    literal = 0
    for ray in BANK:
        for t in range(12):
            literal_check(ray,t)
            literal += 1
    difference_checks = 0
    for ray in BANK:
        values = [ray.literal(t) for t in range(25)]
        for t in range(25):
            for s in range(t):
                difference = values[t]-values[s]
                need(valuation(difference,3) == valuation(t-s,3), 'ternary isometry')
                need(valuation(difference,19) == 1+valuation(t-s,19), '19 scaled isometry')
                difference_checks += 1
    # Enumerate every parameter and every joint address in the declared boxes.
    # The direct parameter language supplies an independent inverse lookup.
    address_checks = 0
    for a in range(4):
        for k in range(1,3):
            period = 3**a*19**(k-1)
            for residue,ray in enumerate(BANK):
                direct = {(ray.residue(t,3**a),ray.residue(t,19**k)):t for t in range(period)}
                need(len(direct) == period, 'direct joint language has no collisions')
                expected = {(x,residue+19*y) for x in range(3**a) for y in range(19**(k-1))}
                need(set(direct) == expected, 'complete joint address language')
                for x,y in expected:
                    result = compile_address(x,a,y,k)
                    need(result.parameter == direct[(x,y)], 'digit compiler equals exhaustive inverse lookup')
                    lifted = compile_address(x,a,y,k,lift=1)
                    need(lifted.parameter == result.parameter+period, 'all-height progression parameter')
                    address_checks += 1
    # Allow formal root padding, so the finite lower-rank control overincludes.
    low_rank = {}
    for length in range(1,4):
        residues = set()
        for word in product(range(1,19),repeat=length):
            value = inverse_literal(word,1,allow_root=True)
            if value is not None:
                residues.add(value % 19)
        need(6 not in residues, 'critical residue cannot have odd rank at most three')
        low_rank[length] = sorted(residues)
    # Three inverse divisions, not one inverse-row parameter, explain freedom.
    critical = BANK[6]
    need(critical.source == 1507 and critical.prefix == (1,7,5), 'critical ray witness')
    need({critical.residue(t,3) for t in range(3)} == {0,1,2}, 'all ternary source rows occur')
    # Same joint residue still does not identify the actual integer.
    sample = compile_address(1,2,6,2)
    sample_source = sample.ray.literal(sample.parameter,bit_cap=20000)
    alias = sample_source+2*3**2*19**2
    need(alias != sample_source and alias % 9 == sample_source % 9 and alias % 361 == sample_source % 361,
         'a constructed certificate is not a certificate for its residue alias')
    huge = compile_address(3**100-1,100,19**80-1,80,lift=1)
    lower,upper = codec.bit_bounds(huge.certificate)
    need(lower > 10**140, 'symbolic source cannot be expanded')
    for args in ((True,0,0,1),(0,-1,0,1),(0,0,0,0),(3,1,0,1),(0,0,19,1)):
        try:
            compile_address(*args)
        except ArithmeticError:
            pass
        else:
            raise ArithmeticError('invalid address accepted')
    print(json.dumps(dict(status='PROVED joint-address rank-four construction; FINITE-EXACT declared controls',
          rays=bank_report,literal_first_hit_checks=literal,two_prime_difference_checks=difference_checks,
          joint_address_universe=dict(ternary_depths=list(range(4)),prime_depths=[1,2],all_addresses=address_checks),
          low_rank_residues=low_rank,
          huge_symbolic=dict(ternary_depth=100,prime_depth=80,parameter_bits=huge.parameter.bit_length(),
                             odd_rank=codec.ranks(huge.certificate)[0],source_bit_lower_bound=lower),
          boundary='This compiler constructs a new certified source in a requested joint address; no arbitrary-source convergence follows.'),indent=2))
    print('PASS: exact arithmetic and independent certificate readers; all checks active under -O')


if __name__ == '__main__':
    main()
