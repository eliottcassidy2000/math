#!/usr/bin/env python3
"""Certified leading binary digits of A000568, without building its full integer.

Davis/Burnside, the non-singleton cycle-core expansion, and Session98's
finite moved-support tail are inherited. Here the rational cutoff form
supports an adaptive exact-head API with an explicit rounding-distance
check. No float or stored sequence value enters the computation.
"""
from fractions import Fraction
from math import factorial, gcd


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def cores(bound):
    """All odd partitions with parts >=3 and total mass <=bound, once each."""
    def rec(parts, mass, length, energy, z, maximum):
        yield parts, mass, length, energy, z
        for ell in range(min(maximum, bound-mass) | 1, 2, -2):
            if ell > bound-mass:
                continue
            cross = sum(gcd(ell, j) * k for j, k in parts)
            zz, ee = z, energy
            for k in range(1, (bound-mass)//ell+1):
                zz *= ell*k
                ee += (ell-1)//2+cross+ell*(k-1)
                yield from rec(parts+((ell,k),), mass+ell*k,
                               length+k, ee, zz, ell-2)
    yield from rec((),0,0,0,1,bound)


def normalized_prefix(n, cutoff):
    falling = [1]
    for j in range(cutoff):
        falling.append(falling[-1]*(n-j))
    total = Fraction(0)
    count = 0
    for _, mass, length, energy, z in cores(cutoff):
        deficit = mass*(mass-1)//2-energy+(mass-length)*(n-mass)
        need(deficit >= 0, 'nonnegative orbit deficit')
        total += Fraction(falling[mass], z << deficit)
        count += 1
    return total, count


def normalized_tail(n, cutoff):
    """Explicit conservative upper bound for omitted normalized weight."""
    if cutoff >= n:
        return Fraction(0)
    c = (n-1)//3
    denominator = 1 << c
    if n >= denominator:
        return None
    r = Fraction(n, denominator)
    first_mass = max(3,cutoff+1)
    if first_mass == 4:
        first_mass = 5  # no odd >=3 core has total mass four
    return r**first_mass/(1-r)


def head(n, bits=128, *, max_core_mass=None):
    """Return exact bit length, leading min(bits,bit_length) bits, and a proof.

    A finite optional max_core_mass caps work and raises if certification has
    not succeeded. With no cap the finite cutoff n gives an exact fallback.
    Small n use this fallback directly. All comparisons are exact rationals.
    """
    if type(n) is not int or n < 0 or type(bits) is not int or bits < 1:
        raise ValueError('n >=0 and bits >=1 must be literal integers')
    if max_core_mass is not None and (type(max_core_mass) is not int or max_core_mass < 0):
        raise ValueError('max_core_mass must be a nonnegative integer')
    cap = n if max_core_mass is None else min(n,max_core_mass)
    fact = factorial(n)
    full_power = n*(n-1)//2
    # Predicted binary magnitude of the identity contribution. A certified
    # quotient determines the true magnitude, including boundary crossings.
    shift = max(0, full_power-fact.bit_length()+1-bits)
    scaled_identity = Fraction(1 << (full_power-shift), fact)
    if n < 13:
        attempts = [n] if cap == n else [cap]
    else:
        attempts = [0]+list(range(3,cap+1))
        if cap in (1,2):
            attempts.append(cap)
    tried = []
    for cutoff in attempts:
        tail = normalized_tail(n,cutoff)
        if tail is None:
            continue
        prefix, count = normalized_prefix(n,cutoff)
        lower = scaled_identity*prefix
        upper = scaled_identity*(prefix+tail)
        lo = lower.numerator//lower.denominator
        hi = upper.numerator//upper.denominator
        tried.append({'cutoff':cutoff,'cores':count,'integer_width':hi-lo})
        if lo != hi:
            continue
        need(lo > 0, 'positive head')
        excess = max(0,lo.bit_length()-bits)
        leading = lo >> excess
        actual_shift = shift+excess
        width = (upper-lower)
        return {
            'n':n,'requested_bits':bits,'bit_length':lo.bit_length()+shift,
            'head_bits':leading.bit_length(),'head_hex':hex(leading),
            'discarded_bits':actual_shift,'cutoff':cutoff,
            'core_types':count,'attempts':tried,
            'factorial_bits':fact.bit_length(),
            'scaled_numerator_bits':full_power-shift+1,
            'lower_numerator_hex':hex(lower.numerator),
            'lower_denominator_hex':hex(lower.denominator),
            'upper_numerator_hex':hex(upper.numerator),
            'upper_denominator_hex':hex(upper.denominator),
            'scaled_floor_hex':hex(lo),
            'width_numerator_bits':width.numerator.bit_length(),
            'width_denominator_bits':width.denominator.bit_length(),
            'tail_numerator_bits':tail.numerator.bit_length(),
            'tail_denominator_bits':tail.denominator.bit_length(),
        }
    raise ArithmeticError(f'no certified head with core mass <= {cap}; increase cap')


def reconstruct(lower, upper, residue, modulus):
    """Unique integer in a rational interval with a supplied modular residue.

    The residue must be established independently. No class-count conjecture
    or approximate bit length may be substituted for it.
    """
    if type(modulus) is not int or modulus < 1:
        raise ValueError('positive integer modulus required')
    if type(residue) is not int:
        raise ValueError('integer residue required')
    lower,upper=Fraction(lower),Fraction(upper)
    if upper < lower or upper-lower >= modulus:
        raise ValueError('ordered interval of width strictly below modulus required')
    residue %= modulus
    offset=(lower-residue)/modulus
    k=-((-offset.numerator)//offset.denominator)
    answer=residue+modulus*k
    if not lower <= answer <= upper:
        raise ValueError('interval and residue are incompatible')
    return answer


if __name__ == '__main__':
    import argparse,json,sys
    if hasattr(sys,'set_int_max_str_digits'):
        sys.set_int_max_str_digits(0)
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n',type=int)
    parser.add_argument('--bits',type=int,default=128)
    parser.add_argument('--max-core-mass',type=int,default=None)
    args=parser.parse_args()
    print(json.dumps(head(args.n,args.bits,max_core_mass=args.max_core_mass),indent=2))
