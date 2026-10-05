"""One finite root proof, a guarded refuel constructor, and its total decoder.

Run with python -B or python -B -O.  No orbit search or file writes occur.
The family is constructed; membership does not assert universal Collatz.
"""
from itertools import product
import collatz_recursive_dependency_kernel_20261004 as kernel

V = (1, 2, 1, 1, 1, 2)
BASE = 3
BASE_WORD = (1, 4)


def need(condition, message):
    if not condition:
        raise ValueError(message)


def nat(value, minimum=0):
    need(type(value) is int and value >= minimum, 'exact integer domain')


def phase_parent(value):
    nat(value, 3)
    need(value % 8 == 3, 'retained parent phase 3 modulo 8')


def address_guard(address):
    need(type(address) is tuple, 'finite tuple address')
    for branch in address:
        nat(branch)


def kappa(parent):
    """Unique k in 1..729 making S^k(parent) equal111 modulo729.

    Six ternary lifts solve 4^k(3*parent+1)=334 modulo2187.
    Only parent modulo729 is observed; no rootedness is inferred.
    """
    nat(parent)
    residue = parent % 729
    k, period = 0, 1
    for precision in range(1, 7):
        modulus = 3**(precision+1)
        choices = [k+d*period for d in range(3)
                   if (pow(4, k+d*period, modulus)*(3*residue+1)-334) % modulus == 0]
        need(len(choices) == 1, 'unique ternary phase lift')
        k = choices[0]
        period *= 3
    return k or 729


def local_child(parent, branch=0):
    """Arithmetic constructor on the phase domain; home needs a parent proof."""
    phase_parent(parent)
    nat(branch)
    k = kappa(parent)+729*branch
    middle = ((1 << (2*k))*(3*parent+1)-1)//3
    numerator = 1024*middle-669
    need(numerator % 729 == 0, 'H inverse integrality')
    child = numerator//729
    need(child % 2048 == 155 and parent < middle < child,
         'native H guard and strict parent-size decrease')
    return child, k


def source(address):
    address_guard(address)
    value = BASE
    for branch in address:
        value, _ = local_child(value, branch)
    return value


def source_mod(address, modulus):
    """Exact unexpanded reader; each generation retains six parent ternary digits."""
    address_guard(address)
    nat(modulus, 1)
    current_modulus = modulus*729**len(address)
    value = BASE % current_modulus
    for branch in address:
        target_modulus = current_modulus//729
        k = kappa(value)+729*branch
        power = pow(4, k, 3*current_modulus)
        shift = (power-1)//3
        numerator = (1024*(power*value+shift)-669) % current_modulus
        need(numerator % 729 == 0, 'reader retains full division precision')
        value = numerator//729
        current_modulus = target_modulus
    return value


def word(address):
    address_guard(address)
    result = BASE_WORD
    for index, branch in enumerate(address):
        k = kappa(source_mod(address[:index], 729))+729*branch
        need(result[0] == 1, 'first exponent retained by the constructor')
        result = V+(result[0]+2*k+2,)+result[1:]
    return result


def ranks(address):
    values = word(address)
    return len(values), sum(values), len(values)+sum(values)


def decode(value):
    """Total family membership for a supplied literal positive odd integer.

    Returns its unique address, or None. Every accepted local step reduces
    value, and only the literal seed3 is a terminal witness.
    """
    nat(value, 1)
    need(value % 2 == 1, 'positive odd supplied source')
    reverse = []
    while value != BASE:
        if value % 2048 != 155:
            return None
        middle = (729*value+669)//1024
        numerator = 3*middle+1
        valuation = (numerator & -numerator).bit_length()-1
        if valuation < 3 or valuation % 2 == 0:
            return None
        k = (valuation-1)//2
        parent_numerator = (numerator >> (2*k))-1
        if parent_numerator % 3:
            return None
        parent = parent_numerator//3
        if parent < 3 or parent % 8 != 3 or parent >= value:
            return None
        least = kappa(parent)
        if k < least or (k-least) % 729:
            return None
        reverse.append((k-least)//729)
        value = parent
    return tuple(reversed(reverse))


def branch_residue(parent, branch, modulus):
    phase_parent(parent)
    nat(branch)
    nat(modulus, 1)
    k = kappa(parent)+729*branch
    numerator = (1024*pow(4, k, 2187*modulus)*(3*parent+1)-3031) % (2187*modulus)
    need(numerator % 2187 == 0, 'branch reader division precision')
    return numerator//2187


def branch_for_address(parent, target, precision):
    """Unique branch t modulo3^precision for a desired ternary source address."""
    phase_parent(parent)
    nat(precision)
    nat(target)
    period, branch = 1, 0
    for _ in range(precision):
        modulus = 3*period
        current = branch_residue(parent, branch, modulus)
        need((target-current) % period == 0, 'previous ternary digits retained')
        digit = ((target-current)//period) % 3
        branch += digit*period
        period = modulus
    need(branch_residue(parent, branch, period) == target % period,
         'target is a residue of this constructed source')
    return branch


def expected_failure(call):
    try:
        call()
    except ValueError:
        return
    raise ValueError('hostile was accepted')


def main():
    kernel.literal(BASE, BASE_WORD, require_root=True)
    for residue in range(729):
        k = kappa(residue)
        solutions = [a for a in range(1, 730)
                     if (pow(4, a, 2187)*(3*residue+1)-334) % 2187 == 0]
        need(solutions == [k], 'independent exhaustive phase clock')
    print('phase clock: all729 parent residues, unique positive exponent1..729')

    addresses = [address for depth in range(5)
                 for address in product((0, 1, 3), repeat=depth)]
    values, total_letters, modular_checks = set(), 0, 0
    for address in addresses:
        value, route = source(address), word(address)
        need(decode(value) == address, 'source membership recovers full address')
        need(value not in values, 'different rooted addresses stay different')
        values.add(value)
        kernel.literal(value, route, require_root=True)
        total_letters += len(route)
        need(len(route) == 6*len(address)+2, 'uniform first-hit odd rank')
        if address:
            current = value
            for index, exponent in enumerate(route[:7]):
                current, actual = kernel.step(current)
                need(actual == exponent, 'independent first-descent replay')
                need((current > value) if index < 6 else (current < value),
                     'exact first descent is already the seventh odd edge')
        for modulus in (1, 8, 19, 729, 2187, 2048, 295*9):
            need(source_mod(address, modulus) == value % modulus,
                 'unexpanded reader agrees with literal source')
            modular_checks += 1
    for address in addresses[:13]:
        cert = kernel.codec_from_word(word(address))
        need(kernel.codec.expand(cert, bit_cap=kernel.codec.bit_bounds(cert)[1]) == source(address),
             'independent inherited ROOT codec')
    print('complete addresses: branches{0,1,3}, depths0..4:', len(addresses),
          'literal odd edges', total_letters, 'modular checks', modular_checks,
          'independent codec exports13; nonseed first descent7')

    path = []
    for depth in range(1, 9):
        address = (0,)*depth
        parent = source(address[:-1])
        path.append((depth, kappa(parent), source(address).bit_length(), ranks(address)[0]))
    print('canonical branch0 path: (depth,k,source bits,odd rank)', path)

    address_checks = derivative_checks = 0
    parents = (BASE, source((0,)), source((1,)))
    for parent in parents:
        for precision in range(1, 6):
            modulus = 3**precision
            images = [branch_residue(parent, t, modulus) for t in range(modulus)]
            need(len(set(images)) == modulus, 'all ternary source addresses realized')
            for target in range(modulus):
                branch = branch_for_address(parent, target, precision)
                need(images[branch] == target, 'digit decoder equals direct branch table')
                address_checks += 1
        for precision in range(6):
            for branch in (0, 1, 7):
                for digit in (1, 2):
                    modulus = 3**(precision+1)
                    difference = branch_residue(parent, branch+digit*3**precision, modulus)-branch_residue(parent, branch, modulus)
                    need((difference-digit*3**precision) % modulus == 0,
                         'normalized ternary derivative is1')
                    derivative_checks += 1
    print('ternary branch bijections:', address_checks, 'targets; derivative controls', derivative_checks)

    # The two generated parents have the same six-digit observation, but their
    # children do not: an observer must retain another six ternary digits.
    p0, p1 = source_mod((0,), 729), source_mod((729,), 729)
    need(p0 == p1 and kappa(p0) == kappa(p1), 'same finite phase observer')
    q0, q1 = source_mod((0, 0), 729), source_mod((729, 0), 729)
    need(q0 != q1, 'finite observer is not a closed transition system')
    print('rooted observer hostile: equal parent residues', p0,
          'and kappa', kappa(p0), 'but child residues', q0, q1)

    # Same source address is not source identity; native H availability is also
    # weaker than family membership. Neither negative is a convergence claim.
    for candidate in (1, 7, 27, 155, 219, 703):
        need(decode(candidate) is None, 'outside this constructed family')
    first = source((0,))
    need(decode(first+2048) is None and (first+2048) % 2048 == first % 2048,
         'source identity hostile within the same native H guard')
    for bad in (True, 3.0, -1, 0, 2):
        expected_failure(lambda bad=bad: decode(bad))
    for bad in ([0], (True,), (-1,), (1.0,)):
        expected_failure(lambda bad=bad: source(bad))
    print('hostiles: native H source155 and same-guard source perturbation rejected; exact types retained')

    giant = (10**100, 7, 10**80)
    residues = tuple(source_mod(giant, m) for m in (8, 19, 729, 2048))
    odd_rank, cost, ordinary_rank = ranks(giant)
    need(odd_rank == 20 and ordinary_rank == cost+20, 'symbolic first-hit rank bookkeeping')
    print('unexpanded address(10^100,7,10^80): residues mod(8,19,729,2048)', residues,
          'odd rank', odd_rank, 'halving-cost digits', len(str(cost)))
    print('PROVED: one checked seed and a decreasing-parent induction close this tree; arbitrary-source coverage is not claimed')


if __name__ == '__main__':
    main()
