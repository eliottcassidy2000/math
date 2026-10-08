"""Marked normalized Gaussian norms transport actual source guards exactly.

No ROOT discovery and no replacement of a supplied source by its norm.
"""
from dataclasses import dataclass, replace
from math import isqrt

import collatz_bounded_inverse_cover_20261007e as bank

CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def integer(n, least=0):
    if type(n) is not int or n < least:
        raise ValueError('exact integer in declared domain required')
    return n


def odd(n):
    integer(n, 1)
    if n % 2 != 1:
        raise ValueError('positive odd source required')
    return n


def norm_coordinate(n):
    odd(n)
    return (n*n+1)//2


def source_from_norm(q):
    integer(q, 1)
    n = isqrt(2*q-1)
    if n*n != 2*q-1 or n % 2 != 1:
        raise ValueError('not the exact norm coordinate of a positive odd source')
    return n


def decode_two(residue, bits, branch):
    integer(bits, 2); integer(residue); integer(branch)
    modulus = 1 << bits
    if branch not in (1, 3) or residue >= modulus or residue % 4 != 1:
        raise ValueError('dyadic norm residue and branch required')
    n = branch
    for j in range(2, bits):
        if ((n*n+1)//2-residue) % (1 << (j+1)):
            n += 1 << j
    return n


def decode_three(residue, depth, branch):
    integer(depth, 1); integer(residue); integer(branch)
    modulus = 3**depth
    if branch not in (1, 2) or residue >= modulus or residue % 3 != 1:
        raise ValueError('ternary unit norm residue and branch required')
    n, power = branch, 3
    for _ in range(1, depth):
        target_modulus = 3*power
        candidates = tuple(n+digit*power for digit in range(3)
                           if ((n+digit*power)**2+1-2*residue) % target_modulus == 0)
        if len(candidates) != 1:
            raise ArithmeticError('unit derivative failed')
        n = candidates[0]
        power = target_modulus
    return n


def decode_mixed(residue, bits, depth, branch_two, branch_three):
    integer(bits, 2); integer(depth, 1); integer(residue)
    a, b = 1 << bits, 3**depth
    if residue >= a*b:
        raise ValueError('canonical mixed residue required')
    x = decode_two(residue % a, bits, branch_two)
    y = decode_three(residue % b, depth, branch_three)
    return x+a*((y-x)*pow(a, -1, b) % b)


@dataclass(frozen=True)
class NormGuard:
    source_residue: int
    bits: int
    depth: int
    branch_two: int
    branch_three: int
    norm_residue: int


def compile_guard(residue, bits, depth):
    integer(bits, 2); integer(depth, 1); integer(residue)
    modulus = (1 << bits)*3**depth
    if residue >= modulus or residue % 2 == 0 or residue % 3 == 0:
        raise ValueError('canonical odd 3-unit source residue required')
    return NormGuard(residue, bits, depth, residue % 4, residue % 3,
                     ((residue*residue+1)//2) % modulus)


def audit_guard(guard):
    if type(guard) is not NormGuard:
        raise ValueError('typed NormGuard required')
    for value in (guard.source_residue, guard.bits, guard.depth,
                  guard.branch_two, guard.branch_three, guard.norm_residue):
        integer(value)
    if guard != compile_guard(guard.source_residue, guard.bits, guard.depth):
        raise ValueError('forged local marker or norm guard')
    return guard


def accepts(guard, source):
    audit_guard(guard); odd(source)
    modulus = (1 << guard.bits)*3**guard.depth
    return (source % 4 == guard.branch_two and source % 3 == guard.branch_three
            and norm_coordinate(source) % modulus == guard.norm_residue)


def native_norm_guard(entry, branch_two):
    bank.audit(entry); integer(branch_two)
    if branch_two not in (1, 3):
        raise ValueError('marked odd dyadic branch required')
    p = entry.P
    r = entry.least_source % p
    r += p*((branch_two-r)*pow(p, -1, 4) % 4)
    return compile_guard(r, 2, len(entry.word))


def norm_receipt(entry, source):
    """The ordinary source is immutable; its norm is only a guard coordinate."""
    bank.audit(entry); odd(source)
    guard = native_norm_guard(entry, source % 4)
    if not accepts(guard, source):
        raise ValueError('supplied source fails the marked norm guard')
    return bank.receipt(entry, source)


def gaussian_image(parameter, x, y):
    integer(parameter)
    if type(x) is not int or type(y) is not int:
        raise ValueError('exact integer coordinates required')
    return parameter*x-y, x+parameter*y


def gaussian_coset(parameter, x, y):
    integer(parameter)
    if type(x) is not int or type(y) is not int:
        raise ValueError('exact integer coordinates required')
    return (x-parameter*y) % (parameter*parameter+1)


def odd_part(n):
    integer(n, 1)
    return n//(n & -n)


def main():
    sequence = tuple(n*n+1 for n in range(12))
    need(sequence == (1, 2, 5, 10, 17, 26, 37, 50, 65, 82, 101, 122),
         'owner sequence with n starting at zero')
    for a in range(20):
        determinant = a*a+1
        need((a+1)**2+1-determinant == 2*a+1, 'successive square shell')
        need(determinant % 3 != 0, 'Gaussian lattice index is a ternary unit')
        need((determinant & -determinant) == (2 if a % 2 else 1),
             'exact dyadic index valuation')
        for x in range(-5, 6):
            for y in range(-5, 6):
                u, v = gaussian_image(a, x, y)
                need(u*u+v*v == determinant*(x*x+y*y), 'Gaussian norm similarity')
                need(gaussian_coset(a, u, v) == 0, 'image lies in exact lattice kernel')
                divisible = (a*x+y) % determinant == 0 and (-x+a*y) % determinant == 0
                need((gaussian_coset(a, x, y) == 0) == divisible,
                     'cyclic quotient detects integral inverse')
        for bits in range(1, 5):
            m = 1 << bits
            kernel = sum(gaussian_image(a, x, y)[0] % m == 0 and
                         gaussian_image(a, x, y)[1] % m == 0
                         for x in range(m) for y in range(m))
            need(kernel == (2 if a % 2 else 1), 'one dyadic bit in full lattice chart')
    for bits in range(2, 10):
        modulus = 1 << bits
        for n in range(1, modulus, 2):
            q = norm_coordinate(n) % modulus
            need(decode_two(q, bits, n % 4) == n, 'dyadic exact branch decoder')
            need(decode_two(q, bits, (-n) % 4) == (-n) % modulus,
                 'other branch is exact negative residue')
    for depth in range(1, 7):
        modulus = 3**depth
        for n in range(modulus):
            if n % 3 == 0:
                continue
            q = (n*n+1)*pow(2, -1, modulus) % modulus
            need(decode_three(q, depth, n % 3) == n, 'ternary unit branch decoder')
    for bits in range(2, 7):
        for depth in range(1, 4):
            modulus = (1 << bits)*3**depth
            for n in range(1, modulus, 2):
                if n % 3 == 0:
                    continue
                guard = compile_guard(n, bits, depth)
                need(decode_mixed(guard.norm_residue, bits, depth, n % 4, n % 3) == n,
                     'joint CRT guard decoder')
                need(accepts(guard, n+2*modulus), 'same source guard at an integral lift')
                need(not accepts(guard, n+2), 'nearby different source cannot replace it')
    for n in range(1, 1000, 2):
        need(source_from_norm(norm_coordinate(n)) == n, 'full positive norm loses no source')
    need(norm_coordinate(3) % 8 == norm_coordinate(5) % 8 == 5,
         'norm-only finite observer collision')
    need(bank.routes.step(3)[1] == 1 and bank.routes.step(5)[1] == 4,
         'collision loses actual valuation')
    need((3*3+1) % 8 == (1*1+1) % 8 and norm_coordinate(3) % 8 != norm_coordinate(1) % 8,
         'raw norm needs one extra dyadic bit before division by2')
    need((3*3+1)*pow(2, -1, 9) % 9 == (6*6+1)*pow(2, -1, 9) % 9,
         'ternary nonunit branch is ramified and excluded from unit decoder')
    for a in range(31):
        for b in range(31):
            if a == b == 0:
                continue
            source = odd_part(a*a+b*b)
            need(source % 4 == 1, 'every odd part of a Gaussian norm is1mod4')
            if source > 1:
                need(bank.routes.step(source)[0] < source, 'known immediate descent at the norm source')
    for e in bank.entries():
        for j in range(3):
            source = bank.ordinary_source(e, j)
            need(norm_receipt(e, source) == bank.receipt(e, source),
                 'all bounded inverse labels retain the exact ordinary source and child')
    for exponent in range(3, 30, 2):
        n = (1 << exponent)-1
        q = norm_coordinate(n)
        need(n % 4 == 3 and n % 3 == 1, 'fixed Mersenne local markers')
        need((q-1) & -(q-1) == 1 << exponent, 'norm coordinate retains the exponent depth')
        need(source_from_norm(q) == n, 'marked norm does not advance the orbit')
    g = compile_guard(11, 3, 1)
    invalid = (
        lambda: norm_coordinate(True), lambda: norm_coordinate(2),
        lambda: source_from_norm(2), lambda: decode_two(5, 3, True),
        lambda: decode_two(3, 3, 1), lambda: decode_three(2, 2, 1),
        lambda: compile_guard(3, 3, 2),
        lambda: audit_guard(replace(g, branch_three=True)),
        lambda: accepts(g, 11.0), lambda: norm_receipt(bank.compile_word((1, 2)), 3),
    )
    for call in invalid:
        try:
            call()
        except (ValueError, TypeError):
            pass
        else:
            raise ValueError('malformed norm/guard accepted')
        need(True, 'typed or non-native hostile rejected')
    print('Square-plus-one sequence:', sequence)
    print('Gaussian lattice quotient: Z/(a^2+1); image test x-a*y=0 modulo index')
    print('Marked normalized odd norm: same-precision dyadic and ternary-unit guard bijections')
    print('Finite-local hostile:3 and5 have q=5 mod8, first valuations1 and4')
    print('Full positive norm is injective; no global information-loss claim')
    print('All252 bounded inverse labels:756 identical source-owned receipts through norm guards')
    print('Odd Gaussian norms already lie1mod4; this gives no descent of a different original source')
    print('No new parameter coverage or ROOT discovery')
    print('Exact checks:', CHECKS)


if __name__ == '__main__':
    main()
