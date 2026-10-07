"""Guarded routes from the symbolic child of the prior eight-bit rule.

Default controls use modular exponentiation and modest ordinary sources.
No Mersenne integer with the enormous displayed exponent is expanded.
"""
from dataclasses import replace
from fractions import Fraction

import collatz_uncovered_join_routes_20261007 as routes
import collatz_completion_anchor_20261007b as anchors
import collatz_eight_bit_completion_20261007c as prior

LEFT = (2,2,4,1,3,3,1,3,3,1,1,1,2,2,3,3,1,2,1,1,2,1,1,2,
        1,3,2,3,2,3,3,3,1,3,4,1,1,3,2)
RIGHT7 = (4,2,3,2,2,1,1,2,2,1,1,2,2,2,1,2,1,1,2,1,1,1,2,1,
          2,1,1,2,3,1,1,2,3,2,2,1,1,1,5,1,1,1,2,2,3,1)
RIGHT8 = (2,2)+RIGHT7[1:]
PHASE = 924745897
PERIOD = 1 << 79
MIN_RUN = 162
CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def integer(n, least=0):
    need(type(n) is int and n >= least, 'exact integer in declared domain')
    return n


def partner(deletion):
    need(type(deletion) is int and deletion in (7, 8), 'declared seven/eight-bit route')
    return RIGHT7 if deletion == 7 else RIGHT8


def audit_rule(deletion=8):
    right = partner(deletion)
    p, q, b = routes.carrier(LEFT)
    pp, qq, bb = routes.carrier(right)
    need(len(LEFT) == 39 and len(right) == 39+deletion, 'exact ordered head lengths')
    need((q, qq) == (1 << 81, 1 << 79), 'costs 81 and 79')
    need(pp == p*3**deletion and q == 4*qq, 'full slope relation')
    need(Fraction(bb-pp, qq) == 4*Fraction(b-p, q)+1, 'full affine carry at minus one')
    need(3*b-3*p+q == 3*bb-3*pp+qq, 'integer collision invariant')
    residue = (q-b)*pow(p, -1, 2*q) % (2*q)
    need(residue % 8 == 1, 'first reset is two')
    e, period = anchors.log_three((residue+1)//2, 81)
    need((e+1, period) == (PHASE, PERIOD), 'minimal native exponent phase')
    need(LEFT[:len(prior.RIGHT)] == prior.RIGHT and LEFT[len(prior.RIGHT)] == 3,
         'retained middle prefix and prior terminal reserve')
    return residue, q


def contains(exponent):
    integer(exponent, 1)
    return exponent-1 >= MIN_RUN and (exponent-PHASE) % PERIOD == 0


def symbolic_source(exponent, modulus, deletion=0):
    integer(modulus, 1)
    need(contains(exponent), 'same guarded exponent')
    need(type(deletion) is int and deletion in (0, 7, 8), 'declared symbolic port')
    return (pow(2, exponent-deletion, modulus)-1) % modulus


def head_endpoint_residue(exponent, modulus, deletion=0):
    """Exact quotient reader, retaining denominator precision."""
    integer(modulus, 1)
    need(contains(exponent), 'same guarded exponent')
    need(type(deletion) is int and deletion in (0, 7, 8), 'declared head endpoint')
    word = LEFT if deletion == 0 else partner(deletion)
    p, q, b = routes.carrier(word)
    x = (2*pow(3, exponent-deletion-1, q*modulus)-1) % (q*modulus)
    numerator = p*x+b
    need(numerator % q == 0, 'native exact modular division')
    return (numerator//q) % modulus


def fixed_run_source(run, parameter):
    integer(run, MIN_RUN)
    integer(parameter)
    residue, q = audit_rule()
    t0 = ((residue+1)//2)*pow(3, -run, q) % q
    need(t0 % 2 == 1, 'source-owned odd cofactor')
    return (1 << (run+1))*(t0+q*parameter)-1


def receipt(source, deletion=8, bit_cap=20000):
    routes.odd(source)
    integer(bit_cap, 1)
    need(source.bit_length() <= bit_cap, 'explicit literal materialization cap')
    right = partner(deletion)
    residue, q = audit_rule(deletion)
    run = (source+1 & -(source+1)).bit_length()-2
    need(run >= MIN_RUN, 'sufficient original-source growth cutoff')
    cofactor = (source+1) >> (run+1)
    need((2*pow(3, run, 2*q)*cofactor-1) % (2*q) == residue, 'native source guard')
    x = 2*3**run*cofactor-1
    p, q, b = routes.carrier(LEFT)
    numerator = p*x+b
    need(numerator % q == 0, 'source head integrality')
    z = numerator//q
    need(z > source and z % 2 == 1, 'growing positive odd source head')
    z3 = 3*z+1
    c = (z3 & -z3).bit_length()-1
    child = ((source+1) >> deletion)-1
    result = routes.Receipt(source, child,
        (1,)*run+LEFT+(c,), (1,)*(run-deletion)+right+(c+2,), z3 >> c)
    return routes.audit(result)


def materialize(exponent, deletion=8, bit_cap=20000):
    integer(bit_cap, 1)
    need(contains(exponent) and exponent <= bit_cap, 'phase and explicit literal cap')
    return receipt((1 << exponent)-1, deletion, bit_cap)


def modular_prefix(exponent, depth, precision):
    """Formal checkpoint reader; lack of bits returns a shorter exact prefix.

    This is not a first-hit ROOT detector or a ROOT certificate.
    """
    integer(exponent, 2)
    integer(depth)
    integer(precision, 2)
    bits = precision
    x = (2*pow(3, exponent-1, 1 << bits)-1) % (1 << bits)
    word, states = [], [(0, -2)]
    cost, invariant = 0, -2
    for _ in range(depth):
        z = (3*x+1) % (1 << bits)
        if z == 0:
            break
        a = (z & -z).bit_length()-1
        if bits-a < 1:
            break
        bits -= a
        x = (z >> a) % (1 << bits)
        word.append(a)
        cost += a
        invariant = 3*invariant+(1 << cost)
        states.append((cost, invariant))
    return tuple(word), tuple(states), bits


def bounded_discovery():
    """Exact declared comparison: D1..128 and source head lengths0..1024."""
    depth, precision = 1024, 8192
    source, states, _ = modular_prefix(PHASE, depth+129, precision)
    need(len(source) == depth+129, 'source precision covers comparison universe')
    hits = []
    for deletion in range(1, 129):
        word, other, _ = modular_prefix(PHASE-deletion, depth+129, precision)
        need(len(word) == depth+129, 'partner precision covers comparison universe')
        for length in range(depth+1):
            if states[length][1] != other[length+deletion][1]:
                continue
            delta = states[length][0]-other[length+deletion][0]
            if delta % 2 == 0 and word[length+deletion] == source[length]+delta:
                hits.append((length, deletion, delta//2))
                break
    need(len(hits) == 24 and min(hits) == (39, 7, 1), 'frozen finite comparison result')
    need((39, 8, 1) in hits, 'same-depth stronger alternative')
    need(source[:39] == LEFT, 'source-authentic discovered head')
    for d in (7, 8):
        word, _, _ = modular_prefix(PHASE-d, 40+d, 256)
        need(word[:39+d] == partner(d), 'source-authentic partner head')
    return tuple(hits)


def reject(function, *args):
    try:
        function(*args)
    except (TypeError, ValueError):
        need(True, 'hostile rejected')
        return
    need(False, 'hostile accepted')


def main():
    for d in (7, 8):
        audit_rule(d)
    hits = bounded_discovery()
    literal = 0
    for run in (162, 163, 200, 256):
        for parameter in (0, 1, 7, 31):
            n = fixed_run_source(run, parameter)
            r7, r8 = receipt(n, 7), receipt(n, 8)
            need(r7.source_word == r8.source_word and r7.endpoint == r8.endpoint,
                 'two authenticated alternatives at the same source and future')
            need(r8.child < r7.child < n, 'strict original-source payment')
            need(len(r8.source_word) == len(r8.child_word), 'equal full odd clocks')
            need(sum(r8.source_word)-sum(r8.child_word) == 8, 'eight-bit actual cost difference')
            p = prior.receipt(((n+1) << 8)-1)
            combined = prior.compose(p, r8)
            need(combined.child == n//256 and combined.source == ((n+1) << 8)-1,
                 'literal parent/child receipt composition keeps exact middle')
            literal += 2
    for t in (0, 1, 2, 9, 10**30):
        E = PHASE+PERIOD*t
        need(contains(E) and prior.contains(E+8), 'nested source phase')
        for m in (2, 3, 19, 729, 1 << 128, 3**20*19):
            z = head_endpoint_residue(E, m)
            for d in (7, 8):
                need(head_endpoint_residue(E, m, d) == (4*z+1) % m,
                     'independent symbolic common-head endpoint')
                need(symbolic_source(E, m, d) == (pow(2, E-d, m)-1) % m,
                     'lossless symbolic source identity')
        need((E-9 & -(E-9)).bit_length()-1 == 5, 'D8 child returns to J3')
    from collatz_general_head_phases_20261007b import finite_bank, phase
    bank = finite_bank()
    need(all((PHASE-8-phase(h, 3).residue) % phase(h, 3).period != 0 for h in bank),
         'new fixed J3 terminal misses all225 existing heads')
    for bad in (True, 924745897.0, -1):
        reject(contains, bad)
    reject(receipt, True)
    reject(receipt, fixed_run_source(162, 0), True)
    reject(materialize, PHASE)
    reject(receipt, fixed_run_source(162, 0)-2)
    r = receipt(fixed_run_source(162, 0))
    reject(routes.audit, replace(r, child=r.child+2))
    need(not contains(PHASE+16), 'same two-twos type does not replace native phase')
    short, _, _ = modular_prefix(PHASE, 100, 16)
    need(len(short) < 100, 'insufficient precision never invents a continuation')
    print('PROVED child phase F=924745897mod2^79, paid deletions7 and8.')
    print('Heads: source39/cost81; child46or47/cost79; positive sibling gap1.')
    print('FINITE-EXACT search D1..128,source-depth0..1024,precision8192:24 matches; minimum depth39.')
    print('FINITE-EXACT modest literal receipts:', literal, '; exact parent composition controls:', literal//2)
    print('PROVED symbolic phase+denominator readers at five parameters; no enormous source expansion.')
    print('OPEN new D8 terminal: exponent924745889mod2^79, J3, no225-bank head.')
    print('No ROOT oracle, no universal coverage, no global minimality claim.')
    print('Checks:', CHECKS+routes.CHECKS+anchors.CHECKS+prior.CHECKS)


if __name__ == '__main__':
    main()
