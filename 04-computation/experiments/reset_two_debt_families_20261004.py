"""Guarded reset-two common-future families, with exact source certificates.

This verifies identities and declared finite controls. It neither assumes a
certificate for an arbitrary smaller source nor claims universal coverage.
Run with ordinary Python or python -O; all checks use explicit exceptions.
"""
from dataclasses import dataclass
from fractions import Fraction

import inverse_ray_ternary_addresses_20261004 as codec


def need(ok, message):
    if not ok:
        raise ValueError(message)


def positive_integer(value):
    return type(value) is int and value > 0


def v2(value):
    need(positive_integer(value), 'positive exact integer valuation')
    return (value & -value).bit_length()-1


def metadata(word):
    word = tuple(word)
    need(bool(word) and all(positive_integer(a) for a in word), 'positive exponent word')
    P, Q, B = 1, 1, 0
    for exponent in word:
        P, Q, B = 3*P, Q*2**exponent, 3*B+Q
    return P, Q, B


def step(source, sign=1):
    need(positive_integer(source) and source % 2, 'positive odd source')
    need(type(sign) is int and sign in (-1, 1), 'retained signed sheet')
    numerator = 3*source+sign
    exponent = v2(numerator)
    return numerator >> exponent, exponent


def replay(source, word, sign=1, first_hit=False):
    metadata(word)
    states = [source]
    for exponent in word:
        if first_hit:
            need(states[-1] != 1, 'root reached before the displayed word ended')
        target, actual = step(states[-1], sign)
        need(actual == exponent, 'actual valuation guard')
        states.append(target)
    return tuple(states)


def source_cylinder(word, sign=1):
    """Exact word sources are residue mod 2Q, not merely the coarse mod Q cell."""
    need(type(sign) is int and sign in (-1, 1), 'retained signed sheet')
    P, Q, B = metadata(word)
    modulus = 2*Q
    residue = (Q-sign*B)*pow(P, -1, modulus) % modulus
    need(residue > 0 and residue % 2, 'positive odd exact-word residue')
    return residue, modulus


@dataclass(frozen=True)
class Portal:
    debt_exponent: int
    word: tuple

    def validate(self):
        need(positive_integer(self.debt_exponent), 'positive debt exponent')
        P, Q, B = metadata(self.word)
        need(len(self.word) == self.debt_exponent+2 and B == Q+7,
             'portal carry and length guard')
        return P, Q, B


PORTALS = (Portal(1, (1, 2, 1)), Portal(3, (1, 6, 1, 3, 1)),
           Portal(4, (2, 2, 1, 3, 3, 1)))


def rewrite_words(run, terminal, portal=PORTALS[0]):
    need(positive_integer(run) and positive_integer(terminal), 'positive run and terminal')
    portal.validate()
    e, S = portal.debt_exponent, sum(portal.word)
    source_word = (1,)*run+(2,)*e+(1, S+2, terminal)
    child_word = (1,)*(run-1)+(2*e+1,)+portal.word+(terminal+2,)
    P, Q, B = metadata(source_word)
    p, q, b = metadata(child_word)
    need(P == p and Q == 2*q and 2*b-B == P, 'exact half-source carry identity')
    return source_word, child_word


def check_rewrite(source, run, terminal, sign=1, portal=PORTALS[0]):
    source_word, child_word = rewrite_words(run, terminal, portal)
    need(positive_integer(source) and source > 1 and source % 2, 'nonroot odd source')
    child = (source-sign)//2
    source_states = replay(source, source_word, sign)
    child_states = replay(child, child_word, sign)
    need(0 < child < source and child % 2, 'strict smaller-source rank')
    need(source_states[-1] == child_states[-1], 'common future')
    return child, source_states, child_states


def applicable_rewrite(source):
    """Exact e1 selector: finite source arithmetic, no root or suffix search."""
    need(positive_integer(source) and source % 2, 'positive odd selector source')
    if source == 1 or source % 4 != 3:
        return None
    run = v2(source+1)-1
    prefix = (1,)*run+(2, 1, 6)
    P, Q, B = metadata(prefix)
    if (P*source+B) % (2*Q) != Q:
        return None
    W = (P*source+B)//Q
    terminal = v2(3*W+1)
    source_word, child_word = rewrite_words(run, terminal)
    residue, modulus = source_cylinder(source_word)
    need(source % modulus == residue, 'exact complete selector cylinder')
    return run, terminal, (source-1)//2, (3*W+1)//2**terminal, source_word, child_word


def transport_certificate(source, run, terminal, child_certificate, portal=PORTALS[0]):
    """Consume an existing first-hit child AST; never search for its suffix."""
    child, source_states, child_states = check_rewrite(source, run, terminal, 1, portal)
    codec.audit_certificate(child_certificate)
    need(codec.expand(child_certificate) == child, 'certificate belongs to the actual smaller source')
    source_word, child_word = rewrite_words(run, terminal, portal)
    suffix = child_certificate
    for exponent in child_word:
        need(suffix != codec.ROOT and codec.exponent(suffix) == exponent,
             'supplied certificate contains the checked child route')
        suffix = suffix.parent
    need(codec.expand(suffix) == child_states[-1], 'retained common endpoint suffix')
    need(1 not in source_states[:-1], 'canonical first-hit source prefix')
    for index in range(len(source_word)-1, -1, -1):
        actual_source, target = source_states[index:index+2]
        row = actual_source % 3
        least = codec.kappa(target % 9, row)
        exponent = source_word[index]
        need(exponent >= least and (exponent-least) % 6 == 0, 'inverse row/block guard')
        suffix = codec.extend(suffix, row, (exponent-least)//6)
    need(codec.expand(suffix) == source, 'transported original source identity')
    return suffix


def normalized_reset_two(source):
    need(positive_integer(source) and source > 1 and source % 4 == 3, 'initial run source')
    run = v2(source+1)-1
    t = (source+1)//2**(run+1)
    b = 1+v2(3**run*t-1)
    need(b >= 3, 'actual first reset must be two')
    M = (3**run*t-1)//2**(b-1)
    child = (source-1)//2
    need(replay(child, (1,)*(run-1)+(b,))[-1] == M, 'retained child-to-hub word')
    e, h = 1, b-2
    Y = replay(source, (1,)*run+(2,))[-1]
    need(Y == 3**e*2**h*M+1, 'initial debt identity')
    while h >= 3:
        Y, exponent = step(Y)
        e, h = e+1, h-2
        need(exponent == 2 and Y == 3**e*2**h*M+1, 'forced debt transition')
    next_value, exponent = step(Y)
    if h == 1:
        need(exponent == 1 and next_value == 3**(e+1)*M+2, 'odd boundary type')
    else:
        raw = 3**(e+1)*M+1
        need(next_value == raw >> v2(raw), 'even boundary type')
    return run, b, M, e, h, Y, next_value


def rational_odd_orbit(source, length):
    source = Fraction(source)
    need(source.denominator % 2 and source.numerator % 2, 'odd rational 2-adic source')
    states, word = [source], []
    for _ in range(length):
        numerator = 3*states[-1]+1
        need(numerator != 0, 'nonzero rational edge')
        exponent = v2(abs(numerator.numerator))
        need(exponent > 0, 'positive rational odd-step exponent')
        word.append(exponent)
        states.append(numerator/2**exponent)
    return tuple(word), tuple(states)


def main():
    print('RESET TWO DEBT FAMILIES: PROVED guarded rewrites; FINITE-EXACT controls; universal closure OPEN')
    print('Portal rows: e, prefix, sum, carry; condition B=2^sum+7')
    for portal in PORTALS:
        _, Q, B = portal.validate()
        print(portal.debt_exponent, portal.word, sum(portal.word), B)

    count = 0
    for portal in PORTALS:
        for run in range(1, 25):
            for terminal in range(1, 9):
                source_word, _ = rewrite_words(run, terminal, portal)
                for sign in (-1, 1):
                    residue, modulus = source_cylinder(source_word, sign)
                    for lift in range(4):
                        n = residue+lift*modulus
                        check_rewrite(n, run, terminal, sign, portal)
                        count += 1
    print('Exact signed word rewrites:', count, '(3 portals; r1..24; a1..8; signs +/-; four lifts)')

    portal_controls = 0
    for portal in PORTALS:
        e, S = portal.debt_exponent, sum(portal.word)
        P = 3**(e+2)
        modulus = 2**(S+3)
        residue = (2**(S+2)-7)*pow(P, -1, modulus) % modulus
        for lift in range(32):
            M = residue+modulus*lift
            Z = 3**(e+1)*M+2
            W, exponent = step(Z)
            need(exponent == S+2, 'exact boundary guard')
            need(replay(M, portal.word)[-1] == 4*W+1, 'boundary reaches a sibling')
            need(step(W)[0] == step(4*W+1)[0], 'sibling common future')
            portal_controls += 1
    print('Boundary-to-sibling controls:', portal_controls, '; first portal guard M=59 mod128')

    growth_controls = 0
    for run in range(8, 129):
        source_word, _ = rewrite_words(run, 1)
        for cut in range(1, len(source_word)+1):
            P, Q, B = metadata(source_word[:cut])
            need(P > Q and B > 0, 'all-height strict source growth at every displayed prefix')
        residue, modulus = source_cylinder(source_word)
        for lift in range(3):
            _, states, _ = check_rewrite(residue+modulus*lift, run, 1)
            need(all(x > states[0] for x in states[1:]), 'literal growing source prefix')
            growth_controls += 1
    print('Growing-family controls:', growth_controls, '(r8..128, three lifts); all prefix multipliers >1')
    word, _ = rewrite_words(8, 1)
    r8, mod8 = source_cylinder(word)
    child, source_states, child_states = check_rewrite(r8, 8, 1)
    print('First r8 cylinder:', r8, 'mod', mod8, '; child', child, '; common endpoint', source_states[-1])
    print('Source r8 states:', source_states)
    print('Child r8 states:', child_states)
    for lift in range(8):
        _, lifted, _ = check_rewrite(r8+mod8*lift, 8, 1)
        need(lifted[-1] == 478505+1062882*lift, 'retained odd-endpoint affine increment')
    need(Fraction(3**11, 2**17) > 1 and Fraction(3**10, 2**16) < 1, 'sharp r8 coefficient boundary')
    w7, _ = rewrite_words(7, 1)
    n7, _ = source_cylinder(w7)
    _, states7, _ = check_rewrite(n7, 7, 1)
    need(states7[-2] < n7, 'r7 all-growth hostile')
    print('r7 hostile:', n7, 'falls to', states7[-2], 'before the join')
    for upper in range(8, 65):
        density = sum((Fraction(1, 2**(r+11)) for r in range(8, upper+1)), Fraction())
        need(density == Fraction(1, 2**18)*(1-Fraction(1, 2**(upper-7))), 'partial cylinder density')
    print('Disjoint growing-family natural density: 1/262144 of all integers; 1/131072 among odds')

    for n in (7, 27, 315):
        print('Retained reset2 debt n', n, ': (r,b,M,e,h,Y,next)=', normalized_reset_two(n))
    need(normalized_reset_two(27)[2] % 128 != 59, '27 does not satisfy the new portal')
    child315, states315, alt315 = check_rewrite(315, 1, 2)
    need(states315[-2] == 25 < 315, '315 is not a growing-prefix efficacy example')
    print('Illustrative 315 join:', states315, '<-', alt315, '; smaller source', child315)

    print('Rational-anchor ansatz controls: length, exact word, endpoint')
    for length in range(3, 7):
        word, states = rational_odd_orbit(Fraction(-7, 3**length), length)
        print(length, word, str(states[-1]))
    word4, states4 = rational_odd_orbit(Fraction(-7, 81), 4)
    need(word4 == (2, 1, 1, 1) and states4[-1] == 3, 'all-height e2 two-step ansatz obstruction')
    need(3 != Fraction(4**1-1, 3) and 3 < Fraction(4**2-1, 3), '3 misses every positive c_s')

    for cycle, exponents in (((1,), (1,)), ((5, 7), (1, 2)),
                              ((17, 25, 37, 55, 41, 61, 91), (1, 1, 1, 2, 1, 1, 4))):
        for i, n in enumerate(cycle):
            need(step(n, -1) == (cycle[(i+1) % len(cycle)], exponents[i]), 'retained minus basin')
        need(6 not in exponents, 'signed cycle cannot trigger e1 portal')
    print('Signed hostile controls: separate minus cycles 1, 5/7, 17/.../91 retained')

    transported = 0
    for run in range(8, 21):
        source_word, _ = rewrite_words(run, 1)
        residue, modulus = source_cylinder(source_word)
        for lift in range(2):
            n = residue+lift*modulus
            child = (n-1)//2
            # This bounded search supplies the premise; transport itself never calls it.
            supplied = codec.encode_source(child, step_cap=5000)
            result = transport_certificate(n, run, 1, supplied)
            need(result == codec.encode_source(n, step_cap=5000), 'post-construction canonical audit')
            codec.literal_certificate_check(result)
            j, ordinary = codec.ranks(supplied)
            need(codec.ranks(result) == (j, ordinary+1), 'ordinary rank rises by one; odd rank stays fixed')
            need(applicable_rewrite(n)[:3] == (run, 1, child), 'growing-cylinder selector applicability')
            transported += 1
    print('Supplied smaller-source AST transports:', transported, '(r8..20, two lifts; premise search explicitly counted)')
    selected = 0
    for n in range(1, 20001, 2):
        result = applicable_rewrite(n)
        current, ones = n, 0
        if n == 1:
            literal = False
        else:
            while True:
                current, exponent = step(current)
                if exponent != 1:
                    break
                ones += 1
            literal = False
            if ones and exponent == 2:
                current, second = step(current)
                current, third = step(current)
                literal = (second, third) == (1, 6)
        need((result is not None) == literal, 'source-only selector versus independent actual-prefix check')
        if result is not None:
            check_rewrite(n, result[0], result[1])
            selected += 1
    print('Selector census: 10000 odd inputs through20000;', selected, 'exact matches; no supplied-home inference')
    for upper in range(1, 65):
        reset_density = sum((Fraction(1, 2**(r+3)) for r in range(1, upper+1)), Fraction())
        need(reset_density == Fraction(1, 8)*(1-Fraction(1, 2**upper)), 'least-counterexample reset2 kernel density')
        portal_density = sum((Fraction(1, 2**(r+10)) for r in range(1, upper+1)), Fraction())
        need(portal_density == Fraction(1, 1024)*(1-Fraction(1, 2**upper)), 'all-terminal portal density')
    print('Necessary least-counterexample kernel: first reset2; density1/8 of integers,1/4 among odds')
    print('All-terminal e1 family density:1/1024 of integers,1/512 among odds; remaining necessary kernel127/512 of odds')
    root_source = (2**168-1279)//243
    need((2**168-1279) % 243 == 0, 'root-to-source ternary guard')
    need(applicable_rewrite(root_source)[1:4] == (158, (root_source-1)//2, 1), 'legitimate first-hit endpoint1')
    terminal_root = transport_certificate(root_source, 1, 158, codec.encode_source((root_source-1)//2))
    need(codec.ranks(terminal_root)[0] == 5, 'no padded root edge in a transported terminal certificate')
    need(applicable_rewrite(68579) is None, 'integral hub alone does not supply the root-to-hub ternary guard')
    print('First-hit root control:(2^168-1279)/243 ->1 with word(1,2,1,6,158); no root padding')
    print('Retained failed-hub lift hostile:68579 does not meet the source cylinder')
    rejection_count = 0
    for thunk in (lambda: rewrite_words(True, 1), lambda: rewrite_words(8, 1.0),
                  lambda: check_rewrite(r8+2, 8, 1),
                  lambda: transport_certificate(r8, 8, 1, codec.ROOT),
                  lambda: source_cylinder((1, 2), True),
                  lambda: Portal(2, (2, 1, 1, 1)).validate()):
        try:
            thunk()
        except ValueError:
            rejection_count += 1
        else:
            raise ValueError('hostile unexpectedly accepted')
    print('Rejected malformed/rank/source/ansatz controls:', rejection_count)
    print('PASS: no unconditional source coverage, no terminal-basin identification, no reset2 debt erasure')


if __name__ == '__main__':
    main()
