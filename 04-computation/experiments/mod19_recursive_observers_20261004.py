"""Exact mod-19 observers and source-certified inverse-address selection.

No finite residue test proves a Collatz route. The certificate constructor
uses the inherited guarded inverse codec; all checks remain active under -O.
"""
from collections import Counter
from functools import lru_cache
from itertools import product
from math import gcd, lcm
from pathlib import Path
import importlib.util
import json
import sys


def need(ok, message):
    if not ok:
        raise ValueError(message)


def vp(n, p):
    need(type(n) is int and n != 0, 'nonzero exact integer valuation')
    n, v = abs(n), 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def precision(depth):
    need(type(depth) is int and depth >= 1, 'positive integer precision')


def word_data(word):
    need(bool(word) and all(type(a) is int and a > 0 for a in word), 'positive valuation word')
    p, q, b = 1, 1, 0
    for a in word:
        p, b, q = 3*p, 3*b+q, q*2**a
    return p, q, b


def affine_mod(word, depth):
    precision(depth)
    p, q, b = word_data(word)
    modulus = 19**depth
    return p*pow(q, -1, modulus) % modulus, b*pow(q, -1, modulus) % modulus


def crt_pair(a, m, b, n):
    """Return the exact common arithmetic progression, or None."""
    need(type(m) is int and type(n) is int and m > 0 and n > 0, 'positive CRT moduli')
    g = gcd(m, n)
    if (b-a) % g:
        return None
    reduced = n//g
    t = 0 if reduced == 1 else ((b-a)//g)*pow(m//g, -1, reduced) % reduced
    return (a+m*t) % lcm(m, n), lcm(m, n)


def endpoint_source_filter(word, endpoint_residue, depth):
    """One exact dyadic word guard refined by its endpoint's 19-adic address."""
    precision(depth)
    p, q, b = word_data(word)
    modulus = 19**depth
    need(type(endpoint_residue) is int and 0 <= endpoint_residue < modulus, 'canonical endpoint address')
    binary = (q-b)*pow(p, -1, 2*q) % (2*q)
    source = (q*endpoint_residue-b)*pow(p, -1, modulus) % modulus
    return crt_pair(binary, 2*q, source, modulus)


def chart112_observer(n, depth):
    """A formal phase annotation, with the one-word integer guard separate."""
    need(type(n) is int, 'exact integer source')
    precision(depth)
    modulus = 19**depth
    coordinate = (11*n+19) % modulus
    shell = depth if coordinate == 0 else vp(coordinate, 19)
    period = 1 if shell == depth else 18*19**(depth-shell-1)
    return dict(coordinate=coordinate, shell=shell, formal_period=period,
                one_word_guard=n > 0 and n % 32 == 7)


def fixed_points(word, depth):
    """Compressed solutions of (Q-P)x=B mod19^depth."""
    precision(depth)
    p, q, b = word_data(word)
    modulus = 19**depth
    g = gcd(q-p, modulus)
    if b % g:
        return None
    reduced = modulus//g
    root = 0 if reduced == 1 else (b//g)*pow((q-p)//g, -1, reduced) % reduced
    return root, reduced, g


def kappa(hub_mod9, row):
    need(type(hub_mod9) is int and hub_mod9 % 3 != 0, 'ternary-unit hub register')
    need(type(row) is int and row in (0, 1, 2), 'source row')
    candidates = [a for a in range(1, 7) if pow(2, a, 9)*hub_mod9 % 9 == 1+3*row]
    need(len(candidates) == 1, 'unique six-block exponent')
    return candidates[0]


def ray19(hub_mod9, hub_mod19, row, block, depth):
    precision(depth)
    need(type(block) is int and block >= 0, 'nonnegative exact block')
    need(type(hub_mod19) is int, 'exact hub register')
    modulus = 19**depth
    exponent = kappa(hub_mod9, row)+6*block
    return (pow(2, exponent, modulus)*hub_mod19-1)*pow(3, -1, modulus) % modulus


def select_ray19(hub_mod9, hub_mod19, row, wanted, depth):
    """All block indices for one requested address; three base tests, then digits.

    Registers need not expand the hub. Return (block, period), or None.
    The value hub_mod19=0 mod19^depth is deliberately a collapsed clock.
    """
    precision(depth)
    need(type(wanted) is int and 0 <= wanted < 19**depth, 'canonical wanted address')
    kappa(hub_mod9, row)
    need(type(hub_mod19) is int, 'exact hub register')
    modulus = 19**depth
    hub = hub_mod19 % modulus
    if hub == 0:
        return (0, 1) if wanted == -pow(3, -1, modulus) % modulus else None
    valuation = vp(hub, 19)
    base_modulus = 19**(valuation+1)
    candidates = [b for b in range(3)
                  if ray19(hub_mod9, hub, row, b, valuation+1) == wanted % base_modulus]
    need(len(candidates) <= 1, 'three distinct normalized base addresses')
    if not candidates:
        return None
    block = candidates[0]
    for a in range(valuation+1, depth):
        current = ray19(hub_mod9, hub, row, block, a+1)
        need((wanted-current) % 19**a == 0, 'retained lower source digits')
        coefficient = ((3*current+1)//19**valuation) % 19
        need(coefficient != 0, 'normalized carry is a unit')
        digit = ((wanted-current)//19**a)*pow(coefficient, -1, 19) % 19
        block += digit*3*19**(a-valuation-1)
    period = 3*19**(depth-valuation-1)
    need(0 <= block < period and ray19(hub_mod9, hub, row, block, depth) == wanted,
         'canonical selected block and exact address')
    return block, period


def select_ray3(hub_register, wanted, depth):
    precision(depth)
    need(type(hub_register) is int and hub_register % 3 != 0, 'ternary-unit register')
    need(type(wanted) is int and 0 <= wanted < 3**depth, 'canonical ternary address')
    row, block = wanted % 3, 0
    least = kappa(hub_register % 9, row)
    for a in range(1, depth):
        modulus = 3**(a+2)
        numerator = (pow(2, least+6*block, modulus)*hub_register-1) % modulus
        need(numerator % 3 == 0, 'ternary division precision')
        current = numerator//3
        need((wanted-current) % 3**a == 0, 'retained ternary address')
        block += ((wanted-current)//3**a % 3)*3**(a-1)
    return row, block, 3**(depth-1)


@lru_cache(None)
def ray_codec():
    # One canonical module identity is essential: certificates passed between
    # sibling scripts must share their dataclass and ROOT, not just field names.
    name = 'inverse_ray_ternary_addresses_20261004'
    if name in sys.modules:
        return sys.modules[name]
    path = Path(__file__).with_name('inverse_ray_ternary_addresses_20261004.py')
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def certificate_mod19(cert, depth):
    """Validate structure, then evaluate without expanding any represented source."""
    precision(depth)
    codec = ray_codec()
    codec.audit_certificate(cert)
    result = 1
    modulus = 19**depth
    for node in reversed(codec.chain(cert)):
        result = (pow(2, codec.exponent(node), modulus)*result-1)*pow(3, -1, modulus) % modulus
    return result


def select_certified_predecessors(parent, wanted3, depth3, wanted19, depth19):
    """Return row, block progression and isolated first-hit exclusion, or None."""
    precision(depth3)
    precision(depth19)
    codec = ray_codec()
    codec.audit_certificate(parent)
    u3 = codec.mod3(parent, depth3+1)
    row, b3, m3 = select_ray3(u3, wanted3, depth3)
    u19 = certificate_mod19(parent, depth19)
    selected = select_ray19(u3 % 9, u19, row, wanted19, depth19)
    if selected is None:
        return None
    joint = crt_pair(b3, m3, *selected)
    if joint is None:
        return None
    block, period = joint
    first = int(parent == codec.ROOT and row == 1 and block == 0)
    return dict(row=row, block=block, period=period, first_parameter=first)


def instantiate_selected(parent, selection, parameter):
    need(type(parameter) is int and parameter >= selection['first_parameter'], 'retained first-hit parameter guard')
    return ray_codec().extend(parent, selection['row'], selection['block']+selection['period']*parameter)


def replay(n, word):
    for expected in word:
        numerator = 3*n+1
        actual = vp(numerator, 2)
        need(actual == expected, 'independent exact valuation replay')
        n = numerator >> actual
    return n


def main():
    report = dict(status='PROVED modular/guard identities; FINITE-EXACT declared universes')
    clocks = []
    for depth in range(1, 7):
        modulus, order = 19**depth, 18*19**(depth-1)
        rho = 27*pow(16, -1, modulus) % modulus
        for base in (2, rho):
            need(pow(base, order, modulus) == 1, 'clock closure')
            need(all(pow(base, order//p, modulus) != 1 for p in (2, 3, 19) if order % p == 0), 'all prime shortening exclusions')
        clocks.append(dict(depth=depth, base2_order=order, rho112_order=order))
    need((2**18-1)//19 % 19 == 3, 'base-two leading lift coefficient')
    need((pow(27*pow(16, -1, 361) % 361, 18, 361)-1)//19 == 4, '112 leading lift coefficient')
    report['clocks'] = clocks

    shells = []
    for depth in range(1, 4):
        modulus = 19**depth
        unseen, histogram = set(range(modulus)), Counter()
        rho, shift = affine_mod((1, 1, 2), depth)
        while unseen:
            start, visited = min(unseen), []
            n = start
            while n in unseen:
                unseen.remove(n)
                visited.append(n)
                n = (rho*n+shift) % modulus
            need(n == start, 'formal permutation cycle closes')
            observations = [chart112_observer(v, depth) for v in visited]
            need(all(o['shell'] == observations[0]['shell'] and o['formal_period'] == len(visited) for o in observations), 'one exact shell orbit')
            histogram[len(visited)] += 1
        need(histogram == Counter({1:1, **{18*19**a:1 for a in range(depth)}}), 'complete shell census')
        shells.append(dict(depth=depth, cycles=dict(sorted(histogram.items()))))
    report['all_residue_shell_census'] = shells

    words = [w for size in range(1, 5) for w in product(range(1, 5), repeat=size)]
    filter_checks = 0
    affine_counts = Counter()
    for word in words:
        p, q, b = word_data(word)
        affine_counts[affine_mod(word, 1)] += 1
        for depth, endpoints in ((1, range(19)), (2, (0, 1, 18, 19, 360))):
            for endpoint in endpoints:
                source, period = endpoint_source_filter(word, endpoint, depth)
                for offset in range(3):
                    n = source+offset*period
                    need(n > 0 and n % 2, 'positive exact source cylinder')
                    result = replay(n, word)
                    need(result == (p*n+b)//q and result % 19**depth == endpoint, 'endpoint filter replay')
                    filter_checks += 1
    report['word_filters'] = dict(words=len(words), depth1_endpoints='all 19', depth2_endpoints=[0,1,18,19,360], offsets=3,
                                  exact_replays=filter_checks, distinct_affine_maps_mod19=len(affine_counts))
    fixed_census = []
    for word in ((4, 4), (1, 7)):
        for depth in range(1, 4):
            multiplier, carry = affine_mod(word, depth)
            actual = [n for n in range(19**depth) if (multiplier*n+carry-n) % 19**depth == 0]
            compressed = fixed_points(word, depth)
            expected = [] if compressed is None else [compressed[0]+j*compressed[1] for j in range(compressed[2])]
            need(actual == expected, 'fixed-axis congruence iff')
            fixed_census.append(dict(word=word, depth=depth, count=len(actual)))
    report['equal_clock_fixed_point_hostile'] = fixed_census
    for depth in range(1, 5):
        long_exponent = 1+18*19**(depth-1)
        need(affine_mod((1,), depth) == affine_mod((long_exponent,), depth), 'full finite-observer opposite-drift alias')
        need(3 > 2 and 3 < 2**long_exponent, 'opposite exact rational drift')
    need(replay(7, (1,1,2)) == 13 and not chart112_observer(13, 2)['one_word_guard'], 'formal chart iteration does not retain source legality')
    report['drift_alias'] = 'For every k, words (1) and (1+18*19^(k-1)) have identical full affine maps mod19^k and opposite drift.'

    hubs, address_tests, lift_tests = (1, 5, 7, 11, 19, 361, 6859), 0, 0
    for hub in hubs:
        v = vp(hub, 19)
        for row in range(3):
            for depth in range(1, 4):
                period = 1 if depth <= v else 3*19**(depth-v-1)
                expected = {ray19(hub % 9, hub, row, b, depth):b for b in range(period)}
                need(len(expected) == period, 'exact inverse channel clock')
                for wanted in range(19**depth):
                    selected = select_ray19(hub % 9, hub % 19**depth, row, wanted, depth)
                    need(selected == (expected[wanted],period) if wanted in expected else selected is None, 'complete address solver iff')
                    address_tests += 1
                if depth > v:
                    for block in range(min(period, 8)):
                        n0 = ray19(hub % 9, hub, row, block, depth+1)
                        coefficient = (3*n0+1)//19**v % 19
                        for digit in range(19):
                            got = ray19(hub % 9, hub, row, block+digit*period, depth+1)
                            need(got == (n0+digit*19**depth*coefficient) % 19**(depth+1), 'nineteen-adic digit carry')
                            lift_tests += 1
    report['ray_address_census'] = dict(hubs=list(hubs), rows=[0,1,2], depths=[1,2,3], all_requested_residues=address_tests, next_digit_checks=lift_tests)

    codec = ray_codec()
    certificates = {u:codec.encode_source(u) for u in hubs}
    joint_cases = joint_hits = 0
    for hub in hubs:
        parent = certificates[hub]
        for wanted3 in range(27):
            for wanted19 in range(19):
                selected = select_certified_predecessors(parent, wanted3, 3, wanted19, 1)
                row = wanted3 % 3
                m19 = 1 if hub % 19 == 0 else 3
                period = lcm(9, m19)
                expected = []
                for b in range(period):
                    literal, _ = codec.ray_source(hub, row, b)
                    if literal % 27 == wanted3 and literal % 19 == wanted19:
                        expected.append(b)
                need(bool(expected) == (selected is not None), 'joint solver complete iff')
                if selected is not None:
                    need(len(expected) == 1 and selected['block'] == expected[0] and selected['period'] == period, 'canonical joint progression')
                    child = instantiate_selected(parent, selected, selected['first_parameter'])
                    literal = codec.literal_certificate_check(child)
                    need(literal % 27 == wanted3 and literal % 19 == wanted19, 'selected first-hit source retains both addresses')
                    joint_hits += 1
                joint_cases += 1
    # An explicit incompatibility: both addresses exist separately in root row2.
    n0, _ = codec.ray_source(1, 2, 0)
    n1, _ = codec.ray_source(1, 2, 1)
    need(select_certified_predecessors(codec.ROOT, n0 % 9, 2, n1 % 19, 1) is None,
         'shared mod-three block digit cannot be discarded')
    root_selection = select_certified_predecessors(codec.ROOT, 1, 2, 1, 1)
    need(root_selection == dict(row=1,block=0,period=3,first_parameter=1), 'isolated root guard, not entire address deletion')
    root_next = codec.expand(instantiate_selected(codec.ROOT, root_selection, 1))
    need(root_next == (2**20-1)//3, 'first retained member beyond self-return')
    report['joint_selector'] = dict(all_cases=joint_cases, compatible_cases=joint_hits,
                                   incompatible_root_row2=dict(ternary=n0%9, mod19=n1%19),
                                   root_exception=root_selection, first_root_address_source=root_next)

    huge = codec.extend(codec.ROOT, 2, 0)
    for i in range(5):
        huge = codec.extend(huge, 1+i%2, 10**100+i)
    wanted3 = 5
    row = wanted3 % 3
    _, block3, _ = select_ray3(codec.mod3(huge, 4), wanted3, 3)
    wanted19 = ray19(codec.mod3(huge, 2), certificate_mod19(huge, 4), row, block3+9*7, 4)
    selection = select_certified_predecessors(huge, wanted3, 3, wanted19, 4)
    need(selection is not None, 'huge certified parent admits compatible refinement')
    child = instantiate_selected(huge, selection, selection['first_parameter'])
    need(codec.mod3(child, 3) == wanted3 and certificate_mod19(child, 4) == wanted19, 'unexpanded child address check')
    lower, upper = codec.bit_bounds(child)
    need(lower > 10**100 and codec.ranks(child)[0] == 7, 'huge certificate keeps exact first-hit rank')
    report['unexpanded_certificate'] = dict(parent_depth=6, child_depth=7, wanted3=wanted3, ternary_depth=3,
        wanted19=wanted19, mod19_depth=4, selection=selection, lower_bits=str(lower), upper_bits=str(upper),
        expanded=False)
    rejects = 0
    for action in (lambda:chart112_observer(True,1), lambda:chart112_observer(7,False),
                   lambda:select_ray19(1,1,2,0,0), lambda:select_ray19(1,1,2,19,1),
                   lambda:instantiate_selected(codec.ROOT,root_selection,0)):
        try:
            action()
        except ValueError:
            rejects += 1
        else:
            raise ValueError('hostile API input accepted')
    need(rejects == 5, 'all declared type/root hostiles rejected')
    report['rejected_hostiles'] = rejects
    report['scope'] = 'Observers filter supplied legal charts and certified predecessor families; no arbitrary-source coverage or residue-based drift theorem.'
    print(json.dumps(report, indent=2))
    print('PASS: exact arithmetic; all checks survive -O')


if __name__ == '__main__':
    main()
