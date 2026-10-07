"""Paid exits from ternary-stripped children, with immutable source comparison.

No ROOT discovery. Conditional receipts consume checked common futures; point
grounding consumes only the preceding package's frozen first-hit certificates.
"""
from dataclasses import dataclass
from hashlib import sha256

import collatz_uncovered_join_routes_20261007 as routes
import collatz_reset2_rules_20261007 as reset
import collatz_child_normal_forms_20261007 as forms


CHECKS = 0
ORPHAN_WORD = (1, 2)+(1,)*5+(1, 2)


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def nat(n, minimum=0):
    need(type(n) is int and n >= minimum, 'exact integer in declared domain')


def valuation(n, prime):
    nat(n, 1)
    need(prime in (2, 3) and type(prime) is int, 'supported exact prime')
    answer = 0
    while n % prime == 0:
        answer += 1
        n //= prime
    return answer


def strip_five(source, repeats):
    routes.odd(source)
    nat(repeats, 0)
    divisor = 9**repeats
    need((source+5) % divisor == 0, 'complete inverse-five guard')
    child = 8**repeats*((source+5)//divisor)-5
    need(child > 0, 'positive inverse-five result')
    return child


def strip_one(source, repeats):
    routes.odd(source)
    nat(repeats, 0)
    divisor = 3**repeats
    need((source+1) % divisor == 0, 'complete inverse-one guard')
    child = 2**repeats*((source+1)//divisor)-1
    need(child > 0, 'positive inverse-one result')
    return child


def inverse_prefix_exit(common, source, word):
    """Prepend a checked actual prefix, then pay against this NEW source."""
    routes.audit(common)
    routes.odd(source)
    need(routes.replay(source, word)[-1] == common.source, 'source-authentic connecting prefix')
    if common.child >= source:
        return None
    return routes.audit(routes.Receipt(source, common.child,
        word+common.source_word, common.child_word, common.endpoint))


def stripped_exit(common, repeats):
    routes.audit(common)
    nat(repeats, 1)
    source = strip_five(common.source, repeats)
    return inverse_prefix_exit(common, source, (1, 2)*repeats)


def payment_threshold(repeats, deletion):
    """Least integer X making h_r > X/2**D-1, or None if impossible."""
    nat(repeats, 1)
    nat(deletion, 1)
    coefficient = 2**deletion*8**repeats-9**repeats
    if coefficient <= 0:
        return None
    toll = 2**(deletion+2)*(9**repeats-8**repeats)
    return toll//coefficient+1


def word_payment_threshold(word, deletion=6):
    """Exact height for any application-order inverse-generator word."""
    nat(deletion, 1)
    record = forms.encode(word)
    p, q = 3**record.B, 2**record.A
    need(p-q <= record.C <= 17*(p-q), 'retained negative-anchor carry envelope')
    coefficient = 2**deletion*q-p
    if coefficient <= 0:
        return None
    return 2**deletion*(q+record.C-p)//coefficient+1


def mixed_word_exit(common, word):
    """Order-preserving mixed inverse prefix followed by an old deletion."""
    routes.audit(common)
    record = forms.encode(word)
    source = forms.apply(record, common.source)
    return inverse_prefix_exit(common, source, forms.forward_word(word))


def word_exit_phase(word, lift=0):
    """K-phase for a paid mixed child of the deep Mersenne rule, or None."""
    nat(lift, 0)
    record = forms.encode(word)
    need(bool(word), 'nonempty inverse word')
    threshold = word_payment_threshold(word)
    phase = forms.mersenne_phase(record)
    if threshold is None or phase is None:
        return None
    exponent, phase_period = phase
    need(exponent % 2 == 1 and phase_period % 2 == 0, 'odd Mersenne phase compatible with deep binary clock')
    odd_period = phase_period//2
    binary = 2**126
    t = (exponent-1457)*pow(binary, -1, odd_period) % odd_period
    first = 1457+binary*t
    period = binary*odd_period
    minimum = max(1457, (threshold-1).bit_length())
    first += period*max(0, (minimum-first+period-1)//period)
    return first+2+period*lift, period, threshold


def mixed_word_residues(word, exponent, modulus):
    nat(exponent, 1459)
    nat(modulus, 1)
    phase = word_exit_phase(word)
    need(phase is not None, 'nonempty paid Mersenne domain')
    first, period, threshold = phase
    need(exponent >= first and (exponent-first) % period == 0, 'supplied source exponent is unchanged and guarded')
    need(exponent-2 >= (threshold-1).bit_length(), 'exact mixed carry height')
    record = forms.encode(word)
    p, q = 3**record.B, 2**record.A
    numerator = (q*pow(2, exponent-2, p*modulus)-q-record.C) % (p*modulus)
    need(numerator % p == 0, 'mixed inverse reader retains denominator precision')
    return numerator//p, (pow(2, exponent-8, modulus)-1) % modulus


def word_native_lift(word, lift=0):
    """Native deep cell intersected with an arbitrary generator-word guard."""
    nat(lift, 0)
    record = forms.encode(word)
    residue, ternary = forms.native_cell(record)
    seed, binary = 2**1457-1, 2**1585
    t = (residue-seed)*pow(binary, -1, ternary) % ternary
    source, period = seed+binary*t, binary*ternary
    threshold = word_payment_threshold(word)
    if threshold is not None:
        source += period*max(0, (threshold-1-source+period-1)//period)
    return source+period*lift


@dataclass(frozen=True)
class MixedPlan:
    exponent: int
    repeats: int


def validate(plan):
    need(type(plan) is MixedPlan, 'typed symbolic mixed plan')
    nat(plan.exponent, 1459)
    nat(plan.repeats, 1)
    need((plan.exponent-1459) % 2**126 == 0, 'deep deletion exponent phase')
    need(1+valuation(plan.exponent-4, 3) >= 2*plan.repeats, 'complete ternary strip precision')
    threshold = payment_threshold(plan.repeats, 6)
    need(threshold is not None and plan.exponent-2 >= (threshold-1).bit_length(), 'immutable mixed-source payment')
    return plan


def mixed_phase(repeats, unit, lift=0):
    """Exact maximal-r Mersenne phase, including unpaid r>=36 for controls."""
    nat(repeats, 1)
    need(type(unit) is int and unit in (1, 2), 'nonzero ternary phase digit')
    nat(lift, 0)
    binary, ternary = 2**126, 3**(2*repeats)
    target = 4+unit*3**(2*repeats-1)
    t = (target-1459)*pow(binary, -1, ternary) % ternary
    return 1459+binary*t+binary*ternary*lift


def plan_residues(plan, modulus):
    """Read both sources without constructing any enormous Mersenne integer."""
    validate(plan)
    nat(modulus, 1)
    divisor = 9**plan.repeats
    value = pow(2, plan.exponent-2, divisor*modulus)+4
    need(value % divisor == 0, 'retained division precision')
    source = (8**plan.repeats*(value//divisor)-5) % modulus
    child = (pow(2, plan.exponent-8, modulus)-1) % modulus
    return source, child


def expand_plan(plan, bit_cap):
    validate(plan)
    nat(bit_cap, 1)
    need(plan.exponent <= bit_cap, 'declared materialization cap')
    middle = 2**(plan.exponent-2)-1
    common = reset.receipt(middle, reset.load_rules()[2])
    need(common is not None, 'deep collision is supplied by its actual guard')
    result = stripped_exit(common, plan.repeats)
    need(result is not None, 'symbolic payment agrees with exact materialization')
    return result


def native_lift(repeats, unit, lift=0):
    """Finite-bit controls in the complete deep native cell, not just Mersennes."""
    nat(repeats, 1)
    need(type(unit) is int and unit in (1, 2), 'nonzero ternary phase digit')
    nat(lift, 0)
    seed, binary, ternary = 2**1457-1, 2**1585, 3**(2*repeats+1)
    target = unit*9**repeats-5
    t = (target-seed)*pow(binary, -1, ternary) % ternary
    return seed+binary*(t+ternary*lift)


def orphan_chain(exponent, bit_cap):
    """The exact G5, J^5, G5 chain; return its genuinely new paid exit."""
    nat(exponent, 1459)
    nat(bit_cap, 1)
    need((exponent-1459) % 2**126 == 0, 'deep binary phase')
    need((exponent-1459) % (2*3**10) == 0, 'complete maximal inverse-chain phase')
    need(exponent <= bit_cap, 'declared materialization cap')
    middle = 2**(exponent-2)-1
    first = strip_five(middle, 1)
    second = strip_one(first, 5)
    source = strip_five(second, 1)
    need(valuation(middle+5, 3) == 2 and valuation(first+1, 3) == 5,
         'first two inverse packets are maximal')
    need(valuation(second+5, 3) == 3 and valuation(source+5, 3) == 1
         and valuation(source+1, 3) == 0, 'last inverse packet and exhausted guards')
    common = reset.receipt(middle, reset.load_rules()[2])
    result = inverse_prefix_exit(common, source, ORPHAN_WORD)
    need(result is not None, 'exhausted inverse child has a paid composed exit')
    return (middle, first, second, source), result


def symbolic_orphan_control(exponent):
    """Evaluate only the precisions needed for the exact inverse-chain guards."""
    nat(exponent, 1459)
    need((exponent-1459) % (2*3**10) == 0, 'ternary chain phase')
    modulus = 3**11
    middle = (pow(2, exponent-2, modulus)-1) % modulus
    need((middle+5) % 9 == 0, 'first inverse packet residue')
    first = (8*((middle+5)//9)-5) % 3**9
    need((first+1) % 3**5 == 0, 'inverse-one packet residue')
    second = (32*((first+1)//3**5)-1) % 3**4
    need((second+5) % 9 == 0, 'second inverse-five residue')
    source = (8*((second+5)//9)-5) % 9
    need(valuation(middle+5, 3) == 2 and valuation(first+1, 3) == 5
         and valuation(second+5, 3) == 3 and source == 1,
         'maximal packet pattern and stopped normalizer')
    return source


def rejects(operation):
    try:
        operation()
    except (ValueError, TypeError):
        return True
    return False


def main():
    print('PROVED: pay after composition against the mixed source, not the larger middle.')
    need(64*8**35 > 9**35 and 64*8**36 < 9**36, 'sharp coefficient crossing')
    need(payment_threshold(35, 6) <= 8192, 'uniform positive height cut')
    for r in range(1, 61):
        bound = payment_threshold(r, 6)
        need((bound is not None) == (r <= 35), 'exact six-deletion budget')
        if bound is not None:
            coefficient = 64*8**r-9**r
            toll = 256*(9**r-8**r)
            need(coefficient*bound > toll and coefficient*(bound-1) <= toll, 'least strict integer height')
    print('D6 strip budget: r1..35 pay for X>=8192; every r>=36 is unpaid for all positive X.')

    deep = reset.load_rules()[2]
    for r in range(1, 37):
        for unit in (1, 2):
            for lift in (0, 2):
                middle = native_lift(r, unit, lift)
                common = reset.receipt(middle, deep)
                need(common is not None and common.child == (middle+1)//64-1, 'complete native deep cell')
                source = strip_five(middle, r)
                need(valuation(middle+5, 3) == 2*r, 'exact maximal strip depth')
                need(valuation(source+5, 2) == 3*r+2 and valuation(source+1, 2) == 2, 'new child leaves the Mersenne one-run state')
                result = stripped_exit(common, r)
                need((result is not None) == (r <= 35), 'literal receipt and sharp source payment agree')
                if result is None:
                    need(common.child > source, 'r36 is a real source-size failure')
    print('FINITE-EXACT native deep cell: r1..36, two ternary units, lifts0/2;140 paid and4 unpaid receipts.')

    for r in range(1, 41):
        for unit in (1, 2):
            for lift in (0, 1, 3):
                k = mixed_phase(r, unit, lift)
                need((k-1459) % 2**126 == 0 and valuation(k-4, 3) == 2*r-1, 'independent exact CRT/maximal-depth guard')
                plan = MixedPlan(k, r)
                if r <= 35:
                    validate(plan)
                    for modulus in (1, 8, 9, 19, 105):
                        hs, _ = plan_residues(plan, modulus)
                        need(type(hs) is int and 0 <= hs < modulus, 'bounded exact modular reader')
                else:
                    need(rejects(lambda: validate(plan)), 'unpaid symbolic plan rejected before expansion')
    concrete = expand_plan(MixedPlan(1459, 1), 2000)
    for modulus in (3, 19, 105, 2187, 2**32):
        need(plan_residues(MixedPlan(1459, 1), modulus) ==
             (concrete.source % modulus, concrete.child % modulus), 'symbolic and materialized source identities')
    print('SYMBOLIC phases: r1..40, two units, lifts0/1/3; no astronomical Mersenne source expansion.')

    need(2187**63 < 64*2048**63 and 2187**64 > 64*2048**64,
         'sharp maximum length of the six-bit mixed budget')
    need(9*2187**62 >= 64*8*2048**62, 'only allG17 survives at length63')
    need(1024*3**441+1 < 2**710, 'every funded mixed height is below the deep phase height')
    cases = ((1,), (1,)*10, (1,)*11, (5,), (5,)*35, (5,)*36,
             (17,), (17,)*63, (17,)*64, (5, 1, 5), (17, 1, 5),
             (5, 17), (17, 5), (17,)+(5,)*34, (5,)*34+(17,),
             (5,)+(1,)*5+(5,))
    paid_words = 0
    for word in cases:
        record = forms.encode(word)
        p, q = 3**record.B, 2**record.A
        threshold = word_payment_threshold(word)
        if threshold is not None:
            coefficient = 64*q-p
            toll = 64*(q+record.C-p)
            need(coefficient*threshold > toll and coefficient*(threshold-1) <= toll,
                 'mixed-word least integer height with ordered carry')
        for lift in (0, 1):
            middle = word_native_lift(word, lift)
            common = reset.receipt(middle, deep)
            result = mixed_word_exit(common, word)
            need((result is not None) == (threshold is not None), 'native mixed word realizes the exact slope/height criterion')
            paid_words += result is not None
        for lift in (0, 1, 3):
            phase = word_exit_phase(word, lift)
            if threshold is None or word[0] == 1:
                need(phase is None, 'unpaid word or empty Mersenne domain remains explicit')
            else:
                k, _, _ = phase
                native = forms.mersenne_phase(record)
                need((k-2-native[0]) % native[1] == 0 and (k-1459) % 2**126 == 0,
                     'independent odd-exponent CRT compatibility')
                for modulus in (1, 9, 19, 105):
                    mixed_word_residues(word, k, modulus)
    need(word_payment_threshold((17,)+(5,)*34) == 2726
         and word_payment_threshold((5,)*34+(17,)) == 3244,
         'same slope does not retain ordered carry or exact height')
    need(word_payment_threshold((17,)*63) == 45594, 'maximal-length mixed positive control')
    for r in range(1, 61):
        need(word_payment_threshold((5,)*r) == payment_threshold(r, 6),
             'general formula specializes to the entire earlier strip budget')
    print('MIXED WORD CONTROLS:16 selected words, two native lifts;', paid_words,
          'paid and', 32-paid_words, 'unpaid. Maximum funded length63, uniquely G17^63 at that length.')
    print('ORDERED HEIGHT HOSTILE: G17 then G5^34 needs X>=2726; reverse order needs X>=3244.')

    states, orphan = orphan_chain(1459, 2000)
    x = 2**1457
    need(states == (x-1, (8*x-13)//9, (256*x-2315)//2187,
                    (2048*x-29455)//19683), 'exact stopped-chain coefficients')
    need(routes.carrier(ORPHAN_WORD) == (19683, 2048, 27407), 'ordered nine-letter return carrier')
    for lift in range(32):
        k = 1459+2**126*3**10*lift
        need(symbolic_orphan_control(k) == 1, 'all-height stopped inverse chain')
    print('NORMALIZER at K1459: G5,J^5,G5; state residues mod2187', tuple(n % 2187 for n in states))
    print('Final z=(2048*2^(K-2)-29455)/19683; z1mod9 stops both inverse guards, but composed D6 exit is paid.')

    words = reset.grounded_named_points()
    for name, common, prefix in (('h1', concrete, (1, 2)), ('z', orphan, ORPHAN_WORD)):
        compiled = routes.discharge(common, words[1451])
        need(compiled == prefix+words[1457], 'independent actual suffix identity, not an extra ROOT search')
        need(routes.replay(common.source, compiled)[-1] == 1, 'fully grounded selected mixed point')
        print('ROOT point', name, 'odd rank', len(compiled), 'valuation cost', sum(compiled),
              'word SHA256', sha256(bytes(compiled)).hexdigest())

    hostiles = (
        lambda: validate(MixedPlan(1459.0, 1)), lambda: validate(MixedPlan(1459, True)),
        lambda: validate(MixedPlan(1461, 1)), lambda: expand_plan(MixedPlan(1459, 1), 100),
        lambda: strip_five(True, 1), lambda: strip_five(3, 1),
        lambda: strip_one(7, 1), lambda: plan_residues(MixedPlan(1459, 1), 0),
        lambda: inverse_prefix_exit(reset.receipt(2**1457-1, deep), states[-1], (1, 2)),
        lambda: routes.discharge(orphan, words[1451][:-1]),
        lambda: word_payment_threshold((True,)),
        lambda: word_exit_phase((5,), True),
        lambda: mixed_word_residues((5,), 1459+2, 19),
        lambda: mixed_word_residues((5,), 1459, 0),
    )
    need(all(rejects(f) for f in hostiles), 'typed, native-guard, source-identity and missing-suffix hostiles')
    print('HOSTILES', len(hostiles), 'rejected. No global closure or unsupplied child ROOT is asserted.')
    print('CHECKS', CHECKS)


if __name__ == '__main__':
    main()
