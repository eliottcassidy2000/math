"""Exact old-controller / actual-step / run bridges, with source-owned payment.

No production orbit discovery.  A paid child is grounded only by a supplied
strict ROOT word.  The complete classification concerns this declared template.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from math import gcd

import adaptive_credit_potential_20261004 as old
import collatz_projective_receipts_20261009 as runs
import collatz_child_normal_forms_20261007 as inv
import translation_phase_decoder_20261005 as five
import paper_local_global_grounding_20261007 as native

CHECKS = 0
NAMES = ('H', 'G', 'A', 'B', 'L', 'LG')
CLASSES = {'H': (52, 54), 'G': (0, 18), 'A': (8, 54),
           'B': (32, 54), 'L': None, 'LG': (36, 54)}
KAPPA = F(2349, 2560)


def need(ok, why):
    if not ok:
        raise ValueError(why)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def integer(n, minimum=0):
    need(type(n) is int and n >= minimum, 'exact integer in range')


def prefix(name):
    need(type(name) is str and name in NAMES, 'declared old prefix')
    p, q, b = five.carrier(name)
    c, modulus = (219, 256) if name == 'LG' else five.LETTERS[name][3:]
    return p, q, b, c, modulus


@dataclass(frozen=True)
class Bridge:
    name: str
    exponent: int
    run: int
    residue: int
    modulus: int


def compile_bridge(name, exponent, run=4):
    """Complete native source ray for the fixed template, or None if empty."""
    p, q, b, c, modulus = prefix(name)
    integer(exponent, 1); integer(run, 4)
    output, increment = (p*c+b)//q, p*modulus//q
    row = runs.row(run)
    target, step = row.base[0], row.period[0]
    denominator = (1 << exponent)*step
    rhs = (1 << exponent)*target - (3*output+1)
    common = gcd(3*increment, denominator)
    if rhs % common:
        return None
    period = denominator//common
    phase = (rhs//common)*pow(3*increment//common, -1, period) % period
    return Bridge(name, exponent, run, c+modulus*phase, modulus*period)


def audit_bridge(spec):
    need(type(spec) is Bridge, 'exact Bridge')
    prefix(spec.name)
    for field in (spec.exponent, spec.run, spec.residue, spec.modulus):
        integer(field, 1)
    need(spec == compile_bridge(spec.name, spec.exponent, spec.run),
         'canonical native progression required')
    return spec


def classified(name, exponent):
    """All-height classification for the R4 connector, not a finite census."""
    prefix(name); integer(exponent, 1)
    entry = CLASSES[name]
    return entry is not None and exponent % entry[1] == entry[0]


def source_word(name, exponent):
    prefix(name); integer(exponent, 1)
    if name in five.ACTUAL:
        return five.ACTUAL[name]+(exponent, 1)
    if name == 'H':
        return old.HEAD+(exponent+2, 1)
    if name == 'LG':
        return (1, 2, 1, 1, exponent+2, 1)
    raise ValueError('L has no one-step bridge')


@dataclass(frozen=True)
class Batch:
    bridge: Bridge
    source: int
    old_output: int
    run_source: int
    child: int
    endpoint: int
    left: tuple
    right: tuple


def receipt(spec, supplied_source):
    """Keep the supplied original source; never replace it by a nearby lift."""
    audit_bridge(spec); old.positive_odd(supplied_source)
    need(supplied_source >= spec.residue and
         (supplied_source-spec.residue) % spec.modulus == 0, 'exact source guard')
    p, q, b, _, _ = prefix(spec.name)
    x = (p*supplied_source+b)//q
    y, a = old.step(x)
    need(a == spec.exponent, 'exact connector valuation')
    run = runs.join(y, spec.run)
    left = source_word(spec.name, a)
    right = run.child_word
    need(run.child < y < supplied_source, 'immutable-source strict payment')
    need(native.replay(supplied_source, left) == run.endpoint and
         native.replay(run.child, right) == run.endpoint, 'strict actual common future')
    return Batch(spec, supplied_source, x, y, run.child, run.endpoint, left, right)


def at(spec, parameter):
    audit_bridge(spec); integer(parameter)
    return receipt(spec, spec.residue+spec.modulus*parameter)


def audit_batch(batch):
    need(type(batch) is Batch, 'exact Batch')
    audit_bridge(batch.bridge)
    for field in (batch.source, batch.old_output, batch.run_source, batch.child, batch.endpoint):
        old.positive_odd(field)
    runs.carrier(batch.left); runs.carrier(batch.right)
    need(batch == receipt(batch.bridge, batch.source), 'authenticated actual batch')
    return batch


def discharge(batch, supplied_child_root):
    """First-hit substitution of a supplied child ROOT proof; no discovery."""
    audit_batch(batch)
    need(native.replay(batch.child, supplied_child_root) == 1, 'supplied child ROOT')
    length = len(batch.right)
    need(supplied_child_root[:length] == batch.right, 'actual common-future prefix')
    word = batch.left+supplied_child_root[length:]
    need(native.replay(batch.source, word) == 1, 'strict original-source ROOT')
    return word


def g_bridge(index):
    integer(index, 1)
    return compile_bridge('G', 18*index)


def inverse_mask(source):
    old.positive_odd(source)
    return tuple(g for g in inv.GENERATORS
                 if source % (3**inv.GENERATORS[g][1]) ==
                 inv.native_cell(inv.encode((g,)))[0])


def energy_band(energy):
    """Incoming exact band b(E); no floating logarithm or terminal inference."""
    need(type(energy) is F and energy >= 116, 'Fraction energy at least 116')
    result = 0
    while energy*KAPPA >= 116:
        energy *= KAPPA
        result += 1
    return result


def credit_rank(source, credits):
    old.positive_odd(source); integer(credits)
    return 3*energy_band(F(source+5)*F(9, 8)**credits)+credits


def g_band_payment(index):
    """Largest number of kappa factors certified by the sharp G-batch ratio."""
    spec = g_bridge(index)
    factor = F(228, spec.residue+5)
    power, count = F(1), 0
    while factor <= power*KAPPA:
        power *= KAPPA
        count += 1
    return factor, count


def reject(call, why):
    try:
        call()
    except (ValueError, TypeError):
        check(True, why)
    else:
        check(False, why)


def main():
    # Independent tiny unit-group reader, rather than trusting the classification.
    for modulus, order in ((27, 18), (81, 54)):
        residues = tuple(pow(2, a, modulus) for a in range(order))
        check(len(set(residues)) == order and pow(2, order, modulus) == 1,
              'complete primitive unit cycle')
    targets = {'H': 61, 'A': 13, 'B': 49, 'LG': 55}
    for name, target in targets.items():
        check([a for a in range(54) if pow(2, a, 81) == target] == [CLASSES[name][0]],
              'independent discrete-log phase')
    controls = 0
    for name in NAMES:
        for a in range(1, 217):
            spec = compile_bridge(name, a)
            check((spec is not None) == classified(name, a), 'full phase vs gcd compiler')
            if spec is None:
                continue
            for t in (0, 1, 2, 7):
                batch = at(spec, t)
                check(native.replay(batch.source, batch.left) == batch.endpoint,
                      'literal original source word')
                check(native.replay(batch.child, batch.right) == batch.endpoint,
                      'literal child word')
                # Independent production library evaluation preserves every native action.
                current = batch.source
                for letter in name:
                    current = five.operate(current, letter)
                check(current == batch.old_output, 'old action guard and carrier')
                check(batch.child < batch.run_source < batch.source, 'paid relative to original')
                controls += 1
    general = 0
    for name in NAMES:
        for a in range(1, 109):
            for run in range(4, 9):
                spec = compile_bridge(name, a, run)
                if spec is not None:
                    for t in (0, 1):
                        batch = at(spec, t)
                        check(batch.child < batch.source and batch.child % 16 == 15,
                              'generic run exact receipt')
                        general += 1
    # Full G sequence, with a closed formula independent of the progression solver.
    for index in range(1, 19):
        a = 18*index; spec = g_bridge(index)
        check(spec.residue == (2264*(1 << a)-23)//27 and
              spec.modulus == 96*(1 << a), 'unbounded G formula, finite controls')
        check(spec.residue % 3 == (2-index) % 3 and spec.residue % 256 == 75,
              'retained ternary and binary source phases')
        for t in (0, 1, 9):
            batch = at(spec, t)
            check((batch.run_source, batch.child, batch.endpoint) ==
                  (283+324*t, 223+256*t, 425+486*t), 'same target ray at every G exponent')
        word = discharge(at(spec, 0), runs.ROOT223)
        check(len(word) == 23 and sum(word) == a+43, 'uniform grounded G base certificate')
    a5 = compile_bridge('A', 8, 5)
    check((a5.residue, a5.modulus) == (81851, 131072), 'A guard unchanged by R5')
    child319 = runs.ROOT223[6:]
    check(native.replay(319, child319) == 1, 'frozen prior proof suffix only')
    word = discharge(at(a5, 0), child319)
    check((len(word), sum(word)) == (18, 45), 'A bridge authenticated ground')
    lg5 = compile_bridge('LG', 36, 5)
    word = discharge(at(lg5, 0), child319)
    check((lg5.residue, lg5.modulus, len(word), sum(word)) ==
          (21972133598683, 35184372088832, 18, 73), 'L barrier crossed by retained G prefix')
    # Exact declared inverse-bank comparison, not a claim about all prior algorithms.
    for index, expected in ((1, F(244, 729)), (2, F(0)), (3, F(1))):
        spec = g_bridge(index)
        hit = sum(bool(inverse_mask(spec.residue+spec.modulus*t)) for t in range(729))
        check(F(hit, 729) == expected, 'source-owned inverse-bank mass')
        for t in range(729):
            n = spec.residue+spec.modulus*t
            check(all(n % five.LETTERS[c][4] != five.LETTERS[c][3] for c in 'HABL'),
                  'no old funded native first action')
    hit = sum(bool(inverse_mask(a5.residue+a5.modulus*t)) for t in range(2187))
    check(F(hit, 2187) == F(973, 2187), 'A inverse-bank comparison')
    # Mod9 loses the carry precision which distinguishes the correct a classes.
    check(compile_bridge('H', 16) is None and compile_bridge('B', 14) is None,
          'coarse output phase is insufficient')
    check(compile_bridge('L', 4) is None and lg5 is not None,
          'one connector impossible does not forbid a longer cascade')
    first = at(g_bridge(1), 0)
    check(first.old_output > first.source and first.child < first.source,
          'atomic payment despite expanding first action')
    check(first.source == 21981259 and not inverse_mask(first.source), 'new paid entry witness')
    check(runs.available(first.source) < 4, 'entry is outside inherited run domain')
    for index in range(1, 7):
        factor, bands = g_band_payment(index)
        check(KAPPA**(bands+1) < factor <= KAPPA**bands, 'exact power budget')
        for t in (0, 1, 7):
            batch = at(g_bridge(index), t)
            ratio = F(batch.child+5, batch.source+5)
            check(ratio <= factor < KAPPA and (ratio == factor) == (t == 0),
                  'sharp carried energy ratio')
            for credits in (0, 3):
                check(credit_rank(batch.source, credits)-credit_rank(batch.child, credits)
                      >= 3*bands, 'incoming integer rank at unchanged credit')
    check(g_band_payment(1)[1] == 133, 'uniform 399-rank decrement')
    check((credit_rank(first.source, 1), credit_rank(first.child, 1)) == (427, 28),
          'sharp actual 399 decrement at retained credit one')
    check((credit_rank(first.source, 0), credit_rank(first.child, 0)) == (423, 21),
          'zero-credit actual decrement need not equal the uniform bound')
    reject(lambda: compile_bridge('G', True), 'bool valuation')
    reject(lambda: compile_bridge('G', 18.0), 'float valuation')
    reject(lambda: compile_bridge('G', 18, 3), 'insufficient run')
    reject(lambda: compile_bridge('GG', 18), 'undeclared template')
    reject(lambda: at(g_bridge(1), -1), 'negative ray parameter')
    reject(lambda: at(g_bridge(1), True), 'bool ray parameter')
    reject(lambda: receipt(g_bridge(1), first.source+2), 'nearby source not substituted')
    reject(lambda: receipt(g_bridge(1), 1), 'ROOT not an entry source')
    reject(lambda: receipt(g_bridge(1), float(first.source)), 'float source')
    reject(lambda: audit_bridge(replace(g_bridge(1), exponent=18.0)), 'forged cached field')
    reject(lambda: audit_bridge(replace(g_bridge(1), residue=first.source+2)), 'forged native row')
    reject(lambda: audit_batch(replace(first, child=223.0)), 'forged child alias')
    reject(lambda: discharge(first, (2,)), 'unrelated terminal ROOT word')
    reject(lambda: discharge(first, runs.ROOT223+(2,)), 'ROOT padding rejected')
    reject(lambda: credit_rank(first.source, True), 'bool credit')
    reject(lambda: energy_band(116.0), 'float energy')
    reject(lambda: energy_band(F(115)), 'energy outside band domain')
    print('Atomic run bridges: exact guarded old-prefix / U_a / R_r batches.')
    print('R4 complete classes: H a52mod54; G a0mod18 positive; A a8mod54; B a32mod54; L empty; LG a36mod54.')
    print(f'Finite controls: {controls} classified R4 lifts; {general} arbitrary-run lifts; 18 G exponent families.')
    print('G: n=(2264*2^(18j)-23)/27+96*2^(18j)*t, j>=1,t>=0; child223+256t.')
    print('G j1 base21981259, period25165824, ROOT rank23/cost61 from supplied ROOT223.')
    print('G inverse-bank complements by j mod3=1,2,0: 485/729, 1, 0; all old HABL entry guards missed.')
    print('A8/R5: 81851+131072t ->319+512t; base ROOT rank18/cost45 from prior ROOT319 suffix.')
    print('LG36/R5: 21972133598683+35184372088832t ->319+512t; base ROOT rank18/cost73.')
    print('Incoming integer-rank bridge: unchanged credit; sharp first ratio228/21981264; at least399 rank units paid.')
    print('No all-parameter ROOT conclusion; no new guard mass over the full ordinary-descent algorithm.')
    print(f'PASS: {CHECKS} exact checks.')


if __name__ == '__main__':
    main()
