"""Complete finite collision search; adaptive depth is not universal arrival.

Run from repository root. Exact integers only; no third-party dependencies.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache
import json

import collatz_uncovered_join_routes_20261007 as routes

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def letters(word):
    need(type(word) is tuple and bool(word), 'nonempty tuple')
    need(all(type(a) is int and a >= 1 for a in word), 'exact positive letters')
    need(word[0] >= 2, 'reduced first letter')


def anchor(word):
    """f_word(-1) = N/2^E in lowest terms."""
    letters(word)
    n, e = -1, word[0]-1
    for a in word[1:]:
        n, e = 3*n+(1 << e), e+a
    need(n % 2 == 1, 'odd numerator retained')
    return n, e


@lru_cache(None)
def _partners(n, e):
    """One lexicographically least witness for EACH attainable word length."""
    found = {1: (e+1,)} if n == -1 else {}
    for f in range(1, e):
        numerator = n-(1 << f)
        if numerator % 3:
            continue
        for length, prefix in _partners(numerator//3, f):
            candidate = prefix+(e-f,)
            key = length+1
            if key not in found or candidate < found[key]:
                found[key] = candidate
    return tuple(sorted(found.items()))


def partners(n, e):
    need(type(n) is int and n % 2 == 1, 'odd exact numerator')
    need(type(e) is int and e >= 1, 'positive exact denominator exponent')
    return dict(_partners(n, e))


def choose_partner(word, run):
    letters(word)
    need(type(run) is int and run >= 1, 'positive exact available run')
    candidates = partners(*anchor(word))
    possible = [p for p in candidates if len(word) < p <= len(word)+run]
    return candidates[max(possible)] if possible else None


def trim_root(source, word):
    """Replay fully, requiring any discarded suffix to be only ROOT loops."""
    routes.odd(source)
    routes.letters(word)
    actual = []
    x = source
    for i, a in enumerate(word):
        if x == 1:
            need(all(b == 2 for b in word[i:]), 'only actual ROOT loops discarded')
            break
        x, b = routes.step(x)
        need(a == b, 'actual first-hit exponent')
        actual.append(a)
    return tuple(actual), x


def collision_receipt(source, core, partner):
    routes.odd(source)
    need(source > 1, 'nonroot source')
    letters(core)
    letters(partner)
    run = routes.v2(source+1)-1
    deletion = len(partner)-len(core)
    need(1 <= deletion <= run, 'strict available trailing-ones deletion')
    need(anchor(core) == anchor(partner), 'exact rational collision')
    child = (source+1)//(1 << deletion)-1
    left, endpoint = trim_root(source, (1,)*run+core)
    right, other = trim_root(child, (1,)*(run-deletion)+partner)
    need(endpoint == other, 'first-hit common endpoint')
    return routes.audit(routes.Receipt(source, child, left, right, endpoint))


def negative_five_lift(receipt, repeats):
    """Prepend (1,2) inverse-cycle receipts; preserve the immutable source."""
    routes.audit(receipt)
    need(type(repeats) is int and repeats >= 0, 'exact nonnegative repetition count')
    divisor = 9**repeats
    need((receipt.child+5) % divisor == 0, 'complete ternary repetition guard')
    child = 8**repeats*((receipt.child+5)//divisor)-5
    need(child > 0, 'positive inverse-cycle child')
    return routes.audit(routes.Receipt(receipt.source, child,
        receipt.source_word, (1, 2)*repeats+receipt.child_word, receipt.endpoint))


def mixed_phase(repeats, unit):
    """An infinite 1459 phase with exactly this many inverse-five repetitions."""
    need(type(repeats) is int and repeats >= 1, 'positive exact repetition depth')
    need(type(unit) is int and unit in (1, 2), 'nonzero ternary digit')
    ternary = 3**(2*repeats)
    target = 4+unit*3**(2*repeats-1)
    binary = 1 << 21
    t = (target-1459)*pow(binary, -1, ternary) % ternary
    return 1459+binary*t, binary*ternary


@dataclass(frozen=True)
class Probe:
    source: int
    core: tuple
    reason: str
    partner: tuple = ()
    receipt: object = None
    root_word: tuple = ()


def probe(source, max_letters=12, max_cost=32):
    """Complete at each inspected prefix; explicit caps retain an OPEN result."""
    routes.odd(source)
    need(type(max_letters) is int and max_letters >= 1, 'positive letter cap')
    need(type(max_cost) is int and max_cost >= 2, 'positive cost cap')
    if source == 1:
        return Probe(source, (), 'ROOT')
    run = routes.v2(source+1)-1
    if run == 0:
        return Probe(source, (), 'NO_TRAILING_ONES')
    t = (source+1)//(1 << (run+1))
    x = 2*3**run*t-1
    core = ()
    cost = 0
    for _ in range(max_letters):
        if x == 1:
            root_word = (1,)*run+core
            need(routes.replay(source, root_word)[-1] == 1, 'first-hit ROOT certificate')
            return Probe(source, core, 'ROOT_CERTIFIED', root_word=root_word)
        next_x, a = routes.step(x)
        if cost+a > max_cost:
            return Probe(source, core, 'COST_FRONTIER')
        core += (a,)
        cost += a
        candidate = choose_partner(core, run)
        if candidate is not None:
            receipt = collision_receipt(source, core, candidate)
            return Probe(source, core, 'PAID', candidate, receipt)
        x = next_x
        if x == 1:
            root_word = (1,)*run+core
            need(routes.replay(source, root_word)[-1] == 1, 'first-hit ROOT at cap')
            return Probe(source, core, 'ROOT_CERTIFIED', root_word=root_word)
    return Probe(source, core, 'LETTER_FRONTIER')


def compositions(total, first=True):
    """Independent exhaustive word universe; first letter at least two."""
    for a in range(2 if first else 1, total+1):
        if a == total:
            yield (a,)
        else:
            for rest in compositions(total-a, False):
                yield (a,)+rest


def independent_anchor(word):
    x = Fraction(-1)
    for a in word:
        x = (3*x+1)/(1 << a)
    return x


def rejects(callback):
    try:
        callback()
    except (ValueError, TypeError):
        return True
    return False


def main():
    brute_words = 0
    brute_values = 0
    # Exhaust all reduced compositions, not a bound on each individual letter.
    for total in range(2, 16):
        census = {}
        for word in compositions(total):
            brute_words += 1
            value = independent_anchor(word)
            by_length = census.setdefault(value, {})
            old = by_length.get(len(word))
            if old is None or word < old:
                by_length[len(word)] = word
            n, e = anchor(word)
            need(value == Fraction(n, 1 << e), 'independent affine anchor')
        for value, expected in census.items():
            brute_values += 1
            e = value.denominator.bit_length()-1
            need(partners(value.numerator, e) == expected,
                 'complete attainable lengths and canonical witnesses')
    print('COMPLETE COMPOSITION CENSUS', json.dumps({
        'total_cost': [2, 15], 'words': brute_words, 'anchor_values': brute_values}))

    all_twos = []
    for d in range(1, 17):
        p = partners(*anchor((2,)*d))
        need(max(p) == d, 'finite all-two longest-partner control')
        all_twos.append([d, max(p)])
    print('ALL-TWO PREFIXES FINITE ONLY', json.dumps(all_twos))

    control = probe(3)
    need(control.reason == 'PAID' and control.receipt.child == 1,
         'empty child ROOT word supported')
    need(control.receipt.child_word == (), 'ROOT suffix removed, not exported')
    cases = []
    for k in (7, 31, 1459):
        item = probe((1 << k)-1, 12, 32)
        if k in (31, 1459):
            need(item.reason == 'PAID', 'known sporadic positive control')
        cases.append({'K': k, 'reason': item.reason, 'core': item.core,
                      'partner': item.partner,
                      'deletion': len(item.partner)-len(item.core)
                      if item.receipt else None})
    print('MERSENNE SELECTOR', json.dumps(cases))
    paid = probe((1 << 1459)-1, 12, 32)
    need(len(paid.core) == 12 and paid.receipt.child == (1 << 1457)-1,
         '1459 direct two-bit deletion')
    need(probe((1 << 1459)-1, 11, 32).reason == 'LETTER_FRONTIER',
         'depth11 remains explicitly unresolved')
    need(probe((1 << 1459)-1, 12, 24).reason == 'COST_FRONTIER',
         'cost24 remains explicitly unresolved')
    mixed = negative_five_lift(paid.receipt, 1)
    need(mixed.child == (2*paid.source-11)//9,
         '1459 binary collision then ternary cycle lift')
    need(rejects(lambda: negative_five_lift(paid.receipt, 2)),
         '1459 pays exactly one ternary inverse cycle')
    import collatz_terminal_lifts_20261007 as lifts
    mixed_controls = []
    for r in range(1, 5):
        exponent = 2+3**(2*r-1)
        source = (1 << (exponent+1))-1
        base = lifts.reset_receipt(source)
        need(base.child == (1 << exponent)-1, 'guarded even-Mersenne reset')
        converted = negative_five_lift(base, r)
        need(converted.child < base.child < source, 'strict mixed guard reduction')
        need(rejects(lambda: negative_five_lift(base, r+1)),
             'exact terminal repetition budget')
        mixed_controls.append([exponent, r, converted.child.bit_length()])
    need(rejects(lambda: negative_five_lift(paid.receipt, True)),
         'typed inverse-cycle repetitions')
    print('MIXED BINARY AND TERNARY CONTROLS', json.dumps(mixed_controls))
    for r in range(1, 13):
        for unit in (1, 2):
            first, period = mixed_phase(r, unit)
            for t in (0, 1, 10**6+7):
                k = first+period*t
                need(k >= 1459 and k % (1 << 21) == 1459,
                     'retained collision phase at arbitrary integer height')
                residue = (pow(2, k-2, 9**(r+1))+4) % 9**(r+1)
                need(residue % 9**r == 0 and residue % 9**(r+1) != 0,
                     'exact maximal ternary guard without Mersenne expansion')
    print('SYMBOLIC MIXED PHASES', 'depths1..12, both nonzero ternary digits, three lifts')

    root_control = probe(127, 9, 32)
    need(root_control.reason == 'ROOT_CERTIFIED', 'ROOT at the final budget letter')
    need(routes.replay(127, root_control.root_word)[-1] == 1, 'saved full ROOT route')
    sample = {'PAID': 0, 'ROOT_CERTIFIED': 0,
              'LETTER_FRONTIER': 0, 'COST_FRONTIER': 0}
    reset_two = {key: 0 for key in sample}
    for source in range(3, 512, 4):
        item = probe(source, 8, 22)
        sample[item.reason] += 1
        k = routes.v2(source+1)
        t = (source+1)//(1 << k)
        if routes.step(2*3**(k-1)*t-1)[1] == 2:
            reset_two[item.reason] += 1
        if item.receipt:
            need(routes.audit(item.receipt) == item.receipt, 'independent receipt replay')
            for position in (0, -1):
                changed = list(item.receipt.source_word)
                changed[position] += 1
                need(rejects(lambda: routes.audit(replace(item.receipt,
                    source_word=tuple(changed)))), 'forged letter rejected')
        if item.reason == 'ROOT_CERTIFIED':
            need(routes.replay(source, item.root_word)[-1] == 1, 'literal supplied ROOT')
    print('FINITE SOURCE UNIVERSE 3 MOD4 THROUGH511', json.dumps(sample, sort_keys=True))
    print('ITS FIRST-RESET-2 SUBSET', json.dumps(reset_two, sort_keys=True))
    for n, e in ((True, 1), (1, True), (2, 3), (1, 0), (1.0, 2)):
        need(rejects(lambda: partners(n, e)), 'typed DP input')
    need(rejects(lambda: choose_partner((2, True), 3)), 'typed core letters')
    need(rejects(lambda: choose_partner((2,), True)), 'typed available run')
    need(rejects(lambda: trim_root(True, (2,))), 'typed ROOT trimming source')
    need(rejects(lambda: collision_receipt(7, (2,), (2, 2))),
         'unequal-anchor rewrite rejected')
    need(rejects(lambda: collision_receipt(7, (8, 1), (4, 1, 1, 3))),
         'valid abstract collision on wrong source rejected')
    print('CACHE', _partners.cache_info())
    print('ALL CHECKS PASSED', CHECKS)


if __name__ == '__main__':
    main()
