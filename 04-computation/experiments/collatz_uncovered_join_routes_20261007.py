"""New guarded smaller dependencies relative to a frozen eight-macro policy.

Exact integer arithmetic; no import-time experiment and no ROOT oracle in any
production receipt/selector API. Run normally and with -O from repository root.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
from math import lcm

import adaptive_boundary_selector_20261004 as prior
import virtual_contraction_ladders_20261004 as virtual
import collatz_complement_routing_20261004 as legacy
import collatz_early_reroute_20261004 as early


CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def odd(n):
    need(type(n) is int and n > 0 and n % 2 == 1, 'exact positive odd source')


def letters(word):
    need(type(word) is tuple and all(type(a) is int and a >= 1 for a in word),
         'exact tuple of positive valuation letters')


def carrier(word):
    letters(word)
    p = q = 1
    b = 0
    for a in word:
        p, q, b = 3*p, q*2**a, 3*b+q
    return p, q, b


def step(n):
    odd(n)
    z, a = 3*n+1, 0
    while z % 2 == 0:
        z //= 2
        a += 1
    return z, a


def replay(n, word):
    odd(n)
    letters(word)
    states = [n]
    for a in word:
        need(n != 1, 'no first-hit ROOT padding')
        n, actual = step(n)
        need(actual == a, 'literal exponent agrees with receipt')
        states.append(n)
    return tuple(states)


@dataclass(frozen=True)
class Family:
    word: tuple
    siblings: int
    p: int
    c: int
    d: int
    least: int
    period: int


def compile_family(word, siblings):
    """All positive sources whose smaller word-predecessor reaches S^k(n)."""
    letters(word)
    need(bool(word), 'nonempty inverse word')
    need(type(siblings) is int and siblings >= 0, 'nonnegative sibling depth')
    p, q, b = carrier(word)
    c, d = q*4**siblings, q*((4**siblings-1)//3)-b
    need(c < p, 'strictly contracting inverse coefficient')
    residue = -d*pow(c, -1, p) % p
    if residue % 2 == 0:
        residue += p
    period = 2*p
    # Strict source>1, child>0, and child<source; equality is not payment.
    minimum = max(1, (-d)//c, d//(p-c))
    least = residue+period*max(0, (minimum-residue)//period+1)
    return Family(word, siblings, p, c, d, least, period)


def canonical(family):
    need(type(family) is Family, 'typed family')
    need(all(type(getattr(family, a)) is int for a in
             ('siblings', 'p', 'c', 'd', 'least', 'period')), 'exact family fields')
    need(family == compile_family(family.word, family.siblings), 'canonical family')


@dataclass(frozen=True)
class Receipt:
    source: int
    child: int
    source_word: tuple
    child_word: tuple
    endpoint: int


def audit(receipt):
    need(type(receipt) is Receipt, 'typed receipt')
    odd(receipt.source)
    odd(receipt.child)
    odd(receipt.endpoint)
    need(receipt.child < receipt.source, 'immutable-source strict size payment')
    for n, w in ((receipt.source, receipt.source_word),
                 (receipt.child, receipt.child_word)):
        p, q, b = carrier(w)
        need(p*n+b == q*receipt.endpoint, 'independent affine common endpoint')
        need(replay(n, w)[-1] == receipt.endpoint, 'literal common endpoint')
    return receipt


def apply_family(source, family):
    odd(source)
    canonical(family)
    if source < family.least or (source-family.least) % family.period:
        return None
    h = (family.c*source+family.d)//family.p
    if family.siblings:
        target, a = step(source)
        left, right = (a,), family.word+(a+2*family.siblings,)
    else:
        target, left, right = source, (), family.word
    return audit(Receipt(source, h, left, right, target))


F91 = compile_family((1, 1, 2, 2), 0)
F111 = compile_family((1, 1, 2, 1, 1, 1, 2), 1)
FAMILIES = (F91, F111)
FRONTIER_WORD = (1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2)


def frontier_receipt(source):
    """Inspect a saved budget frontier; do not discard an unpaid alternative."""
    odd(source)
    if source % 2**30 != 111:
        return None
    p, q, b = carrier(FRONTIER_WORD)
    x = (p*source+b)//q
    h = (x-1)//4
    need(x == 4*h+1 and h % 2 == 1 and 0 < h < source,
         'frontier sibling native guard and immutable payment')
    target, a = step(x)
    need(a >= 3, 'positive child final valuation')
    return audit(Receipt(source, h, FRONTIER_WORD+(a,), (a-2,), target))


def frontier_receipts(selection):
    """Inspect all smaller odd siblings of an authenticated final checkpoint.

    The extra source letter is authenticated by S^k(h)'s exact identity; this
    is one explicit additional receipt edge, not part of the old macro budget.
    """
    need(type(selection) is prior.Selection, 'typed authenticated selection')
    prior.audit_selection(selection)
    need(selection.status == 'PENDING', 'inspect a retained pending frontier')
    n, x = selection.source, selection.endpoint
    siblings = prior.siblings(x)
    if not siblings:
        return ()
    target, a = step(x)
    result = []
    for index, h in enumerate(siblings, 1):
        if h < n:
            need(a-2*index >= 1, 'positive actual child exponent')
            child_word = () if h == 1 else (a-2*index,)
            result.append(audit(Receipt(n, h, selection.word+(a,), child_word, target)))
    return tuple(result)


def new_receipts(source, include_frontier=False, selection=None):
    odd(source)
    need(type(include_frontier) is bool, 'explicit frontier option')
    result = [r for f in FAMILIES if (r := apply_family(source, f)) is not None]
    if include_frontier:
        if selection is None:
            selection = prior.select(source, 8, 'adaptive')
        need(selection.source == source, 'preserved final-frontier source')
        result.extend(frontier_receipts(selection))
    return tuple(sorted(result, key=lambda r: (r.child, len(r.source_word)+len(r.child_word))))


def discharge(receipt, supplied_child_word):
    """Only a supplied checked child ROOT word discharges a dependency."""
    audit(receipt)
    need(replay(receipt.child, supplied_child_word)[-1] == 1,
         'supplied child certificate actually reaches first ROOT')
    r = len(receipt.child_word)
    need(supplied_child_word[:r] == receipt.child_word, 'same actual child prefix')
    word = receipt.source_word+supplied_child_word[r:]
    need(replay(receipt.source, word)[-1] == 1, 'compiled first-hit source certificate')
    return word


def refined_select(source):
    """Preserve the prior result; add a dependency only after its budget ends."""
    odd(source)
    previous = prior.select(source, 8, 'adaptive', word_cap=128, phase_depth=2)
    if previous.status != 'PENDING':
        return previous, None
    options = new_receipts(source)
    found = options[0] if options else None
    return previous, found


@dataclass(frozen=True)
class TraceSeal:
    seed: int
    bits: int
    word: tuple
    # (seed checkpoint, cumulative cost, required current bits, macro kind)
    guards: tuple


def v2(n):
    need(type(n) is int and n > 0, 'positive exact valuation argument')
    return (n & -n).bit_length()-1


@lru_cache(None)
def trace_seal(seed):
    """Sufficient all-height certificate for the two declared pending traces.

    It checks finite branch predicates and the inequalities which must stay
    false. This is deliberately not an unsupported generic symbolic executor.
    """
    odd(seed)
    need(seed in (91, 111), 'two explicitly audited trace seeds')
    result = prior.select(seed, 8, 'adaptive', word_cap=128, phase_depth=2)
    need(result.status == 'PENDING', 'frozen policy exhausted eight macros')
    p = q = 1
    b = A = 0
    precision = sum(result.word)+1  # the terminal oddness bit is retained
    guards = []
    for macro in result.macros:
        x = macro.source
        need(p*seed+b == q*x, 'trace checkpoint affine identity')
        bits = 3
        # Freeze all171 membership flags, not just the ultimately chosen row.
        for row in prior.bank().values():
            diff = x-row['residue']
            depth = row['K'] if diff == 0 else min(row['K'], v2(abs(diff))+1)
            bits = max(bits, depth)
        for z in (x+1, x+5, 11*x+19, 3*x+1):
            bits = max(bits, v2(z)+1)
        need(step(x)[0] > 1, 'no to-ROOT branch at seed; positive lifts stay above it')
        # Stabilize virtual recognition at its native dyadic guard, not by phase.
        need(virtual.recognize(x) is None, 'no recognized virtual row at these checkpoints')
        j = v2(x+1)-1
        if j >= 3:
            s = ((3**j).bit_length()-j+1)//2
            if s >= 1:
                row = virtual.chart(s)
                if row.j == j:
                    need(x % row.modulus != row.r, 'native virtual guard fails')
                    bits = max(bits, row.modulus.bit_length()-1)
        # A non-paid sibling must not become paid at a large source lift.
        siblings = prior.siblings(x)
        bits = max(bits, 2*len(siblings)+3)
        for index, child in enumerate(siblings, 1):
            need(child >= seed and p >= q*4**index,
                 'each sibling minus original source has nonnegative slope and intercept at seed')
        guards.append((x, A, bits, macro.kind))
        precision = max(precision, A+bits)
        for a in macro.word:
            p, q, b = 3*p, q*2**a, 3*b+q
            A += a
            need(p > q, 'every actual source prefix has strictly expanding coefficient')
    need(A == sum(result.word) and len(guards) == 8, 'complete frozen trace')
    return TraceSeal(seed, precision, result.word, tuple(guards))


def pending_progression(family, seed):
    canonical(family)
    seal = trace_seal(seed)
    need(apply_family(seed, family) is not None, 'trace seed lies in new family')
    return seed, lcm(family.period, 2**seal.bits)


def root_word_control(n, cap=1000):
    """FINITE test witness only, never called by receipt production or selection."""
    odd(n)
    word = ()
    while n != 1 and len(word) < cap:
        n, a = step(n)
        word += (a,)
    need(n == 1, 'finite independent ROOT control completed')
    return word


def grounded_closure(seeds, refinements=0, suffix_memory=False, limit=32767):
    """One ascending pass; retain every advertised smaller alternative.

    seeds is an explicit mapping source->supplied first-hit word. The routine
    never obtains an unsupplied ROOT trajectory. Memory stores only suffixes
    of an already validated complete certificate, identically in both policies.
    """
    need(type(seeds) is dict and 1 in seeds and seeds[1] == (), 'explicit ROOT seed bank')
    need(type(refinements) is int and refinements in (0, 1, 2), 'declared rule set')
    need(type(suffix_memory) is bool and type(limit) is int and limit >= 1,
         'finite universe and memory policy')
    certificates, provenance = {}, {}

    def remember(source, word, origin):
        path = replay(source, word)
        need(path[-1] == 1, 'only fully checked ROOT words enter memory')
        cuts = range(len(path)) if suffix_memory else (0,)
        for cut in cuts:
            n, suffix = path[cut], word[cut:]
            if n in certificates:
                need(certificates[n] == suffix, 'unique actual first-hit word')
            else:
                certificates[n] = suffix
                provenance[n] = (origin, source, cut)

    for source, word in sorted(seeds.items()):
        remember(source, word, 'supplied')
    pending, unavailable, direct_new = [], [], []
    for n in range(3, limit+1, 2):
        if n in certificates:
            continue
        result = prior.select(n, 8, 'adaptive')
        if result.status == 'REDUCED':
            choices = [h for h in result.dependencies if h in certificates]
            if not choices:
                unavailable.append(n)
                continue
            h = min(choices)
            tail = certificates[h]
            if result.kind == 'common_future' and h != result.endpoint and h != 1:
                need(replay(h, tail[:1])[-1] == result.endpoint, 'advertised old common future')
                tail = tail[1:]
            remember(n, result.word+tail, 'prior')
        else:
            options = new_receipts(n, include_frontier=refinements == 2,
                                   selection=result) if refinements else ()
            choices = [r for r in options if r.child in certificates]
            if choices:
                receipt = choices[0]
                remember(n, discharge(receipt, certificates[receipt.child]), 'refinement')
                direct_new.append(n)
            elif options:
                unavailable.append(n)
            else:
                pending.append(n)
    grounded = frozenset(n for n in certificates if n <= limit)
    return grounded, tuple(pending), tuple(unavailable), tuple(direct_new), certificates, provenance


def rejects(callback):
    try:
        callback()
    except ValueError:
        return 1
    raise ValueError('hostile unexpectedly accepted')


def main():
    need((F91.p, F91.c, F91.d, F91.least, F91.period) == (81, 64, -73, 91, 162),
         'inverse1122 least algebraic family member')
    need((F111.p, F111.c, F111.d, F111.least, F111.period) ==
         (2187, 2048, -2067, 111, 4374), 'asynchronous sibling family constants')
    # F91's positive paid least member is91, not an assumed seed from a ROOT bank.
    need(F91.least == 91, 'positive paid F91 boundary')
    for family in FAMILIES:
        for t in range(256):
            n = family.least+family.period*t
            receipt = apply_family(n, family)
            need(receipt is not None and receipt.child < n, 'whole-family finite replay')
            need(apply_family(n+2, family) is None, 'adjacent wrong-guard source rejected')
    for t in range(128):
        n = 111+2**30*t
        receipt = frontier_receipt(n)
        need(receipt.child == 81+774840978*t, 'all-height frontier child formula')
        need(prior.select(n, 8, 'adaptive').status == 'PENDING', 'frontier row lies in sealed pending cylinder')
        if t % 2:
            need(receipt.source_word[-1] == 3 and all(y > n for y in replay(n, receipt.source_word)[1:]),
                 'odd-parameter half joins while every source prefix still grows')
    for family, seed in ((F91, 91), (F111, 111)):
        seal = trace_seal(seed)
        n0, period = pending_progression(family, seed)
        for t in (0, 1, 2, 3, 7, 31, 1024, 10**30):
            n = n0+period*t
            result, receipt = refined_select(n)
            need(result.status == 'PENDING' and result.word == seal.word,
                 'independent whole-policy trace replay on positive lifts')
            need(receipt is not None and receipt.child < n, 'new all-height progression payment')
        print('PROVED pending cylinder', seed, 'mod2^'+str(seal.bits),
              '; macro guards', seal.guards)
        print('PROVED added AP:', n0, '+', period, '*t; t>=0.')
    # Explicit finite baseline comparison. No discovery ROOT trajectories are inputs.
    pending, added, still, augmented_pending, frontier_added = [], [], [], [], []
    frontier_probes = 0
    bank = legacy.load_bank()
    for n in range(1, 32768, 2):
        old, new = refined_select(n)
        if old.status == 'PENDING':
            pending.append(n)
            (added if new is not None else still).append(n)
            frontier_probes += bool(prior.siblings(old.endpoint))
            if frontier_receipts(old):
                frontier_added.append(n)
            if (legacy.old_partition(n, bank) is None and not early.candidates(n)
                    and n % 2048 != 155 and n % 256 != 219):
                augmented_pending.append(n)
    need(len(pending) == 243 and added == [91, 111, 8191] and len(still) == 240,
         'frozen odd1..32767 coverage comparison')
    need(all(legacy.old_partition(n, bank) is None and not early.candidates(n)
             for n in added), 'three finite additions also escape old binary/ternary and early banks')
    need(len(augmented_pending) == 116 and all(n in augmented_pending for n in added),
         'stronger finite baseline comparison')
    for n in added:
        receipt = refined_select(n)[1]
        supplied = root_word_control(receipt.child)
        compiled = discharge(receipt, supplied)
        need(compiled == root_word_control(n), 'independent completed finite receipt control')
    # Counter-controls: the user meant332, which halves twice to83, not322.
    need(332 == 4*83 and prior.select(83, 8, 'adaptive').status == 'REDUCED',
         '332 distinguished from historical322')
    need(step(425) == (319, 2), '425 already has a one-step smaller endpoint')
    need(refined_select(27)[1] is None and refined_select(703)[1] is None,
         'declared residuals remain unresolved by these two additions')
    need(all(prior.select(257727+1048576*t, 8, 'adaptive').status == 'REDUCED'
             for t in range(256)), 'concurrent permutation family already covered in this finite control')
    r91 = apply_family(91, F91)
    hostiles = [lambda: apply_family(True, F91), lambda: apply_family(91.0, F91),
                lambda: compile_family((1, 0), 0),
                lambda: apply_family(91, replace(F91, d=-71)),
                lambda: audit(replace(r91, child=73)),
                lambda: audit(replace(r91, source=71, child=91)),
                lambda: discharge(r91, (2,)),
                lambda: replay(1, (2,)), lambda: trace_seal(703)]
    need(sum(rejects(f) for f in hostiles) == 9, 'typed, forged, cyclic and ungrounded hostiles')
    root_boundary = prior.finish(3, 'PENDING', 'budget', 5, (),
                                 (prior.macro('control', 3, (1,)),), ())
    root_receipt = frontier_receipts(root_boundary)[0]
    need(root_receipt.child == 1 and root_receipt.child_word == ()
         and discharge(root_receipt, ()) == (1, 4), 'frontier ROOT is an empty supplied suffix')
    # Directional hostile:71 reaches91 but91 is not smaller than71.
    need(replay(71, (1, 1, 2, 2))[-1] == 91, 'orientation must retain original source')
    print('FINITE-EXACT odd1..32767: priorROOT/reduced', 16384-len(pending),
          '; priorPENDING', len(pending), '; new dependencies', added,
          '; remainingPENDING', len(still), '; first', still[:12])
    print('Finite additions also fail old_partition and the complete early-reroute selector.')
    print('Augmented finite baseline plus old_partition, early.candidates and H/L:116 pending ->113.')
    print('Final-boundary inspection:pending frontiers with siblings', frontier_probes,
          '; paid sibling frontiers', len(frontier_added),
          '; union with two rows', len(set(added)|set(frontier_added)),
          '; paid frontier sources', frontier_added)
    root_only = {1: ()}
    supplied27 = root_word_control(27)
    bank27 = {1: (), 27: supplied27}
    closures = []
    for seeds, memory in ((root_only, False), (bank27, False), (bank27, True)):
        runs = [grounded_closure(seeds, j, memory) for j in (0, 1, 2)]
        need(all(runs[0][0] <= row[0] for row in runs[1:]), 'monotone closure comparison')
        counts = tuple(len(row[0]) for row in runs)
        closures.append(counts)
        print('GROUNDED ascending closure seeds', tuple(seeds), 'suffix_memory', memory,
              '; baseline / two rows / plus frontier', counts,
              '; final new sources', sorted(runs[2][0]-runs[0][0])[:20])
        if memory:
            # These are already authenticated derived words, not new assumed
            # seeds. Preserve the original1/27 bank as their proof provenance.
            for j, first in enumerate(runs):
                second = grounded_closure(dict(first[4]), j, True)
                need(second[4] == first[4], 'second pass preserves the entire source->ROOT-word dictionary')
            sizes = tuple(len(first[4]) for first in runs)
            need(sizes == (22964, 22965, 23418), 'full authenticated cache sizes include outside-universe suffixes')
            print('FINITE fixed point:second pass leaves entire certificate dictionaries unchanged;', sizes)
    print('GROUNDED comparison tuples:', closures)
    need(closures == [(10343, 10343, 10474), (11290, 11290, 11423), (15257, 15258, 15463)],
         'exact fixed-seed grounded comparisons')
    need(frontier_receipt(111).child == 81, 'ROOT-only new completion uses independently grounded81')
    print('Supplied27 certificate used identically:', supplied27)
    print('Family replay512; exact pending-lift controls16; supplied-child discharges3; hostiles9.')
    print('332 has primitive odd part83;425 descends at step1. Neither is promoted as new coverage.')
    print('Permutation-family controls t0..255 all already reduced by this baseline.')
    print('No unconditional ROOT claim for arbitrary child; no universal guard coverage.')
    print('CHECKS', CHECKS)


if __name__ == '__main__':
    main()
