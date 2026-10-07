"""Exact terminal-word lifting and a variable-length reset receipt.

Finite controls are separate from production certificate transport.
Run from the repository root; no third-party dependencies.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import product

import collatz_uncovered_join_routes_20261007 as routes

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def word_guard(word):
    need(type(word) is tuple and bool(word), 'nonempty tuple word')
    need(all(type(a) is int and a >= 1 for a in word), 'positive exact-int letters')


@dataclass(frozen=True)
class DescentRay:
    word: tuple
    p: int
    q: int
    b: int
    residue: int
    period: int
    least: int
    threshold: Fraction


def compile_descent(word):
    """Exact native cylinder, strict original-source payment, no early ROOT."""
    word_guard(word)
    p = q = 1
    b = 0
    cuts = [Fraction(1)]
    for index, a in enumerate(word):
        p, q, b = 3*p, q*2**a, 3*b+q
        if index+1 < len(word):
            cuts.append(Fraction(q-b, p))
    need(q > p, 'positive-carry descent requires contraction')
    cuts.append(Fraction(b, q-p))
    threshold = max(cuts)
    period = 2*q
    residue = (q-b)*pow(p, -1, period) % period
    least = residue
    if least <= threshold:
        least += ((threshold-least)//period+1)*period
    need(type(least) is int and least > threshold and least % 2 == 1,
         'least integer in exact cylinder above every strict cut')
    endpoint = (p*least+b)//q
    need(routes.replay(least, word)[-1] == endpoint < least,
         'independent first-hit boundary replay')
    return DescentRay(word, p, q, b, residue, period, least, threshold)


def audit_ray(ray):
    need(type(ray) is DescentRay, 'typed descent ray')
    need(all(type(getattr(ray, k)) is int for k in ('p', 'q', 'b', 'residue', 'period', 'least')),
         'exact integer metadata')
    need(type(ray.threshold) is Fraction, 'exact rational threshold')
    need(compile_descent(ray.word) == ray, 'canonical recomputation')
    return ray


def apply_ray(ray, source):
    routes.odd(source)
    audit_ray(ray)
    if source < ray.least or source % ray.period != ray.residue:
        return None
    child = (ray.p*source+ray.b)//ray.q
    need(0 < child < source and child % 2 == 1, 'strict paid child')
    return routes.audit(routes.Receipt(source, child, ray.word, (), child))


def lift_terminal(source, supplied_root_word, first_descent=False):
    routes.odd(source)
    need(source > 1 and type(first_descent) is bool, 'nonroot supplied terminal')
    word_guard(supplied_root_word)
    path = routes.replay(source, supplied_root_word)
    need(path[-1] == 1, 'only supplied and replayed complete certificates')
    word = supplied_root_word
    if first_descent:
        cut = next(i for i in range(1, len(path)) if path[i] < source)
        word = word[:cut]
    row = compile_descent(word)
    need(apply_ray(row, source) is not None, 'terminal belongs to its certified ray')
    return row


@dataclass(frozen=True)
class CompletedLift:
    seed: int
    word: tuple
    index: int


def audit_completed(item):
    need(type(item) is CompletedLift, 'typed symbolic completed lift')
    routes.odd(item.seed)
    word_guard(item.word)
    need(item.seed > 1 and type(item.index) is int and item.index >= 0,
         'nonroot seed and exact nonnegative parameter')
    need(routes.replay(item.seed, item.word)[-1] == 1, 'supplied base certificate')
    return routes.carrier(item.word)


def completed_word(item):
    p, _, _ = audit_completed(item)
    return item.word if item.index == 0 else item.word+(2*(1+p*item.index),)


def completed_residue(item, modulus):
    """Evaluate the exact integer source modulo m without expanding it."""
    p, q, _ = audit_completed(item)
    need(type(modulus) is int and modulus >= 1, 'positive exact modulus')
    denominator = 3*p
    numerator = pow(4, p*item.index, denominator*modulus)-1
    need(numerator % denominator == 0, 'retained denominator precision')
    return (item.seed+4*q*(numerator//denominator)) % modulus


def expand_completed(item, bit_cap=10000):
    """Optional bounded realization; symbolic proof does not require expansion."""
    p, q, _ = audit_completed(item)
    need(type(bit_cap) is int and bit_cap >= 1, 'positive explicit bit cap')
    if item.index == 0:
        need(item.seed.bit_length() <= bit_cap, 'base source fits cap')
        return item.seed
    need(2*p*item.index+q.bit_length()+3 <= bit_cap, 'conservative expansion cap')
    numerator = 4**(p*item.index)-1
    need(numerator % (3*p) == 0, 'exact lifted source integer')
    return item.seed+4*q*(numerator//(3*p))


def recognize_completed(seed, supplied_root_word, source):
    """Recognize this explicit family, not an arbitrary ROOT basin."""
    p, q, b = audit_completed(CompletedLift(seed, supplied_root_word, 0))
    routes.odd(source)
    numerator = p*source+b
    if numerator % q:
        return None
    h = numerator//q
    power = 3*h+1
    if power < 4 or power & (power-1):
        return None
    exponent = power.bit_length()-1
    if exponent % 2:
        return None
    j = exponent//2
    if (j-1) % p:
        return None
    return CompletedLift(seed, supplied_root_word, (j-1)//p)


def reset_receipt(source):
    """Known reset switch, with variable-length literal receipts and ROOT guard.

    THM-4555/4556 supply prior mechanism. No terminal oracle is called here.
    """
    routes.odd(source)
    if source == 1:
        return None
    k = routes.v2(source+1)
    if k < 2:
        return None
    t = (source+1)//2**k
    numerator = 3**k*t-1
    r = routes.v2(numerator)
    if r < 2:
        return None
    child = (source-1)//2
    source_word = (1,)*(k-1)+(r+1,)
    child_word = () if child == 1 else (1,)*(k-2)+(2, r-1)
    endpoint = numerator//2**r
    return routes.audit(routes.Receipt(source, child, source_word, child_word, endpoint))


def sporadic_receipt(source):
    """Known collision (2,6,c)~(4,1,1,c+2), retaining variable run length."""
    routes.odd(source)
    if source == 1:
        return None
    k = routes.v2(source+1)
    if k < 2:
        return None
    t = (source+1)//2**k
    before_reset = 2*3**(k-1)*t-1
    x, a = routes.step(before_reset)
    if a != 2 or x == 1:
        return None
    y, b = routes.step(x)
    if b != 6 or y == 1:
        return None
    endpoint, c = routes.step(y)
    return routes.audit(routes.Receipt(source, (source-1)//2,
                        (1,)*(k-1)+(2, 6, c),
                        (1,)*(k-2)+(4, 1, 1, c+2), endpoint))


def rejects(callback):
    try:
        callback()
    except (ValueError, RuntimeError):
        return True
    return False


def first_descent_control(source, cap=1000):
    """Explicit finite discovery control, not used in production compilers."""
    n = source
    word = ()
    for _ in range(cap):
        n, a = routes.step(n)
        word += (a,)
        if n < source:
            return word
    raise ValueError('finite discovery cap exhausted')


def main():
    # Exhaustive declared native-cylinder controls, including ROOT-padded words.
    rows = 0
    for length in range(1, 5):
        for word in product(range(1, 5), repeat=length):
            p, q, b = routes.carrier(word)
            if q <= p:
                continue
            row = compile_descent(word)
            rows += 1
            for source in range(3, 1024, 2):
                native = source >= row.least and source % row.period == row.residue
                try:
                    path = routes.replay(source, word)
                    actual = path[-1] < source
                except ValueError:
                    actual = False
                need(native == actual, 'exact iff between guarded ray and first-hit descent')
            for t in (0, 1, 2, 10, 10**20):
                source = row.least+row.period*t
                receipt = apply_ray(row, source)
                need(receipt.child == (row.p*row.least+row.b)//row.q+2*row.p*t,
                     'all-height affine child lift')
    padded = compile_descent((4, 2))
    need(padded.residue == 5 and padded.least == 133,
         'native source5 excluded by first-hit boundary,133 retained')
    need(apply_ray(padded, 5) is None, 'ROOT padding is not a proof extension')

    terminals = (3, 27, 91, 111, 223, 4591, 32767)
    print('TERMINAL RAYS source,ROOT length/cost,descent length/cost,least')
    for n in terminals:
        root = routes.root_word_control(n)
        full = lift_terminal(n, root)
        short = lift_terminal(n, root, True)
        need(full.least == n and (full.p*n+full.b)//full.q == 1,
             'terminal full-word ray begins at the supplied terminal')
        for t in (0, 1, 7, 10**10):
            receipt = apply_ray(full, n+full.period*t)
            need(receipt.child == 1+2*full.p*t, 'full terminal-word lift formula')
        print(n, (len(full.word), sum(full.word)),
              (len(short.word), sum(short.word)), short.least)

    pending = []
    for source in range(3, 32768, 2):
        receipt = reset_receipt(source)
        k = routes.v2(source+1)
        t = (source+1)//2**k
        expected = k >= 2 and (3**k*t-1) % 4 == 0
        need((receipt is not None) == expected, 'exact reset native guard')
        if receipt is None:
            continue
        # Both the finite child certificate and reconstructed source are explicit controls.
        child_word = routes.root_word_control(receipt.child)
        compiled = routes.discharge(receipt, child_word)
        need(compiled == routes.root_word_control(source), 'first-hit reset substitution')
        selection = routes.prior.select(source, 8, 'adaptive')
        if selection.status == 'PENDING' and not routes.new_receipts(source, True, selection):
            pending.append(source)
    need(len(pending) == 22, 'finite relative coverage delta')
    print('RESET new local reductions beyond the prior31:', pending)

    # A genuinely infinite cap obstruction. Each finite member is only a control;
    # the proof for all multipliers is in the accompanying note.
    for t in (1, 2, 3, 5):
        k = 1458*t
        source = 2**k-1
        selection = routes.prior.select(source, 8, 'adaptive')
        need(selection.status == 'PENDING' and selection.word == (1,)*1024,
             'old8-by128 controller spends its whole budget on expanding ones')
        need(not routes.new_receipts(source, True, selection), 'old inverse and frontier guards fail')
        receipt = reset_receipt(source)
        need(receipt.child == 2**(k-1)-1 and len(receipt.source_word) == k,
             'variable-length reset pays an uncovered source')
        print('CAP ESCAPE exponent', k, 'old word length', len(selection.word),
              'new receipt lengths', len(receipt.source_word), len(receipt.child_word))

    need(pow(2, 1458, 2187) == 1, 'exact congruence underlying all cap-escape rays')
    # A second infinite escape, now in the first-reset2 branch.
    for s in (0,):
        t = 7+32*s
        k = 1458*t+1
        source = 2**k-1
        need(k % 64 == 31 and reset_receipt(source) is None, 'reset2 rather than root switch')
        selection = routes.prior.select(source, 8, 'adaptive')
        need(selection.status == 'PENDING' and not routes.new_receipts(source, True, selection),
             'second infinite family misses the same frozen controller')
        receipt = sporadic_receipt(source)
        need(receipt is not None and receipt.child == 2**(k-1)-1,
             'sporadic collision pays a previously missing reset2 branch')
        need(reset_receipt(receipt.child).child == 2**(k-2)-1,
             'sporadic/root-switch composition decreases the odd exponent by two')
        print('SPORADIC CAP ESCAPE exponent', k, 'receipt length', len(receipt.source_word))
    need(reset_receipt(2**1459-1) is None and sporadic_receipt(2**1459-1) is None,
         'explicit remaining branch, not universal completion')

    hostiles = [
        lambda: compile_descent(()), lambda: compile_descent((True, 4)),
        lambda: compile_descent((0,)), lambda: compile_descent((1,)),
        lambda: lift_terminal(27, (1, 2)), lambda: lift_terminal(1, (2,)),
        lambda: lift_terminal(3, (1, 4, 2)),
        lambda: audit_ray(replace(padded, least=5)),
        lambda: audit_ray(replace(padded, period=64)),
        lambda: reset_receipt(True),
    ]
    need(all(rejects(f) for f in hostiles), 'typed, forged, ungrounded and ROOT-padding controls')
    need(reset_receipt(7) is None and reset_receipt(27) is None,
         'first-reset2 remains an explicit missing branch')
    need(reset_receipt(3).child_word == (), 'ROOT terminal is not padded')

    # Lift only words compiled from the finite terminal bank; no new discovery.
    import collatz_terminal_basis_20261007 as basis
    graph, _, _ = basis.build_observations()
    demand, cycles = basis.frontier(graph)
    need(not cycles, 'inherited finite obstruction has no cycles')
    _, certificates = basis.complete(graph, basis.TERMINAL_WORDS)
    representatives = sorted(min(group) for group in demand.values())
    full_rows = [lift_terminal(n, certificates[n]) for n in representatives]
    short_rows = [lift_terminal(n, certificates[n], True) for n in representatives]
    need(len(representatives) == 37, 'all37 missing components lifted')
    for n, full, short in zip(representatives, full_rows, short_rows):
        need(full.period % short.period == 0 and full.residue % short.period == short.residue,
             'shorter descent guard contains its whole terminal-word cylinder')
        for t in (0, 1, 3, 10**20):
            need(apply_ray(short, short.least+short.period*t) is not None,
                 'lifted missing-component rule at arbitrary parameter scale')
    print('GROUNDED REPRESENTATIVE LIFTS', representatives)
    print('Full/short word costs:', [(n, sum(a.word), sum(b.word))
                                    for n, a, b in zip(representatives, full_rows, short_rows)])
    print('Guard-width exponent savings:', min(sum(a.word)-sum(b.word) for a,b in zip(full_rows,short_rows)),
          max(sum(a.word)-sum(b.word) for a,b in zip(full_rows,short_rows)))

    # Fully grounded subfamilies, not merely arbitrary smaller-child rays.
    for seed in (3, 5, 7, 13):
        base_word = routes.root_word_control(seed)
        for s in range(4):
            item = CompletedLift(seed, base_word, s)
            n = expand_completed(item)
            word = completed_word(item)
            need(routes.replay(n, word)[-1] == 1, 'literal first-hit completed lift')
            need(recognize_completed(seed, base_word, n) == item, 'exact family recognition round trip')
            need(recognize_completed(seed, base_word, n+2) is None, 'adjacent odd nonmember rejected')
            for m in (1, 3, 8, 19, 105, 2187):
                need(completed_residue(item, m) == n % m, 'symbolic/literal source residues agree')
    for seed in representatives:
        item = CompletedLift(seed, certificates[seed], 1)
        p, q, b = audit_completed(item)
        for m in (3, 19, 105, 2187, 2**32, 2*q):
            nmod = completed_residue(item, m)
            hmod = (pow(4, 1+p, 3*m)-1)//3
            need((p*nmod+b-q*hmod) % m == 0, 'symbolic parent/child carrier congruence')
        need(completed_residue(item, 2*q) == seed,
             'completed subfamily lies inside the original exact source guard')
        need(rejects(lambda: expand_completed(item, 10000)), 'large source stays symbolic')
    need(rejects(lambda: completed_word(CompletedLift(3, (1, 4), True))), 'exact symbolic index')
    need(rejects(lambda: completed_residue(CompletedLift(3, (1, 4), 1), True)), 'exact symbolic modulus')
    need(completed_word(CompletedLift(3, (1, 4), 0)) == (1, 4), 's0 never appends a ROOT loop')
    print('COMPLETED SUBFAMILIES all37 representatives; exact modular verification without expansion.')
    print('CONTRACTING WORDS alphabet1..4 length1..4:', rows)
    print('No universal coverage: first-reset2 branch and grounding of arbitrary children remain open.')
    print('CHECKS', CHECKS)


if __name__ == '__main__':
    main()
