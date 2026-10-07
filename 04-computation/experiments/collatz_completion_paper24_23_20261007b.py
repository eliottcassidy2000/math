"""Turn a proposed completed deletion into an exact missing-rule obligation.

Paper mechanisms are inspiration, not imports into the elementary proofs.
Production consumes authenticated actual words; it never discovers ROOT.
"""
from dataclasses import dataclass, replace
from itertools import product
from fractions import Fraction

import collatz_child_normal_forms_20261007 as forms
import collatz_uncovered_join_routes_20261007 as routes


def nonempty(word):
    forms.letters(word)
    if not word:
        raise ValueError("a nonempty generator word is required")
    return word


def strict_depth(numerator, denominator):
    if type(numerator) is not int or type(denominator) is not int or not numerator > denominator > 0:
        raise ValueError("exact positive ratio greater than one required")
    return (numerator // denominator).bit_length()


@dataclass(frozen=True)
class Request:
    parent: int
    word: tuple
    source: int
    minimum: int
    maximum: int


def request(parent, word):
    forms.odd(parent)
    nonempty(word)
    source = forms.apply(forms.encode(word), parent)
    return Request(parent, word, source,
                   strict_depth(parent + 1, source + 1),
                   forms.valuation(parent + 1, 2) - 1)


def audit_request(packet):
    if type(packet) is not Request:
        raise ValueError("exact Request required")
    for x in (packet.parent, packet.source):
        forms.odd(x)
    for x in (packet.minimum, packet.maximum):
        forms.natural(x)
    if request(packet.parent, packet.word) != packet:
        raise ValueError("request does not retain its actual source and word")
    return packet


def deleted_child(parent, depth):
    forms.odd(parent)
    if type(depth) is not int or not 1 <= depth <= forms.valuation(parent + 1, 2) - 1:
        raise ValueError("deletion does not give a positive odd child")
    return (parent + 1) // (1 << depth) - 1


def permitted_depths(packet):
    audit_request(packet)
    return tuple(range(packet.minimum, packet.maximum + 1))


@dataclass(frozen=True)
class UniformBudget:
    word: tuple
    least_parent: int
    least_source: int
    required: int
    length_bound: int
    cost_bound: int


def uniform_budget(word):
    nonempty(word)
    record = forms.encode(word)
    residue, p = forms.native_cell(record)
    parent = residue if residue % 2 else residue + p
    source = forms.apply(record, parent)
    length = len(word)
    dlength = 1
    while (1 << (dlength + length)) <= 3**length:
        dlength += 1
    dcost = 1
    while (1 << (dcost + record.A)) < 5*p:
        dcost += 1
    return UniformBudget(word, parent, source,
                         strict_depth(parent + 1, source + 1), dlength, dcost)


def consume_deletion(packet, supplied):
    """Accept a real parent-to-deleted-child common future, then transport it."""
    audit_request(packet)
    routes.audit(supplied)
    if supplied.source != packet.parent:
        raise ValueError("certificate belongs to another parent")
    quotient, remainder = divmod(packet.parent + 1, supplied.child + 1)
    if remainder or quotient < 2 or quotient & (quotient - 1):
        raise ValueError("supplied smaller child is not a bit deletion")
    depth = quotient.bit_length() - 1
    if depth not in permitted_depths(packet):
        raise ValueError("deletion is unavailable or fails original-source payment")
    if deleted_child(packet.parent, depth) != supplied.child:
        raise ValueError("deleted-child identity failed")
    prefix = forms.forward_word(packet.word)
    if forms.replay(packet.source, prefix) != packet.parent:
        raise ValueError("inverse return word changed the source")
    return routes.audit(routes.Receipt(packet.source, supplied.child,
                                       prefix + supplied.source_word,
                                       supplied.child_word, supplied.endpoint))


def consume_parent_root(packet, parent_word):
    """Stronger premise: an actual parent ROOT proof grounds every legal inverse."""
    audit_request(packet)
    forms.replay(packet.parent, parent_word, root=True)
    result = forms.forward_word(packet.word) + parent_word
    forms.replay(packet.source, result, root=True)
    return result


def main():
    checks = 0

    def check(ok):
        nonlocal checks
        checks += 1
        if not ok:
            raise ArithmeticError("exact control failed")

    def rejects(callback):
        try:
            callback()
        except (ValueError, TypeError, AttributeError):
            check(True)
        else:
            check(False)

    words = [w for ell in range(1, 6) for w in product((1, 5, 17), repeat=ell)]
    native_pairs = 0
    unavailable = feasible = 0
    for word in words:
        budget = uniform_budget(word)
        record = forms.encode(word)
        p, q = 3**record.B, 1 << record.A
        check(p-q <= record.C <= 17*(p-q))
        check(budget.required <= min(budget.length_bound, budget.cost_bound))
        r0 = Fraction(budget.least_parent+1, budget.least_source+1)
        check((1 << budget.required) > r0)
        check((1 << (budget.required-1)) <= r0)
        check(r0 <= Fraction(3, 2)**len(word))
        check(r0 < Fraction(5*p, q))
        previous = r0
        for t in range(17):
            m = budget.least_parent + 2*p*t
            h = budget.least_source + 2*q*t
            packet = request(m, word)
            check(packet.source == h)
            ratio = Fraction(m+1, h+1)
            check(ratio <= previous)
            check(packet.minimum <= budget.required)
            check(forms.replay(h, forms.forward_word(word)) == m)
            previous = ratio
            actual = []
            for D in range(1, packet.maximum+1):
                z = deleted_child(m, D)
                check(z > 0 and z % 2 == 1)
                if z < h:
                    actual.append(D)
                check((z < h) == ((1 << D)*(h+1) > m+1))
            check(tuple(actual) == permitted_depths(packet))
            oddpart = (m+1) >> forms.valuation(m+1, 2)
            check(bool(actual) == (2*oddpart-1 < h))
            if actual:
                feasible += 1
            else:
                unavailable += 1
            native_pairs += 1

    # The former fixed-six-bit failures are numerically restored at depth seven.
    b5, b17 = uniform_budget((5,)*36), uniform_budget((17,)*64)
    check(b5.required == b17.required == 7)
    for budget in (b5, b17):
        check(64*(budget.least_source+1) <= budget.least_parent+1)
        check(128*(budget.least_source+1) > budget.least_parent+1)

    # Every legal Mersenne mixed child has a nonempty arithmetic request interval.
    mersenne_cases = 0
    for word in words:
        phase = forms.mersenne_phase(forms.encode(word))
        if phase is None or phase[0] > 4096:
            continue
        E = phase[0]
        packet = request((1 << E)-1, word)
        check(packet.minimum <= packet.maximum == E-1)
        check(deleted_child(packet.parent, E-1) == 1 < packet.source)
        mersenne_cases += 1

    positive = request(11, (1,))
    check((positive.source, positive.minimum, positive.maximum) == (7, 1, 1))
    parent_deletion = routes.Receipt(11, 5, (1, 2, 3), (), 5)
    transferred = consume_deletion(positive, parent_deletion)
    check((transferred.source, transferred.child, transferred.endpoint) == (7, 5, 5))
    check(transferred.source_word == (1, 1, 2, 3))
    check(forms.replay(7, transferred.source_word) == 5)

    hostile = request(35, (1, 1))
    check((hostile.source, hostile.minimum, hostile.maximum) == (15, 2, 1))
    check(permitted_depths(hostile) == ())
    check(deleted_child(35, 1) == 17 > 15)
    authentic_unpaid = routes.Receipt(35, 17, (1, 5), (2, 3), 5)
    routes.audit(authentic_unpaid)
    rejects(lambda: consume_deletion(hostile, authentic_unpaid))
    # A stronger authentic ROOT premise still succeeds at the identical source.
    check(consume_parent_root(hostile, (1, 5, 4)) == (1, 1, 1, 5, 4))
    check(forms.replay(15, (1, 1, 1, 5, 4), root=True) == 1)

    for bad in (True, 11.0, 0, 2, -1):
        rejects(lambda bad=bad: request(bad, (1,)))
    for bad in ((), (True,), (1.0,), [1], (3,)):
        rejects(lambda bad=bad: uniform_budget(bad))
    rejects(lambda: audit_request(replace(positive, minimum=True)))
    rejects(lambda: audit_request(replace(positive, source=7.0)))
    rejects(lambda: audit_request(replace(positive, maximum=2)))
    rejects(lambda: deleted_child(35, 2))
    rejects(lambda: consume_deletion(positive, True))
    rejects(lambda: consume_deletion(positive, routes.Receipt(11, 5, (1,), (), 5)))
    rejects(lambda: consume_deletion(positive, authentic_unpaid))
    rejects(lambda: consume_parent_root(hostile, (1, 5)))
    rejects(lambda: consume_parent_root(hostile, (1, 5, 4, 2)))

    print("PROVED: exact deletion-request interval and fixed-native-fibre uniform budget.")
    print("Trusted PDF24/23 headlines are premises for mechanisms, not Collatz theorem inputs.")
    print("Word universe: all 363 nonempty generator words through length 5; 17 native pairs each.")
    print("Native pairs", native_pairs, "feasible arithmetic requests", feasible,
          "unavailable", unavailable)
    print("Mersenne literal phase controls", mersenne_cases, "with exponents at most 4096.")
    print("Uniform numerical budget: G5^36 and G17^64 both need 7 bits; 6 fails at the least native pair.")
    print("Positive certificate: parent11->5 gives mixed source7->5, actual word(1,1,2,3).")
    print("Hostile: parent 35 via G1^2 gives 15; required D2 exceeds available D1, target 17.")
    print("Stronger supplied parent ROOT word (1,5,4) still grounds 15; no source replacement or ROOT discovery.")
    print("Completion assumption alone supplies no missing common-future word.")
    print("Exact checks", checks)


if __name__ == '__main__':
    main()
