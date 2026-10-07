"""Guard arbitrary authenticated two-anchor heads on Mersenne exponent phases.

Finite cutoff proves first-hit legality without discovering any ROOT orbit.
The separate decoder is used only to instantiate a declared finite head bank.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product

import collatz_completion_anchor_20261007b as base


CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def integer(n, least=0):
    need(type(n) is int and n >= least, 'exact integer in declared domain')
    return n


@dataclass(frozen=True)
class Head:
    source: tuple
    partner: tuple
    gap: int


def compile_head(source, partner, gap=1):
    base.letters(source)
    base.letters(partner)
    integer(gap, 1)
    need(bool(source), 'nonempty source head')
    need(len(partner) == len(source)+3, 'three additional partner odd edges')
    need(sum(partner) == sum(source)-2*gap, 'exact sibling-gap valuation cost')
    left, right = base.affine(source, 1), base.affine(partner, 1)
    need(right == 4**gap*left+F(4**gap-1, 3), 'full anchored affine identity')
    return Head(source, partner, gap)


def audit_head(head):
    need(type(head) is Head, 'exact Head object')
    need(head == compile_head(head.source, head.partner, head.gap),
         'authenticated head fields')
    return head


def native_cell(head):
    audit_head(head)
    p, q, b = base.carrier(head.source)
    return (q-b)*pow(p, -1, 2*q) % (2*q), 2*q


@dataclass(frozen=True)
class Phase:
    head: Head
    twos: int
    residue: int
    period: int
    first_parameter: int


def phase(head, twos):
    audit_head(head)
    integer(twos, 3)
    need(head.source[0] != 2, 'head begins after the maximal two-run')
    A = sum(head.source)
    x, _ = native_cell(head)
    bits = 2*twos+A
    modulus = 1 << bits
    target = (1+(1 << (2*twos-1))*(x-1)*pow(3, -twos, modulus)) % modulus
    e, period = base.log_three(target, bits)
    k = e+1
    shell = 2*twos-2 if head.source[0] == 1 else 2*twos-1
    need(base.v2(e) == shell, 'native first letter forces exact run shell')
    need(period == 1 << (2*twos+A-2), 'minimal native exponent period')
    minimum_K = 2*A+twos+1
    first = max(0, (minimum_K-k+period-1)//period)
    return Phase(head, twos, k, period, first)


def audit_phase(record):
    need(type(record) is Phase, 'exact Phase object')
    audit_head(record.head)
    for x in (record.twos, record.residue, record.period, record.first_parameter):
        integer(x)
    need(record == phase(record.head, record.twos), 'canonical phase and cutoff')
    return record


def exponent(record, parameter):
    audit_phase(record)
    integer(parameter)
    need(parameter >= record.first_parameter, 'source-growth cutoff retained')
    return record.residue+record.period*parameter


def contains(record, K):
    audit_phase(record)
    integer(K, 1)
    return (K >= 2*sum(record.head.source)+record.twos+1
            and (K-record.residue) % record.period == 0)


def native_guard(head, twos, K):
    """Independent actual-head residue test with enough division precision."""
    audit_head(head)
    integer(twos, 3)
    integer(K, 2)
    need(head.source[0] != 2, 'nonpadding head')
    A = sum(head.source)
    divisor = 1 << (2*twos-1)
    modulus = divisor*(1 << (A+1))
    numerator = (pow(3, twos, modulus)*(pow(3, K-1, modulus)-1)) % modulus
    if numerator % divisor:
        return False
    X = (1+numerator//divisor) % (1 << (A+1))
    return X == native_cell(head)[0]


def source_residue(record, parameter, modulus, deletion=0):
    K = exponent(record, parameter)
    integer(modulus, 1)
    need(type(deletion) is int and deletion in (0, 3, 4), 'source or declared deletion')
    return (pow(2, K-deletion, modulus)-1) % modulus


def materialize(record, parameter, deletion=4, bit_cap=40000):
    K = exponent(record, parameter)
    integer(bit_cap, 1)
    need(K <= bit_cap, 'explicit source-bit materialization cap')
    need(type(deletion) is int and deletion in (3, 4), 'authenticated clearing depth')
    J, head = record.twos, record.head
    n, child = (1 << K)-1, (1 << (K-deletion))-1
    numerator = 3**J*(3**(K-1)-1)
    divisor = 1 << (2*J-1)
    need(numerator % divisor == 0, 'integral run endpoint')
    X = 1+numerator//divisor
    p, q, b = base.carrier(head.source)
    need((p*X+b) % q == 0, 'integral source head')
    Z = (p*X+b)//q
    need(Z % 2 == 1 and Z > n, 'native odd endpoint and immutable growth')
    c = base.v2(3*Z+1)
    endpoint = (3*Z+1) >> c
    left = (1,)*(K-1)+(2,)*J+head.source+(c,)
    clearing = (4, 1, 1) if deletion == 3 else (2, 2, 1, 1)
    right = ((1,)*(K-deletion-1)+clearing+(2,)*(J-3)
             +head.partner+(c+2*head.gap,))
    return base.audit_receipt(base.Receipt(n, child, left, right, endpoint))


discharge = base.discharge


def shell_mass(head):
    audit_head(head)
    need(head.source[0] != 2, 'nonpadding head')
    shift = 1 if head.source[0] == 1 else 2
    return F(1, 1 << (sum(head.source)-shift))


def finite_bank():
    """Declared decoder universe only: source length1..4, alphabet1..12."""
    import collatz_twoanchor_head_decoder_20261007b as decoder
    packets = tuple(packet for length in range(1, 5)
                    for word in product(range(1, 13), repeat=length)
                    for packet in decoder.all_partners(word))
    minimal, _ = decoder.antichain_report(packets)
    by_head = {}
    for p in packets:
        if p.head in minimal and p.head not in by_head:
            by_head[p.head] = compile_head(p.head, p.partner, p.gap)
    return tuple(by_head[w] for w in minimal)


def reject(function, *args):
    try:
        function(*args)
    except ValueError:
        need(True, 'hostile rejected')
        return
    need(False, 'hostile must reject')


def main():
    selected = (
        compile_head((10,), (3, 1, 1, 3)),
        compile_head((1, 10), (1, 1, 2, 2, 3)),
        compile_head((3, 12), (4, 1, 5, 2, 1)),
        compile_head((4, 12), (3, 6, 1, 3, 1)),
        compile_head((5, 8), (3, 1, 3, 3, 1)),
        compile_head((8, 4), (3, 1, 1, 4, 1)),
        compile_head((9, 2), (3, 1, 1, 3, 1)),
        compile_head((1, 2, 9), (1,)*6, 3),
    )
    for head in selected:
        for J in range(3, 19):
            p = phase(head, J)
            for t in (p.first_parameter, p.first_parameter+1, 10**30):
                K = exponent(p, t)
                need(contains(p, K) and native_guard(head, J, K), 'modular native phase')
                need(not native_guard(head, J, K+p.period//2), 'minimal-period hostile')
                need(K-1 >= 2*sum(head.source)+J, 'explicit all-prefix cutoff')
            m = 2*J-2 if head.source[0] == 1 else 2*J-1
            need(F(1 << (m+1), p.period) == shell_mass(head), 'conditional shell mass')
    for J in range(3, 19):
        for index, branch in ((0, 2), (1, 1)):
            old, new = base.phase(J, branch), phase(selected[index], J)
            need((old.residue, old.period) == (new.residue, new.period), 'old compiler recovered exactly')
    print('PROVED general native phase: period2^(2J+A-2), finite cutoffK-1>=2A+J.')
    print('FINITE-EXACT eight heads, J3..18, three parameters including10^30.')
    for A in range(1, 101):
        for J in range(3, 21):
            e = 2*A+J
            need(F(3**(e+J), 2**(e+2*J+A)) == F(9, 8)**(A+J) > 1,
                 'independent uniform all-prefix slope bound')

    literal = []
    for index in (0, 1, 6, 7):
        head = selected[index]
        p = phase(head, 3)
        K = exponent(p, p.first_parameter)
        for D in (3, 4):
            receipt = materialize(p, p.first_parameter, D)
            need(all(x > receipt.source for x in base.replay(receipt.source, receipt.source_word[:-1])[1:]),
                 'all literal prefixes through head grow')
            need(len(receipt.source_word) == len(receipt.child_word), 'equal odd clocks')
            need(sum(receipt.source_word)-sum(receipt.child_word) == D, 'deletion cost sidecar')
            need(source_residue(p, p.first_parameter, 2187, D) == receipt.child % 2187,
                 'symbolic child reader')
        literal.append((head.source, head.gap, K, receipt.source_word[-1]))
    print('Literal controls(head,gap,K,terminal):', literal)
    uncovered = []
    for index in (6, 7):
        h = selected[index]
        p = phase(h, 3)
        t = next(t for t in range(p.first_parameter, p.first_parameter+729)
                 if base.baseline_entry(exponent(p, t)) is None)
        K = exponent(p, t)
        branch = 1 if h.source[0] == 1 else 2
        need(native_guard(h, 3, K) and not base.head_valuation_guard(K, 3, branch),
             'new head guard outside the previous single-head phase')
        need(base.native_baseline_entry(K) is None, 'new point also misses entry bank')
        uncovered.append((h.source, h.gap, K, 729*p.period))
    print('New symbolic rays(head,gap,least-selected-K,period):', uncovered)

    bank = finite_bank()
    import collatz_twoanchor_head_decoder_20261007b as decoder
    need(tuple((h.source, h.partner, h.gap) for h in bank) ==
         tuple((h.head, h.partner, h.gap) for h in decoder.finite_head_bank()),
         'independent bank reconstruction equals the frozen export')
    need(len(bank) == 225, 'declared all-gap prefix-minimal bank')
    need(max(sum(h.source) for h in bank) == 25, 'finite-bank maximum head cost')
    for K in (17, 33, 49):
        for h in bank:
            need(not native_guard(h, 3, K), 'all possible small cutoff exceptions excluded')
    for a in bank:
        for b in bank:
            need(a.source == b.source or b.source[:len(a.source)] != a.source,
                 'actual valuation prefixes are incomparable')
    masses = []
    for first in (1, 3):
        group = [h for h in bank if (h.source[0] == 1 if first == 1 else h.source[0] >= 3)]
        mass = sum(map(shell_mass, group), F())
        masses.append((first, len(group), mass, mass*F(485, 729)))
    need(masses[0][1:3] == (35, F(2233, 524288)), 'even-shell finite bank mass')
    need(masses[1][1:3] == (190, F(142353, 8388608)), 'odd-shell finite bank mass')
    # Cylinder/phase disjointness is checked independently by dyadic congruences.
    for J in (3, 7):
        plans = [phase(h, J) for h in bank]
        need(all(p.first_parameter == 0 for p in plans), 'finite bank retains every positive native member')
        for i, p in enumerate(plans):
            for q in plans[:i]:
                need((p.residue-q.residue) % min(p.period, q.period) != 0,
                     'independent exponent cylinders are disjoint')
        for p in plans:
            counts = {None: 0, 5: 0, 17: 0}
            for t in range(729):
                K = p.residue+p.period*t
                entry = base.baseline_entry(K)
                need(entry == base.native_baseline_entry(K), 'exact old entry predicate')
                counts[entry] += 1
            need(counts == {None: 485, 5: 243, 17: 1}, 'uniform ternary complement')
    print('Finite bank rows(first,count,shellmass,new-relative-to-entry-bank):', masses)
    print('All225 exponent guards disjoint atJ3,7; each has old-entry split243/1/485.')
    print('Finite bank maxcost25; no head matchesK17,33,49. No positive native member is cut.')

    p = phase(selected[0], 3)
    reject(compile_head, (10,), (3, 1, 1, 3), True)
    reject(compile_head, (10,), (3, 1, 2, 2), 1)
    reject(phase, compile_head((2, 10), (2, 3, 1, 1, 3)), 3)
    reject(phase, replace(selected[0], gap=1.0), 3)
    reject(contains, replace(p, first_parameter=False), p.residue)
    reject(contains, replace(p, period=p.period//2), p.residue)
    reject(materialize, phase(selected[2], 3), 0, 4, 40000)
    reject(materialize, p, 0, True)
    reject(source_residue, p, 0, 0)
    need(discharge(base.Receipt(3, 1, (1, 4), (), 1), ()) == (1, 4), 'supplied ROOT transport')
    print('No ROOT search; finite bank mass is paid-dependency coverage, not ROOT completion.')
    print('Checks:', CHECKS+base.CHECKS)


if __name__ == '__main__':
    main()
