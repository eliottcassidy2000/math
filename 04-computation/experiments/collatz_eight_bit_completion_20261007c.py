"""A guarded eight-bit deletion composed through a two-twos child.

Symbolic Mersenne guards never expand their enormous sources. Literal
interfaces retain exact source/child words and consume supplied ROOT proofs.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from math import gcd

import collatz_uncovered_join_routes_20261007 as routes
import collatz_child_normal_forms_20261007 as forms
import collatz_completion_paper24_23_20261007b as budget
import collatz_completion_anchor_20261007b as anchors


LEFT = (2, 2, 2, 1, 2, 9, 16)
J2_LEFT = (1,)*8+(22,)
J2_RIGHT = (4, 1, 3, 3, 1, 3, 3, 1, 1, 1, 2, 2, 3)
RIGHT = (2, 2)+J2_RIGHT
DROP = 8
PHASE = 924745905
PERIOD = 1 << 32
MIN_RUN = 68
CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def integer(value, minimum=0):
    need(type(value) is int and value >= minimum, 'exact integer in domain')
    return value


def affine(word, value):
    p, q, b = routes.carrier(word)
    return F(p*value+b, q)


def audit_rule():
    p, q, b = routes.carrier(LEFT)
    pp, qq, bb = routes.carrier(RIGHT)
    need((len(RIGHT)-len(LEFT), sum(LEFT), sum(RIGHT)) == (8, 34, 32),
         'ordered head lengths and costs')
    need(pp == 3**DROP*p and q == 4*qq, 'slope and clock identity')
    need(affine(RIGHT, -1) == 4*affine(LEFT, -1)+1,
         'full minus-one anchor, including carry')
    residue = (q-b)*pow(p, -1, 2*q) % (2*q)
    need(residue % 4 == 1, 'maximal initial ones has a reset-two head')
    target = (residue+1)//2
    e, period = anchors.log_three(target, 34)
    need((e+1, period) == (PHASE, PERIOD), 'least native Mersenne phase')
    need(PHASE-1 >= MIN_RUN, 'all phase parameters meet growth cutoff')
    return residue, q


def contains(K):
    integer(K, 1)
    audit_rule()
    return K-1 >= MIN_RUN and (K-PHASE) % PERIOD == 0


def symbolic_source(K, modulus, deletion=0):
    integer(modulus, 1)
    need(contains(K), 'same guarded exponent')
    need(type(deletion) is int and deletion in (0, 4, 8), 'declared source port')
    return (pow(2, K-deletion, modulus)-1) % modulus


def fixed_run_source(run, parameter):
    """All native positive sources at a fixed sufficiently long ones-run."""
    integer(run, MIN_RUN)
    integer(parameter)
    residue, q = audit_rule()
    oddpart = (((residue+1)//2)*pow(3, -run, q)) % q
    need(oddpart % 2 == 1, 'source oddpart and exact run')
    return (1 << (run+1))*(oddpart+q*parameter)-1


def receipt(source, bit_cap=20000):
    routes.odd(source)
    integer(bit_cap, 1)
    need(source.bit_length() <= bit_cap, 'explicit literal source cap')
    residue, q = audit_rule()
    run = forms.valuation(source+1, 2)-1
    need(run >= MIN_RUN, 'proved source-growth cutoff')
    oddpart = (source+1) >> (run+1)
    need((2*pow(3, run, 2*q)*oddpart-1) % (2*q) == residue,
         'actual supplied source native guard')
    x = 2*3**run*oddpart-1
    p, q, b = routes.carrier(LEFT)
    need((p*x+b) % q == 0, 'source head integral')
    z = (p*x+b)//q
    need(z % 2 == 1 and z > source, 'source head grows and avoids ROOT')
    c = forms.valuation(3*z+1, 2)
    endpoint = (3*z+1) >> c
    child = ((source+1) >> DROP)-1
    left = (1,)*run+LEFT+(c,)
    right = (1,)*(run-DROP)+RIGHT+(c+2,)
    return routes.audit(routes.Receipt(source, child, left, right, endpoint))


def compose(first, second):
    """Join actual receipts by aligning their immutable middle-source clock."""
    routes.audit(first)
    routes.audit(second)
    need(first.child == second.source, 'same middle source, not another family member')
    middle_a, middle_b = first.child_word, second.source_word
    common = min(len(middle_a), len(middle_b))
    need(middle_a[:common] == middle_b[:common], 'actual middle prefixes agree')
    if len(middle_a) <= len(middle_b):
        left = first.source_word+middle_b[len(middle_a):]
        right, endpoint = second.child_word, second.endpoint
    else:
        left = first.source_word
        right = second.child_word+middle_a[len(middle_b):]
        endpoint = first.endpoint
    return routes.audit(routes.Receipt(first.source, second.child, left, right, endpoint))


def two_stage(source):
    """Independent literal replay of the two four-bit discovery receipts."""
    combined = receipt(source)
    run = forms.valuation(source+1, 2)-1
    n4 = ((source+1) >> 4)-1
    left1 = (1,)*run+(2,)*3+(1, 2, 9, 16)
    right1 = (1,)*(run-4)+(2, 2, 1, 1)+(1,)*6+(22,)
    endpoint1 = routes.replay(source, left1)[-1]
    first = routes.audit(routes.Receipt(source, n4, left1, right1, endpoint1))
    c = combined.source_word[-1]
    left2 = (1,)*(run-4)+(2, 2)+J2_LEFT+(c,)
    second = routes.audit(routes.Receipt(n4, combined.child, left2,
                                         combined.child_word, combined.endpoint))
    need(right1 == left2[:-1], 'the middle source needs exactly one extra edge')
    return first, second


@dataclass(frozen=True)
class MixedPhase:
    word: tuple
    residue: int
    period: int
    numerical_depth: int


def mixed_phase(word):
    """Mersenne parents where the actual D8 receipt pays an inverse child."""
    request = budget.uniform_budget(word)
    need(request.required <= DROP, 'uniform numerical depth fits actual receipt')
    old = forms.mersenne_phase(forms.encode(word))
    need(old is not None, 'inverse word has a Mersenne native phase')
    r, period = old
    divisor = gcd(PERIOD, period)
    need((r-PHASE) % divisor == 0, 'compatible native phases')
    odd_modulus = period//divisor
    lift = ((r-PHASE)//divisor*pow(PERIOD//divisor, -1, odd_modulus)) % odd_modulus
    result = MixedPhase(word, PHASE+PERIOD*lift, PERIOD*odd_modulus, request.required)
    need(contains(result.residue), 'CRT retains deletion source')
    need((result.residue-r) % period == 0, 'CRT retains inverse source')
    return result


def audit_mixed(record):
    need(type(record) is MixedPhase, 'exact MixedPhase object')
    for value in (record.residue, record.period, record.numerical_depth):
        integer(value, 1)
    need(record == mixed_phase(record.word), 'canonical same-source CRT packet')
    return record


def mixed_literal_source(word, run=128, parameter=0):
    """An ordinary, modest-size parent in both native families by exact CRT."""
    integer(run, MIN_RUN)
    integer(parameter)
    b = budget.uniform_budget(word)
    need(b.required <= DROP, 'requested payment depth fits D8')
    record = forms.encode(word)
    residue, modulus = forms.native_cell(record)
    base = fixed_run_source(run, 0)
    step = 1 << (run+1+sum(LEFT))
    lift = ((residue-base)*pow(step, -1, modulus)) % modulus
    return fixed_run_source(run, lift+modulus*parameter)


def transport(word, parent):
    request = budget.request(parent, word)
    return budget.consume_deletion(request, receipt(parent))


def main():
    residue, q = audit_rule()
    # The next child is a new port, not automatically covered by its parent's rule.
    from itertools import product
    import collatz_complement_ladders_20261007c as ladders
    packets = hits = 0
    for length in range(1, 4):
        for head in product(range(1, 13), repeat=length):
            if head[0] == 2:
                continue
            for candidate in ladders.all_ladders(head, ladders.TWO_TWOS):
                ph = ladders.phase(candidate, 2)
                packets += 1
                hits += (PHASE-DROP-ph.residue) % min(PERIOD, ph.period) == 0
    need((packets, hits) == (271, 0), 'declared short J2 bank leaves final child open')
    # Independent component identities, including the changed J2 anchor.
    for y in (F(-1,8), F(1,3), 3, 27, 135):
        x = 81*y+10
        need(affine(J2_RIGHT,y) == 4*affine(J2_LEFT,x)+1, 'J2 shifted anchor')
    for y in (-3, F(1,7), 1, 73):
        x = 3**DROP*(y+1)-1
        need(affine(RIGHT,y) == 4*affine(LEFT,x)+1, 'direct D8 affine identity')
    for run in range(MIN_RUN, MIN_RUN+101):
        need(F(3**run,2**(run+34)) > 1, 'all-prefix growth floor')
    for t in (0,1,2,729,10**30):
        K = PHASE+PERIOD*t
        need(contains(K), 'symbolic phase member')
        need((2*pow(3,K-1,2*q)-1)%(2*q) == residue, 'independent modular source guard')
        need(not contains(K+PERIOD//2), 'least-period hostile')
        need((K-4-177949)%262144 != 0, 'short off-target J2 head is not substituted')
        for m in (19,81,2187):
            need(symbolic_source(K,m,8) == (pow(2,K-8,m)-1)%m, 'same-source symbolic port')
    literal_count = 0
    for run in (68,69,96,128):
        for t in (0,1,17,10**12):
            n = fixed_run_source(run,t)
            direct = receipt(n)
            first, second = two_stage(n)
            need(compose(first,second) == direct, 'two-stage and direct certificates coincide')
            need(sum(direct.source_word)-sum(direct.child_word) == 8, 'paid binary depth')
            need(len(direct.source_word) == len(direct.child_word), 'aligned odd clocks')
            need(all(z>n for z in routes.replay(n,direct.source_word[:-1])[1:]),
                 'literal immutable-source growth at every prefix')
            literal_count += 1
    need(3**13 < (1 << 21), 'every at-most-thirteen-letter inverse word needs at most eight bits')
    mixed=[]
    for word in ((5,)*36,(17,)*64,(5,1,17,5,1)):
        p = audit_mixed(mixed_phase(word))
        for t in (0,1,10**12):
            K = p.residue+p.period*t
            need(contains(K), 'every CRT lift stays in paid deletion guard')
            r,m = forms.mersenne_phase(forms.encode(word))
            need((K-r)%m == 0, 'every CRT lift stays inverse-native')
            n = mixed_literal_source(word,128,t)
            out = transport(word,n)
            need(out.source == forms.apply(forms.encode(word),n), 'transport preserves supplied inverse child')
            need(out.child == ((n+1)>>8)-1 < out.source, 'transport now pays original child')
        mixed.append((word[0],len(word),p.numerical_depth,p.residue.bit_length(),p.period.bit_length()))
    for bad in (True,1.0,-1):
        try: contains(bad)
        except ValueError: need(True,'bad exact type rejected')
        else: need(False,'bad exact type accepted')
    p=mixed_phase((5,)*36)
    for bad in (replace(p,residue=p.residue+1),replace(p,numerical_depth=True)):
        try: audit_mixed(bad)
        except ValueError: need(True,'forged symbolic packet rejected')
        else: need(False,'forged symbolic packet accepted')
    n=fixed_run_source(68,0)
    try: receipt(n,bit_cap=16)
    except ValueError: need(True,'materialization cap retained')
    else: need(False,'materialization cap ignored')
    try: compose(receipt(n),receipt(fixed_run_source(68,1)))
    except ValueError: need(True,'changed middle source rejected')
    else: need(False,'different family member substituted')
    print('PROVED D8 phase K=',PHASE,'mod',PERIOD,'; heads',LEFT,RIGHT)
    print('Discovery: D4 through J3 head129, then D4 through J2 head1^8,22; exact one-edge middle alignment.')
    print('FINITE-EXACT literal ordinary-source controls',literal_count,'; no astronomical Mersenne expansion.')
    print('Mixed inverse transports(label,length,required_depth,CRT_residue_bits,period_bits):',mixed)
    print('PROVED at-most13 inverse letters fit D8; G5^36 andG17^64 now have actual guarded payment.')
    print('FINITE-EXACT next-child probe:271 signed J2 packets, heads1..12,length1..3,nonpadding;0 phase intersections.')
    print('ROOT suffixes remain supplied obligations; no ROOT discovery, no universal coverage claim.')
    print('Exact checks',CHECKS)


if __name__ == '__main__':
    main()
