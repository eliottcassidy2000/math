"""All-J two-anchor deletion phases, with an explicit inverse-entry baseline.

No production orbit search or implicit child ROOT certificate.  Every literal
receipt is replayed; the symbolic phase API does not expand its Mersenne source.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F

import collatz_child_normal_forms_20261007 as forms


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def integer(n, least=0):
    need(type(n) is int and n >= least, 'exact integer in declared domain')
    return n


def odd(n):
    integer(n, 1)
    need(n % 2 == 1, 'positive odd integer')
    return n


def letters(word):
    need(type(word) is tuple and all(type(a) is int and a > 0 for a in word),
         'exact positive valuation tuple')
    return word


def v2(n):
    integer(n, 1)
    return (n & -n).bit_length() - 1


def carrier(word):
    letters(word)
    p = q = 1
    b = 0
    for a in word:
        p, q, b = 3*p, q*2**a, 3*b+q
    return p, q, b


def affine(word, n):
    need(type(n) in (int, F), 'exact affine input')
    p, q, b = carrier(word)
    return F(p*n+b)/q


def replay(n, word):
    odd(n)
    letters(word)
    states = [n]
    for a in word:
        need(n != 1, 'no ROOT padding')
        z = 3*n+1
        need(v2(z) == a, 'exact actual valuation')
        n = z >> a
        states.append(n)
    return tuple(states)


def log_three(unit, bits):
    """Exact exponent in <3> modulo 2**bits; no decimal/log approximations."""
    integer(bits, 3)
    integer(unit)
    modulus = 1 << bits
    unit %= modulus
    need(unit % 8 in (1, 3), 'unit belongs to the exact power-three subgroup')
    exponent = 0 if unit % 8 == 1 else 1
    period = 2
    for precision in range(4, bits+1):
        mod = 1 << precision
        choices = [e for e in (exponent, exponent+period)
                   if pow(3, e, mod) == unit % mod]
        need(len(choices) == 1, 'unique binary exponent lift')
        exponent, period = choices[0], 2*period
    return exponent, period


@dataclass(frozen=True)
class Phase:
    twos: int
    residue: int
    period: int
    branch: int = 2


def phase(twos, branch=2):
    integer(twos, 3)
    need(type(branch) is int and branch in (1, 2), 'exact run-end branch1 or2')
    bits = 2*twos+(11 if branch == 1 else 10)
    modulus = 1 << bits
    if branch == 1:
        target = (1+1017*(1 << (2*twos))*pow(3, -twos-2, modulus)) % modulus
    else:
        target = (1+255*(1 << (2*twos+1))*pow(3, -twos-1, modulus)) % modulus
    exponent, period = log_three(target, bits)
    result = Phase(twos, exponent+1, period, branch)
    need(v2(result.residue-1) == 2*twos+branch-3, 'phase forces the exact two-run shell')
    need(period == 1 << (bits-2), 'minimal exponent period')
    return result


def audit_phase(record):
    need(type(record) is Phase, 'exact Phase object')
    integer(record.twos, 3)
    integer(record.residue, 1)
    integer(record.period, 1)
    need(type(record.branch) is int and record.branch in (1, 2), 'exact branch field')
    need(record == phase(record.twos, record.branch), 'canonical phase fields')
    return record


def exponent(record, parameter):
    audit_phase(record)
    integer(parameter)
    return record.residue+record.period*parameter


def contains(record, K):
    audit_phase(record)
    integer(K, 1)
    return (K-record.residue) % record.period == 0


def source_residue(record, parameter, modulus, child=False, deletion=3):
    """Read the source/child without materializing 2**K."""
    K = exponent(record, parameter)
    integer(modulus, 1)
    need(type(child) is bool, 'exact child flag')
    need(type(deletion) is int and deletion in (3, 4), 'three- or four-bit deletion')
    return (pow(2, K-deletion if child else K, modulus)-1) % modulus


def head_valuation_guard(K, twos, branch=2):
    """Independent numerator test for native head10 or1,10."""
    integer(K, 2)
    integer(twos, 3)
    need(type(branch) is int and branch in (1, 2), 'exact run-end branch')
    bits = 2*twos+(11 if branch == 1 else 10)
    modulus = 1 << bits
    power = twos+(2 if branch == 1 else 1)
    offset = 7*(1 << (2*twos)) if branch == 1 else 1 << (2*twos+1)
    numerator = (pow(3, power, modulus)*(pow(3, K-1, modulus)-1)+offset) % modulus
    return (v2(K-1) == 2*twos+branch-3 and numerator == 1 << (bits-1))


@dataclass(frozen=True)
class Receipt:
    source: int
    child: int
    source_word: tuple
    child_word: tuple
    endpoint: int


def audit_receipt(receipt):
    need(type(receipt) is Receipt, 'exact Receipt object')
    for n in (receipt.source, receipt.child, receipt.endpoint):
        odd(n)
    need(receipt.child < receipt.source, 'strict original-source payment')
    for n, word in ((receipt.source, receipt.source_word),
                    (receipt.child, receipt.child_word)):
        p, q, b = carrier(word)
        need(p*n+b == q*receipt.endpoint, 'independent affine endpoint')
        need(replay(n, word)[-1] == receipt.endpoint, 'strict first-hit replay')
    return receipt


def materialize(record, parameter, bit_cap=20000, deletion=3):
    """Only bounded, selected examples expand their sources and actual words."""
    K = exponent(record, parameter)
    integer(bit_cap, 1)
    need(K <= bit_cap, 'explicit source-bit materialization cap')
    need(type(deletion) is int and deletion in (3, 4), 'three- or four-bit deletion')
    J = record.twos
    n, child = (1 << K)-1, (1 << (K-deletion))-1
    numerator = 3**J*(3**(K-1)-1)
    need(numerator % (1 << (2*J-1)) == 0, 'integral run endpoint')
    X = 1+numerator//(1 << (2*J-1))
    if record.branch == 1:
        need(v2(3*X+1) == 1 and v2(9*X+5) == 11, 'native post-run head1,10')
        source_head, child_head = (1, 10), (1, 1, 2, 2, 3)
        Z = (9*X+5) >> 11
    else:
        need(v2(3*X+1) == 10, 'native post-run head10')
        source_head, child_head = (10,), (3, 1, 1, 3)
        Z = (3*X+1) >> 10
    need(Z > n, 'head still above immutable source, hence above ROOT')
    c = v2(3*Z+1)
    endpoint = (3*Z+1) >> c
    left = (1,)*(K-1)+(2,)*J+source_head+(c,)
    clearing = (4, 1, 1) if deletion == 3 else (2, 2, 1, 1)
    right = (1,)*(K-deletion-1)+clearing+(2,)*(J-3)+child_head+(c+2,)
    return audit_receipt(Receipt(n, child, left, right, endpoint))


def discharge(receipt, supplied_child_word):
    """Transport a supplied first-hit child proof; never discover one."""
    audit_receipt(receipt)
    need(replay(receipt.child, supplied_child_word)[-1] == 1,
         'supplied child proof reaches ROOT')
    cut = len(receipt.child_word)
    need(supplied_child_word[:cut] == receipt.child_word, 'same actual child prefix')
    word = receipt.source_word+supplied_child_word[cut:]
    need(replay(receipt.source, word)[-1] == 1, 'first-hit source ROOT proof')
    return word


def g5_extend(receipt):
    """A native child reduction composes the receipt, preserving its source."""
    audit_receipt(receipt)
    child = forms.apply(forms.encode((5,)), receipt.child)
    need(child < receipt.child, 'additional G5 child is strictly smaller')
    return audit_receipt(Receipt(receipt.source, child, receipt.source_word,
                                 (1, 2)+receipt.child_word, receipt.endpoint))


def baseline_entry(K):
    """First entry of any nonempty G1/G5/G17 word on a Mersenne source."""
    integer(K, 1)
    if K % 6 == 5:
        return 5
    if K % 1458 == 733:
        return 17
    return None


def native_baseline_entry(K):
    """Independent native-source congruences, without using exponent phases."""
    integer(K, 1)
    residue = (pow(2, K, 2187)-1) % 2187
    odd_representative = residue if residue % 2 else residue+2187
    return forms.next_letter(odd_representative)


def phase_counts(record):
    audit_phase(record)
    counts = {5: 0, 17: 0, None: 0}
    for t in range(729):
        K = record.residue+record.period*t
        entry = baseline_entry(K)
        need(entry == native_baseline_entry(K), 'independent baseline guard')
        counts[entry] += 1
    return counts


def reject(function, *args):
    try:
        function(*args)
    except ValueError:
        need(True, 'hostile rejected')
        return
    need(False, 'hostile must be rejected')


def main():
    # Re-derive only the incoming affine template actually consumed here.
    need(carrier((2, 2, 2)) == (27, 64, 37), 'source clearing carrier')
    need(carrier((4, 1, 1)) == (27, 64, 89), 'child clearing carrier')
    need(carrier((2, 2, 1, 1)) == (81, 64, 143), 'four-bit clearing carrier')
    need(carrier((10,)) == (3, 1024, 1), 'source head carrier')
    need(carrier((3, 1, 1, 3)) == (81, 256, 179), 'child head carrier')
    need(carrier((1, 10)) == (9, 2048, 5), 'even-shell source head')
    need(carrier((1, 1, 2, 2, 3)) == (243, 512, 283), 'even-shell child head')
    for x in range(-9, 10):
        x = F(x)
        y = (x+1)/27-1
        need(affine((2, 2, 2), x)-1 == 27*(affine((4, 1, 1), y)-1),
             'clearing polynomial identity at independent points')
        need(affine((2, 2, 1, 1), (x+1)/81-1) == (x+63)/64,
             'four-bit deletion has the identical cleared endpoint')
        Y = (x-1)/27+1
        need(affine((3, 1, 1, 3), Y) == 4*affine((10,), x)+1,
             'head polynomial identity at independent points')
        need(affine((1, 1, 2, 2, 3), Y) == 4*affine((1, 10), x)+1,
             'even-shell head polynomial identity')
    phases = [phase(J) for J in range(3, 33)]
    even_phases = [phase(J, 1) for J in range(3, 33)]
    for record in phases+even_phases:
        J = record.twos
        for t in (0, 1, 2, 7, 10**20):
            K = exponent(record, t)
            need(contains(record, K) and head_valuation_guard(K, J, record.branch), 'symbolic phase guard')
            need(not head_valuation_guard(K+record.period//2, J, record.branch), 'one-less-bit hostile')
        need(phase_counts(record) == {5: 243, 17: 1, None: 485}, 'exact per-phase coverage split')
    need([(p.twos, p.residue, p.period) for p in phases[:3]] ==
         [(3, 1889, 16384), (4, 5249, 65536), (5, 50689, 262144)], 'independent incoming examples and next row')
    need((even_phases[0].residue, even_phases[0].period) == (6129, 32768), 'new even-shell least phase')
    print('PROVED guards: every shell m>=4 has a head phase; child exponent K-D, D=3 or4.')
    print('FINITE-EXACT both branches J=3..32 and five parameters; all baseline splits243/1/485.')
    print('Odd-shell first rows:', [(p.twos, p.residue, p.period) for p in phases[:8]])
    print('Even-shell first rows:', [(p.twos, p.residue, p.period) for p in even_phases[:6]])

    # Count entire fixed shells independently of the specific head phase.
    for m in range(4, 21):
        counts = {5: 0, 17: 0, None: 0}
        for t in range(729):
            K = 1+(1 << m)*(2*t+1)
            need(v2(K-1) == m, 'fixed shell')
            entry = baseline_entry(K)
            need(entry == native_baseline_entry(K), 'shell native baseline')
            counts[entry] += 1
        need(counts == {5: 243, 17: 1, None: 485}, 'shell entry density')
    print('FINITE-EXACT shells m=4..20: inverse-entry bank244/729, complement485/729.')

    selected = ((phases[0], 0), (phases[1], 0), (phases[0], 1), (even_phases[0], 0))
    rows = []
    composed = None
    for p, t in selected:
        K = exponent(p, t)
        receipt = materialize(p, t)
        states = replay(receipt.source, receipt.source_word[:-1])
        need(all(x > receipt.source for x in states[1:]), 'every source prefix through the head grows')
        need(len(receipt.source_word) == len(receipt.child_word) == K+p.twos+3-p.branch,
             'equal actual odd clocks')
        need(sum(receipt.source_word)-sum(receipt.child_word) == 3, 'three-bit actual cost saving')
        for mod in (2, 3, 19, 2187, 65537):
            need(source_residue(p, t, mod) == receipt.source % mod, 'source modular reader')
            need(source_residue(p, t, mod, True) == receipt.child % mod, 'child modular reader')
        rows.append((K, p.twos, p.branch, receipt.source_word[-1], K-3, baseline_entry(K)))
        fourth = materialize(p, t, deletion=4)
        need(fourth.endpoint == receipt.endpoint and fourth.child < receipt.child,
             'same receipt endpoint with four-bit smaller child')
        need(len(fourth.child_word) == len(fourth.source_word)
             and sum(fourth.source_word)-sum(fourth.child_word) == 4,
             'four-bit clock and cost identity')
        need(source_residue(p, t, 2187, True, 4) == fourth.child % 2187,
             'four-bit child modular reader')
        need(v2(K-5) == 2, 'four-bit child has exactly two twos after its ones-run')
        if K in (18273, 6129):
            composed = g5_extend(fourth)
            need(9*composed.child == (1 << (K-1))-13, 'explicit post-deletion G5 formula')
            need(composed.source == fourth.source and composed.endpoint == fourth.endpoint,
                 'composed receipt retains source ownership and common future')
    print('Literal selected receipts (K,J,branch,c,childK,old entry):', rows)
    need(baseline_entry(18273) is None and native_baseline_entry(18273) is None,
         'concrete newly covered parameter')
    for s in range(128):
        K = 18273+729*16384*s
        need(head_valuation_guard(K, 3) and baseline_entry(K) is None,
             'entire excluded ray preserves head and absence of entry')
        need(v2(K-5) == 2 and baseline_entry(K-4) == 5,
             'excluded ray reaches the two-twos type with a legal G5 entry')
        K2 = 6129+729*32768*s
        need(head_valuation_guard(K2, 3, 1) and baseline_entry(K2) is None
             and v2(K2-5) == 2 and baseline_entry(K2-4) == 5,
             'even-shell excluded ray and composable exit')
    print('PROVED excluded ray: K=18273+11943936*s, s>=0; all have the new paid dependency.')
    print('PROVED even-shell excluded ray: K=6129+23887872*s, s>=0.')
    print('D4 normalizes every phase to two twos; excluded ray then has child(2^(K-1)-13)/9 viaG5.')
    need(composed is not None, 'one authenticated composed point control')

    # General transport boundary; these controls do not discover a new ROOT orbit.
    root_receipt = Receipt(3, 1, (1, 4), (), 1)
    need(discharge(root_receipt, ()) == (1, 4), 'supplied empty ROOT suffix')
    reject(discharge, root_receipt, (2,))
    reject(audit_receipt, replace(root_receipt, endpoint=True))
    reject(audit_receipt, replace(root_receipt, child=1.0))
    reject(audit_receipt, replace(root_receipt, source_word=(1, 3)))
    for bad in (True, 3.0, 2, -1):
        reject(phase, bad)
    for bad in (replace(phases[0], twos=3.0), replace(phases[0], residue=1889.0),
                replace(phases[0], period=8192), replace(phases[0], branch=2.0)):
        reject(contains, bad, 1889)
    reject(materialize, phases[2], 0, 20000)
    reject(materialize, phases[0], -1)
    reject(source_residue, phases[0], 0, 0)
    reject(source_residue, phases[0], 0, 3, 1)
    reject(materialize, phases[0], 0, 20000, True)
    reject(materialize, phases[0], 0, 20000, 5)
    reject(log_three, 5, 8)
    reject(phase, 3, True)
    need(affine(forms.GENERATORS[5][3], F(1, 3)) == 1,
         'alternate inverse anchor1/3 really maps to1')
    need(F(2048-2363, 2187) == F(-35, 243), 'signed G17 limit is not an integer-domain claim')
    print('No child ROOT word discovered or assumed; completion requires a supplied certificate.')
    print('Checks:', CHECKS)


if __name__ == '__main__':
    main()
