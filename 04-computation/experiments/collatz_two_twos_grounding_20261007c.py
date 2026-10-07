"""Frozen ROOT certificates and a paid exit from emitted two-twos children.

Default execution only authenticates retained evidence. Bounded orbit discovery
is isolated in discover() and the explicit --rediscover control.
"""
from dataclasses import dataclass, replace
from hashlib import sha256
from pathlib import Path
from fractions import Fraction as F
import base64
import json
import zlib

import collatz_completion_anchor_20261007b as anchor
import collatz_general_head_phases_20261007b as general


DATA = Path(__file__).resolve().parents[2] / '05-knowledge/results/collatz_two_twos_grounding_20261007c.cert.json'
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def integer(n, least=0):
    need(type(n) is int and n >= least, 'exact integer in declared domain')
    return n


def odd(n):
    integer(n, 1)
    need(n % 2 == 1, 'positive odd source')
    return n


def letters(word):
    need(type(word) is tuple and all(type(a) is int and a > 0 for a in word),
         'exact positive valuation tuple')
    return word


def v2(n):
    integer(n, 1)
    return (n & -n).bit_length()-1


@dataclass(frozen=True)
class Stats:
    endpoint: int
    rank: int
    cost: int
    peak_bits: int


def verify_word(source, word, require_root=True):
    """Constant-state-memory, strict first-hit actual replay."""
    odd(source)
    letters(word)
    need(type(require_root) is bool, 'exact ROOT requirement flag')
    n, cost, peak = source, 0, source.bit_length()
    for a in word:
        need(n != 1, 'no first-hit ROOT padding')
        z = 3*n+1
        need(v2(z) == a, 'exact actual valuation')
        n = z >> a
        cost += a
        peak = max(peak, n.bit_length())
    need(not require_root or n == 1, 'supplied word reaches ROOT')
    return Stats(n, len(word), cost, peak)


def reverse_verify(source, word):
    odd(source)
    letters(word)
    n = 1
    for a in reversed(word):
        numerator = (n << a)-1
        need(numerator % 3 == 0, 'independent reverse integrality')
        n = numerator//3
        need(n > 1 and n % 2 == 1, 'reverse word has no earlier ROOT')
    need(n == source, 'independent reverse source identity')
    return True


def frozen_record(E):
    integer(E, 2)
    data = json.loads(DATA.read_text(encoding='utf-8'))
    need(data['format'] == 'leading-ones/zlib-base85-u8-v1', 'declared certificate format')
    need(type(data['tail_discovery_cap']) is int and data['tail_discovery_cap'] == 600000,
         'explicit finite discovery universe')
    matches = [r for r in data['certificates'] if type(r['exponent']) is int and r['exponent'] == E]
    need(len(matches) == 1, 'one frozen source-specific record')
    return matches[0]


def root_word(E):
    """Decode, hash and literally authenticate a frozen source certificate."""
    r = frozen_record(E)
    for name in ('leading_ones', 'rank', 'cost', 'peak_bits', 'tail_length'):
        integer(r[name], 1)
    need(r['leading_ones'] == E-1 and r['rank'] == E-1+r['tail_length'], 'source-owned initial run')
    need(type(r['tail_base85']) is str and type(r['decoded_sha256']) is str, 'exact encoded data fields')
    compressed = base64.b85decode(r['tail_base85'].encode('ascii'))
    inflater = zlib.decompressobj()
    raw = inflater.decompress(compressed, r['tail_length']+1)
    need(inflater.eof and not inflater.unused_data and not inflater.unconsumed_tail
         and len(raw) == r['tail_length'], 'bounded exact payload expansion')
    need(all(a > 0 for a in raw), 'positive byte valuation letters')
    payload = bytes([1])*(E-1)+raw
    need(sha256(payload).hexdigest() == r['decoded_sha256'], 'decoded first-hit word hash')
    word = tuple(payload)
    stats = verify_word((1 << E)-1, word)
    need((stats.rank, stats.cost, stats.peak_bits) == (r['rank'], r['cost'], r['peak_bits']),
         'exact stored certificate metadata')
    return word


def discover(E, tail_cap=600000):
    """Explicit bounded discovery only; never called by production consumers."""
    integer(E, 2)
    integer(tail_cap, 1)
    n = 2*3**(E-1)-1
    tail = []
    while n != 1 and len(tail) < tail_cap:
        z = 3*n+1
        a = (z & -z).bit_length()-1
        n = z >> a
        tail.append(a)
    need(n == 1, 'explicit discovery cap expired before ROOT')
    return (1,)*(E-1)+tuple(tail)


@dataclass(frozen=True)
class Receipt:
    source: int
    child: int
    source_word: tuple
    child_word: tuple
    endpoint: int


def audit_receipt(receipt):
    need(type(receipt) is Receipt, 'exact Receipt')
    for n in (receipt.source, receipt.child, receipt.endpoint):
        odd(n)
    need(receipt.child < receipt.source, 'strict immutable-source payment')
    for source, word in ((receipt.source, receipt.source_word),
                         (receipt.child, receipt.child_word)):
        need(verify_word(source, word, False).endpoint == receipt.endpoint, 'actual common future')
        p, q, b = anchor.carrier(word)
        need(p*source+b == q*receipt.endpoint, 'independent affine common future')
    return receipt


def discharge(receipt, supplied_child_word):
    audit_receipt(receipt)
    verify_word(receipt.child, supplied_child_word)
    cut = len(receipt.child_word)
    need(supplied_child_word[:cut] == receipt.child_word, 'same actual child prefix')
    result = receipt.source_word+supplied_child_word[cut:]
    verify_word(receipt.source, result)
    return result


def reverse_discharge(receipt, supplied_source_word):
    audit_receipt(receipt)
    verify_word(receipt.source, supplied_source_word)
    cut = len(receipt.source_word)
    need(supplied_source_word[:cut] == receipt.source_word, 'same actual source prefix')
    result = receipt.child_word+supplied_source_word[cut:]
    verify_word(receipt.child, result)
    return result


def parent_receipt(K, head, partner, gap, bit_cap=100000):
    """Stream-checked inherited J3/D4 receipt; no giant list of orbit states."""
    integer(K, 5)
    integer(bit_cap, 1)
    need(K <= bit_cap, 'explicit inherited receipt materialization cap')
    packet = general.compile_head(head, partner, gap)
    plan = general.phase(packet, 3)
    need(general.contains(plan, K), 'exact inherited exponent phase and cutoff')
    x = 2*3**(K-1)-1
    X = (27*x+37)//64
    p, q, b = anchor.carrier(head)
    need((p*X+b) % q == 0, 'integral inherited source head')
    Z = (p*X+b)//q
    c = v2(3*Z+1)
    left = (1,)*(K-1)+(2,)*3+head+(c,)
    right = (1,)*(K-5)+(2, 2, 1, 1)+partner+(c+2*gap,)
    return audit_receipt(Receipt((1 << K)-1, (1 << (K-4))-1, left, right, (3*Z+1) >> c))


def load_exit():
    data = json.loads(DATA.read_text(encoding='utf-8'))['paid_exit']
    E = data['exponent']
    integer(E, 2)
    u, v = tuple(data['source_head']), tuple(data['partner_head'])
    letters(u)
    letters(v)
    need(E == 18269 and len(u) == 24 and len(v) == 25, 'declared compact exit')
    p, q, b = anchor.carrier(u)
    pp, qq, bb = anchor.carrier(v)
    need(q == 4*qq and pp == 3*p, 'one-bit deletion and one sibling gap')
    need(F(bb-pp, qq) == 4*F(b-p, q)+1, 'full compensated identity at minus one')
    residue = (q-b)*pow(p, -1, 2*q) % (2*q)
    A = sum(u)
    need(A == 50 and sum(v) == 48, 'exact head costs')
    k, period = anchor.log_three((residue+1)//2, A)
    need(k+1 == E and period == 1 << 48, 'native iff exponent phase')
    return E, u, v, period


def exit_guard(E):
    integer(E, 2)
    seed, _, _, period = load_exit()
    return (E-seed) % period == 0


def exit_receipt(E, bit_cap=100000):
    integer(bit_cap, 1)
    need(exit_guard(E) and E <= bit_cap, 'native exit guard and explicit materialization cap')
    _, u, v, _ = load_exit()
    x = 2*3**(E-1)-1
    p, q, b = anchor.carrier(u)
    need((p*x+b) % q == 0, 'integral source head')
    Z = (p*x+b)//q
    need(Z > (1 << E)-1, 'least-phase growth excludes ROOT padding')
    c = v2(3*Z+1)
    return audit_receipt(Receipt((1 << E)-1, (1 << (E-1))-1,
                                 (1,)*(E-1)+u+(c,), (1,)*(E-2)+v+(c+2,),
                                 (3*Z+1) >> c))


def reject(function, *args):
    try:
        function(*args)
    except (ValueError, TypeError):
        need(True, 'hostile rejected')
        return
    need(False, 'hostile must reject')


def main(rediscover=False):
    parents = {
        6125: ((1, 10), (1, 1, 2, 2, 3), 1),
        18269: ((10,), (3, 1, 1, 3), 1),
        46941: ((9, 2), (3, 1, 1, 3, 1), 1),
        98477: ((1, 2, 9), (1,)*6, 3),
    }
    keep = {}
    for E, packet in parents.items():
        word = root_word(E)
        reverse_verify((1 << E)-1, word)
        if rediscover:
            need(discover(E) == word, 'explicit bounded rediscovery matches frozen evidence')
        receipt = parent_receipt(E+4, *packet)
        lifted = discharge(receipt, word)
        need(len(lifted) == len(word) and sum(lifted) == sum(word)+4, 'grounded parent rank/cost')
        record = frozen_record(E)
        print('GROUNDED childE/parentK/rank/childcost/parentcost/peakbits:',
              E, E+4, len(word), sum(word), sum(lifted), record['peak_bits'], flush=True)
        if E == 18269:
            keep[E] = word

    seed, u, v, period = load_exit()
    receipt = exit_receipt(seed)
    child_word = reverse_discharge(receipt, keep[seed])
    need(len(child_word) == len(keep[seed]) and sum(child_word) == sum(keep[seed])-1,
         'new child grounded by receipt transport, not another orbit search')
    need(receipt.child == (1 << 18268)-1, 'exact new discharged child')
    for t in (0, 1, 9, 10**30):
        E = seed+period*t
        need(exit_guard(E) and E % 16 == 13, 'all-height two-twos exit phase')
        K = E+4
        need(anchor.head_valuation_guard(K, 3), 'exact intersection with old parent phase')
    for t in range(64):
        K = 18273+729*period*t
        need(anchor.baseline_entry(K) is None and exit_guard(K-4), 'excluded two-stage source phase')
    print('PROVED paid exit E=18269mod2^48 -> E-1; head lengths24/25,costs50/48.')
    print('GROUNDED by reverse transport: E18268, rank', len(child_word), 'cost', sum(child_word))
    print('PROVED two-stage excluded ray: K=18273+729*2^48*t -> K-4 -> K-5.')

    # Exact target-mismatch hostiles: short J2 templates need not hit emitted children.
    from collatz_twoanchor_head_decoder_20261007b import decode
    for u, expected in (((1, 12), 27861), ((1, 1, 14), 177949)):
        p, q, b = anchor.carrier(u)
        A = sum(u)
        v = decode(len(u)+4, A-2, b+10*p+(1 << (A-2)))
        need(v is not None, 'authenticated changed-anchor short head')
        native = (q-b)*pow(p, -1, 2*q) % (2*q)
        modulus = 1 << (A+4)
        e, per = anchor.log_three((1+8*(native-1)*pow(3, -2, modulus)) % modulus, A+4)
        need(e+1 == expected, 'exact short-type native phase')
        if expected == 27861:
            need(expected % 16 == 5, 'wrong short head misses every emitted E13mod16 child')
        else:
            bank = general.finite_bank()
            need(all((expected+4-general.phase(h, 3).residue)
                     % min(per, general.phase(h, 3).period) != 0 for h in bank),
                 'all225 parent phases miss this separate J2 head')
    print('HOSTILES: (1,12) has E5mod16; (1,1,14) phase misses all225 emitted parent phases.')
    for bad in (True, 6125.0, 1, -1):
        reject(root_word, bad)
    reject(verify_word, 3, (1, 4, 2))
    reject(verify_word, 3, (1, True))
    reject(audit_receipt, replace(receipt, child=1.0))
    reject(audit_receipt, replace(receipt, endpoint=True))
    reject(exit_guard, 18269.0)
    reject(exit_receipt, 18285)
    reject(exit_receipt, 18269, 1000)
    need(not exit_guard(18285), 'nearby same-shell source is not silently replaced')
    print('No universal coverage: finite ROOT evidence and infinite paid schemas remain distinct.')
    print('Checks:', CHECKS+anchor.CHECKS+general.CHECKS)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--rediscover', action='store_true',
                        help='explicitly rerun the recorded bounded orbit discoveries')
    main(parser.parse_args().rediscover)
