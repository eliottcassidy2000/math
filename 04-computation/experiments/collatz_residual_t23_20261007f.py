"""A source-owned paid exit at the previously missing parameter 23.

Frozen heads use a shared suffix without dropping any child alternative.
The anchored phase compiler avoids solving a large discrete logarithm anew.
No production function searches for, assumes, or manufactures a ROOT word.
"""
from dataclasses import dataclass, replace
from hashlib import sha256
from pathlib import Path
import json
import sys

import collatz_parameter_cover_20261007e as prior
import collatz_eight_child_routes_20261007d as reader
import collatz_context_cancellation_20261007d as context
import collatz_two_twos_grounding_20261007c as streaming

DATA = Path(__file__).resolve().parents[2] / '05-knowledge/results/collatz_residual_t23_20261007f.json'
BASE, STRIDE = prior.BASE, prior.STRIDE
SEED, DEPTH, MAX_DROP, PRECISION = 23, 4096, 256, 16384
DROPS = (3, 4) + tuple(range(7, 29))
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def natural(n, least=0):
    need(type(n) is int and n >= least, 'exact integer in declared domain')


@dataclass(frozen=True)
class Anchored:
    macro: context.Macro
    seed: int
    residue: int
    bits: int
    minimum: int


def compile_anchor(macro, seed):
    """A supplied native source fixes the unique phase by the order of 3."""
    context.audit(macro)
    natural(seed)
    cost = sum(macro.left)
    need(cost >= 34, 'current parameter chart is visible to this head')
    residue, modulus = context.cell(macro)
    E = BASE + STRIDE*seed
    need((2*pow(3, E-1, modulus)-1) % modulus == residue,
         'supplied anchor has the exact native head')
    bits = cost-34
    cutoff = max(0, (max(2*cost, macro.drop+1)+1-BASE+STRIDE-1)//STRIDE)
    need(seed >= cutoff, 'anchor meets all-prefix growth boundary')
    return Anchored(macro, seed, seed % (1 << bits), bits, cutoff)


def audit(packet):
    need(type(packet) is Anchored, 'exact anchored packet')
    for name in ('seed', 'residue', 'bits', 'minimum'):
        natural(getattr(packet, name))
    need(compile_anchor(packet.macro, packet.seed) == packet,
         'recomputed native phase and immutable source cutoff')
    return packet


def contains(packet, parameter):
    audit(packet)
    natural(parameter)
    return parameter >= packet.minimum and parameter % (1 << packet.bits) == packet.residue


def digest(word):
    context.routes.letters(word)
    need(all(a < 256 for a in word), 'declared byte-valuation digest')
    return sha256(bytes(word)).hexdigest()


def load():
    data = json.loads(DATA.read_text(encoding='utf-8'))
    need(data['format'] == 'anchored-mersenne-shared-suffix-v1', 'retained data format')
    need(data['search'] == dict(seed=SEED, max_depth=DEPTH, max_drop=MAX_DROP,
                               precision=PRECISION), 'declared finite proposal box')
    left, suffix = tuple(data['source_head']), tuple(data['partner_suffix'])
    need(digest(left) == data['source_sha256'], 'source head identity')
    answer = []
    for item in data['partners']:
        right = tuple(item['prefix'])+suffix
        need(digest(right) == item['sha256'], 'lossless partner suffix decoding')
        answer.append(compile_anchor(context.Macro(item['drop'], left, right, 1), SEED))
    need(tuple(p.macro.drop for p in answer) == DROPS, 'retain every labelled child')
    need(len(suffix) == 3014, 'declared common suffix compression')
    return tuple(answer)


def discover():
    """Optional bounded proposal replay, separate from the certificate consumer."""
    E = BASE+STRIDE*SEED
    u, source, _ = reader.modular_prefix(E, DEPTH+MAX_DROP+1, PRECISION)
    need(len(u) == DEPTH+MAX_DROP+1, 'sufficient source comparison precision')
    answer = []
    for drop in range(1, MAX_DROP+1):
        v, child, _ = reader.modular_prefix(E-drop, len(u), PRECISION)
        need(len(v) == len(u), 'sufficient child comparison precision')
        for depth in range(DEPTH+1):
            if source[depth][1] != child[depth+drop][1]:
                continue
            delta = source[depth][0]-child[depth+drop][0]
            if delta % 2 == 0 and v[depth+drop] == u[depth]+delta:
                need(delta == 2, 'every retained first collision has positive gap one')
                answer.append(compile_anchor(context.Macro(drop, u[:depth], v[:depth+drop], 1), SEED))
                break
    return tuple(answer)


def endpoint_mod(packet, parameter, modulus, child=False):
    need(contains(packet, parameter), 'same supplied parameter')
    natural(modulus, 1)
    need(type(child) is bool, 'exact side flag')
    m = packet.macro
    E = BASE+STRIDE*parameter-(m.drop if child else 0)
    p, q, b = context.routes.carrier(m.right if child else m.left)
    x = (2*pow(3, E-1, q*modulus)-1) % (q*modulus)
    need((p*x+b) % q == 0, 'division precision retained')
    return ((p*x+b)//q) % modulus


def literal_control(packet):
    """Stream both literal paths for a modest ordinary source, never M_E."""
    audit(packet)
    m = packet.macro
    run = max(2*sum(m.left), m.drop+1)
    source = context.fixed_source(m, run)
    child = ((source+1) >> m.drop)-1
    cofactor = (source+1) >> (run+1)
    x = 2*3**run*cofactor-1
    p, q, b = context.routes.carrier(m.left)
    z = (p*x+b)//q
    terminal = streaming.v2(3*z+1)
    left = (1,)*run+m.left+(terminal,)
    right = (1,)*(run-m.drop)+m.right+(terminal+2,)
    a = streaming.verify_word(source, left, False)
    c = streaming.verify_word(child, right, False)
    need(a.endpoint == c.endpoint and child < source, 'literal paid common future')
    need(a.endpoint > 1, 'literal control leaves its ROOT port open')
    return source.bit_length(), a.rank, c.rank


def reject(function, *args):
    try:
        function(*args)
    except (ValueError, TypeError):
        need(True, 'hostile rejected')
        return
    need(False, 'hostile accepted')


def main(rediscover=False):
    packets = load()
    if rediscover:
        need(discover() == packets, 'complete bounded search reproduces all labels')
    need(all((p.residue, p.bits, p.minimum) == (23, 6518, 0) for p in packets),
         'same exact native guard for all alternatives')
    for p in packets:
        m = p.macro
        need(len(m.left) == 3317 and len(m.right) == 3317+m.drop, 'exact clocks')
        need(sum(m.left) == 6552 and sum(m.right) == 6550, 'exact costs')
        need(contains(p, 23) and not contains(p, 23+(1 << 6517)), 'sharp native precision')
        need(contains(p, 23+(1 << 6518)), 'whole infinite anchored phase')
        for t in (23,):
            for modulus in (3, 19, 64, 729):
                z = endpoint_mod(p, t, modulus)
                y = endpoint_mod(p, t, modulus, True)
                need(y == (4*z+1) % modulus, 'independent same-source modular common future')
        # The old entire residual progression now has a covered subphase;
        # its other binary half remains disjoint from every new word.
        need((13847-p.residue) % (1 << 10) != 0,
             't13847 mod27648 stays outside this new guard')
    # Cheap independent comparison to the inherited discrete-log compiler.
    old = json.loads(prior.DATA.read_text(encoding='utf-8'))['rules'][:12]
    for row in old:
        slow = prior.compile_rule(row)
        m = context.Macro(slow.drop, slow.source_head, slow.partner_head, slow.gap)
        fast = compile_anchor(m, slow.seed_parameter)
        need((fast.residue, fast.bits, fast.minimum) ==
             (slow.parameter_residue, slow.parameter_bits, slow.minimum_parameter),
             'independent discrete-log phase cross-check')
    controls = tuple(literal_control(p) for p in (packets[0], packets[-1]))
    # The inherited eight-bit receipt ends at the exact beginning of this head.
    outer = context.Macro(8, context.old.LEFT, context.old.RIGHT, 1)
    combined = context.cancel(outer, packets[-1].macro)
    need(combined.drop == 36 and sum(combined.left) == 6552, 'paid parent composition')
    # Known local run surgery applies here too; it is not recounted in this packet.
    for k in range(9):
        pumped = prior.extend_ones(packets[-1].macro, k)
        need(pumped.drop == 28 and sum(pumped.left) == 6552+k,
             'all-run structural continuation has the same marked child type')
    reject(compile_anchor, packets[0].macro, True)
    reject(compile_anchor, packets[0].macro, 24)
    reject(audit, replace(packets[0], residue=24))
    reject(compile_anchor, replace(packets[0].macro, right=packets[0].macro.right[:-1]+(2,)), 23)
    reject(contains, packets[0], -1)
    print('PROVED paid exit: E=924745897+2^32*t, t=23 mod2^6518.')
    print('Retained child deletions:', DROPS)
    print('Strongest child exponent: 99708993677; original parent pays36 bits.')
    print('Exact source/partner head costs6552/6550; clocks3317/(3317+D); gap1.')
    print('Shared partner suffix3014 letters; all24 child alternatives retained.')
    print('Ordinary-source stream controls (bits,left rank,right rank):', controls)
    print('Finite proposal box: one seed, D1..256, source depths0..4096,16384bits.')
    print('Anchored phase agrees with12 independent inherited discrete-log compilations.')
    print('OPEN: all original giant child ROOT ports; universal coverage and grounding.')
    print('Checks:', CHECKS)


if __name__ == '__main__':
    main('--rediscover' in sys.argv)
