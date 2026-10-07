"""Guarded reset-two collisions and frozen point certificates.

Production functions consume exact rules or supplied first-hit words. Only the
explicitly labelled discovery control in main performs a bounded ROOT search.
No import-time experiment or default file mutation.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import collatz_uncovered_join_routes_20261007 as routes
import collatz_terminal_lifts_20261007 as old


DATA = Path(__file__).resolve().parents[2] / '05-knowledge/results/collatz_reset2_rules_20261007.json'
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def natural(x, minimum=0):
    need(type(x) is int and x >= minimum, 'exact integer in declared domain')


def v2(x):
    natural(x, 1)
    return (x & -x).bit_length()-1


def fast_step(x):
    routes.odd(x)
    z = 3*x+1
    a = v2(z)
    return z >> a, a


@dataclass(frozen=True)
class Rule:
    name: str
    seed_exponent: int
    drop: int
    source_head: tuple
    partner_head: tuple


def load_rules():
    data = json.loads(DATA.read_text(encoding='utf-8'))
    return tuple(Rule(r['name'], r['seed_exponent'], r['drop'],
                      tuple(r['source_head']), tuple(r['partner_head'])) for r in data['rules'])


def audit_rule(rule):
    need(type(rule) is Rule and type(rule.name) is str and bool(rule.name), 'typed named rule')
    natural(rule.seed_exponent, 2)
    natural(rule.drop, 1)
    routes.letters(rule.source_head)
    routes.letters(rule.partner_head)
    need(bool(rule.source_head) and bool(rule.partner_head), 'nonempty reduced heads')
    need(rule.source_head[0] >= 2 and rule.partner_head[0] >= 2, 'reduced first letters')
    need(len(rule.partner_head)-len(rule.source_head) == rule.drop, 'exact run deletion')
    p, q, b = routes.carrier(rule.source_head)
    pp, qq, bb = routes.carrier(rule.partner_head)
    need(q == 4*qq, 'terminal valuation displacement two')
    need(Fraction(bb-pp, qq) == 4*Fraction(b-p, q)+1, 'exact compensated collision at minus one')
    residue = (q-b)*pow(p, -1, 2*q) % (2*q)
    cost = q.bit_length()-1
    need(cost >= 3, 'Mersenne phase uses the standard order of three')
    period = 2**(cost-2)
    need(rule.seed_exponent-1 >= rule.drop and rule.seed_exponent < period, 'least positive phase representative')
    x = 2*3**(rule.seed_exponent-1)-1
    need(x % (2*q) == residue, 'declared phase seed has the native head')
    need(routes.replay(x, rule.source_head)[-1] > 1, 'least phase member has no ROOT-padded head')
    return p, q, b, residue, period


def mersenne_guard(exponent, rule):
    """Exact all-height exponent guard; no expansion of 2**exponent."""
    natural(exponent, 2)
    _, _, _, _, period = audit_rule(rule)
    return exponent-1 >= rule.drop and (exponent-rule.seed_exponent) % period == 0


def receipt(source, rule):
    """Return a checked strictly-smaller common-future dependency, or None."""
    routes.odd(source)
    p, q, b, residue, _ = audit_rule(rule)
    if source == 1:
        return None
    run = v2(source+1)-1
    if run < rule.drop:
        return None
    cofactor = (source+1) >> run
    x_mod = (pow(3, run, 2*q)*cofactor-1) % (2*q)
    if x_mod != residue:
        return None
    x = 3**run*cofactor-1
    endpoint_head = (p*x+b)//q
    if endpoint_head == 1:
        return None
    endpoint, terminal = fast_step(endpoint_head)
    child = ((source+1) >> rule.drop)-1
    left = (1,)*run+rule.source_head+(terminal,)
    right = (1,)*(run-rule.drop)+rule.partner_head+(terminal+2,)
    return routes.audit(routes.Receipt(source, child, left, right, endpoint))


def fixed_run_source(run, lift, rule):
    """An infinite native source AP at a fixed initial run; exact construction."""
    natural(run, 0)
    natural(lift, 0)
    p, q, b, residue, _ = audit_rule(rule)
    need(run >= rule.drop, 'enough initial ones for the deletion')
    modulus = 3**run
    j0 = -(residue+1)*pow(2*q, -1, modulus) % modulus
    x = residue+2*q*j0
    if (p*x+b)//q == 1:
        j0 += modulus
    x = residue+2*q*(j0+modulus*lift)
    source = 2**run*((x+1)//modulus)-1
    need(receipt(source, rule) is not None, 'constructed native source has the receipt')
    return source


def frozen_seed_word():
    """Decode and authenticate the one stored point certificate; no search."""
    info = json.loads(DATA.read_text(encoding='utf-8'))['terminal']
    for field in ('exponent', 'leading_ones', 'odd_rank', 'valuation_sum', 'peak_bits'):
        natural(info[field], 1)
    need(type(info['tail_hex']) is str and all(c in '123456789abcdef' for c in info['tail_hex']), 'positive hexadecimal valuation data')
    word = (1,)*info['leading_ones']+tuple(int(c, 16) for c in info['tail_hex'])
    states = routes.replay(2**info['exponent']-1, word)
    need(states[-1] == 1, 'stored certificate reaches first ROOT')
    need(len(word) == info['odd_rank'] and sum(word) == info['valuation_sum'], 'stored word rank and cost')
    need(max(x.bit_length() for x in states) == info['peak_bits'], 'stored peak bit length')
    return info['exponent'], word


def reverse_discharge(common, supplied_source_word):
    """Transport an authenticated source ROOT word to its common-future child."""
    routes.audit(common)
    need(routes.replay(common.source, supplied_source_word)[-1] == 1, 'supplied source first-hit proof')
    cut = len(common.source_word)
    need(supplied_source_word[:cut] == common.source_word, 'actual source prefix is retained')
    result = common.child_word+supplied_source_word[cut:]
    need(routes.replay(common.child, result)[-1] == 1, 'transported child first-hit proof')
    return result


def grounded_named_points():
    """Exactly three point certificates, using one frozen proof and two joins."""
    rules = load_rules()
    exponent, word = frozen_seed_word()
    need(exponent == 1457, 'immutable selected seed')
    direct = receipt(2**1459-1, rules[1])
    source_word = routes.discharge(direct, word)
    deeper = receipt(2**1457-1, rules[2])
    child_word = reverse_discharge(deeper, word)
    return {1459: source_word, 1457: word, 1451: child_word}


def discovery_control(exponent, cap):
    """Explicit bounded research control; never called by production routines."""
    natural(exponent, 2)
    natural(cap, 1)
    n = 2**exponent-1
    word = []
    while n != 1 and len(word) < cap:
        n, a = fast_step(n)
        word.append(a)
    need(n == 1, 'declared bounded discovery completed')
    return tuple(word)


def reduced_trace(exponent, cap):
    natural(exponent, 2)
    natural(cap, 0)
    n = 2*3**(exponent-1)-1
    states, word = [n], []
    for _ in range(cap):
        if n == 1:
            break
        n, a = fast_step(n)
        states.append(n)
        word.append(a)
    return tuple(states), tuple(word)


def parallel_join_control(exponent, maximum_drop=64, tail_cap=256):
    """Declared finite collision probe, retaining the first equal-time join."""
    states, word = reduced_trace(exponent, tail_cap)
    found = []
    for drop in range(1, maximum_drop+1):
        partner, letters = reduced_trace(exponent-drop, tail_cap+drop)
        for depth in range(1, min(len(states), len(partner)-drop)):
            if states[depth] == partner[depth+drop]:
                u, v = word[:depth], letters[:depth+drop]
                p, q, b = routes.carrier(u)
                pp, qq, bb = routes.carrier(v)
                need(Fraction(b-p, q) == Fraction(bb-pp, qq), 'integer join authenticates a uniform collision')
                found.append((depth, drop))
                break
    return tuple(sorted(found))


def rejects(operation):
    try:
        operation()
    except (ValueError, TypeError):
        return True
    return False


def main():
    rules = load_rules()
    need(tuple(r.name for r in rules) == ('1459-D1', '1459-D2', '1457-D6'), 'frozen rule inventory')
    print('PROVED conditional collision families; inherited THM4555 minus-one switch mechanism.')
    for rule in rules:
        p, q, b, residue, period = audit_rule(rule)
        for c in range(1, 13):
            pp, qq, bb = routes.carrier(rule.source_head+(c,))
            vp, vq, vb = routes.carrier(rule.partner_head+(c+2,))
            need(Fraction(bb-pp, qq) == Fraction(vb-vp, vq), 'all terminal-letter controls')
        for run in range(rule.drop, rule.drop+6):
            for lift in range(5):
                n = fixed_run_source(run, lift, rule)
                result = receipt(n, rule)
                need(result.child == (n+1)//2**rule.drop-1, 'fixed-run AP child')
        for lift in range(32):
            exponent = rule.seed_exponent+period*lift
            need(mersenne_guard(exponent, rule), 'symbolic all-height exponent membership')
            need((2*pow(3, exponent-1, q)-1) % (2*q) == residue, 'independent modular phase replay')
            need(not mersenne_guard(exponent+2, rule), 'adjacent odd exponent does not satisfy the phase')
        print('RULE', rule.name, 'head lengths/costs', len(rule.source_head), sum(rule.source_head),
              len(rule.partner_head), sum(rule.partner_head), 'native residue', residue,
              'modulus bits', q.bit_length(), 'K phase', rule.seed_exponent,
              'mod 2^'+str(period.bit_length()-1))

    need(rules[1].partner_head == (2, 2)+rules[0].partner_head[1:], 'D2 uses the exact minus-one expansion of first letter4')
    for rule in (rules[1], rules[2]):
        _, _, _, _, period = audit_rule(rule)
        for j in range(20):
            exponent = rule.seed_exponent+729*period*j
            need(exponent % 1458 == rule.seed_exponent % 1458, 'ternary baseline exclusion is retained')
            expected = 1 if rule.seed_exponent == 1459 else 1093
            need((pow(2, exponent, 2187)-1) % 2187 == expected, 'baseline inverse-family exclusion residue')
        source = 2**rule.seed_exponent-1
        selection = routes.prior.select(source, 8, 'adaptive')
        need(selection.status == 'PENDING' and selection.word == (1,)*1024, 'frozen eight-by128 budget is exhausted')
        need(not routes.new_receipts(source, True, selection), 'both inverse rows and final frontier fail')
        need(old.reset_receipt(source) is None and old.sporadic_receipt(source) is None, 'both old uncapped rules fail at the declared source')
        need(receipt(source, rule) is not None, 'new receipt is paid despite frozen baseline failure')
    print('BASELINE: eight adaptive macros, cap128, inverse91/111 and frontier; two old uncapped rules. Both named seeds miss it.')
    print('Baseline-disjoint rays: K=1459+729*2^21*s and K=1457+729*2^126*s, s>=0.')

    joins = parallel_join_control(1457)
    need(joins[:6] == ((62, 3), (62, 4), (62, 5), (62, 6), (81, 1), (81, 2)), 'first six adaptive collision depths')
    print('FINITE-EXACT K1457, D1..64, reduced source depth<=256: first joins', joins[:6], 'total drops with a join', len(joins))
    for m in range(1, 11):
        exponent = 1+2**m
        count = m//2+1
        _, w = reduced_trace(exponent, count+1)
        need(w[:count] == (2,)*count, 'unbounded reduced all-two prefix')
        need(w[count] == 1 if m % 2 == 0 else w[count] >= 3, 'first non-two letter parity split')
    print('PROVED residual sidecar: for m=v2(K-1), initial reduced2-run has floor(m/2)+1 letters.')

    exponent, seed = frozen_seed_word()
    need(discovery_control(exponent, 100000) == seed, 'independent explicitly bounded selected-source discovery reproduces stored data')
    grounded = grounded_named_points()
    expected_costs = {1459: 13104, 1457: 13102, 1451: 13096}
    for exponent, word in grounded.items():
        need(len(word) == 7347 and sum(word) == expected_costs[exponent], 'transported first-hit ranks/costs')
        path = routes.replay(2**exponent-1, word)
        need(path[-1] == 1 and all(n > 1 for n in path[:-1]), 'strict ROOT certificate')
        print('ROOT point exponent', exponent, 'odd steps', len(word), 'valuation cost', sum(word),
              'peak bits', max(n.bit_length() for n in path), 'word SHA256', sha256(bytes(word)).hexdigest())
    need(all(not mersenne_guard(1451, rule) for rule in rules), 'fixed guard bank is not universal, even for the newly certified child')
    print('Frozen JSON SHA256', sha256(DATA.read_bytes()).hexdigest())

    short = rules[1]
    hostile = (
        lambda: receipt(True, short), lambda: receipt(5.0, short), lambda: receipt(2, short),
        lambda: mersenne_guard(True, short), lambda: mersenne_guard(1459.0, short),
        lambda: audit_rule(replace(short, drop=1)),
        lambda: audit_rule(replace(short, drop=2.0)),
        lambda: audit_rule(replace(rules[0], drop=True)),
        lambda: audit_rule(replace(short, seed_exponent=1459.0)),
        lambda: audit_rule(replace(short, source_head=(2.0,)+short.source_head[1:])),
        lambda: audit_rule(replace(short, source_head=(2,)+short.source_head[1:] + (2,))),
        lambda: reverse_discharge(receipt(2**1457-1, rules[2]), seed[:-1]+(seed[-1]+1,)),
        lambda: routes.replay(1, (2,)),
        lambda: fixed_run_source(1, 0, short),
    )
    need(all(rejects(test) for test in hostile), 'type, corruption, unpaid-run and ROOT-padding controls')
    need(receipt(1, short) is None and receipt(3, short) is None, 'ROOT and insufficient-run boundaries')
    print('HOSTILES', len(hostile), 'rejected; arbitrary family children remain supplied obligations.')
    print('CHECKS', CHECKS)


if __name__ == '__main__':
    main()
