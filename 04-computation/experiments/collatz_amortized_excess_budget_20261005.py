"""Finite source recognition with retained excess credit and checked seed leaves.

No target ROOT word enters recognize(). The finite bank is explicit and checked.
Failed policy membership is not a nonconvergence assertion.
"""
from dataclasses import dataclass
from fractions import Fraction as F

from collatz_floor_transport_deadlines_20261005 import step, word_weight
from collatz_backward_measurement_compiler_20261005 import deadline_floor, measurement_plan

CHECKS = 0
HARD_47 = (1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3, 1, 1,
           2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3)
SEED_23 = (1, 1, 5, 4)


def check(ok, witness):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(witness)


def natural(value, least=0):
    if type(value) is not int or value < least:
        raise ValueError('exact integer outside domain')


def source(value):
    natural(value, 1)
    if value % 2 == 0:
        raise ValueError('positive odd source required')


def valuation_word(word):
    if type(word) is not tuple or any(type(a) is not int or a < 1 for a in word):
        raise ValueError('exact positive valuation tuple required')


def replay(n, word):
    source(n)
    valuation_word(word)
    for a in word:
        if n == 1:
            raise ValueError('ROOT padding forbidden')
        n, actual = step(n)
        if actual != a:
            raise ValueError('valuation/source guard mismatch')
    return n


def checked_bank(bank):
    if type(bank) is not dict or not bank:
        raise ValueError('nonempty explicit seed dictionary required')
    result = {}
    for n, word in bank.items():
        if replay(n, word) != 1:
            raise ValueError('seed word does not terminate at ROOT')
        result[n] = word
    return result


def balance(length, cost, q):
    natural(length)
    natural(cost)
    natural(q, 1)
    e = q*cost-length
    if e < 0:
        return F(1, (1 << -e)*3**(q*length))
    return F(1 << e, 3**(q*length))


@dataclass(frozen=True)
class Summary:
    multiplier: F = F(1)
    minimum: F = F(1)


def summary_valid(item):
    if (type(item) is not Summary or type(item.multiplier) is not F or
            type(item.minimum) is not F or not 0 < item.minimum <= 1 or
            item.minimum > item.multiplier):
        raise ValueError('positive prefix-balance summary required')


def concatenate(left, right):
    summary_valid(left)
    summary_valid(right)
    return Summary(left.multiplier*right.multiplier,
                   min(left.minimum, left.multiplier*right.minimum))


def block_summary(word, q):
    valuation_word(word)
    natural(q, 1)
    value = balance(len(word), sum(word), q)
    return Summary(value, min(F(1), value))


def initial_credit(item):
    summary_valid(item)
    b = max(0, item.minimum.denominator.bit_length()-item.minimum.numerator.bit_length())
    while (1 << b)*item.minimum < 1:
        b += 1
    while b and (1 << (b-1))*item.minimum >= 1:
        b -= 1
    return b


def source_cap(n, q, credit=0, lowest_seed=1):
    source(n)
    natural(q, 1)
    natural(credit)
    source(lowest_seed)
    bits = max(0, n.bit_length()-lowest_seed.bit_length())
    if n > lowest_seed << bits:
        bits += 1
    return 207*(q*bits+credit)//197


def recognize(n, q, bank, credit=0, threshold=None, policy='cumulative'):
    """Decide the declared guarded family within a source-only odd-step cap.

    A cumulative receipt checks the running product at every first-descent
    boundary. The individual control checks each factor separately (credit=0).
    All seed receipts are supplied, checked finite objects, not an orbit oracle.
    """
    source(n)
    natural(q, 1)
    natural(credit)
    if threshold is None:
        threshold = 10*q
    natural(threshold, 10*q)
    if type(policy) is not str or policy not in ('cumulative', 'individual'):
        raise ValueError('unknown policy')
    if policy == 'individual' and credit:
        raise ValueError('individual control has no cumulative initial credit')
    seeds = checked_bank(bank)
    cap = source_cap(n, q, credit, min(seeds))
    current, segment_start = n, n
    word, local_word, segments = [], [], []
    item = Summary()

    def answer(status, member=False):
        out = {'status': status, 'member': member, 'source': n, 'endpoint': current,
               'prefix': tuple(word), 'segments': tuple(segments), 'summary': item,
               'initial_credit_required': initial_credit(item), 'search_cap': cap,
               'new_odd_edges': len(word)}
        if member:
            suffix = seeds[current]
            out.update(word=tuple(word)+suffix, seed=current,
                       deadline=cap+len(suffix), seed_suffix_edges=len(suffix))
        return out

    if n in seeds:
        return answer('checked_seed', True)
    if current < threshold:
        return answer('uncovered_small_start')
    for _ in range(cap):
        current, a = step(current)
        word.append(a)
        local_word.append(a)
        if current >= segment_start:
            continue
        piece = tuple(local_word)
        factor = block_summary(piece, q)
        item = concatenate(item, factor)
        segments.append((segment_start, piece, current))
        passes = ((1 << credit)*item.minimum >= 1 if policy == 'cumulative'
                  else factor.multiplier >= 1)
        if not passes:
            return answer('credit_guard_failed')
        if current in seeds:
            return answer('grounded_at_checked_seed', True)
        if current < threshold:
            return answer('uncovered_small_start')
        segment_start, local_word = current, []
    return answer('outside_family_at_proved_deadline')


def funded_measurement(n, q, bank, credit=0, threshold=None):
    receipt = recognize(n, q, bank, credit, threshold)
    if not receipt['member']:
        return receipt
    if n == 1:
        return {**receipt, 'weight_floor': F(1)}
    plan = measurement_plan(n, receipt['deadline'])
    return {**receipt, 'weight_floor': deadline_floor(n, receipt['deadline'])['floor'],
            'measurement': {**plan, 'status': 'grounded_by_cumulative_receipt'}}


def inspect_prefix(n, word, threshold):
    """Validate a supplied finite chain; endpoint need not yet be a ROOT seed."""
    source(n)
    natural(threshold, 1)
    valuation_word(word)
    current, start, pieces, piece = n, n, [], []
    for a in word:
        if start < threshold:
            raise ValueError('segment start below threshold')
        current = replay(current, (a,))
        piece.append(a)
        if current < start:
            pieces.append(tuple(piece))
            start, piece = current, []
    if piece or not pieces:
        raise ValueError('prefix must end at a first-descent boundary')
    return current, tuple(pieces)


def refinance(anchor, word, q, bank, threshold=None):
    """One inverse edge finances a fixed checked excursion, for all phase lifts.

    Output a0 means every a=a0+2t, t>=0, gives an accepted new source.
    This compiles an inherited basin, not a new universal convergence claim.
    """
    source(anchor)
    natural(q, 1)
    if anchor == 1 or anchor % 3 == 0:
        raise ValueError('nonroot 3-unit anchor required')
    if threshold is None:
        threshold = 10*q
    natural(threshold, 10*q)
    seeds = checked_bank(bank)
    endpoint, pieces = inspect_prefix(anchor, word, threshold)
    if endpoint not in seeds:
        raise ValueError('supplied excursion has an undischarged endpoint')
    item = Summary()
    for piece in pieces:
        item = concatenate(item, block_summary(piece, q))
    a = 2 if anchor % 3 == 1 else 3
    while balance(1, a, q)*item.minimum < 1:
        a += 2
    return {'anchor': anchor, 'a0': a, 'period': 2, 'q': q,
            'threshold': threshold, 'word': word, 'seed': endpoint,
            'summary': item, 'odd_rank': 1+len(word)+len(seeds[endpoint])}


def ray47_budget(n):
    """Closed-form source budget for the fixed checked 47-to-23 excursion.

    Exact membership of the inverse ray is tested before returning a budget.
    No iteration or target ROOT certificate is read here.
    """
    source(n)
    numerator = 3*n+1
    if numerator % 47:
        return None
    power = numerator//47
    if power & (power-1):
        return None
    a = power.bit_length()-1
    if a < 3 or a % 2 == 0:
        return None
    return {'source': n, 'a': a, 'q': 3, 'credit': max(0, 37-3*a),
            'anchor': 47, 'seed': 23}


def observed_route(n, cap=256):
    """Finite audit helper only; never called by the production recognizer."""
    source(n)
    natural(cap)
    result = []
    for _ in range(cap):
        if n == 1:
            return tuple(result)
        n, a = step(n)
        result.append(a)
    if n == 1:
        return tuple(result)
    raise ValueError('audit cap exhausted')


def rejects(call, label):
    try:
        call()
    except ValueError:
        check(True, label)
    else:
        check(False, label)


def main():
    bank = {1: (), 23: SEED_23}
    check(replay(23, SEED_23) == 1, 'sole nonroot seed23')
    endpoint, pieces = inspect_prefix(47, HARD_47, 30)
    check(endpoint == 23 and pieces == (HARD_47,), '47 first descent')
    check((len(HARD_47), sum(HARD_47)) == (34, 55), 'hard type')
    check(balance(34, 55, 30) < 1 <= balance(34, 55, 31), 'least q31')
    plan = refinance(47, HARD_47, 3, bank)
    check(plan['a0'] == 13 and plan['odd_rank'] == 39, 'least financing phase')
    print('Exact funding family: n=(47*2^(13+2t)-1)/3, t>=0; q3, credit0, bank{1,23}.')
    print('  47 -> 23: length34/cost55, least individual rate31; source ROOT rank39.')
    for a in range(13, 142, 2):
        n = (47*(1 << a)-1)//3
        result = recognize(n, 3, bank)
        direct = recognize(n, 3, bank, policy='individual')
        check(result['member'] and not direct['member'], ('strict extension', a))
        check(result['word'] == (a,)+HARD_47+SEED_23, ('word', a))
        check(replay(n, result['word']) == 1, ('independent replay', a))
        check(len(result['word']) <= result['deadline'], ('deadline', a))
        check((1 << (3*a+130)) >= 3**105, ('all-height base inequality', a))
    near = recognize((47*(1 << 11)-1)//3, 3, bank)
    rescued = recognize((47*(1 << 11)-1)//3, 3, bank, credit=4)
    check(not near['member'] and rescued['member'], 'declared initial debt changes family')
    check(near['initial_credit_required'] == 4, 'exact missing credit')
    check(replay(rescued['source'], rescued['word']) == 1, 'initial-credit route')
    witness = funded_measurement(128341, 3, bank)
    check(witness['deadline'] == 57, 'source-only deadline57')
    check(witness['weight_floor'] <= word_weight(witness['word']), 'grounded source floor')
    leaf = witness['measurement']
    leaf_word = ((leaf['section_valuation'],) if leaf['section_valuation'] else ())+witness['word']
    check(replay(leaf['leaf'], leaf_word) == 1, 'selected-leaf transport')
    check(leaf['floor'] <= word_weight(leaf_word), 'selected-leaf floor')
    check(leaf['measurement_lower'] > 0, 'positive actual readout guarantee')
    print('  least source128341: cumulative deadline57; q31 bit-bound1058; exact structural rank39.')
    targeted = funded_measurement(128341, 3, {23: SEED_23})
    check(targeted['member'] and targeted['deadline'] == 44,
          'known bank minimum gives a sharper source-only deadline')
    check(targeted['weight_floor'] <= word_weight(targeted['word']), 'targeted floor')
    print('  singleton bank{23} gives deadline44, using its actual minimum endpoint; threshold30 cannot be substituted.')
    print('  selected leaf=%d, degree=%d; source floor=%s.' %
          (leaf['leaf'], leaf['degree'], witness['weight_floor']))
    print('  a11 source32085 fails credit0, requires exactly4 initial binary credit units.')
    check(2**166 < 3**105 < 2**167, 'exact logarithmic rounding for source credit')
    for a in range(3, 34, 2):
        n = (47*(1 << a)-1)//3
        budget = ray47_budget(n)
        check(budget is not None and budget['a'] == a and
              budget['credit'] == max(0, 37-3*a), 'source-only phase and credit')
        item = recognize(n, 3, {23: SEED_23}, credit=budget['credit'])
        check(item['member'] and replay(n, item['word']) == 1, 'adaptive source budget')
        check(item['initial_credit_required'] == budget['credit'], 'least credit equality')
        if a <= 11:
            check(item['deadline'] == 42, 'small-phase deadline42')
        if budget['credit']:
            insufficient = recognize(n, 3, {23: SEED_23}, credit=budget['credit']-1)
            check(not insufficient['member'], 'one less credit fails')
        check(ray47_budget(n+2) is None, 'nearby source is not the same inverse ray')
    check(ray47_budget(31) is None, 'a1 grows to47 and is not the decreasing ray')
    rejects(lambda: ray47_budget(True), 'source-budget exact type')
    print('  source-only ray test:3n+1=47*2^a, odd a>=3; least q3 credit=max(0,37-3a).')
    print('  phases a3,5,7,9,11 have singleton-bank deadline42; fixed shared excursion/seed costs retained.')
    fixed_one = refinance(47, HARD_47, 1, bank)
    check(fixed_one['a0'] == 37, 'fixed rate1 financing')
    for a in range(37, 54, 2):
        n = (47*(1 << a)-1)//3
        one = recognize(n, 1, bank)
        check(one['member'] and replay(n, one['word']) == 1, 'fixed rate1 family')
    print('  fixed q1 also works on the phase a=37+2t, with the same finite excursion and seed.')

    # Summary composition is exact; all actual-source guards remain separate.
    left, right = block_summary((13,), 3), block_summary(HARD_47, 3)
    check(right.multiplier < 1 and concatenate(left, right).minimum == 1,
          'history reset destroys valid funding')
    pieces0 = ((1, 2), (4,), HARD_47, (13,))
    associative = 0
    for u in pieces0:
        for v in pieces0:
            for w in pieces0:
                x, y, z = (block_summary(t, 3) for t in (u, v, w))
                check(concatenate(concatenate(x, y), z) == concatenate(x, concatenate(y, z)),
                      'prefix-summary associativity')
                associative += 1
    check(block_summary((1, 3), 3) == block_summary((3, 1), 3), 'balance is not a source guard')
    check(replay(19, (1, 3)) == 11 and replay(29, (3, 1)) == 17,
          'distinct lawful sources with the same balance')
    rejects(lambda: replay(19, (3, 1)), 'same totals cannot forge a source word')

    # Same finite bank, source universe and rate; policy cost is counted separately.
    small_bank = {n: observed_route(n, 128) for n in range(1, 30, 2)}
    counts, queries = {}, {}
    configurations = [('individual', 0), ('cumulative', 0), ('cumulative', 4), ('cumulative', 16)]
    accepted = {}
    for policy, credit in configurations:
        tag = '%s/B%d' % (policy, credit)
        members, edges = set(), 0
        for n in range(31, 2048, 2):
            item = recognize(n, 3, small_bank, credit, policy=policy)
            edges += item['new_odd_edges']
            if item['member']:
                members.add(n)
                check(replay(n, item['word']) == 1, ('finite independent replay', tag, n))
                check(len(item['word']) <= item['deadline'], ('finite deadline', tag, n))
        accepted[tag] = members
        counts[tag], queries[tag] = len(members), edges
    check(accepted['individual/B0'] <= accepted['cumulative/B0'] <=
          accepted['cumulative/B4'] <= accepted['cumulative/B16'], 'nested controls')
    print('Finite universe:1009 odd sources31..2047; explicit15-seed bank1..29, max seed rank=%d.' %
          max(map(len, small_bank.values())))
    print('  membership counts=%s' % counts)
    print('  newly queried odd edges by policy=%s; seed verification separate.' % queries)
    extra = sorted(accepted['cumulative/B0']-accepted['individual/B0'])
    check(extra[0] == 203, 'least finite strict-extension source')
    example = recognize(203, 3, small_bank)
    check(example['segments'] == ((203, (1, 2, 4), 43), (43, (1, 2, 2), 37),
                                  (37, (4,), 7)), 'retained funding example203')
    print('  least extra203:124 ->43,122 ->37,4 -> checked seed7; first ten extras=%s.' % extra[:10])
    print('  summary associativity:%d triples; generic summaries do not supply source legality.' % associative)

    # Bounded failure is a policy verdict, even for inputs known independently to converge.
    for n in (27, 55, 32085):
        item = recognize(n, 3, bank)
        check(not item['member'], ('policy miss', n))
        check(replay(n, observed_route(n)) == 1, ('policy miss still ROOT', n))
    for invalid in (True, 1.0, 0, 2):
        rejects(lambda invalid=invalid: recognize(invalid, 3, bank), 'source type')
    rejects(lambda: recognize(128341, True, bank), 'rate type')
    rejects(lambda: recognize(128341, 3, bank, credit=-1), 'negative credit')
    rejects(lambda: recognize(128341, 3, bank, threshold=29), 'low threshold')
    rejects(lambda: recognize(128341, 3, {23: (1, 1, 5)}), 'ungrounded bank')
    rejects(lambda: recognize(128341, 3, {1: (2,)}), 'root self proof')
    rejects(lambda: refinance(27, observed_route(27), 3, bank), '3-divisible anchor')
    rejects(lambda: refinance(47, HARD_47[:-1], 3, bank), 'unfinished excursion')
    rejects(lambda: concatenate(Summary(F(1), F(2)), Summary()), 'forged summary')
    check(F(207, 197) < F(1051, 1000), 'repaired decimal deadline')
    print('PASS: %d exact checks; no assertion statements; no target ROOT words supplied to recognition.' % CHECKS)


if __name__ == '__main__':
    main()
