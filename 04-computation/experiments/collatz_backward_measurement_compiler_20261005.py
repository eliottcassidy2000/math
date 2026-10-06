"""Compile a proved odd-step deadline into a quantitative measurement plan.

The compiler is conditional on the supplied deadline. Its arithmetic does not
prove that deadline. Literal replay is an optional separate finite checker.
"""
from fractions import Fraction as F
from math import comb
import json
from collatz_floor_transport_deadlines_20261005 import step, word_weight, control_route


def natural(n, least=0):
    if type(n) is not int or n < least:
        raise ValueError('exact natural integer required')


def source(n):
    natural(n, 3)
    if n % 2 == 0:
        raise ValueError('nonroot positive odd source required')


def deadline_floor(n, time):
    """Conditional rational W-floor; no route or stored certificate is read."""
    source(n)
    natural(time, 1)
    halving_ceiling = (n * 10**time // 3**time).bit_length()-1
    counter_ceiling = (halving_ceiling+time-1)//2-1
    if counter_ceiling < 1:
        raise ValueError('deadline incompatible with necessary ROOT counter bounds')
    c = counter_ceiling
    floor = F(2, (c+2)*comb(c+1, (c+1)//2))
    return {'status': 'conditional_on_deadline', 'source': n, 'odd_step_deadline': time,
            'halving_ceiling': halving_ceiling, 'counter_ceiling': c, 'floor': floor}


def degree_for_floor(floor):
    if type(floor) is not F or not 0 < floor <= 1:
        raise ValueError('positive exact probability floor required')
    degree, error = 0, F(8)
    while error > floor/2:
        degree += 1
        error *= F(16, 25)
    return degree, error


def selected_leaf(n):
    source(n)
    if n % 3 == 0:
        return n, 0, None
    a = {1:6, 2:5, 4:4, 5:1, 7:2, 8:3}[n % 9]
    leaf = ((1 << a)*n-1)//3
    if leaf % 6 != 3 or step(leaf) != (n, a):
        raise ValueError('section guard failed')
    return leaf, 1, a


def measurement_plan(n, time):
    """One canonical plan for this particular supplied deadline, still conditional."""
    source(n)
    natural(time, 1)
    leaf, extra, a = selected_leaf(n)
    bound = deadline_floor(leaf, time+extra)
    degree, error = degree_for_floor(bound['floor'])
    return {**bound, 'original_source': n, 'original_deadline': time,
            'leaf': leaf, 'index': (leaf-3)//6, 'section_valuation': a,
            'degree': degree, 'selector_error': error,
            'measurement_lower': bound['floor']-error}


def check_deadline(n, time):
    """Finite literal validation of a candidate; no claim about later ROOT visits."""
    source(n)
    natural(time, 1)
    current, word = n, []
    for used in range(time+1):
        if current == 1:
            return {'met': True, 'word': tuple(word)}
        if used == time:
            break
        current, a = step(current)
        word.append(a)
    return {'met': False, 'prefix': tuple(word), 'endpoint': current}


def hybrid_measurement_plan(n):
    """Discharge the deadline premise by a proved proper-family recognizer.

    No stored ROOT word or measured moment is supplied. Membership is a
    structural convergence certificate, and is not asserted for all sources.
    """
    source(n)
    from collatz_inductive_floor_receipts_20261005 import recognize_hybrid
    if not recognize_hybrid(n)['member']:
        return {'status': 'outside_proved_family', 'source': n}
    plan = measurement_plan(n, 6*n.bit_length()+1)
    return {**plan, 'status': 'proved_on_hybrid_family',
            'deadline_provenance': 'canonical mixed-tree rank theorem',
            'certified_readout_floor': plan['floor']/2}


def main():
    checks = 0

    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            raise ValueError(label)

    # Independent finite ROOT routes test the conditional algebra, not its premise.
    examples = []
    for n in range(3, 1024, 2):
        word = control_route(n, cap=512)
        tau, total = len(word), sum(word)
        even_count = sum(a % 2 == 0 for a in word)
        l, k = tau-1, sum((a-1)//2 for a in word)
        exact = word_weight(word)
        check(even_count >= 1, 'final valuation is even')
        check(2*(l+k+1) == total+tau-even_count, 'exact counter conversion')
        check((1 << total)*3**tau <= n*10**tau, 'telescoped halving budget')
        old = None
        for extra in (0, 1, 5):
            bound = deadline_floor(n, tau+extra)
            check(total <= bound['halving_ceiling'], 'halving ceiling')
            check(l+k <= bound['counter_ceiling'], 'total counter ceiling')
            check(0 < bound['floor'] <= exact, 'independent deadline gives valid floor')
            if old is not None:
                check(bound['floor'] <= old, 'weaker deadline cannot improve floor')
            old = bound['floor']
        if n in (3, 5, 7, 27, 85):
            plan = measurement_plan(n, tau)
            leaf_word = control_route(plan['leaf'], cap=513)
            check(len(leaf_word) <= plan['odd_step_deadline'], 'section transports deadline')
            check(word_weight(leaf_word) >= plan['floor'], 'actual leaf floor control')
            check(plan['measurement_lower'] >= plan['floor']/2 > 0, 'positive kernel plan')
            check(8*F(16,25)**plan['degree'] == plan['selector_error'], 'declared error')
            if plan['degree']:
                check(8*F(16,25)**(plan['degree']-1) > plan['floor']/2,
                      'unique least sufficient degree')
            examples.append({'source': n, 'premise_odd_steps': tau,
                             'leaf': plan['leaf'], 'total_counter_ceiling': plan['counter_ceiling'],
                             'conditional_floor': str(plan['floor']), 'degree': plan['degree'],
                             'scope': 'known route used only to validate this finite control'})

    # The central binomial denominator is the exact formal minimum at fixed N.
    for c in range(1, 51):
        lower = F(2, (c+2)*comb(c+1, (c+1)//2))
        attained = False
        for n in range(1, c+1):
            for k in range(1, n+1):
                value = F(2, (n+2)*comb(n+1, k))
                check(value >= lower, 'finite counter simplex minimum')
                attained |= value == lower
        check(attained, 'formal bound attained; integer realizability not asserted')

    hostiles = []
    for label, rule in (('bitlength', lambda b:b),
                        ('twice bitlength', lambda b:2*b),
                        ('four bitlength plus one', lambda b:4*b+1)):
        for n in range(3, 1024, 2):
            time = rule(n.bit_length())
            verdict = check_deadline(n, time)
            if not verdict['met']:
                actual = len(control_route(n, cap=512))
                check(actual > time, 'candidate deadline refuted')
                hostiles.append({'candidate': label, 'first_failed_source_in_universe': n,
                                 'proposed_deadline': time, 'actual_odd_steps': actual})
                break
        else:
            raise ValueError('hostile unexpectedly absent')
    for call in (lambda:deadline_floor(True, 2), lambda:deadline_floor(4, 2),
                 lambda:deadline_floor(3, 0), lambda:deadline_floor(3, 1),
                 lambda:deadline_floor(3, 2.0), lambda:degree_for_floor(F(0))):
        try:
            call()
        except ValueError:
            check(True, 'invalid premise or type rejected')
        else:
            raise ValueError('bad input accepted')
    mixed = hybrid_measurement_plan(739)
    check(mixed['status'] == 'proved_on_hybrid_family', 'independent family discharges premise')
    check(mixed['leaf'] == 15765 and mixed['degree'] == 151, 'composed finite measurement plan')
    check(mixed['floor'] == F(1, 9448295337167360434640071920), 'telescoped family floor')
    check(mixed['certified_readout_floor'] <= mixed['measurement_lower'], 'retained exact margin')
    check(hybrid_measurement_plan(7)['status'] == 'outside_proved_family', 'proper domain retained')
    print(json.dumps({'status':'PROVED conditional reverse compiler; all-source premise OPEN',
                      'checks':checks, 'controls':examples, 'missing_guard_hostiles':hostiles,
                      'proved_hybrid_control': {'source':739, 'leaf':mixed['leaf'],
                          'independent_family_deadline':mixed['original_deadline'],
                          'degree':mixed['degree'],
                          'positive_readout_floor':str(mixed['certified_readout_floor']),
                          'scope':'structural family proof; no target ROOT record or moment input'},
                      'unique': 'least degree for one supplied deadline; not unique existence of a deadline'},
                     indent=2))


if __name__ == '__main__':
    main()
