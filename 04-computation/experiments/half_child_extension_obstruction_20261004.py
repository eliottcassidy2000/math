"""Exact obstruction to indefinitely extending the affine half-child decoder.

This is a finite certificate audit, not a termination algorithm. All 5000
sources n=3 mod4 below20000 and their children are explicitly replayed home.
Checks remain active under python -O. The all-length conclusions have proofs
in the matching note; no finite padding cutoff is used as that proof.
"""
from collections import Counter
from fractions import Fraction
import json


def need(condition, message):
    if not condition:
        raise ValueError(message)


def positive_odd(n):
    need(type(n) is int and n > 0 and n % 2 == 1, 'exact positive odd source')


def step(n):
    positive_odd(n)
    t = 3*n+1
    a = (t & -t).bit_length()-1
    return t >> a, a


def word_data(word):
    p = q = 1
    b = 0
    for a in word:
        need(type(a) is int and a >= 1, 'exact positive valuation')
        p, q, b = 3*p, (1 << a)*q, 3*b+q
    return p, q, b


def replay(n, word, allow_root_padding=False):
    positive_odd(n)
    path = [n]
    for a in word:
        need(type(a) is int and a >= 1, 'exact positive valuation')
        need(allow_root_padding or n != 1, 'first root is terminal')
        nxt, actual = step(n)
        need(actual == a, 'supplied word has the exact valuations')
        n = nxt
        path.append(n)
    return path


def home_word(n, cap=5000):
    """Bounded control only: obtain an explicit first-hit certificate."""
    positive_odd(n)
    out = []
    while n != 1:
        need(len(out) < cap, 'finite control failed to reach root within cap')
        n, a = step(n)
        out.append(a)
    return tuple(out)


def validate_home(n, word):
    need(replay(n, word)[-1] == 1, 'supplied first-hit home certificate')
    return len(word), sum(word)


def charge(word):
    return sum(word)-2*len(word)


def half_child_identity(left, right):
    """Equality as affine functions of X, not equality only at a seed."""
    p, q, b = word_data(left)
    pp, qq, bb = word_data(right)
    return 2*p*qq == pp*q and 2*b*qq == q*(2*bb-pp)


def classify_certified_source(n, source_word, child_word):
    """Consume existing certificates; never claim an unobserved source homes."""
    positive_odd(n)
    need(n >= 3 and n % 4 == 3, 'half-child is positive odd')
    m = (n-1)//2
    rn, an = validate_home(n, source_word)
    rm, am = validate_home(m, child_word)
    delta = (an-2*rn)-(am-2*rm)
    return dict(source=n, child=m, source_rank=rn, child_rank=rm,
                source_halvings=an, child_halvings=am, charge_difference=delta,
                padded_possible=(delta == 1),
                strict_possible=(rn == rm and an-am == 1))


def decide_strict_with_supplied_child(n, child_word):
    """Decide this exact ansatz in at most the supplied child's odd rank.

    A rejection says nothing against other children or unequal word lengths.
    The actual source/frontier/prefix remain explicit for a caller to retain.
    """
    positive_odd(n)
    need(n >= 3 and n % 4 == 3, 'half-child is positive odd')
    rank, cost = validate_home((n-1)//2, child_word)
    frontier = n
    prefix = []
    while len(prefix) < rank and frontier != 1:
        frontier, a = step(frontier)
        prefix.append(a)
    accepted = frontier == 1 and len(prefix) == rank and sum(prefix) == cost+1
    if accepted:
        need(half_child_identity(prefix, child_word), 'accepted complete affine join')
    return dict(source=n, child=(n-1)//2, child_rank=rank,
                observed_prefix=tuple(prefix), frontier=frontier,
                strict_possible=accepted,
                scope='only the exact affine half-child ansatz without root padding')


def lift_join(n, m, left, right):
    """Inherited least simultaneous exact-cylinder periods for a checked join."""
    positive_odd(n)
    positive_odd(m)
    need(0 < m < n, 'strictly smaller actual child')
    endpoint = replay(n, left)[-1]
    need(replay(m, right)[-1] == endpoint, 'actual common endpoint')
    r, s = len(left), len(right)
    a, d = sum(left), sum(right)
    period_n = (1 << (a+1))*3**max(s-r, 0)
    period_m = (1 << (d+1))*3**max(r-s, 0)
    return period_n, period_m, endpoint


def independent_route(n):
    """Independent ordinary-step reader, avoiding bit valuations."""
    path, word = [n], []
    while n != 1:
        need(len(word) < 5000, 'independent finite cap')
        n = 3*n+1
        a = 0
        while n % 2 == 0:
            n //= 2
            a += 1
        path.append(n)
        word.append(a)
    return path, tuple(word)


def main():
    counts = Counter()
    differences = Counter()
    independent_prefixes = slope_controls = repair_controls = bounded_source_steps = 0
    least_blocked = least_padding_only = None
    examples = {}
    for n in range(3, 20000, 4):
        m = (n-1)//2
        wn, wm = home_word(n), home_word(m)
        pn, wn2 = independent_route(n)
        pm, wm2 = independent_route(m)
        need(wn == wn2 and wm == wm2, 'independent ordinary route reader')
        rec = classify_certified_source(n, wn, wm)
        bounded = decide_strict_with_supplied_child(n, wm)
        need(bounded['strict_possible'] == rec['strict_possible'],
             'supplied child rank bounds a complete strict-ansatz decision')
        need(len(bounded['observed_prefix']) <= len(wm), 'promised finite source bound')
        bounded_source_steps += len(bounded['observed_prefix'])
        delta = rec['charge_difference']
        differences[delta] += 1
        counts['sources'] += 1
        counts['padded_possible' if rec['padded_possible'] else 'all_length_blocked'] += 1
        counts['strict_possible' if rec['strict_possible'] else 'strict_blocked'] += 1
        if not rec['padded_possible'] and least_blocked is None:
            least_blocked = n
        if rec['padded_possible'] and not rec['strict_possible']:
            counts['padding_only'] += 1
            if least_padding_only is None:
                least_padding_only = n

        # Read all possible pre-root synchronous joins literally. Beyond both
        # roots the cost difference is constant, so there is no cutoff guess.
        length = max(len(wn), len(wm))
        pn += [1]*(length-len(wn))
        pm += [1]*(length-len(wm))
        padded_n = wn+(2,)*(length-len(wn))
        padded_m = wm+(2,)*(length-len(wm))
        an = am = 0
        saw_strict = False
        saw_padded = False
        for j in range(length+1):
            same = pn[j] == pm[j]
            permitted = same and an-am == 1
            saw_padded |= permitted
            saw_strict |= permitted and j <= min(len(wn), len(wm))
            need(half_child_identity(padded_n[:j], padded_m[:j]) == permitted,
                 'pointwise join plus exact slope iff affine identity')
            independent_prefixes += 1
            if j < length:
                an += padded_n[j]
                am += padded_m[j]
        need(saw_padded == rec['padded_possible'] and saw_strict == rec['strict_possible'],
             'literal join test agrees with rank-charge classification')
        for extra in (0, 1, 2, 7):
            left = padded_n+(2,)*extra
            right = padded_m+(2,)*extra
            need(charge(left)-charge(right) == delta, 'root padding retains charge')
            need(half_child_identity(left, right) == rec['padded_possible'],
                 'padded affine coefficients independently agree')
        if n in (3, 7, 27, 315, 391, 703, 13483):
            examples[str(n)] = rec

        # General slope orbit checked at several different odd-time offsets.
        if n < 404:
            for k in range(-4, 5):
                rn, rm = len(wn), len(wm)
                extra_n = max(0, k-rn+rm)
                extra_m = rn+extra_n-rm-k
                left = wn+(2,)*extra_n
                right = wm+(2,)*extra_m
                p, q, b = word_data(left)
                pp, qq, bb = word_data(right)
                lam = Fraction(p*qq, q*pp)
                need(lam == Fraction(2)**(-delta)*Fraction(3, 4)**k,
                     'complete root-padded slope orbit')
                beta = m-lam*n
                need(Fraction(p, q) == Fraction(pp, qq)*lam and
                     Fraction(b, q) == Fraction(pp, qq)*beta+Fraction(bb, qq),
                     'general child slope and seed fix the full affine map')
                slope_controls += 1

    need(least_blocked == 7 and least_padding_only == 3, 'least typed hostiles')
    need(examples['7']['charge_difference'] == 0, '7/3 invariant obstruction')
    seven_cut = decide_strict_with_supplied_child(7, (1, 4))
    need(seven_cut['observed_prefix'] == (1, 1) and seven_cut['frontier'] == 17
         and not seven_cut['strict_possible'], '7 strict ansatz blocked after two source steps')
    # Explicit all-height repaired source and child families.
    period_n, period_m, endpoint = lift_join(7, 3, (1, 1, 2, 3), (1,))
    need((period_n, period_m, endpoint) == (256, 108, 5), 'unequal-rank repair periods')
    for t in list(range(128))+[10**20, 10**100]:
        n, m = 7+256*t, 3+108*t
        need(0 < m < n and Fraction(27*n+3, 64) == m, 'repaired decreasing child map')
        need(replay(n, (1, 1, 2, 3))[-1] == replay(m, (1,))[-1] == 5+162*t,
             'repaired family exact guards and common endpoint')
        # A different intercept also gives a valid affine identity. At t=0
        # its child root padding must be trimmed from a certificate.
        n, m = 7+4096*t, 1+2048*t
        need(replay(n, (1, 1, 2, 3, 4))[-1] ==
             replay(m, (2,)*5, allow_root_padding=True)[-1] == 1+486*t,
             'changed-intercept repair retains source and integer guard')
        need(half_child_identity((1, 1, 2, 3, 4), (2,)*5) is False,
             'changed intercept is not the forbidden half-child identity')
        repair_controls += 2

    malformed = 0
    for n in (True, 1.0, 0, -1, 2, 5):
        try:
            classify_certified_source(n, (), ())
        except ValueError:
            malformed += 1
        else:
            raise ValueError('malformed source accepted')
    try:
        classify_certified_source(7, (1, 4), (1, 4))
    except ValueError:
        pass
    else:
        raise ValueError('wrong-source certificate accepted')
    # An alleged home certificate may not quietly include the root loop.
    try:
        classify_certified_source(3, (1, 4, 2), ())
    except ValueError:
        pass
    else:
        raise ValueError('root-padded first-hit certificate accepted')

    result = dict(
        status='FINITE-EXACT controls; all-length statements proved in matching note',
        universe='all5000 n=3 mod4 in3..19999; supplied home words checked with cap5000',
        counts=dict(sorted(counts.items())), least_all_length_obstruction=least_blocked,
        least_padding_only=least_padding_only, charge_difference_counts=dict(sorted(differences.items())),
        literal_synchronous_prefix_controls=independent_prefixes,
        supplied_child_bounded_decisions=counts['sources'], bounded_source_steps=bounded_source_steps,
        seven_bounded_decision=seven_cut,
        general_slope_orbit_controls=slope_controls, all_height_repair_controls=repair_controls,
        malformed_sources_rejected=malformed, wrong_source_and_root_padding_rejected=True,
        examples=examples,
        unequal_length_repair=dict(source='7+256t', child='3+108t', endpoint='5+162t',
                                   source_word=[1, 1, 2, 3], child_word=[1], t='all integers>=0'),
        limitation='This classifies supplied home certificates, not arbitrary unfinished states.')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
