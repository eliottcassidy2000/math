"""Exact first-return cut sets, symbolic final exponents, and retained debt.

Stdlib only. All checks execute normally and under -O. The declared complete
universes are (modulus,budget)=(223,8),(233,10), all cut cells, all positive
compositions through each budget, and the finite controls printed below.
"""
from dataclasses import dataclass
from fractions import Fraction
from math import gcd, lcm
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def natural(n, minimum=0):
    need(type(n) is int and n >= minimum, 'exact integer in domain')


def source(n):
    natural(n, 1)
    need(n % 2 == 1, 'positive odd source')


def modulus_guard(modulus):
    natural(modulus, 3)
    need(modulus % 2 == 1, 'odd return modulus')


def word_guard(word):
    need(type(word) is tuple, 'immutable valuation word')
    for a in word:
        natural(a, 1)


def valuation2(n):
    natural(n, 1)
    return (n & -n).bit_length()-1


def data(word):
    word_guard(word)
    p = q = 1
    b = 0
    for a in word:
        p, q, b = 3*p, q*(1 << a), 3*b+q
    return p, q, b


def direct_carry(word):
    return sum(3**(len(word)-1-i)*2**sum(word[:i])
               for i in range(len(word)))


def step(n):
    source(n)
    a = valuation2(3*n+1)
    return (3*n+1) >> a, a


def replay(n, word):
    source(n)
    word_guard(word)
    path = [n]
    for a in word:
        need(n != 1, 'a first-hit word cannot continue beyond root')
        n, actual = step(n)
        need(actual == a, 'exact valuation guard')
        path.append(n)
    return tuple(path)


def crt(a, m, b, n):
    need(gcd(m, n) == 1, 'coprime CRT')
    return (a+m*((b-a)*pow(m, -1, n) % n)) % (m*n)


def word_residue(word):
    p, q, b = data(word)
    return ((q-b)*pow(p, -1, 2*q)) % (2*q), 2*q


@dataclass(frozen=True)
class Tail:
    prefix: tuple
    minimum_last: int

    def __post_init__(self):
        word_guard(self.prefix)
        natural(self.minimum_last, 1)

    @property
    def coefficients(self):
        p, q, b = data(self.prefix)
        return 3*p, q, 3*b+q

    @property
    def cylinder(self):
        p, q, c = self.coefficients
        period = q*(1 << self.minimum_last)
        return (-c*pow(p, -1, period)) % period, period

    @property
    def relative_mass(self):
        return Fraction(2, self.cylinder[1])

    def matches(self, n):
        source(n)
        r, period = self.cylinder
        return n % period == r


def first_return_tail(tail, modulus):
    modulus_guard(modulus)
    p = q = 1
    b = 0
    for a in tail.prefix:
        p, q, b = 3*p, q*(1 << a), 3*b+q
        if b % modulus == 0:
            return False
    return (3*b+q) % modulus == 0


def uniform_descent_tail(tail, modulus):
    """Sufficient all-height cutoff using n>=modulus, never a necessity claim."""
    need(first_return_tail(tail, modulus), 'first-return tail required')
    p, q, c = tail.coefficients
    k = tail.minimum_last
    while modulus*(q*(1 << k)-p) <= c:
        k += 1
    return Tail(tail.prefix, k)


def apply_return_tail(n, tail, modulus):
    """Original-source guard first; output contains the complete checked word."""
    source(n)
    modulus_guard(modulus)
    need(n % modulus == 0, 'source lies in the marked return section')
    need(first_return_tail(tail, modulus), 'symbolic first-return property')
    need(tail.matches(n), 'original-source tail congruence')
    p, q, c = tail.coefficients
    k = valuation2(p*n+c)-(q.bit_length()-1)
    need(k >= tail.minimum_last, 'actual final exponent meets tail guard')
    word = tail.prefix+(k,)
    path = replay(n, word)
    need(path[-1] == (p*n+c)//(q*(1 << k)), 'affine endpoint')
    need(path[-1] % modulus == 0 and all(x % modulus for x in path[1:-1]),
         'first marked return')
    return dict(source=n, word=word, path=path, endpoint=path[-1],
                descends=path[-1] < n, last=k)


def compile_cut(modulus, budget):
    """Partition the formal cost tree: bounded first returns or first overflow."""
    modulus_guard(modulus)
    natural(budget)
    returns, debts = [], []

    def visit(word, p, q, b):
        if word and b % modulus == 0:
            returns.append(word)
            return
        minimum = budget-(q.bit_length()-1)+1
        debts.append(Tail(word, minimum))
        for a in range(1, minimum):
            visit(word+(a,), 3*p, q*(1 << a), 3*b+q)

    visit((), 1, 1, 0)
    return tuple(returns), tuple(debts)


def observe_to_budget(n, modulus, budget):
    """Retain the first crossing edge as debt; do not discard observed work.

    Budget counts halvings, not calls to this function or odd edges. The
    crossing edge is evaluated exactly once and reported separately.
    """
    source(n)
    modulus_guard(modulus)
    natural(budget)
    need(n % modulus == 0, 'marked source section')
    original, current, used = n, n, 0
    word, path = [], [n]
    while True:
        if current == 1:
            return dict(kind='ROOT', source=original, word=tuple(word),
                        path=tuple(path), used=used, frontier=1, exit=None)
        nxt, a = step(current)
        if used+a > budget:
            tail = Tail(tuple(word), budget-used+1)
            need(tail.matches(original), 'debt retains original cylinder')
            return dict(kind='DEBT', source=original, word=tuple(word),
                        path=tuple(path), used=used, frontier=current,
                        exit=(a, nxt), tail=tail,
                        exit_root=nxt == 1, exit_return=nxt % modulus == 0)
        word.append(a)
        path.append(nxt)
        used += a
        current = nxt
        if current % modulus == 0:
            return dict(kind='RETURN', source=original, word=tuple(word),
                        path=tuple(path), used=used, frontier=current, exit=None)


def compositions(total):
    if total == 0:
        yield ()
    else:
        for first in range(1, total+1):
            for rest in compositions(total-first):
                yield (first,)+rest


def cut_compositions(total):
    for mask in range(1 << (total-1)):
        ends = [0]+[i for i in range(1, total) if mask >> (i-1) & 1]+[total]
        yield tuple(b-a for a, b in zip(ends, ends[1:]))


def order2(modulus):
    need(modulus >= 3 and modulus % 2 == 1, 'odd observation modulus')
    x = 2 % modulus
    d = 1
    while x != 1:
        x = x*2 % modulus
        d += 1
    return d


def least_multiple(r, period, modulus):
    n = modulus*((r*pow(modulus, -1, period)) % period)
    need(n > 0 and n % 2, 'least positive odd multiple exists')
    return n


def audit_cut(modulus, budget):
    returns, debts = compile_cut(modulus, budget)
    period = 1 << (budget+1)
    cells = {}
    for word in returns:
        r, q = word_residue(word)
        for x in range(r, period, q):
            need(x not in cells, 'bounded return cylinders disjoint')
            cells[x] = ('return', word)
    for tail in debts:
        r, q = tail.cylinder
        need(q == period and r not in cells, 'one new overflow cell')
        cells[r] = ('debt', tail)
    need(set(cells) == set(range(1, period, 2)), 'exact complete odd cut set')

    independent_returns = []
    for total in range(1, budget+1):
        recursive = set(compositions(total))
        need(recursive == set(cut_compositions(total)), 'independent word generation')
        for word in recursive:
            p, q, b = data(word)
            need(b == direct_carry(word), 'independent carry sum')
            if b % modulus == 0 and all(data(word[:i])[2] % modulus
                                       for i in range(1, len(word))):
                independent_returns.append(word)
    need(set(returns) == set(independent_returns), 'independent first-return census')
    baseline = sum((Fraction(1, data(w)[1]) for w in returns), Fraction())
    overflow_mass = sum((t.relative_mass for t in debts), Fraction())
    need(baseline+overflow_mass == 1, 'finite cylinder masses sum to one')

    terminal = tuple(t for t in debts if first_return_tail(t, modulus))
    refined = tuple(uniform_descent_tail(t, modulus) for t in terminal)
    added = sum((t.relative_mass for t in refined), Fraction())
    final_period = max([period]+[t.cylinder[1] for t in refined])
    expected = (baseline+added)*final_period/2
    need(expected.denominator == 1, 'integer finite-cell census')
    hit, root, crossing_root, literals = 0, 0, 0, 0
    for j in range(final_period//2):
        n = modulus*(2*j+1)
        observation = observe_to_budget(n, modulus, budget)
        memberships = [t for t in refined if t.matches(n)]
        need(len(memberships) <= 1, 'refined tails remain disjoint')
        if observation['kind'] == 'ROOT':
            root += 1
            need(not memberships, 'root cannot lie in a marked-return tail')
        elif observation['kind'] == 'RETURN':
            need(observation['word'] in returns and not memberships,
                 'bounded direct return agrees with cut')
            need(observation['frontier'] < n, 'inherited bounded return descends')
            hit += 1
        else:
            need(cells[n % period] == ('debt', observation['tail']),
                 'literal first overflow agrees with symbolic cut')
            crossing_root += observation['exit_root']
            if memberships:
                row = apply_return_tail(n, memberships[0], modulus)
                need(row['descends'] and row['word'][:-1] == observation['word']
                     and row['path'][-2:] == (observation['frontier'], observation['exit'][1]),
                     'symbolic debt completion reuses the exact crossing edge')
                hit += 1
        literals += 1
    need(hit == expected, 'complete direct source census equals exact mass')

    controls = 0
    table = []
    for original, tail in zip(terminal, refined):
        p, q, c = tail.coefficients
        r, mod = tail.cylinder
        least = least_multiple(r, mod, modulus)
        least_row = apply_return_tail(least, tail, modulus)
        early_descent = None
        for length in range(1, len(tail.prefix)+1):
            ep, eq, eb = data(tail.prefix[:length])
            if modulus*(eq-ep) > eb:
                early_descent = tail.prefix[:length]
                break
        all_proper_grow = all(data(tail.prefix[:i])[0] > data(tail.prefix[:i])[1]
                              for i in range(1, len(tail.prefix)+1))
        for k in range(tail.minimum_last, tail.minimum_last+6):
            word = tail.prefix+(k,)
            wr, wmod = word_residue(word)
            seed = least_multiple(wr, wmod, modulus)
            for height in (0, 1, 17):
                row = apply_return_tail(seed+height*modulus*wmod, tail, modulus)
                need(row['descends'] and row['last'] == k, 'all declared tail controls')
                if early_descent is not None:
                    need(row['path'][len(early_descent)] < row['source'],
                         'return coverage is not new ordinary-descent coverage')
                if all_proper_grow:
                    need(all(x > row['source'] for x in row['path'][1:-1]),
                         'all-growing prefixes give a genuine first descent at return')
                controls += 1
        # Every discarded low exponent is genuinely expanding in these two cases.
        for k in range(original.minimum_last, tail.minimum_last):
            need(q*(1 << k) < p, 'discarded exponent has uniformly expanding slope')
        table.append(dict(prefix=tail.prefix, first_overflow_last=original.minimum_last,
                          certified_last=tail.minimum_last, P=p, Q_prefix=q, carry=c,
                          residue=r, dyadic_period=mod, mass=str(tail.relative_mass),
                          least_source=least, actual_last=least_row['last'],
                          actual_endpoint=least_row['endpoint'],
                          shorter_descent_prefix=early_descent,
                          all_proper_prefixes_grow=all_proper_grow))
    return dict(modulus=modulus, budget=budget, bounded_words=returns,
                bounded_mass=str(baseline), debt_cells=len(debts),
                terminal_debt_cells=len(terminal), completion_rows=table,
                added_mass=str(added), total_return_descent_mass=str(baseline+added),
                residual_modular_rule_mass=str(1-baseline-added),
                source_period=modulus*final_period, direct_odd_multiples=literals,
                direct_certified_count=hit, in_budget_root_points=root,
                crossing_root_points=crossing_root, independent_tail_controls=controls)


def phase_hostiles():
    v = (1, 3, 1, 1, 1, 3)
    p, q, c = Tail(v, 1).coefficients
    need((p, q, c) == (2187, 1024, 5359), 'incoming expanding carry')
    rows = []
    for a in (1, 2, 3, 4):
        for depth in (1, 2):
            m = 3**a*19**depth
            clock = order2(m)
            need(clock == lcm(2*3**(a-1), 18*19**(depth-1)), 'joined exact clock')
            pair = []
            for k in (1, 1+clock):
                word = v+(k,)
                r, period = word_residue(word)
                mp, oddperiod = crt(0, 233, 1, m), 233*m
                n = crt(r, period, mp, oddperiod)
                row = apply_return_tail(n, Tail(v, k), 233)
                need(n % m == 1 and row['last'] == k, 'same observed source, separate guards')
                pair.append(row)
            need(not pair[0]['descends'] and pair[1]['descends'], 'same finite observer, opposite drift')
            need(pair[0]['endpoint'] % m == pair[1]['endpoint'] % m, 'same observed endpoint')
            need((p*pow(q*2, -1, m)) % m ==
                 (p*pow(q*(1 << (1+clock)), -1, m)) % m,
                 'same entire affine multiplier, not only tested point')
            need((c*pow(q*2, -1, m)) % m ==
                 (c*pow(q*(1 << (1+clock)), -1, m)) % m,
                 'same entire affine carry observation')
            rows.append(dict(ternary_depth=a, observer19_depth=depth, clock=clock,
                             growing_source_bits=pair[0]['source'].bit_length(),
                             descending_source_bits=pair[1]['source'].bit_length(),
                             source_residue=1, endpoint_residue=pair[0]['endpoint'] % m))
    return rows


def unsectioned_families():
    """The two genuine first descents do not require the 223 return filter."""
    rows = []
    for prefix, minimum in (((1, 1, 2, 1, 2, 1), 4),
                            ((1, 1, 1, 1, 2, 1, 1), 5)):
        tail = Tail(prefix, minimum)
        p, q, c = tail.coefficients
        r, period = tail.cylinder
        e = (p*r+c)//period
        need((p*r+c) % period == 0 and r*(period-p) > c,
             'least-source bound proves unconditional descent on the whole cylinder')
        need(all(data(prefix[:i])[0] > data(prefix[:i])[1]
                 for i in range(1, len(prefix)+1)), 'every proper prefix grows')
        for parameter in range(64):
            n = r+period*parameter
            extra = valuation2(e+p*parameter)
            path = replay(n, prefix+(minimum+extra,))
            need(path[-1] == (e+p*parameter) >> extra and path[-1] < n and
                 all(x > n for x in path[1:-1]), 'unfiltered exact first-descent family')
        rows.append(dict(nominal_word=prefix+(minimum,), residue=r, period=period,
                         endpoint_intercept=e, endpoint_slope=p, first_parameter=0,
                         replayed_parameters=64,
                         scope='whole dyadic first-descent family; 223 only selects marked returns'))
    return rows


def main():
    results = [audit_cut(223, 8), audit_cut(233, 10)]
    need(results[0]['added_mass'] == '19/4096' and
         results[0]['total_return_descent_mass'] == '67/4096', '223 exact extension')
    need(results[1]['added_mass'] == '9/2048' and
         results[1]['total_return_descent_mass'] == '25/2048', '233 exact extension')
    need(all(row['shorter_descent_prefix'] is not None and
             len(row['shorter_descent_prefix']) <= 2
             for row in results[1]['completion_rows']),
         'every 233 completion was already a one/two-step ordinary descent')
    need(sum(row['all_proper_prefixes_grow'] for row in results[0]['completion_rows']) == 2,
         'two 223 exit tails really first descend on their returning edge')
    hostile = apply_return_tail(707155, Tail((1, 3, 1, 1, 1, 3), 1), 233)
    need(hostile['endpoint'] == 755153 and not hostile['descends'], 'incoming expanding first return')
    growing = observe_to_budget(707155, 233, 10)
    need(growing['used'] == 10 and growing['frontier'] == 503435 and
         growing['exit'] == (1, 755153), 'retained first beyond-budget edge')
    root = observe_to_budget(3, 3, 10)
    need(root['kind'] == 'ROOT' and root['word'] == (1, 4), 'first root, no artificial (2) suffix')
    crossing = observe_to_budget(5, 5, 0)
    need(crossing['kind'] == 'DEBT' and crossing['exit_root'] and
         crossing['exit'] == (4, 1), 'beyond-budget root edge is retained')
    w1, w2 = (5, 2, 1), (2, 5, 1)
    d1, d2 = data(w1), data(w2)
    need(sorted(w1) == sorted(w2) and d1[:2] == d2[:2] and
         d1[2] == 233 and d2[2] == 149, 'ordered carry versus commutative counts')
    need(first_return_tail(Tail(w1[:-1], 1), 233) and
         not first_return_tail(Tail(w2[:-1], 1), 233), 'order changes modular return')
    invalid = [(True, 233, 10), (0, 233, 10), (2, 233, 10),
               (233, 233, -1), (233, 233, True), (233, 2, 10)]
    for args in invalid:
        try:
            observe_to_budget(*args)
        except ValueError:
            pass
        else:
            raise ValueError('invalid source/budget/modulus accepted')
    print(json.dumps(dict(status='PROVED cut/tail/observer theorems; FINITE-EXACT declared audits',
                          cuts=results, growing_debt=dict(source=707155, frontier=503435,
                          exact_exit=(1, 755153), within_cost=10, actual_return_cost=11),
                          ordered_carry_hostile=(d1, d2),
                          finite_observer_drift_aliases=phase_hostiles(),
                          frontier_family_crosslinks=unsectioned_families(),
                          root_controls=('3->5->1 stops', '5->1 crossing edge retained'),
                          invalid_inputs_rejected=len(invalid)), indent=2))
    print('PASS: source guards and first-hit root checks remain active under -O')


if __name__ == '__main__':
    main()
