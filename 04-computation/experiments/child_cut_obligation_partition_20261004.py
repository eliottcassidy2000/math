"""Exact parameter obligations for fixed-word child-cut repair.

The algebraic APIs allow U(1)=1; strict first-hit callers must retain the
root sidecar. The finite route audit uses only supplied first-hit routes.
No API decides the universal child-cut conjecture or searches for a proof
of an arbitrary source.
"""
from dataclasses import dataclass
from fractions import Fraction
from collections import Counter

import supplied_state_join_extensions_20261004 as inherited


def need(ok, message):
    if not ok:
        raise ValueError(message)


def replay(n, word, strict=False):
    need(type(n) is int and n > 0 and n % 2, 'positive odd integer')
    values = [n]
    for a in word:
        need(type(a) is int and a > 0, 'positive exact valuation')
        need(not strict or n != 1, 'strict first-hit root sidecar')
        z = 3*n+1
        need(z % 2**a == 0 and (z//2**a) % 2, 'exact word guard')
        n = z//2**a
        values.append(n)
    return tuple(values)


@dataclass(frozen=True)
class PairCell:
    source_word: tuple
    child_word: tuple
    endpoint: int
    endpoint_period: int
    source: int
    child: int
    source_period: int
    child_period: int

    def at(self, t):
        need(type(t) is int and t >= 0, 'nonnegative exact parameter')
        return (self.source+self.source_period*t,
                self.child+self.child_period*t,
                self.endpoint+self.endpoint_period*t)


def normalize_pair(source_word, child_word):
    """Return all positive compatible realizations, or None for disjoint guards."""
    w, v = tuple(source_word), tuple(child_word)
    p, q, b = inherited.data(w)
    pp, qq, bb = inherited.data(v)
    left, right = (b*pow(q, -1, p)) % p, (bb*pow(qq, -1, pp)) % pp
    if (left-right) % min(p, pp):
        return None
    H = max(p, pp)
    J = left if p >= pp else right
    if J % 2 == 0:
        J += H
    n, m = (q*J-b)//p, (qq*J-bb)//pp
    P, Q = 2*q*(H//p), 2*qq*(H//pp)
    need(0 < n < P and 0 < m < Q, 'canonical positive seed bounds')
    need(replay(n, w)[-1] == J == replay(m, v)[-1], 'canonical exact realization')
    return PairCell(w, v, J, 2*H, n, m, P, Q)


def partition(cell):
    """Finite smaller-child interval and upper-ray repair for an expanding pair.

    Each cut row is (cut index, seed value, period, first repair parameter).
    None for the last coordinate means that this cut never repairs the slope
    and strict source-size comparison together. The endpoints of intervals
    can be negative: [0,-1] denotes empty.
    """
    P, Q = cell.source_period, cell.child_period
    need(Q > P, 'expanding fixed pair')
    maximum = (cell.source-cell.child-1)//(Q-P)
    states = replay(cell.child, cell.child_word)
    period = Q
    rows = []
    for i, x in enumerate(states):
        lower = None
        if period < P:
            lower = max(0, (x-cell.source)//(P-period)+1)
        elif period == P and x < cell.source:
            lower = 0
        rows.append((i, x, period, lower))
        if i < len(cell.child_word):
            a = cell.child_word[i]
            need((3*period) % 2**a == 0, 'integral cut period')
            period = 3*period//2**a
    bounds = [row[3] for row in rows if row[3] is not None]
    first = min(bounds) if bounds else None
    last_obligation = maximum if first is None else min(maximum, first-1)
    count = max(0, maximum+1)
    need(count <= (P+(Q-P)-1)//(Q-P), 'finite initial seed count bound')
    return dict(maximum=maximum, first_repair=first,
                last_obligation=last_obligation, cuts=tuple(rows))


def residue_product(states):
    product = Fraction(1)
    for x in states[:-1]:
        product *= 1+Fraction(1, 3*x)
    return product


def compositions(total):
    if total == 0:
        yield ()
        return
    for first in range(1, total+1):
        for suffix in compositions(total-first):
            yield (first,)+suffix


def rational_replay(n, word):
    values = [Fraction(n)]
    for a in word:
        n = (3*values[-1]+1)/2**a
        need(n > 0 and n.numerator % 2 and n.denominator % 2,
             'positive odd rational exact valuation')
        values.append(n)
    return tuple(values)


def main():
    print('CHILD-CUT OBLIGATION PARTITION: exact fixed-word intervals; universal integer repair OPEN')
    # Strong rational hostile: both routes have finite rational first-hit tails.
    hostile = rational_replay(Fraction(77, 27), (1, 1, 3))
    need(hostile == (Fraction(77, 27), Fraction(43, 9), Fraction(23, 3), Fraction(3)),
         'positive rational hostile')
    need(all(x > 1 for x in hostile) and hostile[0] < 3 < min(hostile[1:-1]),
         'all-above-one hostile and failed smaller cuts')
    need(rational_replay(3, (1, 4))[-1] == 1, 'rational hostile has a completed common tail')
    cell = normalize_pair((), (1, 1, 3))
    need((cell.endpoint, cell.source, cell.child) == (47, 47, 55), 'integer guard misses bad interval')
    need(Fraction(23, 16) < 3 < Fraction(19, 5) < 47, 'exact rational interval versus integer lattice')
    scaled = [77]
    for a in (1, 1, 3):
        z = 3*scaled[-1]+27
        need(z % 2**a == 0 and (z//2**a) % 2, 'scaled 3x+27 actual valuation')
        scaled.append(z//2**a)
    need(scaled == [77, 129, 207, 81], 'integer different-map hostile')
    small = rational_replay(Fraction(11, 81), (1, 2))
    need(small == (Fraction(11, 81), Fraction(19, 27), Fraction(7, 9)), 'small rational control')
    need(rational_replay(Fraction(5, 27), (1,))[-1] == small[-1], 'rational common endpoint')
    print('Rational hostile: 77/27 ->43/9 ->23/3 ->3; source3; slope32/27; all values>1')
    print('Repair obstruction interval23/16<J<19/5; integral odd endpoint J=47+54t excludes it')
    print('Scaled different-map control3x+27: 77->129->207->81; no conclusion for3x+1 integers')

    # Independent backwards-rational endpoint enumeration checks the full
    # compatibility iff, including incompatibility and root-padding boundaries.
    small_words = [()] + [w for total in range(1, 5) for w in compositions(total)]
    endpoint_checks = 0
    for w in small_words:
        for v in small_words:
            cell = normalize_pair(w, v)
            for J in range(1, 82, 2):
                starts = []
                for word in (w, v):
                    x = Fraction(J)
                    for a in reversed(word):
                        x = (2**a*x-1)/3
                    starts.append(x)
                valid = all(x.denominator == 1 and x > 0 and x.numerator % 2 for x in starts)
                represented = cell is not None and J >= cell.endpoint and (J-cell.endpoint) % cell.endpoint_period == 0
                need(valid == represented, 'independent exact endpoint compatibility iff')
                endpoint_checks += 1
    print('Independent backwards endpoint checks:', endpoint_checks)

    # Direct word replay in a complete modest symbolic universe.
    source_words = [()] + [w for total in range(1, 7) for w in compositions(total) if len(w) <= 2]
    child_words = [()] + [w for total in range(1, 13) for w in compositions(total)]
    counts = Counter()
    largest_seed_count = 0
    for w in source_words:
        p, q, b = inherited.data(w)
        for v in child_words:
            pp, qq, bb = inherited.data(v)
            if p*qq <= pp*q:
                continue
            counts['expanding'] += 1
            cell = normalize_pair(w, v)
            if cell is None:
                continue
            counts['compatible'] += 1
            out = partition(cell)
            for t in (0, 1):
                n, m, J = cell.at(t)
                ns, ms = replay(n, w), replay(m, v)
                need(ns[-1] == ms[-1] == J, 'independent actual lift replay')
                predicted = out['first_repair'] is not None and t >= out['first_repair']
                actual = any(x < n and period <= cell.source_period
                             for (_, _, period, _), x in zip(out['cuts'], ms))
                need(predicted == actual, 'cut threshold iff')
                need((m < n) == (t <= out['maximum']), 'smaller-child interval iff')
            if out['maximum'] >= 0:
                counts['smaller_seed_shapes'] += 1
                largest_seed_count = max(largest_seed_count, out['maximum']+1)
                need(out['last_obligation'] < 0, 'finite symbolic probe has no remaining seed obligations')
    need((len(source_words), len(child_words)) == (22, 4096), 'declared symbolic universe')
    need(dict(counts) == {'expanding':44600, 'compatible':15116, 'smaller_seed_shapes':11},
         'finite symbolic counts')
    print('Symbolic shapes: source length<=2/cost<=6 plus empty22; child cost<=12 plus empty4096')
    print('Expanding44600; compatible15116; smaller-seed shapes11; unfilled obligations0; maxseedcount', largest_seed_count)

    # Reuse the exact previous universe; this is a changed obligation analysis,
    # not an enlarged source convergence census.
    routes = {n: inherited.explicit_route(n) for n in range(1, 1002, 2)}
    indexed = {n: {x:i for i,x in enumerate(states)} for n,(states,_) in routes.items()}
    counts = Counter()
    thresholds = Counter()
    product_checks = 0
    examples = []
    for n in range(3, 1002, 2):
        ns, nw = routes[n]
        for m in range(1, n, 2):
            ms, mw = routes[m]
            r = next(i for i,x in enumerate(ns) if x in indexed[m])
            s = indexed[m][ns[r]]
            w, v = nw[:r], mw[:s]
            p, q, _ = inherited.data(w)
            pp, qq, _ = inherited.data(v)
            if p*qq <= pp*q:
                continue
            counts['expanding'] += 1
            if ns[r] < n:
                counts['endpoint_below'] += 1
                continue
            counts['nontrivial'] += 1
            cell = normalize_pair(w, v)
            need(cell is not None, 'observed compatible pair')
            t, remainder = divmod(n-cell.source, cell.source_period)
            need(remainder == 0 and t >= 0 and cell.at(t) == (n, m, ns[r]), 'original identity in normalized family')
            out = partition(cell)
            need(out['last_obligation'] < 0, 'all positive smaller-child seeds of this pair shape are repaired')
            thresholds[(out['maximum'], out['first_repair'])] += 1
            left_product = residue_product(ns[:r+1])
            for cut in range(s+1):
                x = ms[cut]
                suffix_product = residue_product(ms[cut:s+1])
                ps, qs, _ = inherited.data(v[cut:])
                lam = Fraction(p*qs, q*ps)
                need(lam == Fraction(x, n)*suffix_product/left_product, 'exact two-path residue-product criterion')
                product_checks += 1
            if (n,m) == (233,231):
                need((out['maximum'], out['first_repair']) == (0,0), 'canonical inherited repair parameter')
                examples = [(i,x,str(Fraction(period,cell.source_period)),lower)
                            for i,x,period,lower in out['cuts'] if lower is not None]
    need(dict(counts) == {'expanding':2914, 'endpoint_below':2822, 'nontrivial':92}, 'inherited universe counts')
    print('Inherited source universe odds<=1001: expanding2914; directendpoint2822; nontrivial92')
    print('All92 entire fixed-word smaller-seed intervals discharged; (lastseed,firstrepair) frequencies:', sorted(thresholds.items()))
    print('Exact two-path residue-product checks:', product_checks)
    print('Pair233/231 repairable cuts (index,seed,slope,firstparameter):', examples)

    # Export type/guard controls survive Python -O.
    for thunk in (lambda: normalize_pair((True,),()), lambda: replay(3,(1.0,)),
                  lambda: cell.at(True), lambda: replay(1,(2,),strict=True)):
        try:
            thunk()
        except ValueError:
            pass
        else:
            raise ValueError('hostile type or first-hit guard accepted')
    need(normalize_pair((1,), (2,)) is None, 'disjoint ternary endpoint guards')
    need(replay(1,(2,2)) == (1,1,1), 'formal root padding is explicitly a separate domain')
    print('Status: no integer counterexample found in inherited universe; universal cut existence remains OPEN')


if __name__ == '__main__':
    main()
