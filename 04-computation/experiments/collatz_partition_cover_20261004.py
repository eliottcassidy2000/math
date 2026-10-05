"""Exact blind cells for bounded-depth smaller-child Collatz join portfolios.

No convergence oracle is used. The inverse search has an exact size bound,
not a valuation-exponent cutoff. The all-depth theorem is in the proof note.
"""
from fractions import Fraction
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def exact_nonnegative(n):
    need(type(n) is int and n >= 0, 'exact nonnegative depth')


def odd(n):
    need(type(n) is int and n > 0 and n % 2, 'positive exact odd integer')


def step(n):
    odd(n)
    z = 3*n+1
    a = (z & -z).bit_length()-1
    return z >> a, a


def replay(n, word):
    odd(n)
    states = [n]
    for a in word:
        need(type(a) is int and a >= 1, 'exact positive valuation')
        n, actual = step(n)
        need(actual == a, 'actual word, including any checked root loops')
        states.append(n)
    return states


def word_data(word):
    p = q = 1
    b = 0
    for a in word:
        need(type(a) is int and a >= 1, 'exact positive valuation')
        p, q, b = 3*p, q*2**a, 3*b+q
    return p, q, b


def exact_cylinder(word):
    p, q, b = word_data(word)
    return (q-b)*pow(p, -1, 2*q) % (2*q), 2*q


def blind_cell(R, S):
    """All positive members n>1 miss every smaller join with depths <=R,S."""
    exact_nonnegative(R)
    exact_nonnegative(S)
    dyadic, ternary = 2**(R+1), 3**S
    residue = ternary*((-pow(ternary, -1, dyadic)) % dyadic)
    modulus = dyadic*ternary
    return residue, modulus


def checked_dependency(n, m, left, right):
    odd(n)
    odd(m)
    need(m < n, 'strict original-source comparison')
    a, b = replay(n, left), replay(m, right)
    need(a[-1] == b[-1], 'retained actual common future')
    return a[-1]


def smaller_ancestors(n, endpoint, S, counters):
    """Exhaust every ancestor m<n through depth S, with unrestricted exponents.

    If h inverse steps remain, any eventual m<n requires
    2^h*(current+1)<3^h*(n+1). This bounds each next exponent exactly.
    """
    odd(n)
    odd(endpoint)
    exact_nonnegative(S)
    hits = []

    def visit(x, word):
        counters['inverse_nodes'] += 1
        depth = len(word)
        if x < n:
            need(replay(x, word)[-1] == endpoint, 'independent actual ancestor word')
            hits.append((x, word))
        remaining = S-depth
        if not remaining or x % 3 == 0:
            return
        if 2**remaining*(x+1) >= 3**remaining*(n+1):
            return
        h = remaining-1
        a = 1
        while 2**h*((x << a)+2) < 3**(h+1)*(n+1):
            counters['inverse_exponents'] += 1
            numerator = (x << a)-1
            if numerator % 3 == 0:
                y = numerator//3
                need(y > 0 and y % 2, 'positive exact inverse branch')
                visit(y, (a,)+word)
            a += 1

    visit(endpoint, ())
    return tuple(hits)


def bounded_joins(n, R, S, counters):
    odd(n)
    exact_nonnegative(R)
    exact_nonnegative(S)
    endpoint = n
    left = ()
    hits = []
    for r in range(R+1):
        counters['forward_endpoints'] += 1
        for m, right in smaller_ancestors(n, endpoint, S, counters):
            need(checked_dependency(n, m, left, right) == endpoint, 'checked point join')
            hits.append((r, len(right), m, left, right))
        if r < R:
            endpoint, a = step(endpoint)
            left += (a,)
    return tuple(hits)


def compositions(total):
    if total == 0:
        yield ()
        return
    for a in range(1, total+1):
        for tail in compositions(total-a):
            yield (a,)+tail


def escape_family(R, S):
    """A reset>=3 subfamily inside the blind cell, with source depth R+1."""
    need(type(R) is int and R >= 1, 'positive forward bound for this repair')
    exact_nonnegative(S)
    left, right = (1,)*R+(3,), (1,)*(R-1)+(2, 1)
    residue, dyadic = exact_cylinder(left)
    ternary = 3**S
    base = residue+dyadic*((-residue)*pow(dyadic, -1, ternary) % ternary)
    return base, dyadic*ternary, left, right


def main():
    counters = dict(forward_endpoints=0, inverse_nodes=0, inverse_exponents=0)
    examples = []
    cell_count = sources = 0
    for R in range(9):
        for S in range(9):
            residue, modulus = blind_cell(R, S)
            cell_count += 1
            need(residue % 2**(R+1) == 2**(R+1)-1 and residue % 3**S == 0,
                 'independent CRT residue checks')
            for t in (0, 1, 7):
                n = residue+modulus*t
                if n == 1:
                    n += modulus
                sources += 1
                need(replay(n, (1,)*R)[-1] == (3**R*(n+1)//2**R)-1,
                     'actual all-ones prefix is fixed by dyadic guard')
                need(not bounded_joins(n, R, S, counters),
                     'no smaller common-future source in the full bounded-depth search')
            if R == S:
                examples.append(dict(R=R, S=S, residue=residue, modulus=modulus,
                                     odd_relative_density=str(Fraction(1, 2**R*3**S))))

    # Independent finite check: enumerate all smaller odd labels, rather than
    # inverse words, for small sources and compare the exact sets of diagrams.
    brute_controls = 0
    for n in range(3, 152, 2):
        for R, S in ((0, 3), (1, 2), (2, 3), (3, 4)):
            fast = {(r, s, m) for r, s, m, _, _ in bounded_joins(n, R, S, counters)}
            forward = [n]
            for _ in range(R):
                forward.append(step(forward[-1])[0])
            slow = set()
            for m in range(1, n, 2):
                z = m
                for s in range(S+1):
                    for r, endpoint in enumerate(forward):
                        if z == endpoint:
                            slow.add((r, s, m))
                    z = step(z)[0]
            need(fast == slow, 'inverse enumeration agrees with all-smaller-label replay')
            brute_controls += 1

    # All words of cost<=16; test the formal-source-zero sign barrier using
    # a separate closed carry sum and explicit rational inverse replay.
    words = inverse_sign_controls = contracting = 0
    for cost in range(1, 17):
        for v in compositions(cost):
            words += 1
            s = len(v)
            p, q, b = word_data(v)
            direct_b = sum(3**(s-j-1)*2**sum(v[:j]) for j in range(s))
            need(b == direct_b, 'independent carry formula')
            for r in range(s):
                lam = Fraction(3**r*q, 2**r*p)
                intercept = lam-Fraction(q+b, p)
                if lam < 1:
                    contracting += 1
                    need(intercept < 0 and intercept.denominator != 1,
                         'contracting map has no integral formal-zero intercept')
                if intercept.denominator == 1:
                    need(intercept > 0, 'longer integral inverse of z_r is positive')
                if cost <= 8:
                    x = Fraction(3**r, 2**r)-1
                    for a in reversed(v):
                        x = (2**a*x-1)/3
                    need(x == intercept, 'separate rational inverse chain')
                    inverse_sign_controls += 1
    need(words == 65535, 'complete positive-composition universe')

    positive_controls = 0
    for R in range(1, 9):
        for S in range(9):
            base, period, left, right = escape_family(R, S)
            blind_residue, blind_modulus = blind_cell(R, S)
            for t in (0, 1, 17, 10**30):
                n = base+period*t
                need((n-blind_residue) % blind_modulus == 0, 'repair stays inside blind class')
                need(checked_dependency(n, (n-1)//2, left, right) > 1,
                     'one greater forward depth supplies a guarded smaller child')
                need(len(left) == R+1, 'repair honestly exceeds the excluded forward bound')
                positive_controls += 1

    # Omitting the ternary condition would be false: even arbitrarily long
    # all-one prefixes can have a smaller one-edge ancestor on another row.
    dyadic_hostiles = []
    for R in range(9):
        n, m = 3*2**(R+1)-1, 2**(R+2)-1
        need(n % 2**(R+1) == 2**(R+1)-1 and n % 3 == 2,
             'same deep dyadic class, different ternary row')
        need(checked_dependency(n, m, (), (1,)) == n,
             'mixed-guard inverse rule refutes a dyadic-only point obstruction')
        dyadic_hostiles.append((R, n, m))
    need(replay(581, (4, 3))[-1] == replay(27, (1,))[-1] == 41,
         'a larger actual ancestor exists; only smaller ones are excluded')
    need(Fraction(64, 3)*27+5 == 581, 'positive integral formal-zero intercept control')
    mersenne_controls = []
    for p in range(3, 62, 2):
        n = 2**p-1
        if p % 6 == 5:
            need(n % 9 == 4, 'inherited two-edge inverse guard')
            m = (8*n-5)//9
            need(checked_dependency(n, m, (), (1, 2)) == n,
                 'old predecessor rule handles arbitrarily long binary-one prefixes')
            mersenne_controls.append(p)
        else:
            need(n % 9 in (1, 7), 'remaining odd exponent classes are outside that guard')
    need(replay(29, (3,))[-1] == 11 and 29 % 8 == 5 and 11 % 8 != 5,
         'a descending family need not contain its own smaller obligation')
    try:
        checked_dependency(27, 31, (1, 2), ())
    except ValueError:
        pass
    else:
        raise ValueError('descent from a later frontier confused with payment of original source')

    malformed = 0
    for args in ((True, 2), (-1, 2), (1, 2.0), (1, -1)):
        try:
            blind_cell(*args)
        except ValueError:
            malformed += 1
        else:
            raise ValueError('malformed depth accepted')
    result = dict(
        status='PROVED all-depth blind-cell theorem; FINITE-EXACT independent controls',
        theorem='n=-1 mod2^(R+1), n=0 mod3^S, n>1 => no m<n with U^r(n)=U^s(m), r<=R,s<=S',
        universe='R,S=0..8; each cell at parameters0,1,7 (replace root1 by next member)',
        cells=cell_count, blind_sources_checked=sources, diagonal_cells=examples,
        independent_all_smaller_label_controls=brute_controls, search_counters=counters,
        child_words_cost_at_most16=words, contracting_word_pairs=contracting,
        rational_inverse_controls=inverse_sign_controls,
        depth_growing_repair_controls=positive_controls,
        dyadic_only_hostiles=dyadic_hostiles,
        larger_ancestor_control=dict(source=27, child=581, endpoint=41),
        family_not_closed_control=dict(source=29, smaller=11, source_family='5mod8'),
        inherited_mod9_mersenne_exponents=mersenne_controls,
        unpaid_original_source_rejected=True,
        malformed_depths_rejected=malformed,
        scope='No convergence conclusion. No exponent cutoff in the bounded inverse search.')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
