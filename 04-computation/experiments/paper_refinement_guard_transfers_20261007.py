"""Exact coarse-tail sibling guards and trusted-paper transfer controls.

No ROOT lookup in the production compiler. Run normally and with -O.
The finite posterior controls illustrate the supplied paper's mechanism;
they do not attempt to verify its reconstruction or Ramsey headline.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product, permutations


CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def odd(n):
    need(type(n) is int and n > 0 and n % 2 == 1, 'exact positive odd integer')


def natural(n):
    need(type(n) is int and n >= 0, 'exact nonnegative integer')


def word_guard(word):
    need(type(word) is tuple and all(type(a) is int and a > 0 for a in word),
         'exact positive valuation tuple')


def valuation(n):
    need(type(n) is int and n > 0, 'positive valuation argument')
    return (n & -n).bit_length()-1


def step(n):
    odd(n)
    a = valuation(3*n+1)
    return (3*n+1) >> a, a


def carrier(word):
    word_guard(word)
    p = q = 1
    b = 0
    for a in word:
        p, q, b = 3*p, q << a, 3*b+q
    return p, q, b


def replay(n, word):
    odd(n)
    word_guard(word)
    path = [n]
    for a in word:
        need(n != 1, 'first-hit ROOT is not padded')
        n, actual = step(n)
        need(actual == a, 'exact actual valuation')
        path.append(n)
    return tuple(path)


@dataclass(frozen=True)
class TailRow:
    word: tuple
    depth: int
    p: int
    q: int
    b: int
    residue: int
    modulus: int
    least: int
    child0: int


def compile_tail(word, depth):
    """All-height paid S^depth ancestor, retaining arbitrary final valuation."""
    word_guard(word)
    need(type(depth) is int and depth >= 1, 'positive sibling depth')
    p, q, b = carrier(word)
    power = 4**depth
    c = (power-1)//3
    delta = q*power-p
    need(delta > 0, 'contracting child coefficient')
    modulus = 2*q*power
    residue = (q*(c+power)-b)*pow(p, -1, modulus) % modulus
    # All inequalities are strict. The floor formulation includes negatives.
    lower = max(1, (q*c-b)//p, (b-q*c)//delta)
    least = residue+modulus*max(0, (lower-residue)//modulus+1)
    child0 = (p*least+b-q*c)//(q*power)
    need(0 < child0 < least and child0 % 2, 'nonempty paid odd row')
    return TailRow(word, depth, p, q, b, residue, modulus, least, child0)


def canonical(row):
    need(type(row) is TailRow, 'typed symbolic row')
    need(all(type(getattr(row, key)) is int for key in
             ('depth', 'p', 'q', 'b', 'residue', 'modulus', 'least', 'child0')),
         'exact row fields')
    need(row == compile_tail(row.word, row.depth), 'authenticated canonical row')


@dataclass(frozen=True)
class Receipt:
    source: int
    child: int
    source_word: tuple
    child_word: tuple
    endpoint: int


def audit(receipt):
    need(type(receipt) is Receipt, 'typed receipt')
    odd(receipt.source)
    odd(receipt.child)
    odd(receipt.endpoint)
    need(receipt.child < receipt.source, 'original-source payment')
    for n, word in ((receipt.source, receipt.source_word),
                    (receipt.child, receipt.child_word)):
        p, q, b = carrier(word)
        need(p*n+b == q*receipt.endpoint, 'independent affine endpoint check')
        need(replay(n, word)[-1] == receipt.endpoint, 'actual common endpoint')
    return receipt


def apply_tail(row, n):
    canonical(row)
    odd(n)
    if n < row.least or (n-row.residue) % row.modulus:
        return None
    t = (n-row.least)//row.modulus
    h = row.child0+2*row.p*t
    target, residual = step(h)
    right = () if h == 1 else (residual,)
    return audit(Receipt(n, h, row.word+(2*row.depth+residual,), right, target))


def discharge(receipt, supplied_child_word):
    audit(receipt)
    need(replay(receipt.child, supplied_child_word)[-1] == 1,
         'supplied child certificate actually reaches ROOT')
    d = len(receipt.child_word)
    need(supplied_child_word[:d] == receipt.child_word, 'matching child prefix')
    result = receipt.source_word+supplied_child_word[d:]
    need(replay(receipt.source, result)[-1] == 1, 'compiled first-hit ROOT')
    return result


def reserve_phase(row, extra):
    """Parameter t class with at least extra further odd S-ancestors."""
    canonical(row)
    natural(extra)
    modulus = 4**extra
    z = (3*row.child0+1)//2
    residue = -z*pow(3*row.p, -1, modulus) % modulus if modulus > 1 else 0
    return residue, modulus


def spread(posterior, grid):
    """Three-state barycentric marginal spread, not the conditional kernel."""
    need(type(posterior) is tuple and len(posterior) == 3
         and all(type(x) is F and x >= 0 for x in posterior)
         and sum(posterior) == 1, 'exact three-state posterior')
    need(type(grid) is int and grid > 0, 'positive exact grid')
    scaled = tuple(grid*x for x in posterior)
    low = tuple(x.numerator//x.denominator for x in scaled)
    rem = tuple(scaled[i]-low[i] for i in range(3))
    missing = grid-sum(low)
    choices = {}
    if missing == 0:
        choices[tuple(F(x, grid) for x in low)] = F(1)
    else:
        for i in range(3):
            if missing == 1:
                z = tuple(F(low[j]+(i == j), grid) for j in range(3))
                mass = rem[i]
            else:
                need(missing == 2, 'three-state rounding remainder')
                z = tuple(F(low[j]+(i != j), grid) for j in range(3))
                mass = 1-rem[i]
            if mass:
                choices[z] = choices.get(z, F(0))+mass
    need(sum(choices.values()) == 1, 'spread total mass')
    need(all(sum(m*z[i] for z, m in choices.items()) == posterior[i]
             for i in range(3)), 'barycenter is the old labeled posterior')
    return choices


def independent_membership(word, depth, n):
    x = n
    for a in word:
        if x == 1:
            return False
        x, actual = step(x)
        if actual != a:
            return False
    for _ in range(depth):
        if x % 8 != 5:
            return False
        x = (x-1)//4
    return 0 < x < n


def reject(f):
    try:
        f()
    except ValueError:
        return
    raise ValueError('hostile accepted')


def main():
    rows = instances = 0
    for length in range(6):
        for word in product(range(1, 5), repeat=length):
            p, q, _ = carrier(word)
            depth = 1
            while q*4**depth <= p:
                depth += 1
            for k in (depth, depth+1):
                row = compile_tail(word, k)
                rows += 1
                for t in (0, 1, 2, 7, 31, 10**12):
                    n = row.least+row.modulus*t
                    receipt = apply_tail(row, n)
                    need(receipt is not None and receipt.child == row.child0+2*p*t,
                         'positive parameter lift')
                    need(apply_tail(row, n+2) is None, 'adjacent nonmember')
                    instances += 1
    iff_controls = 0
    for length in range(4):
        for word in product(range(1, 4), repeat=length):
            p, q, _ = carrier(word)
            for k in (1, 2):
                if q*4**k <= p:
                    continue
                row = compile_tail(word, k)
                for n in range(1, 514, 2):
                    got = apply_tail(row, n) is not None
                    need(got == independent_membership(word, k, n),
                         'independent actual-word/sibling iff')
                    iff_controls += 1
    frontier = (1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2)
    named = compile_tail(frontier, 1)
    need((named.least, named.modulus, named.child0, 2*named.p) ==
         (111, 2**30, 81, 774840978), 'inherited frontier row recovered exactly')
    need(all(3**i > 2**sum(frontier[:i]) for i in range(1, 19)),
         'all retained frontier prefixes grow')
    need(3**19 > 2**30, 'residual-one common endpoint expands')
    short = compile_tail((1, 1, 1), 1)
    need((short.least, short.modulus, short.child0, 2*short.p) == (15, 64, 13, 54),
         'short symbolic row with a genuinely growing terminal half')
    for t in range(1, 64, 2):
        n = 15+64*t
        rec = apply_tail(short, n)
        need(rec.source_word == (1, 1, 1, 3) and rec.endpoint == 20+81*t
             and all(x > n for x in replay(n, rec.source_word)[1:]),
             'short reroute pays child while every source edge grows')
    for t in range(1, 256, 2):
        n = named.least+named.modulus*t
        rec = apply_tail(named, n)
        need(rec.source_word[-1] == 3 and all(x > n for x in replay(n, rec.source_word)[1:]),
             'paid sibling with every actual source checkpoint still above original')
    phase_controls = 0
    for row in (named, compile_tail((), 1), compile_tail((1, 2), 1)):
        horizon = 10
        counts = [0]*(horizon+2)
        for t in range(2**horizon):
            h = row.child0+2*row.p*t
            b = valuation(3*h+1)
            counts[min(b, horizon+1)] += 1
            for j in range(5):
                residue, modulus = reserve_phase(row, j)
                need((t % modulus == residue) == (b >= 2*j+1),
                     'exact extra-reserve parameter congruence')
                phase_controls += 1
        need(counts[1:horizon+1] == [2**(horizon-b) for b in range(1, horizon+1)]
             and counts[horizon+1] == 1, 'geometric finite-prefix law including unresolved tail')
    posterior_controls = 0
    for total in range(1, 13):
        for a in range(total+1):
            for b in range(total-a+1):
                p = (F(a, total), F(b, total), F(total-a-b, total))
                for grid in range(1, 6):
                    law = spread(p, grid)
                    for i in range(3):
                        if p[i]:
                            conditional = {z: mass*z[i]/p[i] for z, mass in law.items()}
                            need(sum(conditional.values()) == 1, 'spin-conditioned refinement kernel')
                            for z, mass in law.items():
                                need(p[i]*conditional[z] == mass*z[i], 'calibrated joint posterior')
                    posterior_controls += 1
    # Marginal spreading independent of S is not the refinement kernel.
    p = (F(1, 2), F(1, 3), F(1, 6))
    pure = spread(p, 1)
    need(pure[(F(1), F(0), F(0))] == F(1, 2), 'marginal pure weight')
    need(p[0] != 1, 'independent marginal draw has old posterior, not pure posterior')
    squared = tuple(x*x/sum(y*y for y in p) for x in p)
    need(squared == (F(9, 14), F(2, 7), F(1, 14)) and squared != p,
         'a duplicated datum is not two conditionally independent branch observations')
    # Perfectly labeled pure posteriors become no information if only their orbit is retained.
    labels = tuple(set(permutations((F(1), F(0), F(0)))))
    need(len(labels) == 3 and tuple(sum(z[i] for z in labels)/3 for i in range(3)) ==
         (F(1, 3),)*3, 'unlabeled permutation orbit loses spin information')
    # A proof suffix has unlimited reuse; injective matching would reject this lawful pair.
    need(replay(13, (3, 4))[-1] == replay(53, (5, 4))[-1] == 1,
         'two source obligations reuse the same certified child5')
    need(step(13)[0] == step(53)[0] == 5, 'two roles, one allowed reusable label')
    root_row = compile_tail((), 1)
    root_receipt = apply_tail(root_row, 5)
    need(root_receipt.child == 1 and root_receipt.child_word == ()
         and discharge(root_receipt, ()) == (4,), 'ROOT has empty suffix')
    receipt = apply_tail(root_row, 13)
    need(receipt is not None and discharge(receipt, (1, 4)) == (3, 4),
         'authenticate child3 and splice common endpoint5')
    # Dropping the exact word/source sidecar is unsound even for the old golden alias.
    need(step(223)[1] == 1 and step(233)[1] == 2, 'same golden observation, incompatible source guard')
    hostiles = [lambda: compile_tail((1, True), 1), lambda: compile_tail((1,), True),
                lambda: compile_tail((1,)*6, 1), lambda: apply_tail(root_row, True),
                lambda: apply_tail(root_row, 5.0), lambda: apply_tail(root_row, 4),
                lambda: apply_tail(replace(root_row, child0=3), 5),
                lambda: audit(replace(root_receipt, child=3)),
                lambda: discharge(root_receipt, (2,)), lambda: reserve_phase(root_row, True),
                lambda: replay(1, (2,)), lambda: spread((F(1), F(0)), 2)]
    for f in hostiles:
        reject(f)
    print('PROVED symbolic rows: original source, exact prefix, sibling depth, strict cuts retained.')
    print('FINITE-EXACT words length0..5 over1..4:', rows, 'rows;', instances, 'positive lifts.')
    print('Independent literal membership iff controls:', iff_controls)
    print('Recovered frontier:', named.least, '+', named.modulus, '*t ->',
          named.child0, '+', 2*named.p, '*t; t>=0.')
    print('Short row15+64t ->13+54t; odd t has actual1113 and all source prefixes growing.')
    print('Residual b>=1: exact frequency2^-b; extra j siblings frequency4^-j;', phase_controls, 'controls.')
    print('Three-state rational posterior/grid controls:', posterior_controls)
    print('Hostiles:', len(hostiles), '; no ROOT oracle or independent-child assumption in production.')
    print('Paper headlines are trusted premises; their global theorems are not re-audited here.')
    print('No new global Collatz coverage; named111 row is inherited, now generically compiled.')
    print('CHECKS', CHECKS)


if __name__ == '__main__':
    main()
