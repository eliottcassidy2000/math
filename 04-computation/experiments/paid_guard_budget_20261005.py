"""Four-credit paid dependency on219mod256; exact guard-preserving composition.

Run normally and with -O. Imported prior checks have no import-time experiment.
All coverage statements are restricted to explicit native guards.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
import adaptive_credit_potential_20261004 as prior


def need(ok, message):
    if not ok:
        raise ValueError(message)


@dataclass(frozen=True)
class Guarded:
    p: int
    q: int
    b: int
    residue: int
    modulus: int


I = Guarded(1, 1, 0, 1, 2)
H = Guarded(729, 1024, 669, 155, 2048)
G = Guarded(9, 8, 5, 11, 16)
K = Guarded(81, 128, 53, 219, 256)
L = Guarded(9, 16, -3, 219, 256)
OPS = {'H': H, 'G': G, 'L': L}
CREDIT = {'H': 2, 'G': -1, 'L': 4}
RHO = F(6561, 7168)
ETA = F(243, 256)


def compose(u, v):
    """Chronological composition, retaining exact dyadic source intersections."""
    a = ((u.q*v.residue-u.b)*pow(u.p, -1, u.q*v.modulus)) % (u.q*v.modulus)
    m = u.q*v.modulus
    need((a-u.residue) % min(m, u.modulus) == 0, 'empty native guard intersection')
    r, modulus = (a, m) if m >= u.modulus else (u.residue, u.modulus)
    return Guarded(v.p*u.p, v.q*u.q, v.p*u.b+v.b*u.q, r, modulus)


def carrier(word):
    need(type(word) is str and all(c in OPS for c in word), 'declared alphabet')
    out = I
    for c in word:
        out = compose(out, OPS[c])
    return out


def value(op, n):
    prior.positive_odd(n)
    need(n % op.modulus == op.residue, 'retained native source guard')
    p = op.p*n+op.b
    need(p > 0 and p % op.q == 0, 'positive integral endpoint')
    y = p//op.q
    prior.positive_odd(y)
    return y


def local_l_receipt(n):
    child = value(L, n)
    k = value(K, n)
    need(prior.actual_word(child, (1, 2)) == k, 'inverse-G predecessor is actual')
    checkpoint = prior.actual_word(n, (1, 2, 1, 1))
    need(checkpoint == 4*k+1, 'source quarter-child checkpoint')
    target, a = prior.step(k)
    source_word = (1, 2, 1, 1, a+2)
    child_word = (1, 2, a)
    need(prior.actual_word(n, source_word) == target, 'source native join')
    need(prior.actual_word(child, child_word) == target, 'child native join')
    need(child < k < n, 'two strict smaller dependencies')
    need(prior.root_rank(child) < prior.root_rank(n), 'original-source rank')
    need(prior.potential(child, 4) <= RHO*prior.potential(n, 0), 'four-credit mint')
    need(prior.potential(k, 3) == prior.potential(child, 4), 'inverse-G credit returned')
    need((k+5) % 9 == 0 and (k+5) % 27 != 0 and (child+5) % 3 != 0,
         'exactly one integral inverse-G step')
    return child, source_word, child_word


def execute(n, word, funded=True):
    op = carrier(word)
    expected = value(op, n)
    x, credit, funders = n, 0, 0
    for c in word:
        x = value(OPS[c], x)
        credit += CREDIT[c]
        funders += c != 'G'
        if funded:
            need(credit >= 0, 'nonnegative actual prefix credit')
            need(prior.potential(x, credit) <= (n+5)*ETA**funders, 'global potential')
            need(x < n and prior.root_rank(x) < prior.root_rank(n), 'immutable payment')
    need(x == expected, 'guarded composition endpoint')
    return x


def compile_join(n, word):
    """Actual common-future words; extend the shared tail only when needed."""
    states = [n]
    for c in word:
        states.append(value(OPS[c], states[-1]))
    terminal = states[-1]
    join_value = terminal
    route = ()
    terminal_route = ()
    for c in reversed(word):
        required = {'H': 1, 'G': 0, 'L': 3}[c]
        while len(route) < required:
            need(join_value != 1, 'native dependency forbids premature root padding')
            join_value, a = prior.step(join_value)
            route += (a,)
            terminal_route += (a,)
        if c == 'G':
            route = (1, 2)+route
        elif c == 'H':
            route = prior.HEAD+(route[0]+2,)+route[1:]
        else:
            need(route[:2] == (1, 2), 'L child retains forced12 prefix')
            route = (1, 2, 1, 1, route[2]+2)+route[3:]
    need(prior.actual_word(n, route) == join_value, 'compiled source word')
    need(prior.actual_word(terminal, terminal_route) == join_value, 'compiled terminal word')
    need(len(route)-len(terminal_route) == 6*word.count('H')+2*word.count('G')+2*word.count('L'),
         'relative odd-depth identity')
    need(sum(route)-sum(terminal_route) == 10*word.count('H')+3*word.count('G')+4*word.count('L'),
         'relative halving-cost identity')
    return len(route), len(terminal_route)


def rejected(f):
    try:
        f()
    except ValueError:
        return 1
    raise ValueError('hostile accepted')


def run():
    need(compose(L, G) == K, 'exact guarded LG equals K')
    need(prior.interface_bound((L.p, L.q, L.b), 4, 219) == RHO, 'sharp L interface')
    need(prior.interface_bound((K.p, K.q, K.b), 3, 219) == RHO, 'sharp K interface')
    need(RHO < ETA, 'import into previous finite controller')
    for t in (0, 1, 2, 7, 31, 127, 1024):
        local_l_receipt(219+256*t)
    source_checks = funded_checks = receipt_checks = 0
    # All three-letter words, including unpaid controls, and three positive lifts.
    for length in range(1, 7):
        for letters in product('HGL', repeat=length):
            word = ''.join(letters)
            op = carrier(word)
            balance = 0
            funded = True
            for c in word:
                balance += CREDIT[c]
                funded &= balance >= 0
            for lift in (0, 1, 7):
                n = op.residue+lift*op.modulus
                execute(n, word, funded=funded)
                compile_join(n, word)
                source_checks += 1
                receipt_checks += 1
                if funded:
                    funded_checks += 1
    # Native precision, not just final oddness: L's formula admits27mod32.
    need((L.q-L.b)*pow(L.p, -1, 2*L.q) % (2*L.q) == 27, 'coarse affine cylinder')
    need((9*27-3)//16 == 15 and 27 % 256 != 219, 'lost-three-bit witness')
    need(prior.actual_word(27, (1, 2)) == 31 and (8*31-5)//9 == 27,
         'unpaid G/inverse-G round trip')
    # Five growth credits really overspend on a fully native cylinder.
    overspend = carrier('L'+'G'*5)
    need((overspend.p, overspend.q, overspend.b) == (531441, 524288, 1925333),
         'fifth-credit formal map')
    n = overspend.residue
    y = execute(n, 'L'+'G'*5, funded=False)
    need(y > n, 'legal fifth-credit size failure')
    compile_join(n, 'L'+'G'*5)
    rank_failure = None
    for lift in range(4):
        a = n+lift*overspend.modulus
        b = execute(a, 'L'+'G'*5, funded=False)
        if prior.root_rank(b) > prior.root_rank(a):
            rank_failure = (a, b)
            break
    need(rank_failure is not None, 'native fifth-credit rank hostile')
    bad = [lambda: value(L, 27), lambda: execute(n, 'L'+'G'*5),
           lambda: carrier('X'), lambda: value(L, 219.0)]
    need(sum(rejected(f) for f in bad) == 4, 'malformed and budget hostiles')
    print('PROVED scoped dependency219+256t ->123+144t; intermediate K139+162t.')
    print('Exact guarded L then G equals K; exactly one integral inverse-G step is available.')
    print('L earns4 credits, K earns3; sharp potential ratio6561/7168 at source219.')
    print('Imported controller factor243/256 survives; G count<=2H+4L+2P bounds all legal episodes.')
    print('All H/G/L words length1..6, lifts0,1,7:', source_checks,
          'native/composition controls;', funded_checks, 'funded prefix controls;', receipt_checks, 'compiled joins.')
    print('Affine-only L forgets3 guard bits: formula27mod32, native219mod256;27->15 rejected.')
    print('Fifth-credit legal cylinder:', overspend.residue, 'mod', overspend.modulus, ';', n, '->', y)
    print('Fifth-credit rank failure:', rank_failure, '; four malformed/budget controls rejected.')
    print('No universal guard coverage or supplied ROOT proof for arbitrary child is asserted.')


if __name__ == '__main__':
    run()
