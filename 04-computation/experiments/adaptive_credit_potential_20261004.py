"""Exact guarded adaptive payment; no arbitrary-source or ROOT coverage claim.

Run with python -X utf8 -B, and repeat with -O. No import-time experiments.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations, product
from math import comb

H = (729, 1024, 669)
G = (9, 8, 5)
IDENTITY = (1, 1, 0)
HEAD = (1, 2, 1, 1, 1, 2)
KAPPA = F(2349, 2560)
ETA = F(15, 16)
ADAPTIVE_ETA = F(243, 256)


def need(condition, message):
    if not condition:
        raise ValueError(message)


def integer(x, low=0):
    need(type(x) is int and x >= low, 'exact integer in declared range')


def positive_odd(x):
    integer(x, 1)
    need(x % 2 == 1, 'positive odd source')


def symbols(word):
    need(type(word) is str and all(c in 'HG' for c in word), 'H/G word')


def compose(first, second):
    """Chronological composition; triples represent (P*x+B)/Q."""
    p, q, b = first
    r, s, c = second
    return r*p, s*q, r*b+c*q


@dataclass(frozen=True)
class Summary:
    affine: tuple = IDENTITY
    h: int = 0
    g: int = 0
    delta: int = 0
    minimum: int = 0


def join(left, right):
    need(type(left) is Summary and type(right) is Summary, 'summary type')
    return Summary(compose(left.affine, right.affine), left.h+right.h,
                   left.g+right.g, left.delta+right.delta,
                   min(left.minimum, left.delta+right.minimum))


def summary(word):
    symbols(word)
    result = Summary()
    for letter in word:
        item = Summary(H, 1, 0, 2, 0) if letter == 'H' else Summary(G, 0, 1, -1, -1)
        result = join(result, item)
    return result


def repeat_summary(item, count):
    need(type(item) is Summary, 'summary type')
    integer(count)
    result, operations = Summary(), 0
    while count:
        if count & 1:
            result = join(result, item)
            operations += 1
        count >>= 1
        if count:
            item = join(item, item)
            operations += 1
    return result, operations


def cylinder(item):
    need(type(item) is Summary, 'summary type')
    p, q, b = item.affine
    modulus = 2*q
    return ((q-b)*pow(p, -1, modulus)) % modulus, modulus


def affine_value(data, x):
    p, q, b = data
    return F(p*x+b, q)


def step(x):
    positive_odd(x)
    y = 3*x+1
    a = (y & -y).bit_length()-1
    return y >> a, a


def actual_word(x, word):
    positive_odd(x)
    need(type(word) is tuple and all(type(a) is int and a >= 1 for a in word),
         'positive exact valuation tuple')
    for a in word:
        need(x != 1, 'first-hit padding forbidden')
        x, seen = step(x)
        need(seen == a, 'actual valuation mismatch')
    return x


def action(x, letter, audit=True):
    positive_odd(x)
    need(letter in ('H', 'G'), 'primitive operation')
    if letter == 'H':
        need(x % 2048 == 155, 'H native guard')
        y = (729*x+669)//1024
        need(111 <= y < x, 'H output floor and strict payment')
        if audit:
            checkpoint = actual_word(x, HEAD)
            need(checkpoint == 4*y+1, 'H quarter-child identity')
            need(step(checkpoint)[0] == step(y)[0], 'H common future')
        return y
    need(x % 16 == 11, 'G native guard')
    y = (9*x+5)//8
    need(y > x, 'G strict growth')
    if audit:
        need(actual_word(x, (1, 2)) == y, 'G actual prefix')
    return y


def potential(x, credit):
    positive_odd(x)
    integer(credit)
    return (x+5)*F(9, 8)**credit


def interface_bound(data, credit_change, minimum):
    """Exact supremum of the potential ratio on the real half-line.

    An arithmetic/native guard is a separate required premise.
    """
    need(type(data) is tuple and len(data) == 3
         and all(type(v) is int for v in data), 'exact affine triple')
    p, q, b = data
    integer(p, 1)
    integer(q, 1)
    need(type(credit_change) is int, 'exact signed credit change')
    integer(minimum, 1)
    need(affine_value(data, minimum) > 0, 'positive half-line image')
    return F(9, 8)**credit_change * max(F(p, q),
           (affine_value(data, minimum)+5)/(minimum+5))


def root_rank(x):
    positive_odd(x)
    if x == 1:
        return (0, 0)
    delta = x-1
    k = (delta & -delta).bit_length()-1
    return (3**k*(delta >> k)**2, k)


def verify_episode(source, word):
    positive_odd(source)
    symbols(word)
    item = summary(word)
    need(item.minimum >= 0, 'unfunded prefix')
    residue, modulus = cylinder(item)
    need(source % modulus == residue, 'source cylinder')
    x, credit, h, g = source, 0, 0, 0
    trace = [x]
    for letter in word:
        old = potential(x, credit)
        if letter == 'G':
            need(credit >= 1, 'G credit guard')
            credit -= 1
            g += 1
        else:
            credit += 2
            h += 1
        x = action(x, letter)
        new = potential(x, credit)
        if letter == 'G':
            need(new == old, 'exact neutral growth debit')
        else:
            need(new <= KAPPA*old, 'uniform paid credit minting')
        need(116 <= new <= (source+5)*KAPPA**h, 'potential bounds')
        need(x < source and root_rank(x) < root_rank(source), 'immutable source payment')
        need(g <= 2*h and credit == 2*h-g, 'prefix token account')
        trace.append(x)
    need(F(x) == affine_value(item.affine, source), 'summary endpoint')
    return x, credit, tuple(trace)


def common_future_certificate(source, word):
    """Compile a paid dependency into two first-hit-safe actual prefixes."""
    endpoint, credit, trace = verify_episode(source, word)
    target, last = step(endpoint)
    route = (last,)
    for letter in reversed(word):
        route = HEAD+(route[0]+2,)+route[1:] if letter == 'H' else (1, 2)+route
    need(actual_word(source, route) == target, 'compiled original-source path')
    need(actual_word(endpoint, (last,)) == target, 'terminal common-future path')
    item = summary(word)
    need(len(route) == 6*item.h+2*item.g+1, 'compiled odd depth')
    need(sum(route) == 10*item.h+3*item.g+last, 'compiled halving cost')
    return endpoint, route, (last,)


def budget_bound(source):
    positive_odd(source)
    need(source >= 155, 'nonempty episode source floor')
    # Largest h allowed by 116 <= (source+5)*KAPPA**h.
    h, left, right = 0, 116, source+5
    while left*2560 <= right*2349:
        h += 1
        left *= 2560
        right *= 2349
    return h


@lru_cache(None)
def trees(size):
    integer(size)
    if size == 0:
        return (None,)
    result = []
    for a in range(size):
        for b in range(size-a):
            c = size-1-a-b
            for left, middle, right in product(trees(a), trees(b), trees(c)):
                result.append((left, middle, right))
    return tuple(result)


def tree_word(tree):
    if tree is None:
        return ''
    need(type(tree) is tuple and len(tree) == 3, 'ordered ternary node')
    a, b, c = tree
    return 'H'+tree_word(a)+'G'+tree_word(b)+'G'+tree_word(c)


def parse_tree(word):
    symbols(word)
    need(summary(word).minimum >= 0 and summary(word).delta == 0, 'balanced tree word')

    def parse(position):
        if position == len(word) or word[position] == 'G':
            return None, position
        left, position = parse(position+1)
        need(position < len(word) and word[position] == 'G', 'first close')
        middle, position = parse(position+1)
        need(position < len(word) and word[position] == 'G', 'second close')
        right, position = parse(position+1)
        return (left, middle, right), position

    tree, end = parse(0)
    need(end == len(word), 'tree parser consumed source')
    return tree


def rejected(callback):
    try:
        callback()
    except (ValueError, TypeError):
        return 1
    raise ValueError('hostile was incorrectly accepted')


def four_slot_controls():
    """Incoming paid permutations, with native guards and first-hit truncation."""
    count = rooted = 0
    for a in tuple(range(3, 15))+(31, 64):
        for tail in sorted(set(permutations((1, 2, a)))):
            word = (1,)+tail
            data = IDENTITY
            for valuation in word:
                data = compose(data, (3, 1 << valuation, 1))
            p, q, b = data
            need(p == 81 and q == 1 << (a+4), 'four-slot clocks')
            need(b <= 45+14*(1 << a), 'ordered carry upper bound')
            need(interface_bound(data, 1, 7) <= F(1023, 1024),
                 'one-credit uniform four-slot interface')
            modulus = 2*q
            residue = ((q-b)*pow(p, -1, modulus)) % modulus
            need(interface_bound(data, 1, residue) <= ETA,
                 'sharp family bound uses native source floor')
            for lift in (0, 1, 7):
                source = residue+lift*modulus
                need(source >= 7 and source % 4 == 3, 'four-slot source floor')
                x = source
                for index, valuation in enumerate(word):
                    if x == 1:
                        need(all(v == 2 for v in word[index:]), 'only root padding remains')
                        rooted += 1
                        break
                    x, seen = step(x)
                    need(seen == valuation, 'four-slot native valuation')
                need(F(x) == affine_value(data, source), 'four-slot endpoint')
                need(potential(x, 1) <= ETA*potential(source, 0),
                     'four-slot credit payment')
                award = 1 if source == 7 else 2
                need(potential(x, award) <= ADAPTIVE_ETA*potential(source, 0),
                     'source-sensitive four-slot credit payment')
                need(x < source and root_rank(x) < root_rank(source),
                     'four-slot immutable-source rank')
                count += 1
    # These exact incoming cylinders can be imported at zero credit too.
    need(interface_bound((81, 128, 85), 0, 187) == F(31, 48), '1213 zero-credit bound')
    need(interface_bound((81, 128, 73), 0, 7) == F(5, 6), '1123 zero-credit bound')
    need(actual_word(7, (1, 1, 2, 3)) == 5, 'two-credit hostile is native')
    need(potential(5, 1)/potential(7, 0) == ETA, 'sharp bank contraction')
    need(actual_word(19, (1, 3, 1, 2)) == 13, 'sharp adaptive bank witness')
    need(potential(13, 2)/potential(19, 0) == ADAPTIVE_ETA, 'sharp adaptive contraction')
    need(potential(5, 2) > potential(7, 0), 'two credits fail for the full bank')
    return count, rooted


def extended_controls(source_sensitive=False):
    """Adaptive compositions of H, G, and the two concrete imported patterns."""
    patterns = {'H': H, 'G': G, 'A': (81, 128, 85), 'B': (81, 128, 73)}
    credits = {'H': 2, 'G': -1, 'A': 1, 'B': 1}
    if source_sensitive:
        credits.update(A=2, B=2)
    factor = ADAPTIVE_ETA if source_sensitive else ETA
    actual = {'A': (1, 2, 1, 3), 'B': (1, 1, 2, 3)}
    words = [''.join(w) for size in range(1, 5) for w in product('HGAB', repeat=size)]
    words += ['HGGAGBG', 'AHGGBG', 'BGAHGGBG']
    checked = 0
    for word in words:
        balance = 0
        funded = True
        data = IDENTITY
        for letter in word:
            balance += credits[letter]
            funded &= balance >= 0
            data = compose(data, patterns[letter])
        if not funded:
            continue
        p, q, b = data
        modulus = 2*q
        residue = ((q-b)*pow(p, -1, modulus)) % modulus
        for lift in (0, 1, 7):
            source = residue+lift*modulus
            need(source % 4 == 3, 'funded extended entry has rank K1')
            x, balance, funding = source, 0, 0
            for letter in word:
                award = credits[letter]
                if source_sensitive and letter in 'AB' and x == 7:
                    award = 1
                balance += award
                need(balance >= 0, 'actual source-sensitive prefix credit')
                funding += letter != 'G'
                x = action(x, letter) if letter in 'HG' else actual_word(x, actual[letter])
                energy = potential(x, balance)
                need(6 <= energy <= (source+5)*factor**funding,
                     'extended prefix potential and finite bound')
                need(x < source and root_rank(x) < root_rank(source),
                     'extended prefix original-source payment')
            need(F(x) == affine_value(data, source), 'extended exact cylinder endpoint')
            checked += 1
    return checked


def experiment():
    cylinder_checks = ballot_checks = compiled_checks = 0
    # All words, including unfunded and G-first controls; all positive lifts named.
    for length in range(1, 9):
        for letters in product('HG', repeat=length):
            word = ''.join(letters)
            item = summary(word)
            residue, modulus = cylinder(item)
            for lift in (0, 1, 7):
                source = residue+lift*modulus
                x = source
                for letter in word:
                    x = action(x, letter)
                need(F(x) == affine_value(item.affine, source), 'all-word cylinder endpoint')
                cylinder_checks += 1
                if item.minimum >= 0:
                    endpoint, credit, trace = verify_episode(source, word)
                    need(endpoint == x, 'independent transition agreement')
                    need(item.h <= budget_bound(source), 'finite episode bound')
                    common_future_certificate(source, word)
                    need(interface_bound(item.affine, item.delta, residue)
                         <= KAPPA**item.h, 'composite pattern interface bound')
                    ballot_checks += 1
                    compiled_checks += 1

    counts = []
    tree_checks = 0
    for m in range(6):
        bank = trees(m)
        need(len(bank) == comb(3*m, m)//(2*m+1), 'ternary count')
        counts.append(len(bank))
        for tree in bank:
            word = tree_word(tree)
            need(parse_tree(word) == tree, 'lossless ternary round trip')
            item = summary(word)
            need((item.h, item.g, item.delta, item.minimum) == (m, 2*m, 0, 0), 'tree credits')
            if m:
                residue, modulus = cylinder(item)
                need(modulus == 1 << (16*m+1), 'balanced source precision')
                verify_episode(residue, word)
            tree_checks += 1

    short = [''.join(w) for size in range(4) for w in product('HG', repeat=size)]
    joins = 0
    for a, b in product(short, repeat=2):
        need(join(summary(a), summary(b)) == summary(a+b), 'summary composition')
        joins += 1
    powers = 0
    for word in ('HGG', 'HHGGGG', 'HGHGGG', 'HGGHGG'):
        for count in (0, 1, 2, 3, 8, 32):
            result, operations = repeat_summary(summary(word), count)
            need(result == summary(word*count), 'binary recursive summary')
            need(operations <= 2*max(1, count.bit_length()), 'matrix operation count')
            powers += 1

    need(potential(111, 2)/potential(155, 0) == KAPPA, 'sharp H credit factor')
    need(interface_bound(H, 2, 155) == KAPPA, 'H finite interface contract')
    need(interface_bound(G, -1, 11) == 1, 'G neutral interface contract')
    need(interface_bound(H, 3, 155) > 1, 'third token cannot be minted freely')
    hostile = summary('HGGG')
    need(hostile.affine == (531441, 524288, 1598741), 'third-credit affine witness')
    source, modulus = cylinder(hostile)
    x = source
    for letter in 'HGGG':
        x = action(x, letter)
    need(x > source, 'third-credit size-payment failure')
    rank_source = source+modulus
    rank_endpoint = rank_source
    for letter in 'HGGG':
        rank_endpoint = action(rank_endpoint, letter)
    need((rank_source, rank_endpoint) == (1159323, 1175143), 'rank hostile source identity')
    need(root_rank(rank_endpoint) > root_rank(rank_source), 'third-credit rank failure')
    need(summary('HGG').delta == summary('GGH').delta == 0, 'same net credit')
    need(summary('HGG').minimum == 0 and summary('GGH').minimum == -2, 'order loss hostile')
    bad = [lambda: verify_episode(11, 'G'),
           lambda: verify_episode(155, 'HG'),
           lambda: verify_episode(source, 'HGGG'),
           lambda: verify_episode(155.0, 'H'),
           lambda: verify_episode(155, ('H',)),
           lambda: parse_tree('HGGG'),
           lambda: parse_tree('GHG'),
           lambda: repeat_summary(summary('H'), -1)]
    rejected_count = sum(rejected(f) for f in bad)
    four_slot_checks, first_hit_truncations = four_slot_controls()
    extended_checks = extended_controls()
    adaptive_checks = extended_controls(source_sensitive=True)
    print('Native types: H is a common-future dependency; G is actual word12.')
    print('E=(x+5)*(9/8)^credit; G preserves E; H ratio<=2349/2560, equality at155.')
    print('Every nonempty legal ballot prefix pays its immutable source and inherited root rank.')
    print('Every legal episode is finite:116*2560^h<=(N+5)*2349^h and g<=2h.')
    print('All H/G words length1..8, lifts0,1,7:', cylinder_checks, 'exact-cylinder/native checks.')
    print('Ballot prefix payment and common-future certificate replays:', ballot_checks, compiled_checks)
    print('Ordered ternary tree counts m0..5:', counts, '; exact round trips:', tree_checks)
    print('Summary concatenations:', joins, '; compressed repeat controls:', powers)
    print('Exact half-line interface bound audits each funded composite with its original credit ledger.')
    print('Incoming four-slot bank earns one credit, sharp uniform ratio<=15/16:',
          four_slot_checks, 'native cases;', first_hit_truncations, 'first-hit truncations.')
    print('Extended finite controller: E>=6, G-count<=2*H-count+four-slot-count; two-credit bank hostile7->5.')
    print('Mixed H/G/1213/1123 funded programs through length4 plus three deeper switches:',
          extended_checks, 'exact source and prefix replays.')
    print('Source-sensitive award P: one credit at7, two elsewhere; sharp factor243/256;',
          adaptive_checks, 'mixed source replays.')
    print('Three-credit native witness:', source, '->', x, '; source modulus', modulus)
    print('Three-credit native rank failure:', rank_source, '->', rank_endpoint)
    print('Same net credit does not preserve prefix legality: HGG versus GGH.')
    print('Malformed/type/guard/budget controls rejected:', rejected_count)
    print('Core H/G entry is exactly155mod2048; extended entry adds the guarded four-slot bank, not arbitrary sources.')
    print('No automatic external refill or universal ROOT claim.')


if __name__ == '__main__':
    experiment()
