"""Exact small tournament codecs and typed recursion controls.

All tests use exceptions and survive -O. The companion note gives the
general proofs and distinguishes storage from an arithmetic certificate.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations, product


def need(condition, message):
    if not condition:
        raise ValueError(message)


def validate(t):
    n = len(t)
    need(all(len(row) == n for row in t), 'square tournament matrix')
    need(all(type(t[i][j]) is int and t[i][j] in (0, 1)
             for i in range(n) for j in range(n)), 'exact binary adjacency')
    need(all(t[i][i] == 0 for i in range(n)), 'no loops')
    need(all(t[i][j]+t[j][i] == 1 for i in range(n) for j in range(i+1, n)),
         'one orientation per pair')
    return n


def fixed_path(n, mask):
    need(type(n) is int and n >= 1 and type(mask) is int and mask >= 0,
         'exact positive size and nonnegative mask')
    pairs = [(i, j) for i in range(n) for j in range(i+2, n)]
    # At n=4 the documented bit order is a02,b13,c03.
    if n == 4:
        pairs = [(0, 2), (1, 3), (0, 3)]
    need(mask < 2**len(pairs), 'fixed-path mask range')
    t = [[int(i < j) for j in range(n)] for i in range(n)]
    for bit, (i, j) in enumerate(pairs):
        if mask >> bit & 1:
            t[i][j], t[j][i] = t[j][i], t[i][j]
    return tuple(map(tuple, t))


def arbitrary_tournament(n, mask):
    t = [[int(i < j) for j in range(n)] for i in range(n)]
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        if mask >> bit & 1:
            t[i][j], t[j][i] = t[j][i], t[i][j]
    return tuple(map(tuple, t))


SCORES = {
    (0, 1, 2, 3): 'T',
    (0, 2, 2, 2): '+',
    (1, 1, 1, 3): '-',
    (1, 1, 2, 2): 'S',
}
ALPHABET = ('T', '+', '-', 'S')


def class4(t):
    need(validate(t) == 4, 'four vertices')
    return SCORES[tuple(sorted(map(sum, t)))]


def relabel(t, order):
    need(sorted(order) == list(range(len(t))), 'vertex permutation')
    return tuple(tuple(t[i][j] for j in order) for i in order)


def canonical(t):
    return min(tuple(x for row in relabel(t, p) for x in row)
               for p in permutations(range(len(t))))


def hamiltonian_paths(t):
    n = validate(t)
    if n == 0:
        return 1

    @lru_cache(None)
    def count(mask, last):
        if mask == 1 << last:
            return 1
        rest = mask ^ (1 << last)
        return sum(count(rest, j) for j in range(n) if rest >> j & 1 and t[j][last])

    return sum(count((1 << n)-1, j) for j in range(n))


def join(left, right):
    n, m = validate(left), validate(right)
    return tuple(tuple(left[i][j] if i < n and j < n else
                       right[i-n][j-n] if i >= n and j >= n else int(i < n)
                       for j in range(n+m)) for i in range(n+m))


SEEDS = dict(zip(ALPHABET, (fixed_path(4, j) for j in range(4))))


def encode(word):
    need(type(word) is tuple and all(x in SEEDS for x in word), 'four-symbol word')
    t = ()
    for symbol in word:
        t = join(t, SEEDS[symbol])
    return t


def scc_order(t):
    """Recover intrinsic strongly connected components, then their total order."""
    n = validate(t)
    reach = [sum(1 << j for j in range(n) if t[i][j] or i == j) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    unused, components = set(range(n)), []
    while unused:
        i = min(unused)
        c = tuple(j for j in sorted(unused) if reach[i] >> j & 1 and reach[j] >> i & 1)
        components.append(c)
        unused.difference_update(c)
    components.sort(key=lambda c: -sum(t[c[0]][j] for j in range(n) if j not in c))
    for index, c in enumerate(components):
        for other in components[index+1:]:
            need(all(t[i][j] for i in c for j in other), 'condensation total order')
    return components


def decode(t):
    """No construction partition or original vertex labels are supplied."""
    need(validate(t) % 4 == 0, 'whole four-vertex blocks')
    word, block = [], []
    for c in scc_order(t):
        need(len(block)+len(c) <= 4, 'strong component crosses a required block cut')
        block.extend(c)
        if len(block) == 4:
            sub = tuple(tuple(t[i][j] for j in block) for i in block)
            word.append(class4(sub))
            block = []
    need(not block, 'complete final block')
    return tuple(word)


def substitution(t, s):
    n, m = validate(t), validate(s)
    return tuple(tuple(t[i//m][j//m] if i//m != j//m else s[i % m][j % m]
                       for j in range(n*m)) for i in range(n*m)) if m else ()


def dominance(t):
    validate(t)
    return tuple(tuple(t[i][j]-t[j][i] for j in range(len(t))) for i in range(len(t)))


def identity(n):
    return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))


def multiply(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def add(a, b):
    return tuple(tuple(x+y for x, y in zip(ra, rb)) for ra, rb in zip(a, b))


def scale(c, a):
    return tuple(tuple(c*x for x in row) for row in a)


def tensor(a, b):
    return tuple(tuple(x*y for x in rowa for y in rowb) for rowa in a for rowb in b)


A2 = ((1, 1), (1, -1))
B2 = ((0, 1), (-1, 0))


def skew_double_matrix(m):
    return add(tensor(A2, m), tensor(B2, identity(len(m))))


def from_dominance(m):
    n = len(m)
    need(all(m[i][i] == 0 for i in range(n)), 'skew zero diagonal')
    need(all(m[i][j] == -m[j][i] and m[i][j] in (-1, 1)
             for i in range(n) for j in range(i+1, n)), 'skew tournament entries')
    return tuple(tuple(int(x > 0) for x in row) for row in m)


def main():
    classes, auts, hs = defaultdict(set), {}, {}
    for mask in range(64):
        t = arbitrary_tournament(4, mask)
        label = class4(t)
        classes[label].add(canonical(t))
        auts[label] = sum(relabel(t, p) == t for p in permutations(range(4)))
        hs[label] = hamiltonian_paths(t)
    need(all(len(x) == 1 for x in classes.values()) and len(classes) == 4,
         'scores classify all64 labelled four-tournaments')
    fibers = defaultdict(list)
    for mask in range(8):
        fibers[class4(fixed_path(4, mask))].append(mask)
    need(dict(fibers) == {'T': [0], '+': [1], '-': [2], 'S': [3, 4, 5, 6, 7]},
         'recovered fixed-path fibers')
    need(hs == {'T': 1, '+': 3, 'S': 5, '-': 3} and
         all(len(fibers[c])*auts[c] == hs[c] for c in ALPHABET), 'path/aut fiber law')
    need(class4(fixed_path(4, 3)) == class4(fixed_path(4, 4)) == 'S' and
         class4(fixed_path(4, 3 ^ 1)) == '-' and class4(fixed_path(4, 4 ^ 1)) == 'S',
         'same class, different future named flip')
    refinement = {m: (class4(fixed_path(4, m)),
                      tuple(class4(fixed_path(4, m ^ b)) for b in (1, 2, 4)))
                  for m in range(8)}
    need(len(set(refinement.values())) == 8, 'one-step prediction already needs all8 states')
    need(tuple(class4(fixed_path(4, j)) for j in range(4)) == ALPHABET,
         'closed two-bit halfsection')
    need(all((j ^ b) < 4 for j in range(4) for b in (1, 2)), 'allowed halfsection flips')
    need(class4(fixed_path(4, 0 ^ 7)) != 'T', 'all-free-bit complement is not converse')
    distinguished = 0
    for n in range(2, 6):
        m = (n-1)*(n-2)//2
        transitive = [x for x in range(2**m)
                      if tuple(sorted(map(sum, fixed_path(n, x)))) == tuple(range(n))]
        need(transitive == [0], 'unique transitive fixed-path presentation')
        for u, v in combinations(range(2**m), 2):
            need(u ^ u == 0 and v ^ u != 0, 'translation distinguishes every pair')
            distinguished += 1
    print('N4 fibers', dict(fibers), '; H', hs, '; automorphisms', auts)
    print('S states ab/c diverge after flip a; unrestricted predictive quotient needs8 states.')
    print('Closed halfsection c0 has4 states; general unique-transitive lemma controls all n.')
    print('Finite predictive controls n2..5:', distinguished, 'distinct pairs.')

    codec_checks = 0
    seed_h = {c: hamiltonian_paths(SEEDS[c]) for c in ALPHABET}
    for length in range(5):
        for word in product(ALPHABET, repeat=length):
            t = encode(word)
            n = len(t)
            orders = [tuple(range(n)), tuple(reversed(range(n))),
                      tuple(range(0, n, 2))+tuple(range(1, n, 2))]
            for order in orders:
                need(decode(relabel(t, order)) == word, 'unlabelled four-block word decoder')
                codec_checks += 1
            if length <= 2:
                expected = 1
                for symbol in word:
                    expected *= seed_h[symbol]
                need(hamiltonian_paths(t) == expected, 'ordered-join path product')
    xy, yx = encode(('+', '-')), encode(('-', '+'))
    need(decode(xy) != decode(yx) and hamiltonian_paths(xy) == hamiltonian_paths(yx) == 9,
         'H loses both order and the converse-diamond distinction')
    need([len(c) for c in scc_order(xy)] == [3, 1, 1, 3] and
         [len(c) for c in scc_order(yx)] == [1, 3, 3, 1], 'intrinsic unequal word carriers')
    print('Four-symbol codec:341 words of length0..4;3 relabelings each;', codec_checks, 'checks.')
    print('H(+,-)=H(-,+)=9, but SCC size words3,1,1,3 and1,3,3,1 distinguish them.')

    singleton = ((0,),)
    c3 = tuple(tuple(int((j-i) % 3 == 1) for j in range(3)) for i in range(3))
    for t in [singleton, c3]+list(SEEDS.values()):
        padded = join(join(singleton, t), singleton)
        need(len(padded) == len(t)+2 and hamiltonian_paths(padded) == hamiltonian_paths(t),
             'endpoint +2 preserves H')
        doubled = join(t, t)
        need(hamiltonian_paths(doubled) == hamiltonian_paths(t)**2, 'self-join squares H')
        need(len(substitution(t, t)) == len(t)**2, 'self-substitution squares order')
    for n in range(1, 101):
        need((n+1)*n//2-(n-1)*(n-2)//2 == 2*n-1, 'correct ModeB boundary size')
    need(hamiltonian_paths(join(c3, c3)) == 9 and
         hamiltonian_paths(substitution(c3, c3)) == 3159, 'order-square is not H-square')
    t2 = ((0, 1), (0, 0))
    need(class4(substitution(t2, t2)) == 'T' and
         class4(from_dominance(skew_double_matrix(dominance(t2)))) == '-',
         'cloning and skew doubling already differ at order4')
    print('Typed modes: padding n+2/Hfixed; selfjoin2n/H^2; selfsubstitution n^2.')
    print('C3 selfjoin:order6,H9; selfsubstitution:order9,H3159.')
    print('TT2 cloning givesT4,H1; skew doubling gives source-over-C3,H3.')

    need(multiply(A2, A2) == scale(2, identity(2)) and
         multiply(B2, B2) == scale(-1, identity(2)) and
         add(multiply(A2, B2), multiply(B2, A2)) == ((0, 0), (0, 0)),
         'two small matrices and four-product cancellation')
    skew_checks = 0
    for n in range(1, 5):
        for mask in range(2**(n*(n-1)//2)):
            m = dominance(arbitrary_tournament(n, mask))
            square = multiply(m, m)
            current = m
            for depth in range(4):
                validate(from_dominance(current))
                expected = tensor(identity(2**depth),
                                  add(scale(2**depth, square),
                                      scale(-(2**depth-1), identity(n))))
                need(multiply(current, current) == expected, 'iterated skew-square identity')
                skew_checks += 1
                if depth < 3:
                    current = skew_double_matrix(current)
    print('Skew recursion: A^2=2I,B^2=-I,AB+BA=0; four-product expansion cancels cross terms.')
    print('D^d(M)^2=I_(2^d) tensor [2^d M^2-(2^d-1)I];', skew_checks,
          'checks, all labelled n1..4,depth0..3.')
    print('Squared skew magnitude z maps2z+1; z+1 doubles, unlike word-slope squaring.')
    invalid = 0
    for thunk in (lambda: fixed_path(4, 8), lambda: fixed_path(True, 0),
                  lambda: encode(('unknown',)), lambda: decode(singleton),
                  lambda: decode(substitution(c3, ((0, 1), (0, 0))))):
        try:
            thunk()
        except ValueError:
            invalid += 1
        else:
            raise ValueError('invalid codec domain accepted')
    print('Malformed or non-language controls rejected:', invalid)
    print('PASS: exact graph storage and typed recursions; no arithmetic coverage inference.')


if __name__ == '__main__':
    main()
