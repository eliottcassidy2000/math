"""Complete dyadic-response censuses and exact characteristic-two decoders.

Only modular determinants are observed, not exact singularity or marked
cokernel classes. Active gates survive python -O. No external libraries.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import json

GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def norm(a, b):
    return a*a + a*b + b*b


def signature(m, z, a, b, full=False):
    n = m if full else min(m, 5)
    return tuple((z-norm(a+u, b+v))**2 % m
                 for u, v in product(range(n), repeat=2))


def key(k, z, a, b):
    m = 2**k
    if k == 0:
        return (0, 0, 0)
    if k <= 2:
        return (z % 2, a % 2, b % 2)
    if k == 3:
        p = 4 if z % 2 else 2
    else:
        p = m//4 if z % 2 else m//2
    return (z % (m//2), a % p, b % p)


def dyadic_count(k):
    return 1 if k == 0 else 8 if k <= 2 else 5*2**(3*k-6)


def count(m):
    a = b = 0
    while m % 2 == 0:
        m //= 2
        a += 1
    while m % 3 == 0:
        m //= 3
        b += 1
    return dyadic_count(a)*(1 if b == 0 else 7*27**(b-1))*m**3


def census(k):
    m = 2**k
    by_key, by_response = {}, {}
    histogram = Counter()
    fixed_nine = set()
    for a, b in product(range(m), repeat=2):
        q = [norm(a+u, b+v) % m
             for u, v in product(range(min(m, 5)), repeat=2)]
        for z in range(m):
            response = bytes((z-x)**2 % m for x in q)
            predicted = key(k, z, a, b)
            need(by_key.setdefault(predicted, response) == response,
                 'same predicted class has the same complete function')
            need(by_response.setdefault(response, predicted) == predicted,
                 'complete functions separate all predicted classes')
            histogram[response] += 1
            if z == 9 % m:
                fixed_nine.add(response)
            if m <= 8:
                need(tuple(response) == signature(m, z, a, b), 'grid path')
                # Direct full response is determined by the predicted key;
                # compare against its representative, without interpolation.
                need(signature(m, z, a, b, True) ==
                     signature(m, *predicted, True), 'every residue query')
    need(len(by_response) == dyadic_count(k), 'dyadic count formula')
    fibers = dict(sorted(Counter(histogram.values()).items()))
    need(sum(size*n for size, n in fibers.items()) == m**3, 'all states counted')
    if k >= 3:
        need(fibers == {8: 4*2**(3*k-6), 32: 2**(3*k-6)}, 'fiber histogram')
    children = Counter()
    if k:
        for child in by_key:
            children[key(k-1, *child)] += 1
        for parent, number in children.items():
            expected = (8 if k == 1 else 1 if k == 2 else
                        (8 if parent[0] % 2 else 2) if k == 3 else
                        (2 if parent[0] % 2 else 32) if k == 4 else 8)
            need(number == expected, 'exact precision-transition branching')
    fixed = (1 if k == 0 else 4 if k <= 2 else 16 if k == 3
             else 2**(2*k-4))
    need(len(fixed_nine) == fixed, 'fixed source z=9 phase count')
    return {'modulus': m, 'states': m**3, 'classes': len(by_response),
            'fibers': fibers, 'fixed_z9_classes': len(fixed_nine),
            'children_per_parent_histogram': dict(sorted(Counter(children.values()).items())),
            'all_residue_queries': m <= 8}


def polyrem(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a


def polygcd(a, b):
    while b:
        a, b = b, polyrem(a, b)
    return a


def mul(a, b, modulus):
    out = 0
    top = 1 << (modulus.bit_length()-1)
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        if a & top:
            a ^= modulus
        b >>= 1
    return out


def carryless(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        b >>= 1
    return out


def field_test(r, polynomial):
    q = 2**r
    need(polynomial.bit_length() == r+1, 'specified field polynomial degree')
    # Rabin's elementary irreducibility criterion, explicitly evaluated.
    powers = [2]
    for _ in range(r):
        powers.append(mul(powers[-1], powers[-1], polynomial))
    divisors = [p for p in range(2, r+1) if r % p == 0 and
                all(p % d for d in range(2, p))]
    need(powers[r] == 2, 'x^(2^r)=x modulo the field polynomial')
    for p in divisors:
        need(polygcd(powers[r//p] ^ 2, polynomial) == 1,
             'no irreducible factor of a proper divisor degree')
    table = [[mul(a, b, polynomial) for b in range(q)] for a in range(q)]
    for a, b in product(range(q), repeat=2):
        need(table[a][b] == polyrem(carryless(a, b), polynomial),
             'independent multiplication by long division')
    squares = [table[a][a] for a in range(q)]
    need(len(set(squares)) == q, 'squaring is a permutation')
    root = [0]*q
    for a, aa in enumerate(squares):
        root[aa] = a

    def Q(a, b):
        return squares[a] ^ table[a][b] ^ squares[b]

    cube_roots = [x for x in range(q) if Q(x, 1) == 0]
    need(len(cube_roots) == 2, 'the even-degree field splits the ternary clock')
    need({squares[x] for x in cube_roots} == set(cube_roots) and
         all(squares[x] != x for x in cube_roots), 'Frobenius exchanges the two clock roots')
    states = 0
    for a, b in product(range(q), repeat=2):
        q00, q10, q01 = Q(a, b), Q(a ^ 1, b), Q(a, b ^ 1)
        need(q00 ^ q10 == 1 ^ b and q00 ^ q01 == 1 ^ a,
             'polar identity from direct field arithmetic')
        for z in range(q):
            # Three direct outputs, then the stated decoder. Exhaustive,
            # including singular outputs; no random or generic-rank filter.
            w0 = root[squares[z ^ q00]]
            w1 = root[squares[z ^ q10]]
            w2 = root[squares[z ^ q01]]
            aa, bb = 1 ^ w0 ^ w2, 1 ^ w0 ^ w1
            need(aa == a and bb == b and (w0 ^ Q(aa, bb)) == z,
                 'three-query finite-field decoder')
            states += 1
    return {'field_size': q, 'polynomial_binary': hex(polynomial),
            'nontrivial_cube_roots': cube_roots, 'Frobenius_swaps_clock_roots': True,
            'complete_states_decoded': states, 'classes': q**3,
            'ring_of_same_size_classes': count(q)}


def lift_tests():
    tested = 0
    for k in range(4, 19):
        m = 2**k
        for z, a, b in product(range(4), repeat=3):
            step = m//4 if z % 2 else m//2
            baseline = signature(m, z, a, b)
            for dz, da, db in product(range(2), range(m//step), range(m//step)):
                need(signature(m, z+dz*(m//2), a+da*step, b+db*step) == baseline,
                     'every claimed lift in selected high-precision fibers')
                tested += 1
            need(signature(m, z+m//4, a, b) != baseline, 'visible source bit')
            need(signature(m, z, a+step//2, b) != baseline, 'visible phase bit')
            # k output bits require at most k+2 observation bits; k+1
            # fails for this odd-source witness when k >= 3.
        need(signature(2**(k+1), 1, 0, 0) ==
             signature(2**(k+1), 1, 2**(k-1), 0), 'one extra bit is insufficient')
        need(signature(2**(k+2), 1, 0, 0) !=
             signature(2**(k+2), 1, 2**(k-1), 0), 'two extra bits distinguish')
    # At modulus8 the even and odd source strata swap their fiber dimensions.
    need(signature(8, 0, 0, 0) == signature(8, 0, 2, 0), 'even k3 collapse')
    need(signature(8, 1, 0, 0) != signature(8, 1, 2, 0), 'odd k3 separation')
    need(signature(16, 0, 0, 0) != signature(16, 0, 4, 0), 'even k4 separation')
    need(signature(16, 1, 0, 0) == signature(16, 1, 4, 0), 'odd k4 collapse')
    return tested


def main():
    rows = [census(k) for k in range(7)]
    lifts = lift_tests()
    fields = [field_test(2, 0b111), field_test(6, 0x43), field_test(8, 0x11b)]
    crt = {m: count(m) for m in (5,10,11,12,13,25,50,76,100,228,256,768)}
    need(crt[50] == crt[100] == 125000, 'odd-part doubling plateau')
    need(crt[76] == 54872 and crt[228] == 384104, 'previous CRT controls')
    need(crt[256] == 1310720, 'ring256 theorem value')
    out = {'scope': 'modular determinant functions, no exact singularity flag',
           'gates': GATES, 'dyadic_censuses': rows, 'selected_lift_checks': lifts,
           'fields': fields, 'crt_theorem_values_not_full_censuses': crt,
           'dyadic_all_order_status': 'PROVED in DYADIC_AND_FIELD_PORTS.md'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
