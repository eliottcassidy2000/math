"""Exact 223/Fermat/Collatz certificates; stdlib only, no sampled probabilities.

Run: python 04-computation/experiments/crossroads223_20260926_geometry.py
Every finite universe and external theorem boundary is printed explicitly.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import sys


def primes_through(n):
    return [p for p in range(2, n + 1)
            if all(p % d for d in range(2, isqrt(p) + 1))]


def order(a, p):
    x = 1
    for k in range(1, p):
        x = x * a % p
        if x == 1:
            return k
    raise AssertionError("not a unit")


def fermat_via_image(p):
    h = {pow(x, 6, p) for x in range(1, p)}
    d = gcd(6, p - 1)
    pairs = sum((-1 - u) % p in h for u in h)
    torus = d * d * pairs
    boundary = 3 * d if p - 1 in h else 0
    return torus + boundary, torus, len(h), pairs


def fermat_direct(p):
    # Disjoint normalized charts: z=1, then z=0,y=1.
    sixth = [pow(x, 6, p) for x in range(p)]
    affine = sum((sixth[x] + sixth[y] + 1) % p == 0
                 for x in range(p) for y in range(p))
    torus = sum((sixth[x] + sixth[y] + 1) % p == 0
                for x in range(1, p) for y in range(1, p))
    infinity = sum((sixth[x] + 1) % p == 0 for x in range(p))
    return affine + infinity, torus


def compositions(total):
    if total == 0:
        yield ()
    for first in range(1, total + 1):
        for tail in compositions(total - first):
            yield (first,) + tail


def affine(word):
    a, b, c = 1, 1, 0
    for k in word:
        a, c, b = 3 * a, 3 * c + b, b * 2 ** k
    return a, b, c


def least_source(word, p=1):
    a, b, c = affine(word)
    # Requiring the final numerator to have valuation exactly S fixes all
    # intermediate valuations too, since all letters are positive.
    modulus = 2 * b
    residue = ((b - c) * pow(a, -1, modulus)) % modulus
    n = p * ((residue * pow(p, -1, modulus)) % modulus)
    return n if n else p * modulus


def orbit_word(n, word, shift=1):
    nodes = [n]
    for k in word:
        value = 3 * n + shift
        assert value % 2 ** k == 0 and value // 2 ** k % 2 == 1
        n = value // 2 ** k
        nodes.append(n)
    return nodes


def decode_carry(c, r):
    answer = []
    while r > 1:
        diff = c - 3 ** (r - 1)
        if diff <= 0:
            return None
        k = (diff & -diff).bit_length() - 1
        if k == 0:
            return None
        answer.append(k)
        c = diff // 2 ** k
        r -= 1
    return tuple(answer) if c == 1 else None


def gcd_free_rows(values):
    """Exact valuation row space without factoring; bases end pairwise coprime.

    Invariant: value[j] = product(base**row[j]). Splitting b,c at g=gcd(b,c)
    replaces their rows by rows at g,b/g,c/g; equal bases are merged.
    A final base may be composite, but its distinct prime valuation rows
    are nonzero scalar multiples of this row, so rational rank is exact.
    """
    bases = {}

    def add(base, row):
        if base > 1:
            bases.setdefault(base, Counter()).update(row)

    for j, value in enumerate(values):
        assert value >= 1
        add(value, Counter({j: 1}))
    while True:
        keys = list(bases)
        pair = None
        for i, b in enumerate(keys):
            for c in keys[i + 1:]:
                g = gcd(b, c)
                if g > 1:
                    pair = b, c, g
                    break
            if pair:
                break
        if pair is None:
            break
        b, c, g = pair
        rb, rc = bases.pop(b), bases.pop(c)
        add(g, rb + rc)
        add(b // g, rb)
        add(c // g, rc)
    for j, value in enumerate(values):
        product = 1
        for base, row in bases.items():
            product *= base ** row[j]
        assert product == value
    return list(bases.values())


def rank_mod(rows, width, prime=1000000007):
    pivots = {}
    for raw in rows:
        row = [x % prime for x in raw]
        for k, pivot in pivots.items():
            factor = row[k]
            if factor:
                row = [(x - factor * y) % prime for x, y in zip(row, pivot)]
        leading = next((j for j, x in enumerate(row) if x), None)
        if leading is not None:
            inverse = pow(row[leading], -1, prime)
            pivots[leading] = [(x * inverse) % prime for x in row]
            # Insertion order is pivot discovery order; each earlier pivot
            # remains zero in every subsequently inserted row.
            if len(pivots) == width:
                return width
    return len(pivots)


def paired_rank(sources):
    n = len(sources)
    values = [3 * m for m in sources] + [3 * m + 1 for m in sources]
    rows = gcd_free_rows(values)
    first = [[row[j] for j in range(n)] for row in rows]
    second = [[row[n + j] for j in range(n)] for row in rows]
    return rank_mod(first + second, n), rank_mod(first, n), len(rows)


def trial_factor_paired_rank(sources):
    """Independent small-control path via actual prime factorizations."""
    n = len(sources)
    rows = {}
    for j, m in enumerate(sources):
        for slot, value in enumerate([3 * m, 3 * m + 1]):
            p = 2
            while p * p <= value:
                exponent = 0
                while value % p == 0:
                    value //= p
                    exponent += 1
                if exponent:
                    rows.setdefault((slot, p), [0] * n)[j] = exponent
                p += 1
            if value > 1:
                rows.setdefault((slot, value), [0] * n)[j] = 1
    return rank_mod(rows.values(), n)


def run_rank_round():
    print("G2 HIGH-HEIGHT ROUND: gcd-free exact row space; rank modulo1000000007")
    print("Both coordinate slots retained; first-slot sign torsion is sum(c) mod2.")
    print("Full modular rank certifies rational independence regardless of torsion.")
    pset = [11, 43, 65, 253]
    nset = [13, 23, 121, 215]
    hostile = pset + nset
    assert paired_rank(hostile)[0] < len(hostile)
    assert paired_rank(hostile)[0] == trial_factor_paired_rank(hostile) == 7
    print("GENERIC_NONCHRONOLOGICAL_HOSTILE", hostile, paired_rank(hostile))
    assert paired_rank([1])[:2] == (1, 1)
    assert paired_rank([1, 1])[:2] == (1, 1)
    print("TRIVIAL_CYCLE_CONTROLS one_edge", paired_rank([1]),
          "two_laps", paired_rank([1, 1]))
    control = []
    n = 27
    while n != 1:
        control.append(n)
        value = 3 * n + 1
        n = value // (value & -value)
    assert paired_rank(control)[:2] == (41, 36)
    assert trial_factor_paired_rank(control) == 41
    print("INJECTIVE_ORBIT27_CONTROL_NODES", len(control) + 1,
          "TRANSITIONS", len(control), "RANKS", paired_rank(control))
    certificate_count = 0
    maximum_bits = 0
    # These are prescribed exact valuation cylinders, not more small starts.
    families = [
        ("critical223return", (2, 3, 1, 1), [4, 8, 16, 32]),
        ("expanding223return", (1, 1, 1, 1, 2, 1, 1, 1), [2, 4, 8, 16]),
        ("one_run", (1,), [32, 64, 128, 222]),
    ]
    for name, block, repeats in families:
        for count in repeats:
            word = block * count
            total = sum(word)
            least = least_source(word, 223)
            step = 223 * 2 ** (total + 1)
            for parameter in [0, 1, 223, 2 ** 64 + 223, 2 ** 256 + 223]:
                source = least + step * parameter
                nodes = orbit_word(source, word)
                assert len(nodes) == len(set(nodes))
                if name.endswith("223return"):
                    assert all(nodes[j] % 223 == 0
                               for j in range(0, len(nodes), len(block)))
                rank, first_rank, bases = paired_rank(nodes[:-1])
                assert rank == len(word), (name, count, parameter, rank)
                maximum_bits = max(maximum_bits, max(x.bit_length() for x in nodes))
                certificate_count += 1
                print("RANK", name, "r", len(word), "parameter", parameter,
                      "source_bits", source.bit_length(), "pair", rank,
                      "first", first_rank, "coprime_bases", bases)
    print("FULL_RANK_CYLINDERS", certificate_count, "MAX_NODE_BITS", maximum_bits)
    print("G2 remains OPEN; these cylinders do not exhaust possible trajectories.")


def main():
    print("SCRIPT_SHA256", sha256(Path(__file__).read_bytes()).hexdigest())
    print("FINITE UNIVERSE: all primes <=433, every nonzero sixth-power residue.")
    projective_empty = []
    torus_empty = []
    for p in primes_through(433):
        projective, torus, _, _ = fermat_via_image(p)
        if projective == 0:
            projective_empty.append(p)
        if torus == 0:
            torus_empty.append(p)
    assert projective_empty == [7, 31, 67, 79, 139, 223]
    assert torus_empty == [2, 5, 7, 13, 31, 61, 67, 79, 97, 139, 157, 223, 277]
    print("PROJECTIVE_EMPTY", projective_empty)
    print("TORUS_EMPTY", torus_empty)
    print("INDEPENDENT DIRECT CHART CONTROLS (p, projective points, torus points)")
    for p in [2, 3, 5, 7, 31, 67, 79, 139, 223, 277, 397, 401, 433, 439]:
        direct = fermat_direct(p)
        assert direct == fermat_via_image(p)[:2]
        print(p, *direct)
    print("CITED Hasse-Weil genus10 excludes projective emptiness for p>=401;")
    print("exact squared margin at401:", 402 ** 2 - 400 * 401)
    print("boundary <=18; excludes torus emptiness for p>=439;")
    print("exact squared margin at439:", (439 - 17) ** 2 - 400 * 439)
    assert all(not all(p % d for d in range(2, isqrt(p) + 1))
               for p in range(434, 439))

    p = 223
    h = {pow(2, k, p) for k in range(37)}
    assert h == {pow(x, 6, p) for x in range(1, p)}
    assert p - 1 not in h and p - 3 not in h
    print("ORDERS_223", {a: order(a, p) for a in [2, 3, 5, 7]})
    print("FACTOR_2^37-1", (2 ** 37 - 1) // p)
    assert (2 ** 37 - 1) % p == 0 and (2 ** 37 - 1) % p ** 2 != 0
    assert pow(3, 6, p) == pow(2, 21, p)
    print("H_223", sorted(h))
    classes = {pow(3, a, p) * x % p: a for a in range(6) for x in h}
    counts = [[0] * 6 for _ in range(6)]
    for n in range(1, p):
        y = (3 * n + 1) % p
        if y:
            counts[classes[n]][classes[y]] += 1
    print("SEXTIC_TRANSITION_COUNTS_223")
    for row in counts:
        print(row)
    zero_edges = [(a, b) for a in range(6) for b in range(6) if not counts[a][b]]
    assert zero_edges == [(5, 3)]
    print("ONLY_FORBIDDEN_NONZERO_EDGE", zero_edges)
    print("TWO_STEP_223_RETURN: impossible, since -3 notin H")
    pairs = [(a, b) for a in range(1, 38) for b in range(1, 38)
             if (9 + 3 * pow(2, a, p) + pow(2, a + b, p)) % p == 0]
    print("THREE_STEP_RETURN_INTERNAL_EXPONENTS_MOD37", sorted(pairs, key=lambda ab: (sum(ab), ab)))

    print("EXHAUSTIVE positive compositions total halvings S<=9:")
    by_sum = {}
    for total in range(1, 10):
        survivors = []
        for word in compositions(total):
            a, b, c = affine(word)
            if c % p == 0:
                n = least_source(word, p)
                nodes = orbit_word(n, word)
                assert nodes[-1] % p == 0
                survivors.append((word, a, b, c, n, nodes[-1]))
        by_sum[total] = survivors
        print("S", total, "return_words", survivors)
    assert not any(by_sum[s] for s in range(1, 7))
    assert [x[0] for x in by_sum[7]] == [(2, 3, 1, 1)]
    assert [x[0] for x in by_sum[8]] == [(2, 3, 1, 2)]
    assert all(last < first for total in [7, 8]
               for _, _, _, _, first, last in by_sum[total])
    print("FIRST EXPANDING RETURN CONTROLS")
    for word, a, b, c, n, last in by_sum[9]:
        if a > b:
            nodes = orbit_word(n, word)
            assert all(m % 223 for m in nodes[1:-1])
            print(word, nodes)
    word = (2, 3, 1, 1)
    assert affine(word) == (81, 128, 223)
    rat = [Fraction(223, 47)]
    for k in word:
        rat.append((3 * rat[-1] + 1) / 2 ** k)
    assert rat[0] == rat[-1]
    print("RATIONAL_PLUS_CYCLE", [str(x) for x in rat])
    print("INTEGER_3n+47_CYCLE", orbit_word(223, word, 47))
    print("CARRY223_DECODER_r1_to5", [(r, decode_carry(223, r)) for r in range(1, 6)])
    print("SCOPE: local arithmetic-geometric classification and bounded return gate;")
    print("neither orbit termination nor temporal equidistribution is proved.")
    if "--rank" in sys.argv:
        run_rank_round()


if __name__ == "__main__":
    main()
