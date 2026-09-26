"""Exact common-tail, dyadic-cylinder, and persistent-prime-motif families.

Run with python or python -O. Standard library only, no disabled asserts.
"""
from fractions import Fraction
from itertools import product
from math import gcd


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def half(n):
    return (3 * n + 1) // 2 if n % 2 else n // 2


def odd_step(n):
    x, k = 3 * n + 1, 0
    while x % 2 == 0:
        x //= 2
        k += 1
    return x, k


def odd_trace(n, count):
    nodes, sums, word = [n], [0], []
    for _ in range(count - 1):
        n, k = odd_step(n)
        nodes.append(n)
        sums.append(sums[-1] + k)
        word.extend([1] + [0] * (k - 1))
    return nodes, sums, tuple(word)


def affine(word):
    a, b, c, e = 1, 1, 0, 0
    for bit in word:
        if bit:
            a, c, e = 3 * a, 3 * c + b, e + 1
        b *= 2
    return a, b, c, e


def replay(n, word):
    start, a, b, c, e = n, 1, 1, 0, 0
    for bit in word:
        check(n % 2 == bit, "actual parity guard")
        if bit:
            a, c, e = 3 * a, 3 * c + b, e + 1
        b *= 2
        n = half(n)
    check(b * n == a * start + c, "affine replay")
    return n


def length_to_one(n, classical=False):
    count = 0
    while n != 1:
        n = (3 * n + 1 if n % 2 else n // 2) if classical else half(n)
        count += 1
        check(count < 1000, "finite length control")
    return count


def motif_check(nodes, support):
    for p, expected in support.items():
        actual = tuple(i for i, n in enumerate(nodes) if n % p == 0)
        check(actual == expected, f"exclusive prime support {p}")


def first_u(a, modulus, k):
    residue = ((a + 1) * pow(pow(3, k, modulus), -1, modulus)) % modulus
    return residue + max(0, (a + 1 - residue + modulus - 1) // modulus) * modulus


def family_count(a, modulus, k, x):
    least = 2 ** k * first_u(a, modulus, k) - 1
    return max(0, (x - least) // (2 ** k * modulus) + 1)


def union_count(a, modulus, minimum, x):
    return sum(family_count(a, modulus, k, x)
               for k in range(minimum, (x + 1).bit_length()))


def main():
    check(length_to_one(27) == 70 and length_to_one(27, True) == 111, "27 time convention")
    comb = [(82 * 4 ** r - 1) // 3 for r in range(12)]
    for r, n in enumerate(comb):
        check(odd_step(n) == (41, 2 * r + 1), "41 common-tail comb")
        check(length_to_one(n) == 70 + 2 * r, "shortcut total stopping time")
        check(length_to_one(n, True) == 111 + 2 * r, "classical total stopping time")
        check(n % 3 == r % 3, "one-third dead inverse branches")
        if r:
            check(half(half(n)) < n, "common tail is not long first-descent time")
        if r + 1 < len(comb):
            check(comb[r + 1] == 4 * n + 1, "comb recursion")
    print(f"direct odd predecessors of41: {comb}; total shortcut lengths={[70+2*r for r in range(12)]}")
    for x in (26, 27, 1000, 1000000):
        explicit = [27 * 2 ** j for j in range(x.bit_length() + 1) if 27 * 2 ** j <= x]
        expected = max(0, (x // 27).bit_length())
        check(len(explicit) == expected, "exact hits27 count")
    check(all((2 ** k * 27 - 1) % 3 for k in range(1, 101)), "27 has no odd predecessors")
    print("hits27 exactly27*2^r; odd predecessors absent; count floor(log2(X/27))+1 forX>=27")

    # Construct a positive basin41 point in every parity cylinder through depth10.
    logarithms = {}
    for e in range(1, 11):
        mod, period, x, table = 3 ** e, 2 * 3 ** (e - 1), 1, {}
        for r in range(period):
            check(x not in table, "primitive order modulo3^e")
            table[x] = r
            x = (2 * x) % mod
        check(x == 1 and len(table) == period, "unit group enumeration")
        logarithms[e] = table
    dense_cases = 0
    for length in range(11):
        for word in product((0, 1), repeat=length):
            a, b, c, e = affine(word)
            if e:
                r = logarithms[e][c * pow(b * 41, -1, a) % a]
                period = 2 * 3 ** (e - 1)
            else:
                r, period = 0, 1
            while b * 2 ** r * 41 <= c:
                r += period
            source = (b * 2 ** r * 41 - c) // a
            end = replay(source, word)
            check(end == 2 ** r * 41 and end >> r == 41, "dense positive basin construction")
            dense_cases += 1
    print(f"basin41 dyadic density construction: {dense_cases} parity cylinders, every length0..10")

    fixed_depth = []
    for depth in range(1, 5):
        count = 0
        for source in range(1, 10001, 2):
            n = source
            for _ in range(depth):
                n = odd_step(n)[0]
            count += n == 41
        cap = (3 ** depth * 10001).bit_length() ** depth
        check(count <= cap, "fixed odd-depth polylog bound")
        fixed_depth.append((depth, count, cap))
    print(f"fixed odd-depth U^d(n)=41, oddn<=10000, (d,count,proven cap): {fixed_depth}")

    bases = [
        (27, 36, (8, 13, 35), 455, {5: (1, 2), 7: (0, 1), 13: (0, 2)}, 3 * 2 ** 49),
        (4347, 10, tuple(range(10)), 4669, {7: (0, 6), 23: (0, 3), 29: (3, 6)}, 59136),
    ]
    motif_cases = 0
    for a, node_count, selected, odd_modulus, support, period in bases:
        nodes, sums, word = odd_trace(a, node_count)
        length, ell = len(word), len(word) + 1
        modulus = odd_modulus * 2 ** ell
        motif_check([nodes[i] for i in selected], support)
        e, minimum = 0, Fraction(100)
        for j, bit in enumerate(word, 1):
            e += bit
            slope = Fraction(3 ** e, 2 ** j)
            minimum = min(minimum, slope)
            check(slope > 1, "all base prefix slopes expanding")
        check(minimum == Fraction(9, 8), "exact minimum base slope")
        effective = modulus // gcd(a + 1, modulus)
        check(pow(3, period, effective) == 1, "residue clock period")
        divisors = (2, 3) if a == 27 else (2, 3, 7, 11)
        for p in divisors:
            check(pow(3, period // p, effective) != 1, "residue clock minimality")
        for k in range(31):
            u0 = first_u(a, modulus, k)
            check((u0 & -u0) == 4, "disjoint shell v2(u)=2")
            check(first_u(a, modulus, k + period) == u0, "exact phase recurrence")
            for t in (0, 1, 7, 2 ** 20):
                u = u0 + modulus * t
                source = 2 ** k * u - 1
                landing = replay(source, (1,) * k)
                check(landing == 3 ** k * u - 1 and landing >= a and (landing - a) % modulus == 0,
                      "prepended family landing")
                lifted, lifted_sums, lifted_word = odd_trace(landing, node_count)
                check(lifted_word == word and lifted_sums == sums, "same exact valuation word")
                motif_check([lifted[i] for i in selected], support)
                x = source
                for bit in (1,) * k + word:
                    check(x % 2 == bit, "actual prepended parity")
                    x = half(x)
                    check(x > source, "actual uniform no-dip prefix")
                motif_cases += 1
        print(f"base{a}: L={length}, ell={ell}, M={modulus}, min slope={minimum}, residue period={period}")
        for kmin in (0, 1, 5, 20):
            density = Fraction(2, modulus * 2 ** kmin)
            for x in (2 ** 40, 2 ** 80, 2 ** 120):
                count = union_count(a, modulus, kmin, x)
                # One O(1) floor error per active scale and a bounded geometric tail.
                error = abs(Fraction(count) - density * x)
                bound = (2 + Fraction(a + 1, modulus)) * (x + 1).bit_length() + 3
                check(error <= bound, "uniform logarithmic discrepancy control")
            print(f"  tailk>={kmin}: natural density={density}; Z2 closure Haar={Fraction(2,2**(ell+kmin))}")
    print(f"persistent motifs and no-dip: {motif_cases} exact lifted family realizations")

    # The longer 27 prefix DOES descend: local motif persistence is not endless no-dip.
    nodes, sums, word = odd_trace(27, 38)
    a, b, c, _ = affine(word)
    check(nodes[-1] == 23 and sums[-1] == 59 and a < b, "longer contracting 27 prefix")
    cutoff = Fraction(c, b - a)
    for t in (0, 1, 17):
        source = 27 + 2 ** 60 * t
        end = replay(source, word)
        check(source > cutoff and end < source, "longer cylinder descends")
    print(f"59-step27 cylinder: slope={Fraction(a,b)}; descent cutoff={cutoff}; source27 and all positive lifts descend")
    print("PASS: exact finite families/frequencies; natural density differs from Z2 closure Haar measure; no convergence inference")


if __name__ == "__main__":
    main()
