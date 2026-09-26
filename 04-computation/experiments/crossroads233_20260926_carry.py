"""Exact depth-233 carry controls and elementary height/rank certificates.

Run: python 04-computation/experiments/crossroads233_20260926_carry.py
All proof checks remain enabled under python -O. No third-party libraries.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import factorial, gcd, log2, prod, sqrt


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def shortcut(n):
    return (3 * n + 1) // 2 if n % 2 else n // 2


def replay(n, length):
    word, nodes = [], []
    a, b, carry = 1, 1, 0
    for j in range(length):
        odd = n % 2
        word.append(odd)
        if odd:
            nodes.append(n)
            a, carry = 3 * a, 3 * carry + b
        b *= 2
        n = shortcut(n)
    return tuple(word), n, a, b, carry, nodes


def gcd_products(values):
    return [prod(gcd(a, b) for j, b in enumerate(values) if i != j)
            for i, a in enumerate(values)]


def private_cofactors(values):
    out = []
    for i, a in enumerate(values):
        r = a
        for j, b in enumerate(values):
            if i == j:
                continue
            while (g := gcd(r, b)) > 1:
                r //= g
        out.append(r)
    return out


def factors(n):
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def valuation_rank(values):
    fs = [factors(n) for n in values]
    primes = sorted(set().union(*(f.keys() for f in fs)))
    rows = [[Fraction(f.get(p, 0)) for f in fs] for p in primes]
    rank = 0
    for col in range(len(values)):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        v = rows[rank][col]
        rows[rank] = [x / v for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][col]:
                v = rows[i][col]
                rows[i] = [x - v * y for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def odd_nodes(n, count):
    nodes, sums = [n], [0]
    for _ in range(count - 1):
        x, k = 3 * nodes[-1] + 1, 0
        while x % 2 == 0:
            x //= 2
            k += 1
        nodes.append(x)
        sums.append(sums[-1] + k)
    return nodes, sums


def incident_carries(nodes, sums):
    count = len(nodes)
    qs = [1] * count
    for i in range(count):
        for j in range(i + 1, count):
            c = 2 ** (sums[j] - sums[i]) * nodes[j] - 3 ** (j - i) * nodes[i]
            independent = sum(3 ** (j - 1 - ell) * 2 ** (sums[ell] - sums[i])
                              for ell in range(i, j))
            check(c == independent and c > 0, "two carry formulas")
            check(c % gcd(nodes[i], nodes[j]) == 0, "pair gcd divides carry")
            qs[i] *= c
            qs[j] *= c
    return qs


def r_bound(count, i):
    e3 = (count - 1) * (count - 2) // 2 + i * (i - 1) // 2
    e2 = i * (count - 1) - i * (i + 1) // 2
    return Fraction(2 ** (count - 1) * factorial(i) * factorial(count - 1 - i)
                    * 3 ** e3, 2 ** e2)


def height_bound(count):
    return 3 ** (count - 2) * r_bound(count, count - 1)


def main():
    a, b = 3 ** 153, 2 ** 233
    u = pow(a, -1, 2 ** 80)
    high, low = 2 ** 153 * u - 1, 2 ** 153 * u - 5
    hub = (a * u - 1) // 2 ** 80
    tail = "01110111100111110111011110110011111111110010011110000100101111010011001101010101"
    words = [(1,) * 153 + (0,) * 80, tuple(map(int, "110" * 51 + tail))]
    carries = [a - 2 ** 153, 5 * a - 2 ** 153]
    for name, source, word, carry in zip(("N", "M"), (high, low), words, carries):
        actual, end, aa, bb, cc, nodes = replay(source, 233)
        check((actual, end, aa, bb, cc) == (word, hub, a, b, carry), "233 collision")
        odd = 0
        for j, digit in enumerate(word, 1):
            odd += digit
            check(3 ** odd > 2 ** j, "strict prefix slope")
        check(len(nodes) == 153 and all(x % 3 for x in nodes), "odd node types")
        private = private_cofactors(nodes)
        check(all(r > 1 for r in private), "153 independent private-prime columns")
        gp = gcd_products(nodes)
        print(f"collision {name}: source={source}; hub={end}; private certificates=153/153; "
              f"raw gcd dominance={sum(x > y for x,y in zip(nodes,gp))}/153")
        for lift in (1, 17, 2 ** 100):
            w, h, *_ = replay(source + lift * b, 233)
            check(w == word and h == hub + lift * a, "arbitrarily high collision lifts")

    block_words, realizations = 0, 0
    for t in range(1, 11):
        mod, bins = a ** t, Counter()
        for digits in product((0, 1), repeat=t):
            carry = 0
            for j, digit in enumerate(digits):
                carry = a * carry + carries[digit] * b ** j
            bins[carry % mod] += 1
            block_words += 1
            if t <= 6:
                denominator = b ** t
                source = (-carry * pow(a ** t, -1, denominator)) % denominator
                if source == 0:
                    source += denominator
                actual, end, aa, bb, cc, _ = replay(source, 233 * t)
                check(actual == sum((words[d] for d in digits), ()), "concatenation guards")
                check(aa == a ** t and bb == denominator and cc == carry, "block affine lift")
                realizations += 1
        check(len(bins) == 2 ** (t - 1) and set(bins.values()) == {2}, "suffix decoding")
    print(f"block digit test: {block_words} words through 10 blocks; exactly 2 per carry class; "
          f"{realizations} actual positive replays")

    tuples, certs, dependent = 0, 0, 0
    for length in (2, 3, 4):
        for values in combinations_with_replacement(range(2, 21), length):
            tuples += 1
            gp, rank = gcd_products(values), valuation_rank(values)
            certified = all(a > q for a, q in zip(values, gp))
            if certified:
                certs += 1
                check(rank == length, "gcd dominance implies independence")
            if rank < length:
                dependent += 1
                check(any(q % a == 0 for a, q in zip(values, gp)), "dependent divisibility witness")
    check(valuation_rank((6, 10, 15)) == 3, "independent without private primes")
    check(not any(r > 1 for r in private_cofactors((6, 10, 15))), "private-prime non-equivalence")
    check(valuation_rank((2, 3, 6)) == 2, "one relation is still a relation")
    check(valuation_rank((2, 4)) == 1, "gcd kernel positivity does not imply rank")
    print(f"gcd lemma universe: {tuples} tuples, {certs} sufficient certificates, "
          f"{dependent} dependent tuples all have divisibility witness")

    bound_checks = 0
    for count in range(2, 101):
        rs = [r_bound(count, i) for i in range(count)]
        check(max(rs) == rs[-1], "endpoint maximum")
        for i in range(count - 1):
            expected = Fraction((i + 1) * 3 ** i * 2 ** i,
                                (count - 1 - i) * 2 ** (count - 2))
            check(rs[i + 1] / rs[i] == expected, "ratio formula")
            bound_checks += 1
        for i in range(count - 2):
            check(rs[i + 1] / rs[i] < rs[i + 2] / rs[i + 1], "strict log convexity")
    for count, expected in ((2, 2), (3, 108), (4, 39366), (5, 86093442)):
        check(height_bound(count) == expected, "small H_N")
        print(f"H_{count}={expected}")
    print(f"incident bound formulas: {bound_checks} exact ratios through N=100")

    fixed_word_cases = 0
    for count in range(2, 7):
        for ks in product(range(1, 5), repeat=count - 1):
            sums, carry = [0], 0
            for k in ks:
                carry = 3 * carry + 2 ** sums[-1]
                sums.append(sums[-1] + k)
            modulus = 2 ** (sums[-1] + 1)
            residue = ((2 ** sums[-1] - carry) * pow(3 ** (count - 1), -1, modulus)) % modulus
            low_nodes, low_sums = odd_nodes(residue, count)
            check(low_sums == sums, "exact valuation cylinder")
            qs = incident_carries(low_nodes, sums)
            cutoff = max(Fraction(3 ** (count - 2) * q * 2 ** s, 3 ** i)
                         for i, (q, s) in enumerate(zip(qs, sums)))
            lower = cutoff.numerator // cutoff.denominator + 1
            source = residue + max(0, (lower - residue + modulus - 1) // modulus) * modulus
            nodes, actual_sums = odd_nodes(source, count)
            check(actual_sums == sums and source > cutoff, "fixed-word high specialization")
            check(all(3 * m > g for m, g in zip(nodes, gcd_products([3 * m for m in nodes]))),
                  "fixed-word rank cutoff")
            fixed_word_cases += 1
    print(f"fixed-word cutoff: {fixed_word_cases} exact valuation words, N=2..6, k_i=1..4")

    nodip, proved_cases, horizon_cases = 0, 0, 0
    for source in range(3, 100001, 2):
        nodes, sums = odd_nodes(source, 5)
        x, horizon_nodip = source, True
        for _ in range(source.bit_length() - 1):
            x = shortcut(x)
            if x < source:
                horizon_nodip = False
                break
        for count in (2, 3, 4, 5):
            ns, ss = nodes[:count], sums[:count]
            if source >= count and min(ns) >= source:
                nodip += 1
                qs = incident_carries(ns, ss)
                check(all(q < r_bound(count, i) for i, q in enumerate(qs)), "no-dip incident bound")
                if source > height_bound(count):
                    check(all(3 * m > g for m, g in zip(ns, gcd_products([3 * m for m in ns]))),
                          "pointwise no-dip rank certificate")
                    proved_cases += 1
            if count <= 4 and horizon_nodip:
                ell = (3 ** count).bit_length() + 1
                if source >= max(count, 2 ** ell):
                    check(sums[count] <= ell, "finite no-descent horizon sidecar")
                    horizon_cases += 1
    for count in range(2, 31):
        threshold = height_bound(count)
        bits = max(count + 1, threshold.numerator.bit_length() + 2)
        source = 2 ** bits - 1
        nodes, sums = odd_nodes(source, count)
        qs = incident_carries(nodes, sums)
        check(min(nodes) >= source > threshold, "high no-dip control")
        check(all(q < r_bound(count, i) for i, q in enumerate(qs)), "high incident bounds")
        check(all(3 * m > g for m, g in zip(nodes, gcd_products([3 * m for m in nodes]))),
              "high all-node gcd certificate")
        proved_cases += 1
    print(f"no-dip universe: odd sources 3..100000, N=2..5; {nodip} prefixes; "
          f"{proved_cases} above-threshold certificates including 29 high controls through N=30")
    print(f"no-descent horizon sidecar: {horizon_cases} exact guarded cases, N=2..4")
    for exponent in (2, 5, 20, 100):
        source = (4 ** exponent - 1) // 3
        nodes, _ = odd_nodes(source, 3)
        check(nodes[1:] == [1, 1], "height without no-dip hostile")
    print("hostiles: arbitrarily high immediate descent repeats 1; (6,10,15) defeats necessity; "
          "(2,3,6) defeats unique-relation/Hamiltonian inference; (2,4) defeats gcd-kernel rank")
    print(f"asymptotic all-source no-dip coefficient: {1/sqrt(log2(3)-0.5):.15f} (strictly smaller c required)")
    print("PASS: all checks exact except final displayed asymptotic decimal; Collatz and G2 remain open")


if __name__ == "__main__":
    main()
