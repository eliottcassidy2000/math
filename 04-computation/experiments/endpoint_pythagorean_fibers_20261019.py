#!/usr/bin/env python3
"""Exact endpoint-tournament / Gaussian-norm fiber bridge.

Standard library; no repository writes. All checks survive python -O.
Outputs a deterministic JSON audit. The proofs are in the companion note.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd, isqrt, prod
import json

CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def integer(n, low, name):
    if type(n) is not int or n < low:
        raise ValueError(name)
    return n


def arc(n, a, b):
    return (a < b and (a, b) != (0, n - 1)) or (a, b) == (n - 1, 0)


def path_from_mask(r, mask):
    integer(r, 1, "positive number of interior coordinates")
    integer(mask, 0, "integer mask")
    if mask >= 1 << r:
        raise ValueError("mask outside coordinate range")
    return tuple(j + 1 for j in range(r) if mask >> j & 1) + (r + 1, 0) + tuple(
        j + 1 for j in range(r) if not (mask >> j & 1)
    )


def mask_from_path(r, path):
    integer(r, 1, "positive number of interior coordinates")
    n = r + 2
    if not isinstance(path, tuple) or any(type(v) is not int for v in path):
        raise ValueError("path must be an integer tuple")
    if len(path) != n or set(path) != set(range(n)):
        raise ValueError("path must use every vertex once")
    if not all(arc(n, a, b) for a, b in zip(path, path[1:])):
        raise ValueError("not a directed Hamiltonian path of E_n")
    if path == tuple(range(n)):
        return None
    k = path.index(n - 1)
    if k + 1 >= n or path[k + 1] != 0:
        raise ValueError("nonexceptional path must use the reversed edge")
    mask = sum(1 << (v - 1) for v in path[:k])
    if path_from_mask(r, mask) != path:
        raise ValueError("incorrect increasing-block order")
    return mask


def mul(z, w):
    x, y = z
    a, b = w
    return x * a - y * b, x * b + y * a


def power(z, e):
    w = (1, 0)
    for _ in range(e):
        w = mul(w, z)
    return w


def conjugate(z):
    return z[0], -z[1]


def quotient(z, w):
    num = mul(z, conjugate(w))
    den = w[0] * w[0] + w[1] * w[1]
    return (num[0] // den, num[1] // den) if num[0] % den == num[1] % den == 0 else None


def trial_factor(n):
    ans = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                e += 1
                n //= p
            ans.append((p, e))
        p = 3 if p == 2 else p + 2
    if n > 1:
        ans.append((n, 1))
    return ans


@lru_cache(maxsize=None)
def gaussian_prime(p):
    if p % 4 != 1 or trial_factor(p) != [(p, 1)]:
        raise ValueError("requires a prime congruent to one modulo four")
    for b in range(1, isqrt(p // 2) + 1):
        a = isqrt(p - b * b)
        if a * a + b * b == p and a > b:
            return a, b
    raise RuntimeError("split prime has no selected Gaussian factor")


def data_from_factors(factors):
    if not factors or any(p % 4 != 1 for p, _ in factors):
        raise ValueError("nondegenerate primitive hypotenuse needs only split primes")
    return tuple((p, e, gaussian_prime(p)) for p, e in factors)


def fiber_data(c):
    integer(c, 2, "hypotenuse greater than one")
    return data_from_factors(trial_factor(c))


def alpha_from_mask(data, mask):
    r = len(data)
    integer(mask, 0, "integer mask")
    if mask >= 1 << r:
        raise ValueError("mask outside coordinate range")
    z = (1, 0)
    for j, (_, e, pi) in enumerate(data):
        z = mul(z, power(pi if mask >> j & 1 else conjugate(pi), e))
    return z


def ppt_from_alpha(z):
    n, m = sorted(map(abs, z), reverse=True)
    if m == 0 or n == m or gcd(n, m) != 1 or (n + m) % 2 != 1:
        raise ValueError("requires primitive, odd, nonaxis Gaussian norm")
    return n * n - m * m, 2 * n * m, n * n + m * m


def canonical_mask_from_ppt(data, ppt):
    if not isinstance(ppt, tuple) or len(ppt) != 3 or any(type(x) is not int or x <= 0 for x in ppt):
        raise ValueError("positive integer triple required")
    u, v, c = ppt
    if u * u + v * v != c * c or gcd(u, v) != 1 or (u + v) % 2 != 1:
        raise ValueError("primitive Pythagorean triple required")
    if c != prod(p ** e for p, e, _ in data):
        raise ValueError("wrong marked norm fiber")
    a = u if u % 2 else v
    n2, m2 = (c + a) // 2, (c - a) // 2
    n, m = isqrt(n2), isqrt(m2)
    if n * n != n2 or m * m != m2 or not n > m > 0:
        raise ValueError("Euclid inverse failed")
    z = (n, m)
    mask = 0
    for j, (_, e, pi) in enumerate(data):
        valuations = []
        for q in (pi, conjugate(pi)):
            w, value = z, 0
            while True:
                w = quotient(w, q)
                if w is None:
                    break
                value += 1
            valuations.append(value)
        if tuple(valuations) == (e, 0):
            mask |= 1 << j
        elif tuple(valuations) != (0, e):
            raise ValueError("nonprimitive or wrong Gaussian allocation")
    return mask if mask & 1 else mask ^ ((1 << len(data)) - 1)


def spf_sieve(limit):
    spf = list(range(limit + 1))
    for p in range(2, isqrt(limit) + 1):
        if spf[p] == p:
            for n in range(p * p, limit + 1, p):
                if spf[n] == n:
                    spf[n] = p
    return spf


def sieve_factor(n, spf):
    result = []
    while n > 1:
        p, e = spf[n], 0
        while n % p == 0:
            n //= p
            e += 1
        result.append((p, e))
    return result


def euclid_fibers(limit):
    ans = defaultdict(set)
    for n in range(2, isqrt(limit) + 1):
        for m in range(1, min(n, isqrt(limit - n * n) + 1)):
            if gcd(n, m) == 1 and (n + m) % 2:
                ans[n * n + m * m].add((n * n - m * m, 2 * n * m, n * n + m * m))
    return ans


def hp_count_dp(n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v]:
                for w in range(n):
                    if not (mask >> w & 1) and arc(n, v, w):
                        dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[-1])


def rejected(func, *args):
    try:
        func(*args)
    except ValueError:
        check(True, "expected rejection")
    else:
        check(False, "invalid input accepted")


def run():
    path_counts = {}
    for n in range(3, 11):
        r = n - 2
        expected = {path_from_mask(r, m) for m in range(1 << r)} | {tuple(range(n))}
        for m in range(1 << r):
            p = path_from_mask(r, m)
            check(all(arc(n, a, b) for a, b in zip(p, p[1:])), "address gives a path")
            check(mask_from_path(r, p) == m, "path/address inverse")
            check(path_from_mask(r, m ^ ((1 << r) - 1)) != p, "free address complement")
        check(mask_from_path(r, tuple(range(n))) is None, "exceptional increasing path")
        count = hp_count_dp(n)
        check(count == 1 + (1 << r) == len(expected), "DP count")
        path_counts[n] = count
        if n <= 8:
            actual = {p for p in permutations(range(n)) if all(arc(n, a, b) for a, b in zip(p, p[1:]))}
            check(actual == expected, "literal exhaustive path census")

    automorphism_counts = {}
    for n in range(3, 8):
        automorphisms = [p for p in permutations(range(n)) if all(
            arc(n, a, b) == arc(n, p[a], p[b]) for a in range(n) for b in range(n) if a != b)]
        check(len(automorphisms) == (3 if n == 3 else 1), "automorphism count")
        r = n - 2
        for p in automorphisms:
            acts_as_complement = all(tuple(p[v] for v in path_from_mask(r, m)) == path_from_mask(r, m ^ ((1 << r) - 1)) for m in range(1 << r))
            check(not acts_as_complement, "no automorphism induces address complement")
        automorphism_counts[n] = len(automorphisms)

    limit = 100000
    spf = spf_sieve(limit)
    independent = euclid_fibers(limit)
    by_support = Counter()
    eligible_norms = 0
    raw_addresses = 0
    units = ((1, 0), (-1, 0), (0, 1), (0, -1))
    for c in range(1, limit + 1):
        factors = sieve_factor(c, spf)
        if not factors or any(p % 4 != 1 for p, e in factors):
            check(not independent[c], "excluded norm has no primitive triple")
            continue
        data = data_from_factors(factors)
        r = len(data)
        eligible_norms += 1
        raw_addresses += 1 << r
        fibers = defaultdict(list)
        for mask in range(1 << r):
            z = alpha_from_mask(data, mask)
            check(z[0] * z[0] + z[1] * z[1] == c, "Gaussian norm")
            check(alpha_from_mask(data, mask ^ ((1 << r) - 1)) == conjugate(z), "complement is conjugation")
            triple = ppt_from_alpha(z)
            check(triple[0] % 2 == 1 and triple[1] % 4 == 0, "leg parity")
            check(triple in independent[c], "Gaussian-to-Euclid census")
            canonical = canonical_mask_from_ppt(data, triple)
            expected = mask if mask & 1 else mask ^ ((1 << r) - 1)
            check(canonical == expected, "full inverse to complement class")
            check(canonical_mask_from_ppt(data, (triple[1], triple[0], c)) == canonical, "size/parity leg order")
            for unit in units:
                check(ppt_from_alpha(mul(unit, z)) == triple, "unit quotient")
                check(ppt_from_alpha(mul(unit, conjugate(z))) == triple, "conjugate-unit quotient")
            fibers[triple].append(mask)
        check(set(fibers) == independent[c], "all and only primitive triples in fiber")
        check(len(fibers) == 1 << (r - 1), "fiber quotient count")
        for masks in fibers.values():
            check(len(masks) == 2 and masks[0] ^ masks[1] == (1 << r) - 1, "exact complement pairs")
        by_support[r] += len(fibers)

    angular_norms = [c for c in range(2, 1001) if independent[c]] + [1105, 32045, 1185665]
    for c in angular_norms:
        data = fiber_data(c)
        r = len(data)
        for k in range(4):
            phases = []
            for p, e, pi in data:
                numerator = power(pi, 8 * k * e)
                denominator = p ** (4 * k * e)
                phases.append(tuple(Fraction(v, denominator) for v in numerator))
            weighted_sum = (Fraction(0), Fraction(0))
            for mask in range(1 << r):
                weight = (Fraction(1), Fraction(0))
                for j, phase in enumerate(phases):
                    weight = mul(weight, phase if mask >> j & 1 else conjugate(phase))
                weighted_sum = tuple(a + b for a, b in zip(weighted_sum, weight))
            factored = prod(2 * z[0] for z in phases)
            check(weighted_sum == (factored, 0), "Laurent path sum factors exactly")
            angular_sum = Fraction(0)
            for mask in range(1, 1 << r, 2):
                a, b, norm = ppt_from_alpha(alpha_from_mask(data, mask))
                angular_sum += Fraction(power((a, b), 4 * k)[0], c ** (4 * k))
            check(weighted_sum[0] / 2 == angular_sum, "weighted path orbit/PPT angular identity")

    prime_power_counts = {}
    for p in (5, 13, 17, 29, 41):
        z, c, accepted = (1, 0), 1, 0
        for e in range(1, 1001):
            z, c = mul(z, gaussian_prime(p)), c * p
            a, b, norm = ppt_from_alpha(z)
            check(norm == c, "exact prime-power norm and primitive decoder")
            check(c * c != 4 * a * b, "no angular boundary exponent")
            accepted += c * c > 4 * a * b
        check(0 < accepted < 1000, "both prime-power outcomes in finite sample")
        prime_power_counts[p] = {"exponents": 1000, "positive_ap": accepted}

    # General formal z-coefficients, constructed by independent finite Euler multiplication.
    coefficient_limit = 10000
    formal = {1: {0: 1}}
    for p in range(5, coefficient_limit + 1, 4):
        if spf[p] != p:
            continue
        old = list(formal.items())
        for base, poly in old:
            q = p
            while base * q <= coefficient_limit:
                target = formal.setdefault(base * q, {})
                for degree, coefficient in poly.items():
                    target[degree + 1] = target.get(degree + 1, 0) + coefficient
                q *= p
    for c in range(1, coefficient_limit + 1):
        factors = sieve_factor(c, spf)
        expected = {len(factors): 1} if all(p % 4 == 1 for p, _ in factors) else {}
        check(formal.get(c, {}) == expected, "formal coefficient is z^omega")
        specialized = sum(value * 2 ** degree for degree, value in formal.get(c, {}).items())
        check(specialized == (1 if c == 1 else 0) + 2 * len(independent[c]), "F(2) coefficient")

    # Exact elementary symmetric sums for four specified primes at s=2.
    test_primes = (5, 13, 17, 29)
    xs = [Fraction(1, p * p - 1) for p in test_primes]
    es = [Fraction(1)] + [Fraction(0)] * len(xs)
    for x in xs:
        for degree in range(len(xs), 0, -1):
            es[degree] += x * es[degree - 1]
    for degree in range(5):
        direct = sum((prod(subset) for subset in combinations(xs, degree)), Fraction(0))
        check(es[degree] == direct, "elementary symmetric coefficient")
    check(sum(2 ** r * es[r] for r in range(5)) == prod(1 + 2 * x for x in xs), "F(2) elementary symmetric specialization")

    # Exponents are independent of support size; enumerate all four-prime exponent vectors 0..3.
    truncated_xs = [sum((Fraction(1, p ** (2 * e)) for e in range(1, 4)), Fraction(0)) for p in test_primes]
    direct_by_support = Counter()
    for exponents in product(range(4), repeat=4):
        c = prod(p ** e for p, e in zip(test_primes, exponents))
        r = sum(e > 0 for e in exponents)
        direct_by_support[r] += Fraction(1, c * c)
    for r in range(5):
        check(direct_by_support[r] == sum((prod(subset) for subset in combinations(truncated_xs, r)), Fraction(0)), "bounded exponents/support refinement")

    # Primitive Gaussian-lattice coefficients test the exact Eisenstein normalization.
    lattice_limit = 1000
    primitive_lattice = Counter()
    for x in range(-isqrt(lattice_limit), isqrt(lattice_limit) + 1):
        for y in range(-isqrt(lattice_limit), isqrt(lattice_limit) + 1):
            c = x * x + y * y
            if 0 < c <= lattice_limit and gcd(x, y) == 1:
                primitive_lattice[c] += 1
    for c in range(1, lattice_limit + 1):
        f_c = (1 if c == 1 else 0) + 2 * len(independent[c])
        f_half = ((1 if c == 2 else 0) + 2 * len(independent[c // 2])) if c % 2 == 0 else 0
        check(primitive_lattice[c] == 4 * (f_c + f_half), "E(i,s)=2(1+2^-s)F(2,s)")

    examples = []
    for c in (5, 25, 41, 65, 1105, 3125, 32045, 1185665):
        data = fiber_data(c)
        r = len(data)
        rows = []
        for mask in range(1, 1 << r, 2):
            z = alpha_from_mask(data, mask)
            a, b, c = ppt_from_alpha(z)
            comp = mask ^ ((1 << r) - 1)
            check(canonical_mask_from_ppt(data, (a, b, c)) == mask, "displayed example inverse")
            rows.append({"S": [j + 1 for j in range(r) if mask >> j & 1], "alpha": z,
                         "path": path_from_mask(r, mask), "complement_path": path_from_mask(r, comp),
                         "parity_ordered_ppt": (a, b, c), "positive_ap": c * c > 4 * a * b,
                         "signed_ap": tuple(2 * t for t in (c * c - 4 * a * b, c * c, c * c + 4 * a * b))})
        check(len(rows) == 1 << (r - 1), "example support quotient count")
        examples.append({"C": c, "r": r, "rows": rows})

    for bad in (0, 1, 2, 3, 10, 15, 45, 75, True, 5.0):
        rejected(fiber_data, bad)
    for args in ((0, 0), (True, 0), (1, -1), (2, 4), (2, True), (2, 1.0)):
        rejected(path_from_mask, *args)
    for path in ((0, 2, 1), (0, 1, 1), (0, 1), (0, 1, True), [0, 1, 2]):
        rejected(mask_from_path, 1, path)
    rejected(canonical_mask_from_ppt, fiber_data(5), (6, 8, 10))
    rejected(canonical_mask_from_ppt, fiber_data(65), (3, 4, 5))
    rejected(ppt_from_alpha, (1, 0))
    rejected(ppt_from_alpha, (1, 1))
    rejected(ppt_from_alpha, (4, 2))

    return {"status": "PASS", "checks": CHECKS,
            "hypotenuse_limit": limit, "primitive_triples": sum(map(len, independent.values())),
            "eligible_norms": eligible_norms, "raw_addresses": raw_addresses,
            "primitive_triples_by_support_size": dict(sorted(by_support.items())),
            "prime_power_counts": prime_power_counts,
            "weighted_angular_fiber_count": len(angular_norms), "weighted_angular_modes": list(range(4)),
            "hamiltonian_path_counts": path_counts, "automorphism_counts": automorphism_counts,
            "formal_coefficient_limit": coefficient_limit, "primitive_gaussian_norm_limit": lattice_limit,
            "elementary_symmetric_test_primes": test_primes, "examples": examples}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
