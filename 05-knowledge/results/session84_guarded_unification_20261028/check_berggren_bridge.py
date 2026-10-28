#!/usr/bin/env python3
"""Exact, dependency-free checks for the guarded Berggren/semiprime bridge."""
from bisect import bisect_right
from hashlib import sha256
import argparse
import json
from math import gcd, isqrt


def require(test, message):
    if not test:
        raise RuntimeError(message)


def sieve(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p * p:limit + 1:p] = b"\x00" * ((limit - p * p) // p + 1)
    return flags, [p for p in range(2, limit + 1) if flags[p]]


def prime_trial(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def sigma_direct(n):
    result = 0
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            result += d
            if d * d != n:
                result += n // d
    return result


def phi(m, n):
    return m * m - n * n, 2 * m * n, m * m + n * n


def child(m, n, branch):
    if branch == "U":
        return 2 * m - n, m
    if branch == "A":
        return 2 * m + n, m
    if branch == "D":
        return m + 2 * n, n
    raise ValueError(branch)


def inverse_word(m, n):
    result = []
    while (m, n) != (2, 1):
        old_m = m
        if m > 3 * n:
            m, n = m - 2 * n, n
            result.append("D")
        elif m > 2 * n:
            m, n = n, m - 2 * n
            result.append("A")
        else:
            m, n = n, 2 * n - m
            result.append("U")
        require(m > n > 0 and m < old_m, "inverse descent")
    return "".join(reversed(result))


def odd_pair(m, n, scale, gap_leg):
    a, b, c = phi(m, n)
    r = a if gap_leg == "a" else b
    return scale * (c - r) - 1, scale * (c + r) - 1


def pair_children(p, q):
    z = isqrt((p + 1) * (q + 1))
    return [((q, p + 4 * q + 4 - 4 * z), 2 * (q + 1) - z),
            ((q, p + 4 * q + 4 + 4 * z), 2 * (q + 1) + z),
            ((p, 4 * p + q + 4 + 4 * z), 2 * (p + 1) + z)]


def encode_odd(p, q):
    r, z, c = (q - p) // 2, isqrt((p + 1) * (q + 1)), (p + q + 2) // 2
    scale = gcd(gcd(r, z), c)
    r0, z0, c0 = r // scale, z // scale, c // scale
    gap_leg = "a" if r0 % 2 else "b"
    a, b = (r0, z0) if gap_leg == "a" else (z0, r0)
    m, n = isqrt((c0 + a) // 2), isqrt((c0 - a) // 2)
    require(phi(m, n) == (a, b, c0), "inverse Euclid lift")
    require(odd_pair(m, n, scale, gap_leg) == (p, q), "odd inverse")
    require(gap_leg == "a" or scale % 2 == 0, "even-gap parity")
    return inverse_word(m, n), m, n, scale, gap_leg


def main(limit):
    require(limit >= 1000, "use a source bound at least 1000")
    flags, primes = sieve(limit)
    direct = {}
    tested_pairs = 0
    for i, p in enumerate(primes):
        if p * p > limit:
            break
        for q in primes[i + 1:bisect_right(primes, limit // p)]:
            tested_pairs += 1
            total = (p + 1) * (q + 1)
            z = isqrt(total)
            if z * z == total:
                require(prime_trial(p) and prime_trial(q), "independent prime oracle")
                require(sigma_direct(p * q) == total, "independent divisor sum")
                require(p != 3 and q != 3, "three exclusion")
                direct[p * q] = (p, q, p + q + 1, z)

    for p, q, _, _ in direct.values():
        for (pp, qq), zz in pair_children(p, q):
            require(pp < qq and zz > 0 and (pp + 1) * (qq + 1) == zz * zz,
                    "universal square-product child")
        if p == 2:
            first, middle, last = pair_children(p, q)
            require(first[0][1] > 2 and first[0][1] % 2 == 0, "even U exclusion")
            require(middle[0][1] > 2 and middle[0][1] % 2 == 0, "even A exclusion")
            kk = isqrt((q + 1) // 12)
            require(last[0] == (2, 12 * (kk + 1) ** 2 - 1), "even D advance")
        else:
            _, m, n, scale, gap_leg = encode_odd(p, q)
            order = "UAD" if gap_leg == "a" else "DAU"
            for ((pp, qq), _), branch in zip(pair_children(p, q), order):
                mm, nn = child(m, n, branch)
                require(odd_pair(mm, nn, scale, gap_leg) == (pp, qq), "child chart intertwining")

    # All odd admissible pairs have p>=5, hence their centered c is at most
    # floor((floor(limit/5)+7)/2); p+limit/p decreases on [5,sqrt(limit)].
    cmax = (limit // 5 + 7) // 2
    tree = {}
    stack = [(2, 1, "")]
    while stack:
        m, n, word = stack.pop()
        if m * m + n * n > cmax:
            continue
        require((m, n) not in tree, "duplicate Berggren address")
        require(inverse_word(m, n) == word, "word roundtrip")
        tree[m, n] = word
        for branch in "UAD":
            mm, nn = child(m, n, branch)
            stack.append((mm, nn, word + branch))
    euclid = {(m, n) for m in range(2, isqrt(cmax) + 1) for n in range(1, m)
              if m * m + n * n <= cmax and gcd(m, n) == 1 and (m + n) % 2 == 1}
    require(set(tree) == euclid, "independent Euclid/tree census")

    generated = {}
    chart_counts = {"a": 0, "b": 0, "even_D": 0}
    for (m, n), word in tree.items():
        for gap_leg in "ab":
            scale = 1 if gap_leg == "a" else 2
            step = 1 if gap_leg == "a" else 2
            while True:
                p, q = odd_pair(m, n, scale, gap_leg)
                source = p * q
                if source > limit:
                    break
                if p >= 2 and flags[p] and flags[q]:
                    require(source not in generated, "noninjective decorated map")
                    require(encode_odd(p, q) == (word, m, n, scale, gap_leg), "decorated roundtrip")
                    generated[source] = (p, q, p + q + 1, isqrt((p + 1) * (q + 1)))
                    chart_counts[gap_leg] += 1
                scale += step

    k = 1
    while 2 * (12 * k * k - 1) <= limit:
        q = 12 * k * k - 1
        if flags[q]:
            a, b, c = phi(2 * k, 1)
            require((q - 2, 12 * k, q + 4) == (3 * a, 3 * b, 3 * c), "three-times D ray")
            require(inverse_word(2 * k, 1) == "D" * (k - 1), "D word")
            generated[2 * q] = (2, q, q + 3, 6 * k)
            chart_counts["even_D"] += 1
        k += 1
    require(generated == direct, "complete prime-pair/decorated-tree comparison")

    # Six exact branch failures, in both chart orientations.
    controls = [
        ("a", "U", 4, 3, 1, (17, 31), (31, 49)),
        ("a", "A", 2, 1, 6, (11, 47), (47, 299)),
        ("a", "D", 2, 1, 3, (5, 23), (5, 95)),
        ("b", "U", 2, 1, 12, (11, 107), (11, 299)),
        ("b", "A", 2, 1, 8, (7, 71), (71, 391)),
        ("b", "D", 2, 1, 12, (11, 107), (107, 299)),
    ]
    for leg, branch, m, n, scale, before, after in controls:
        require(odd_pair(m, n, scale, leg) == before, "control source")
        require(all(map(prime_trial, before)), "control prime source")
        mm, nn = child(m, n, branch)
        require(odd_pair(mm, nn, scale, leg) == after, "control child")
        require(not all(map(prime_trial, after)), "guard failure control")
    require([prime_trial(12 * k * k - 1) for k in range(1, 6)] == [True] * 4 + [False], "D ray first failure")
    require(sigma_direct(22) == 36 and sigma_direct(66) == 144, "2-3-11 anchor")
    require(sigma_direct(70) == 144, "2-5-7 distinct mechanism")
    require(sigma_direct(115) == 144 and sigma_direct(517) == 576, "valid scaling")
    require(sigma_direct(2185) == 2880 != 2304, "scaling prime-guard hostile")
    require(odd_pair(2, 1, 6, "a") == (11, 47), "gap sidecar a")
    require(odd_pair(2, 1, 6, "b") == (5, 53), "gap sidecar b")
    require(odd_pair(2, 1, 1, "a") == (1, 7) and not prime_trial(1),
            "invalid root guard")
    require(odd_pair(3, 2, 1, "a") == (7, 17) and sigma_direct(119) == 144,
            "valid child of invalid parent forbids subtree pruning")
    require(sigma_direct(517) != sigma_direct(265), "same target different squares")
    for kk in range(1, 100):
        q = 12 * kk * kk - 1
        if prime_trial(q):
            target = 72 * kk * kk + 6
            require(sigma_direct(6 * q) - 6 * q == target, "lifted D target")
            require(target - 3 == 3 * (24 * kk * kk + 1) and not prime_trial(target - 3),
                    "entire lifted target has no class-II source")
    require(prime_trial(293) and sigma_direct(293 ** 2) - 293 ** 2 == 294,
            "lifted fibers need not be singletons")
    fixed_gap_counts = {}
    for r in range(1, 101):
        factor_candidates = set()
        for aa in range(1, r):
            if r * r % aa:
                continue
            bb = r * r // aa
            if (aa + bb) % 2:
                continue
            cc = (aa + bb) // 2
            p, q = cc - r - 1, cc + r - 1
            if p % 2 == 1 and prime_trial(p) and prime_trial(q):
                factor_candidates.add((p, q))
        pmax = (r * r + 1) // 2 - r - 1
        scan_candidates = set()
        for p in range(3, pmax + 1, 2):
            q = p + 2 * r
            ss = (p + 1) * (q + 1)
            if isqrt(ss) ** 2 == ss and prime_trial(p) and prime_trial(q):
                scan_candidates.add((p, q))
        require(factor_candidates == scan_candidates, "complete fixed-gap factor/scan comparison")
        fixed_gap_counts[r] = len(factor_candidates)
    require(all(fixed_gap_counts[r] == 0 for r in range(1, 5)) and fixed_gap_counts[5] == 1,
            "least odd-prime gap is ten")
    receipt = {
        "status": "PASS", "source_limit": limit, "distinct_prime_pairs_checked": tested_pairs,
        "square_sigma_sources": len(direct), "chart_counts": chart_counts,
        "primitive_hypotenuse_bound": cmax, "primitive_triples_checked": len(tree),
        "first_sources": [[source, *direct[source]] for source in sorted(direct)[:20]],
        "source_records_sha256": sha256(json.dumps(sorted(direct.items()), separators=(",", ":")).encode()).hexdigest(),
        "branch_failure_controls": controls,
        "fixed_gap_half_gaps_checked": 100,
        "named_checks": ["complete distinct-prime-pair census", "independent divisor sums",
                         "independent prime trial division for accepted pairs", "Euclid/tree census",
                         "decorated inverse and injectivity", "odd/even complete partition",
                         "three-times D ray", "six branch failures", "first even-D prime failure",
                         "2-3-11 and 2-5-7 anchors", "scaling and gap-sidecar hostiles",
                         "universal pair-child formulas", "lifted target class-II exclusion",
                         "k=2 nonsingleton fiber hostile", "invalid-parent valid-child no-pruning",
                         "fixed-gap factor pairs versus bounded-p scan", "least odd-prime gap ten"],
    }
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=1_000_000)
    main(parser.parse_args().limit)
