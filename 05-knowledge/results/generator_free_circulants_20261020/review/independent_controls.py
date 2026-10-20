#!/usr/bin/env python3
"""Independent audit controls; no imports from the producer's implementation.

Source fixed points use GF(2) elimination on m variables. Target exact
period/parity counts use direct Mobius inversion, not divisor recursion.
Exhaustive controls use integer subset masks and image actions.
"""
import argparse
import json
from math import gcd
from pathlib import Path

CHECKS = 0


def require(value, description):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(description)


def divisors(n):
    return [x for x in range(1, n + 1) if n % x == 0]


def mobius(n):
    count, prime = 0, 2
    while prime * prime <= n:
        if n % prime == 0:
            n //= prime
            count += 1
            if n % prime == 0:
                return 0
        prime += 1
    return (-1) ** (count + (n > 1))


def phi(n):
    return sum(gcd(i, n) == 1 for i in range(1, n + 1))


def rotation(mask, n, k):
    k %= n
    return ((mask << k) | (mask >> (n - k))) & ((1 << n) - 1)


def representative(mask, n):
    return min(rotation(mask, n, k) for k in range(n))


def action(mask, n, a, k=0):
    result = 0
    while mask:
        low = mask & -mask
        position = low.bit_length() - 1
        result |= 1 << ((a * position + k) % n)
        mask ^= low
    return result


def source_words(m):
    lowmask = (1 << m) - 1
    return [low | ((low ^ lowmask) << m) for low in range(1 << m)]


def primitive_length(mask, m):
    return next(d for d in divisors(m) if rotation(mask, m, d) == mask)


def target_words(m):
    result = []
    for mask in range(1 << m):
        d = primitive_length(mask, m)
        if (mask & ((1 << d) - 1)).bit_count() % 2 == 0:
            result.append(mask)
    return result


def canonical_projection(mask, m):
    h = sum(1 << (2 * j) for j in range(m))
    selected = mask & h
    if selected.bit_count() % 2:
        selected ^= h
    return selected | rotation(selected, 2 * m, m)


def source_linear_fixed(m, a, k):
    """f_i=f_(a*i+k), with f_(i+m)=1+f_i, by row elimination."""
    pivots = {}
    variable_mask = (1 << m) - 1
    for i in range(m):
        j = (a * i + k) % (2 * m)
        row = (1 << i) ^ (1 << (j % m)) ^ ((j // m) << m)
        while row & variable_mask:
            pivot = ((row & variable_mask) & -(row & variable_mask)).bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
        else:
            if row:
                return 0
    return 1 << (m - len(pivots))


def affine_cycles(e, a, k):
    unused = set(range(e))
    lengths = []
    while unused:
        start = min(unused)
        cycle = {start}
        cur = (a * start + k) % e
        while cur != start:
            cycle.add(cur)
            cur = (a * cur + k) % e
        unused -= cycle
        lengths.append(len(cycle))
    return lengths


def target_mobius_fixed(m, a, k):
    result = 0
    for d in divisors(m):
        exact_even = 0
        for e in divisors(d):
            mu = mobius(d // e)
            if not mu:
                continue
            lengths = affine_cycles(e, a, k)
            needs_parity = (d // e) % 2 and any(x % 2 for x in lengths)
            exact_even += mu * (1 << (len(lengths) - bool(needs_parity)))
        require(exact_even >= 0, "Mobius exact-period count nonnegative")
        result += exact_even
    return result


def fixed_type_counts(m, a):
    source_sum = sum(source_linear_fixed(m, a, k) for k in range(2 * m))
    target_sum = sum(target_mobius_fixed(m, a, k) for k in range(m))
    require(source_sum % (2 * m) == target_sum % m == 0, "coset Burnside integrality")
    return [source_sum // (2 * m), target_sum // m]


def common_type_total(m):
    numerator = sum(phi(d) * (1 << (m // d)) for d in divisors(m) if d % 2)
    require(numerator % (2 * m) == 0, "closed common-type formula integral")
    return numerator // (2 * m)


def small_controls():
    rows = []
    for m in range(1, 14):
        n = 2 * m
        xs, ys = source_words(m), target_words(m)
        xorbits = {representative(x, n) for x in xs}
        yorbits = {representative(y, m) for y in ys}
        units = [a for a in range(1, n) if gcd(a, n) == 1]
        require(len(xorbits) == len(yorbits) == common_type_total(m), "common total m=%d" % m)
        fixed = []
        for a in units:
            brute = [sum(representative(action(x, n, a), n) == x for x in xorbits),
                     sum(representative(action(y, m, a), m) == y for y in yorbits)]
            require(brute == fixed_type_counts(m, a), "type fixed points m=%d,a=%d" % (m, a))
            # Direct labelled affine controls, stronger than the quotient sum.
            if m <= 8:
                for k in range(n):
                    require(sum(action(x, n, a, k) == x for x in xs) == source_linear_fixed(m, a, k),
                            "GF2 individual affine fixed count")
                for k in range(m):
                    require(sum(action(y, m, a, k) == y for y in ys) == target_mobius_fixed(m, a, k),
                            "Mobius individual affine fixed count")
            fixed.append([a, *brute])
        if m % 2:
            fibers = {}
            type_map = {}
            for x in xs:
                y = canonical_projection(x, m)
                require(rotation(y, n, m) == y, "projection symmetric")
                low = y & ((1 << m) - 1)
                require(low in ys, "projection target membership")
                fibers.setdefault(y, []).append(x)
                type_map.setdefault(representative(x, n), set()).add(representative(low, m))
                for a in units:
                    require(canonical_projection(action(x, n, a), m) == action(y, n, a), "Aut(C) equivariance")
                require(canonical_projection(rotation(x, n, 1), m) == rotation(y, n, 1), "C equivariance")
                # The output at identity is parity of all H coordinates except identity.
                other_h = sum(1 << (2 * j) for j in range(1, m))
                require((y & 1) == ((x & other_h).bit_count() % 2), "FO parity reduction bit")
            require(len(fibers) == len(ys), "onto labelled target")
            for pair in fibers.values():
                require(len(pair) == 2 and pair[0] ^ pair[1] == (1 << n) - 1,
                        "exact complementary labelled fiber")
            require(all(len(images) == 1 for images in type_map.values()), "well-defined type map")
            images = {next(iter(values)) for values in type_map.values()}
            require(images == yorbits and len(images) == len(xorbits), "type bijection")
        elif m >= 4:
            witness = ((1 << m) - 1) << m
            image = action(witness, n, m + 1)
            require((witness ^ rotation(witness, n, 1)).bit_count() == 2, "witness transitions")
            require((image ^ rotation(image, n, 1)).bit_count() == n - 2, "moved witness transitions")
            require(representative(image, n) != representative(witness, n), "nonfixed source type")
            require(all(action(y | (y << m), n, m + 1) == y | (y << m) for y in ys),
                    "target pointwise fixed by alpha")
            if m % 4 == 0:
                require(fixed_type_counts(m, m + 1)[0] == 0, "strong no-fixed-type result")
        else:
            require(len(xorbits) == len(yorbits) == 1, "m=2 exception")
        rows.append({"m": m, "types": len(xorbits), "fixed": fixed})
    return rows


def prime_power_controls():
    expected = {9: [3, 2, 2, 0, 2], 25: [5, 2, 172, 0, 46],
                27: [3, 3, 316, 16, 16], 49: [7, 2, 349536, 0, 3372],
                81: [3, 4, 13743895360, 0, 2522]}
    rows = []
    for q, (p, r, total, xf, yf) in expected.items():
        m = (q - 1) // 2
        require(fixed_type_counts(m, 1) == [total, total], "prime-power totals")
        require(common_type_total(m) == total, "closed total formula")
        require(fixed_type_counts(m, p) == [xf, yf], "prime-power Frobenius fixed types")
        rows.append({"q": q, "p": p, "degree": r, "types": total, "fixed": [xf, yf]})
    return rows


def f9_hostile():
    # Integers a+3b encode a+b*t with t^2=2. No producer field code reused.
    def plus(x, y): return (x % 3 + y % 3) % 3 + 3 * ((x // 3 + y // 3) % 3)
    def minus(x): return (-x % 3) % 3 + 3 * ((-(x // 3)) % 3)
    def times(x, y):
        a, b, c, d = x % 3, x // 3, y % 3, y // 3
        return (a*c + 2*b*d) % 3 + 3 * ((a*d + b*c) % 3)
    powers, cur = [], 1
    for _ in range(8):
        powers.append(cur)
        cur = times(cur, 4)
    require(set(powers) == set(range(1, 9)) and cur == 1, "F9 primitive element")
    connection = {powers[i] for i in [2, 3, 6, 7]}
    require(connection == {3, minus(3), 7, minus(7)}, "F9 stated connection set")
    encoding = {(u, v): plus(times(u, 3), times(v, 7)) for u in range(3) for v in range(3)}
    require(len(set(encoding.values())) == 9, "F9 stated basis")
    edges = {(i, j) for i in range(9) for j in range(i + 1, 9) if plus(j, minus(i)) in connection}
    rook_edges = {tuple(sorted((encoding[u], encoding[v]))) for u in encoding for v in encoding
                  if u != v and (u[0] == v[0] or u[1] == v[1])}
    require(edges == rook_edges and len(edges) == 18, "F9 rook identification")
    perm = {encoding[(u, v)]: encoding[((-u) % 3, v)] for u, v in encoding}
    require({tuple(sorted((perm[i], perm[j]))) for i, j in edges} == edges, "rook involution automorphism")
    reversed_edges = sum(perm[i] > perm[j] for i, j in edges)
    swapped_edges = sum(perm[i] == j and perm[j] == i for i, j in edges)
    require(swapped_edges == 3 and reversed_edges % 2 == 1, "rook odd edge sign")
    xs = {representative(x, 8) for x in source_words(4)}
    require(len(xs) == 2 and all(representative(action(x, 8, 3), 8) != x for x in xs),
            "Frobenius swaps the two F9 multiplier tournament types")
    return {"vertices": 9, "edges": 18, "swapped_edges": swapped_edges,
            "reference_orientation_reversed_edges": reversed_edges,
            "source_multiplier_types": 2, "source_semilinear_types": 1,
            "target_even_root_types": 2, "rook_fully_sign_even": False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = {"status": "PASS", "small": small_controls(),
              "prime_powers": prime_power_controls(), "f9": f9_hostile(), "checks": CHECKS}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    Path(args.output).write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
