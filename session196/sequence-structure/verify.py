#!/usr/bin/env python3
"""Exact, dependency-free controls for the S196 Hamiltonian pair-mass theorem."""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import comb, factorial
from pathlib import Path
import hashlib
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def pair_mass(n):
    if n == 0:
        return 1
    return factorial(n) + sum(
        2**r * comb(j - 1, r - 1) * comb(n - 2*j, r) * factorial(n - 2*j)
        for j in range(1, (n - 1)//2 + 1)
        for r in range(1, min(j, n - 2*j) + 1)
    )


def perm_mass(n):
    total = 0
    for p in permutations(range(n)):
        if any(a == b + 1 for a, b in zip(p, p[1:])):
            continue
        total += 1 << sum(a + 1 == b for a, b in zip(p, p[1:]))
    return total


def edges(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def rows_from_mask(n, mask):
    rows = [0]*n
    for b, (i, j) in enumerate(edges(n)):
        if mask >> b & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return rows


def hamilton(rows):
    n = len(rows)
    if not n:
        return 1
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full + 1):
        for v in range(n):
            count = dp[mask][v]
            if not count:
                continue
            nxt = rows[v] & (full ^ mask)
            while nxt:
                bit = nxt & -nxt
                u = bit.bit_length() - 1
                dp[mask | bit][u] += count
                nxt -= bit
    return sum(dp[full])


def literal_hamilton(rows):
    return sum(all(rows[a] >> b & 1 for a, b in zip(p, p[1:]))
               for p in permutations(range(len(rows))))


def relabel_mask(rows, order, es):
    return sum(1 << b for b, (i, j) in enumerate(es)
               if rows[order[i]] >> order[j] & 1)


def exhaustive_moments(max_n=6):
    census = []
    expected_classes = [1, 1, 1, 2, 4, 12, 56]
    for n in range(1, max_n + 1):
        es = edges(n)
        count = 1 << len(es)
        hs = []
        for mask in range(count):
            rows = rows_from_mask(n, mask)
            h = hamilton(rows)
            require(h % 2 == 1, ('parity', n, mask))
            if n <= 5:
                require(h == literal_hamilton(rows), ('literal count', n, mask))
            hs.append(h)
        w = pair_mass(n)
        require(sum(hs) == Fraction(count * factorial(n), 2**(n-1)), ('first moment', n))
        require(sum(h*h for h in hs) == Fraction(count * factorial(n) * w, 4**(n-1)), ('second moment', n))

        seen = bytearray(count)
        raw_spectrum = Counter()
        weighted1 = Fraction(0)
        weighted2 = Fraction(0)
        all_perms = list(permutations(range(n)))
        orbit_sizes = Counter()
        for mask in range(count):
            if seen[mask]:
                continue
            rows = rows_from_mask(n, mask)
            orbit = {relabel_mask(rows, p, es) for p in all_perms}
            for other in orbit:
                require(hs[other] == hs[mask], ('isomorphism HP invariant', n, mask, other))
                seen[other] = 1
            require(factorial(n) % len(orbit) == 0, 'orbit divisor')
            aut = factorial(n)//len(orbit)
            require(hs[mask] % aut == 0, ('free Hamilton action', n, mask))
            orbit_sizes[aut] += 1
            raw_spectrum[hs[mask]] += 1
            weighted1 += Fraction(hs[mask], aut)
            weighted2 += Fraction(hs[mask]**2, aut)
        require(sum(raw_spectrum.values()) == expected_classes[n], ('A000568', n))
        require(weighted1 == 2**((n-1)*(n-2)//2), ('orbit first moment', n))
        require(weighted2 * 4**(n-1) == count*w, ('orbit second moment', n))
        census.append({
            'n': n, 'labelled_tournaments': count,
            'ordinary_classes_A000568': sum(raw_spectrum.values()),
            'sum_H': sum(hs), 'sum_H_squared': sum(h*h for h in hs),
            'W_n': w, 'class_H_multiplicities': dict(sorted(raw_spectrum.items())),
            'automorphism_order_multiplicities': dict(sorted(orbit_sizes.items())),
            'orbit_first_mass': str(weighted1), 'orbit_second_mass': str(weighted2),
        })
    return census


def multiply(a, b, size):
    out = [Fraction(0)]*size
    for i, x in enumerate(a[:size]):
        for j, y in enumerate(b[:size-i]):
            out[i+j] += x*y
    return out


def series_mass(max_n):
    size = max_n + 1
    block = [Fraction(0)]*size
    block[1] = 1
    for k in range(3, size, 2):
        block[k] = 2
    power = [Fraction(1)] + [Fraction(0)]*max_n
    answer = [Fraction(0)]*size
    for m in range(size):
        for k in range(size):
            answer[k] += factorial(m)*power[k]
        power = multiply(power, block, size)
    require(all(value.denominator == 1 for value in answer), 'series integrality')
    return [int(value) for value in answer]


def asymptotic_coefficients(degree):
    size = degree + 1
    answer = [Fraction(0)]*size
    for j in range(1, size):
        for r in range(1, j + 1):
            shift = 2*j - r
            if shift >= size:
                continue
            product = [Fraction(1)] + [Fraction(0)]*degree
            for i in range(r):
                product = multiply(product, [1, -(2*j+i)], size)
            for i in range(2*j):
                product = multiply(product, [Fraction(i**k) for k in range(size)], size)
            factor = Fraction(2**r * comb(j-1, r-1), factorial(r))
            for k in range(size-shift):
                answer[k+shift] += factor*product[k]
    return answer


def check_tail_bounds():
    tested = 0
    for n in range(12, 151):
        layers = [Fraction(0)]*((n-1)//2 + 1)
        for j in range(1, len(layers)):
            layers[j] = sum((Fraction(2**r * comb(j-1, r-1) * comb(n-2*j, r) * factorial(n-2*j), factorial(n))
                             for r in range(1, min(j, n-2*j)+1)), Fraction(0))
        require(sum(layers) == Fraction(pair_mass(n), factorial(n))-1, ('tail mass', n))
        for J in range(min(8, len(layers))):
            bound = Fraction(2, 3)*Fraction(6, n)**(J+1)/(1-Fraction(6, n))
            require(sum(layers[J+1:]) <= bound, ('tail bound', n, J))
            tested += 1
    return tested


def weighted_hamilton(weights):
    n = len(weights)
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full+1):
        outside = full ^ mask
        for v in range(n):
            value = dp[mask][v]
            if not value:
                continue
            nxt = outside
            while nxt:
                bit = nxt & -nxt
                u = bit.bit_length()-1
                dp[mask | bit][u] += value*weights[v][u]
                nxt -= bit
    return sum(dp[full])


def deterministic_witness(n):
    weights = [[int(u != v) for v in range(n)] for u in range(n)]
    for u in range(n-1):
        weights[u][u+1], weights[u+1][u] = 2, 0
    current = weighted_hamilton(weights)
    require(current == pair_mass(n), ('planted mean', n))
    trace = []
    for u, v in edges(n):
        if v == u + 1:
            continue
        weights[u][v], weights[v][u] = 2, 0
        forward = weighted_hamilton(weights)
        weights[u][v], weights[v][u] = 0, 2
        backward = weighted_hamilton(weights)
        require(forward + backward == 2*current, ('conditional mean', n, u, v))
        chosen_forward = forward >= backward
        if chosen_forward:
            weights[u][v], weights[v][u] = 2, 0
        new = max(forward, backward)
        require(new >= current, ('mean decreased', n, u, v))
        trace.append({'pair': [u, v], 'before_scaled_mean': current,
                      'forward_scaled_mean': forward, 'backward_scaled_mean': backward,
                      'chosen_forward': chosen_forward})
        current = new
    rows = [sum(1 << v for v in range(n) if weights[u][v] == 2) for u in range(n)]
    h = hamilton(rows)
    require(current == 2**(n-1)*h, ('final DP', n))
    require(current >= pair_mass(n), ('constructive bound', n))
    if n <= 8:
        require(literal_hamilton(rows) == h, ('witness literal', n))
    return {'n': n, 'outgoing_masks': rows, 'H': h,
            'guaranteed_lower_bound': str(Fraction(pair_mass(n), 2**(n-1))),
            'decision_trace': trace}


def main():
    brute = []
    for n in range(1, 10):
        value = perm_mass(n)
        require(value == pair_mass(n), ('independent permutation mass', n))
        brute.append(value)
    series = series_mass(30)
    require(series == [pair_mass(n) for n in range(31)], 'formal composition identity')
    coefficients = asymptotic_coefficients(7)
    expected = [Fraction(0), Fraction(2), Fraction(0), Fraction(-14, 3),
                Fraction(-58, 3), Fraction(-916, 15), Fraction(-6986, 45), Fraction(-13010, 63)]
    require(coefficients == expected, 'asymptotic coefficients')
    receipt = {
        'status': 'PASS', 'scope': 'All-order proof in PROOF.md; explicitly finite implementation controls here.',
        'labelled_and_orbit_census': exhaustive_moments(),
        'permutation_W_values_n1_to_n9': brute,
        'formal_composition_verified_through_n': 30,
        'variance_asymptotic_coefficients_x0_to_x7': [str(c) for c in coefficients],
        'uniform_tail_checks': check_tail_bounds(),
        'deterministic_witnesses': [deterministic_witness(n) for n in range(1, 11)],
        'W_n_values_n0_to_n30': series,
    }
    encoded = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    Path(__file__).with_name('receipt.json').write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
