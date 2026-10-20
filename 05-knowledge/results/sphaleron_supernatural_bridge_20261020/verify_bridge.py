#!/usr/bin/env python3
"""Exact finite checks for the typed sphaleron/Steinitz reconstruction.

Dependency: SymPy (only for the symbolic degree check).
No computation asserts a physical/arithmetic dynamical equivalence.
"""

import json
import math
from collections import Counter

import sympy as sp


def require(condition, explanation):
    if not condition:
        raise RuntimeError(explanation)


def degree_check():
    mu, theta, phi = sp.symbols("mu theta phi", real=True)
    y = sp.Matrix([
        -sp.sin(mu) * sp.sin(theta) * sp.sin(phi),
        -sp.sin(mu) * sp.sin(theta) * sp.cos(phi),
        sp.sin(mu) * sp.cos(mu) * (sp.cos(theta) - 1),
        sp.cos(mu)**2 + sp.sin(mu)**2 * sp.cos(theta),
    ])
    require(sp.trigsimp(y.dot(y)) == 1, "Boundary map must land in S^3")
    basepoint = sp.Matrix([0, 0, 0, 1])
    for variable, value in [(mu, 0), (mu, sp.pi), (theta, 0)]:
        difference = y.subs(variable, value) - basepoint
        require(all(sp.trigsimp(z) == 0 for z in difference), "Smash basepoint condition")
    jac = sp.trigsimp(sp.det(sp.Matrix.hstack(y, y.diff(mu), y.diff(theta), y.diff(phi))), method="fu")
    expected = (1 - sp.cos(theta)) * sp.sin(mu)**2 * sp.sin(theta)
    require(sp.trigsimp(jac - expected) == 0, "Signed pullback volume")
    integral = sp.integrate(jac, (phi, 0, 2*sp.pi), (theta, 0, sp.pi), (mu, 0, sp.pi))
    degree = sp.simplify(integral / (2*sp.pi**2))
    require(degree == 1, "Degree-one control")
    require(sp.simplify((-integral)/(2*sp.pi**2)) == -1, "Orientation-reversal hostile")
    return {"unit_norm": True, "smash_basepoints": True,
            "pullback_volume_density": str(expected), "degree": int(degree),
            "reversed_orientation_degree": -1}


def poly_mod(a, f):
    while a and a.bit_length() >= f.bit_length():
        a ^= f << (a.bit_length() - f.bit_length())
    return a


def irreducible(f):
    n = f.bit_length() - 1
    for d in range(1, n//2 + 1):
        for g in range(1 << d, 1 << (d+1)):
            if poly_mod(f, g) == 0:
                return False
    return True


def multiply(a, b, f):
    value = 0
    while b:
        if b & 1:
            value ^= a
        a <<= 1
        b >>= 1
    return poly_mod(value, f)


def power(a, e, f):
    value = 1
    while e:
        if e & 1:
            value = multiply(value, a, f)
        a = multiply(a, a, f)
        e >>= 1
    return value


def binary_rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            bit = value.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = value
                break
            value ^= pivots[bit]
    return len(pivots)


def cycles(permutation):
    require(sorted(permutation) == list(range(len(permutation))), "Must be a permutation")
    seen = set()
    counts = Counter()
    for start in range(len(permutation)):
        if start in seen:
            continue
        x, length = start, 0
        while x not in seen:
            seen.add(x)
            x = permutation[x]
            length += 1
        require(x == start, "Cycle closure")
        counts[length] += 1
    return dict(sorted(counts.items()))


def finite_field_checks():
    rows = []
    for n in range(1, 9):
        f = next(f for f in range((1 << n) + 1, 1 << (n+1), 2) if irreducible(f))
        q = 1 << n
        frobenius = [multiply(a, a, f) for a in range(q)]
        fc = cycles(frobenius)
        require(math.lcm(*fc) == n, "Frobenius order equals extension degree")
        necklace_counts = {int(d): int(sum(sp.mobius(d//e) * (2**e) for e in sp.divisors(d)) // d)
                           for d in sp.divisors(n)}
        require(fc == necklace_counts, "Independent binary-necklace cycle count")
        for k in range(n+1):
            fixed = sum(power(a, 2**k, f) == a for a in range(q))
            require(fixed == 2**math.gcd(n, k), "Independent Frobenius fixed-point formula")
        for a in range(q):
            require(power(a, q, f) == a, "Finite-field q-power identity")
            for b in range(q):
                require(frobenius[a ^ b] == frobenius[a] ^ frobenius[b], "Frobenius additivity")
        normal = None
        for a in range(1, q):
            orbit, current = [], a
            for _ in range(n):
                orbit.append(current)
                current = frobenius[current]
            if binary_rank(orbit) == n:
                normal = orbit
                break
        require(normal is not None, "A normal basis was found")
        require(all(frobenius[normal[i]] == normal[(i+1) % n] for i in range(n)), "Normal-basis cyclic shift")
        alpha = next(a for a in range(1, q) if all(power(a, (q-1)//p, f) != 1 for p in sp.factorint(q-1)))
        mult = [multiply(alpha, a, f) for a in range(q)]
        mc = cycles(mult)
        require(math.lcm(*mc) == q-1, "Primitive multiplication order")
        require((q-1) % 2 == 1, "Multiplication cycle has no factor two")
        rows.append({"n": n, "field_modulus_bits": bin(f), "normal_basis": normal,
                     "frobenius_order": n, "frobenius_cycle_counts": fc,
                     "necklace_and_fixed_point_crosschecks": True,
                     "primitive_multiplier": alpha, "multiplication_order": q-1,
                     "multiplication_cycle_counts": mc})
    return rows


def modular_hostiles():
    odd_prime_power_witnesses = []
    for p in [3, 5, 7, 11, 13]:
        for a in range(1, 5):
            modulus = p**a
            order = int(sp.n_order(2, modulus))
            require(pow(2, order, modulus) == 1, "Odd prime-power divisor witness")
            odd_prime_power_witnesses.append({"prime": p, "exponent": a, "degree": order})
    require(not any((2*x) % 2 == 1 for x in range(2)), "No half of one modulo two")
    for m in range(1, 41):
        for n in range(1, m+1):
            if m % n == 0:
                for k in range(-m, m+1):
                    require((k % m) % n == k % n, "Winding residues are coherent")
    return {"half_of_one_mod_2_exists": False,
            "coherent_winding_residues_checked_through_modulus": 40,
            "odd_prime_power_witnesses": odd_prime_power_witnesses}


def main():
    report = {"status": "FINITE-EXACT PASS; not a physical isomorphism",
              "symbolic_degree": degree_check(),
              "finite_fields": finite_field_checks(),
              "modular_controls": modular_hostiles()}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
