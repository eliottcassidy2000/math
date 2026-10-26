#!/usr/bin/env python3
"""Pepin's test as the valuation-two Collatz branch A(x) = (3x+1)/4 (opus-2026-10-26-S58).

klein's THM-4793 (P) reads Pepin through the odd branch x -> 3x+1 on Z/F_m.  Here the carrier is the valuation-two
branch A(x) = (3x+1)/4, the Syracuse step on x = 1 (mod 8):
  (i)   on integers, the Fermat line (THM-4793 (F)) is the A-orbit of F_m: S^k(F_m) = A^k(F_m) = 3^k 2^(2^m-2k) + 1;
  (ii)  A fixes 1 (the trivial cycle) and is conjugate to y -> (3/4) y, y = x - 1;
  (iii) mod F_m (m >= 1): F_m is prime iff A has exactly two cycles on Z/F_m, {1} and one cycle of length 2^(2^m)
        through 0 = F_m; equivalently A^((F_m-1)/2) is the reflection x -> 2 - x (half-turn about the trivial cycle);
  (iv)  the integer Fermat line is the initial arc of that cycle (A^k(F_m) = A^k(0) mod F_m);
  (v)   reflection x -> 2 - x commutes with A; it sends F_m to 1 - 2^(2^m) and the run of 1 - 2^K mirrors the line.
Mersenne counterpart: B(x) = (3x+1)/2 (valuation one) fixes -1 (the negative trivial cycle), the Mersenne line of
M_p = 2^p - 1 (THM-4556) is the B-orbit of M_p, and mod a Mersenne prime B^((M_p-1)/2) is the reflection about -1
(Euler's criterion, 3/2 a non-residue) -- necessary, not sufficient; the cycle type of B is computed.
Usage: python collatz_pepin_valuation_two_20261026.py
"""
import json, os, sys, math, time
from sympy import factorint, isprime, n_order

HERE = os.path.dirname(os.path.abspath(__file__))

# complete factorisations of composite Fermat numbers used here (classical: Euler 1732; Landry 1880;
# Morrison-Brillhart 1970)
FERMAT_FACTORS = {
    5: [641, 6700417],
    6: [274177, 67280421310721],
    7: [59649589127497217, 5704689200685129054721],
}


def order_mod(a, q, qm1_fact):
    """multiplicative order of a mod prime q, given the factorisation of q - 1."""
    o = q - 1
    for p, e in qm1_fact.items():
        for _ in range(e):
            if pow(a, o // p, q) == 1:
                o //= p
            else:
                break
    return o


def affine_cycle_count_squarefree(primes, a_num, a_den):
    """cycles of y -> a y on Z/N, N = prod primes (squarefree): sum over subsets T of prod_{q in T}(q-1)/lcm orders."""
    ords = []
    for q in primes:
        a = a_num * pow(a_den, -1, q) % q
        ords.append(order_mod(a, q, factorint(q - 1)))
    total = 0
    k = len(primes)
    for mask in range(1 << k):
        num = 1
        L = 1
        for i in range(k):
            if mask >> i & 1:
                num *= primes[i] - 1
                L = L * ords[i] // math.gcd(L, ords[i])
        total += num // L
    return total, ords


def affine_power(a, b, e, N):
    """(a, b) represents x -> a x + b mod N; return its e-th iterate."""
    ra, rb = 1, 0
    pa, pb = a % N, b % N
    while e:
        if e & 1:
            ra, rb = ra * pa % N, (ra * pb + rb) % N
        pa, pb = pa * pa % N, (pa * pb + pb) % N
        e >>= 1
    return ra, rb


if __name__ == '__main__':
    out = {}
    t0 = time.time()
    # (i), (iv), (v): integer checks on the line and its mirror
    rows = []
    for m in range(2, 12):
        K = 2 ** m
        F = 2 ** K + 1
        x = F
        mirror = 1 - 2 ** K
        ok = True
        for k in range(2 ** (m - 1)):          # k = 0 .. 2^(m-1)-1 valuation-two steps
            if (3 * x + 1) % 4 or ((3 * x + 1) // 4) % 2 == 0 and k < 2 ** (m - 1) - 1:
                pass
            assert x == 3 ** k * 2 ** (K - 2 * k) + 1
            assert mirror == 2 - x                 # the mirror orbit
            v = 0
            y = 3 * x + 1
            while y % 2 == 0:
                y //= 2
                v += 1
            if k < 2 ** (m - 1) - 1:
                assert v == 2, (m, k, v)
                x = (3 * x + 1) // 4
                ym = 3 * mirror + 1
                assert ym % 4 == 0
                mirror = ym // 4
            else:
                assert v == 3 and y == (3 ** (2 ** (m - 1)) + 1) // 2
        rows.append(m)
    out['line_and_mirror_checked_m'] = rows
    # (ii)-(iv) mod F_m
    res = {}
    for m in range(1, 8):
        F = 2 ** (2 ** m) + 1
        inv4 = pow(4, -1, F)
        a, b = 3 * inv4 % F, inv4 % F                # A(x) = a x + b
        assert (a * 1 + b) % F == 1                   # fixed point 1
        ha, hb = affine_power(a, b, (F - 1) // 2, F)  # half-turn
        half_is_reflection = (ha == F - 1 and hb == 2)
        prime = isprime(F)
        row = {'F_m_prime': prime, 'half_period_iterate_is_x->2-x': half_is_reflection}
        if prime:
            o = n_order(a, F)
            row['ord_3/4'] = o
            row['cycles'] = 1 + (F - 1) // o
            # the integer line is the initial arc of the cycle through 0
            x, z = F, 0
            for k in range(2 ** (m - 1)):
                assert x % F == z
                x = (3 * x + 1) // 4 if k < 2 ** (m - 1) - 1 else x
                z = (a * z + b) % F
        else:
            prs = FERMAT_FACTORS[m]
            assert math.prod(prs) == F
            cyc, ords = affine_cycle_count_squarefree(prs, 3, 4)
            L = 1
            for o in ords:
                L = L * o // math.gcd(L, o)
            row['factor_orders_of_3/4'] = ords
            row['ord_3/4_mod_F'] = L
            row['cycles'] = cyc
            row['cycle_through_F_m_length'] = L
            row['pepin_deficit_log2'] = round(math.log2(F - 1) - math.log2(L), 3)
        res[m] = row
        print(f'm={m}: F_m prime {prime}; half-period iterate of A is x -> 2-x: {half_is_reflection}; '
              f'{ {k: v for k, v in row.items() if k not in ("F_m_prime",)} }', flush=True)
    out['A_mod_F_m'] = res
    # (vi) antipodal arcs: F_(2^n) = 4^M + 1 = 2 (mod F_n), its line is 2 - (line of F_n) mod F_n, and
    #      kind-pasteur's Theorem L (F_n | 3^M + 1) is the half-turn
    anti = {}
    for n in range(1, 6):
        F = 2 ** (2 ** n) + 1
        M = (F - 1) // 2
        b2 = (pow(4, M, F) + 1) % F
        inv4 = pow(4, -1, F)
        a = 3 * inv4 % F
        x0, x2, ok = 0, b2, True
        for j in range(min(M, 4000)):
            ok &= (x0 + x2) % F == 2
            x0, x2 = (a * x0 + inv4) % F, (a * x2 + inv4) % F
        anti[n] = {'F_(2^n) mod F_n': b2, 'antipodal_arcs': ok, 'F_n | 3^M + 1': (pow(3, M, F) + 1) % F == 0}
        print(f'antipodal n={n}: {anti[n]}', flush=True)
    out['antipodal_arcs'] = anti
    # Mersenne counterpart: B(x) = (3x+1)/2 mod M_p, fixed point -1
    mer = {}
    for p in [3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127]:
        M = 2 ** p - 1
        a, b = 3 * pow(2, -1, M) % M, pow(2, -1, M)
        assert (a * (M - 1) + b) % M == M - 1          # fixed point -1
        ha, hb = affine_power(a, b, (M - 1) // 2, M)
        half_reflection = (ha == M - 1 and hb == M - 2)   # x -> -2 - x
        fac = factorint(M - 1)
        o = order_mod(a, M, fac)
        mer[p] = {'half_period_iterate_is_x->-2-x': half_reflection, 'ord_3/2': o, 'M-1': M - 1,
                  'cycles': 1 + (M - 1) // o, '3/2_primitive_root': o == M - 1}
        # line check: S^k(M_p) = 3^k 2^(p-k) - 1 = B^k(M_p), mod M_p equal to B^k(0)
        x, z = M, 0
        for k in range(p - 1):
            assert x == 3 ** k * 2 ** (p - k) - 1 and x % M == z
            x = (3 * x + 1) // 2
            z = (a * z + b) % M
        print(f'M_{p}: B^((M-1)/2) = reflection about -1: {half_reflection}; ord(3/2) = {o}; '
              f'cycles {1 + (M - 1) // o}; primitive {o == M - 1}', flush=True)
    out['B_mod_M_p'] = mer
    print(f'({time.time() - t0:.0f}s)')
    json.dump(out, open(os.path.join(HERE, 'collatz_pepin_valuation_two.json'), 'w'), indent=1, default=str)
