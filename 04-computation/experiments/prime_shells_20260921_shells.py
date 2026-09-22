"""Exact controls: triangle shell cosets, signed orders and Wieferich lifts.

No external packages. Every check remains active under python -O.
"""
from collections import Counter
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path
import json

CHECKS = 0


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def factors(n):
    out = {}
    q = 2
    while q*q <= n:
        while n % q == 0:
            out[q] = out.get(q, 0)+1
            n //= q
        q += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def phi(n):
    ans = n
    for p in factors(n):
        ans = ans//p*(p-1)
    return ans


def order2(n):
    ans = phi(n)
    for p in factors(ans):
        while ans % p == 0 and pow(2, ans//p, n) == 1:
            ans //= p
    return ans


def signed_order(n):
    o = order2(n)
    return o//2 if o % 2 == 0 and pow(2, o//2, n) == n-1 else o


def valuation(n, p):
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def oddrep(x, u):
    x %= u
    if x == 0:
        raise ValueError('zero is outside the unit shell')
    return x if x % 2 else u-x


def cycles(u):
    unseen = {v for v in range(1, u, 2) if gcd(u, v) == 1}
    result = []
    while unseen:
        start = min(unseen)
        v = start
        cyc = []
        while v in unseen:
            unseen.remove(v)
            cyc.append(v)
            v = abs(u-2*v)
        require(v == start, 'permutation cycle closes at its own start')
        result.append(cyc)
    return result


def triangle(u, v):
    return (u*v, (u*u-v*v)//2, (u*u+v*v)//2)


# Exact arithmetic in Z[z]/(1+z+...+z^16), with 16 coefficient slots.
def ring(v):
    v = list(v)+[0]*max(0, 17-len(v))
    w = [0]*17
    for i, c in enumerate(v):
        w[i % 17] += c
    return tuple(w[i]-w[16] for i in range(16))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def neg(a):
    return tuple(-x for x in a)


def mul(a, b):
    v = [0]*31
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            v[i+j] += x*y
    return ring(v)


ONE = ring([1])
ZERO = ring([])


def zpower(k):
    v = [0]*17
    v[k % 17] = 1
    return ring(v)


def root(v):
    return add(zpower(v), zpower(-v))


def orbit_polynomial(orbit):
    a = [ONE]
    for v in orbit:
        b = [ZERO]*(len(a)+1)
        for i, coefficient in enumerate(a):
            b[i] = add(b[i], neg(mul(coefficient, root(v))))
            b[i+1] = add(b[i+1], coefficient)
        a = b
    return a


def run():
    c17 = cycles(17)
    require(c17 == [[1, 15, 13, 9], [3, 11, 5, 7]], 'specified 17 cycles')
    for cyc in c17:
        character = {pow(v, 8, 17) for v in cyc}
        require(len(character) == 1, 'Legendre label is constant on each cycle')
        for v in cyc:
            require(oddrep(2*v, 17) == abs(17-2*v), 'fold equals modular doubling')
            require(pow(oddrep(3*v, 17), 8, 17) != pow(v, 8, 17), '3 switches labels')
            require(oddrep(9*v, 17) == oddrep(8*v, 17), 'M3 squared equals D cubed')
            a, b, c = triangle(17, v)
            require(a*a+b*b == c*c and gcd(a, b) == 1, 'primitive shell triple')
    scalar_involutions = [a for a in range(1, 17, 2) if oddrep(a*a, 17) == 1]
    require(scalar_involutions == [1, 13], 'only scalar involutions preserve the two sectors')
    require(all(pow(a, 8, 17) == 1 for a in scalar_involutions), 'no scalar XOR swap')

    eta = ZERO
    for v in c17[0]:
        eta = add(eta, root(v))
    require(add(add(mul(eta, eta), eta), ring([-4])) == ZERO, 'period eta²+eta−4')
    polys = []
    for cyc in c17:
        coefficients = orbit_polynomial(cyc)
        encoded = []
        for c in coefficients:
            beta = -c[3]
            alpha = c[0]+beta
            require(c == add(ring([alpha]), tuple(beta*x for x in eta)), 'quartic over Q(eta)')
            encoded.append([alpha, beta])
        polys.append(encoded)
    product = [ZERO]*9
    for i, a in enumerate(orbit_polynomial(c17[0])):
        for j, b in enumerate(orbit_polynomial(c17[1])):
            product[i+j] = add(product[i+j], mul(a, b))
    full = [1, -4, -10, 10, 15, -6, -7, 1, 1]
    require(product == [ring([v]) for v in full], 'real cyclotomic polynomial product')
    for v in range(1, 17):
        require(add(mul(root(v), root(v)), ring([-2])) == root(2*v), 'exact Chebyshev conjugacy')
    for cyc in c17:
        multiplier = ONE
        for v in cyc:
            multiplier = mul(multiplier, add(root(v), root(v)))
        require(multiplier == ring([-16]), 'four-cycle derivative multiplier is minus16')

    counts = {}
    for u in range(3, 502, 2):
        cs = cycles(u)
        d = signed_order(u)
        require(all(len(c) == d for c in cs), 'uniform signed-order cycle length')
        require(len(cs)*d == phi(u)//2, 'shell partition count')
        local = [signed_order(p**e) for p, e in factors(u).items()]
        L = lcm(*local)
        signs = [pow(2, L, p**e) == 1 for p, e in factors(u).items()]
        require(d == (L if len(set(signs)) == 1 else 2*L), 'CRT signs synchronize')
        counts[u] = [len(cs), d]
    require(signed_order(15) == 4 and lcm(signed_order(3), signed_order(5)) == 2,
            'hostile to discarding local return signs')
    require(counts[63] == [3, 6] and counts[65] == [4, 6], '63 and65 sign mechanisms')

    prime_lifts = []
    primes = [p for p in range(3, 502, 2) if factors(p) == {p: 1}]+[1093, 3511]
    for p in primes:
        d = signed_order(p)
        sign = 1 if pow(2, d, p) == 1 else -1
        r = valuation(2**d-sign, p)
        require(r == valuation(2**(p-1)-1, p), 'same Wieferich valuation')
        row = []
        for e in range(1, 5):
            modulus = p**e
            period = d*p**max(0, e-r)
            number = (p-1)//(2*d)*p**min(e-1, r-1)
            require(pow(2, period, modulus) in (1, modulus-1), 'predicted period returns')
            for q in factors(period):
                require(pow(2, period//q, modulus) not in (1, modulus-1), 'predicted period minimal')
            require(period*number == (p-1)*p**(e-1)//2, 'predicted count partitions units')
            row.append({'exponent': e, 'period': period, 'cycles': number})
            if modulus <= 40000:
                cs = cycles(modulus)
                require(len(cs) == number and all(len(c) == period for c in cs), 'direct lift census')
        if p in (3, 5, 7, 17, 1093, 3511):
            prime_lifts.append({'p': p, 'd': d, 'sign': sign, 'valuation': r, 'levels': row})

    # Compare inverse-affine braids on odd residues modulo 2*p^s.
    braid_rows = []
    for p in (3, 5, 7, 17, 1093, 3511):
        d = order2(p)
        q = 2**d
        c = (q-1)//p
        r = valuation(q-1, p)
        for s in range(1, 4):
            period = p**max(0, s-r+1)
            modulus = 2*p**s
            # Binary exponentiation of the affine map; no huge q**period.
            A, B, aa, bb, n = 1, 0, q % modulus, c % modulus, period
            while n:
                if n & 1:
                    A, B = aa*A % modulus, (aa*B+bb) % modulus
                aa, bb = aa*aa % modulus, (aa*bb+bb) % modulus
                n >>= 1
            for x in (1, 3, 5, 19):
                require((A*x+B) % modulus == x % modulus, 'inverse affine predicted return')
            braid_rows.append([p, s, period, p**min(s, r-1)])

    return {'status': 'FINITE-EXACT PASS', 'checks_passed': CHECKS,
            'universe': {'all_odd_shells': [3, 501], 'prime_lifts': len(primes),
                         'lift_exponents': [1, 4], 'direct_census_max_modulus': 40000},
            'cycles17': c17, 'triangles17': [[triangle(17, v) for v in c] for c in c17],
            'quartic_coefficients_ascending_alpha_beta_eta': polys,
            'real17_minimal_polynomial_ascending': full,
            'prime_lifts': prime_lifts, 'inverse_braid_p_s_period_count': braid_rows,
            'source_sha256_lf': sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest()}


if __name__ == '__main__':
    output = Path(__file__).with_suffix('.json')
    output.write_text('{"status":"RUNNING"}\n', encoding='utf-8', newline='\n')
    result = run()
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status': result['status'], 'checks_passed': CHECKS}))
