#!/usr/bin/env python3
"""S36 referee audit, part E: Pisano periods, Wall-Sun-Sun primes, THM-4657 orbit law, Fibonacci Paley primes.

Independent of the session's fibonacci_pisano_wall_sun_sun_20261014.py:
  * pi(n) is computed as the multiplicative order of Q = [[1,1],[1,0]] in GL_2(Z/n) by testing divisors of the
    a-priori bound 6n with integer matrix powers (not by iterating the sequence), and cross-checked by iteration;
  * Wall-Sun-Sun test via the Lucas sequence: p is WSS iff L_p == 1 (mod p^2) (equivalent classical form, p != 2, 5),
    checked against F_(p - (5/p)) == 0 (mod p^2) computed by matrix powers;
  * the golden ratio mod p^2 is found by brute force over Z/p^2 (no Hensel step); orbits of <g> are counted on Z/p^2;
  * the explosion branch of the orbit law (which no known prime reaches) is tested with units g of odd order whose
    p-th power structure forces ord_(p^2)(g) = ord_p(g), i.e. the 'Wieferich' situation, for every p <= 61.
"""
from math import gcd
from sympy import primerange, factorint, divisors, isprime


def matmul(A, B, m):
    return [[(A[0][0] * B[0][0] + A[0][1] * B[1][0]) % m, (A[0][0] * B[0][1] + A[0][1] * B[1][1]) % m],
            [(A[1][0] * B[0][0] + A[1][1] * B[1][0]) % m, (A[1][0] * B[0][1] + A[1][1] * B[1][1]) % m]]


def matpow(A, e, m):
    R = [[1 % m, 0], [0, 1 % m]]
    while e:
        if e & 1:
            R = matmul(R, A, m)
        A = matmul(A, A, m)
        e >>= 1
    return R


Q = [[1, 1], [1, 0]]


def pisano_order(n):
    """order of Q mod n: first verify Q^M = I for a multiple M built from p^(k-1) * {p-1, 2(p+1), 20, 3}, then take the
    least divisor d of M with Q^d = I (the verification makes this independent of the theory being audited)."""
    if n == 1:
        return 1
    I = [[1 % n, 0], [0, 1 % n]]
    M = 1
    for q, k in factorint(n).items():
        base = 3 if q == 2 else (20 if q == 5 else (q - 1 if q % 5 in (1, 4) else 2 * (q + 1)))
        t = q ** (k - 1) * base
        M = M * t // gcd(M, t)
    if matpow(Q, M, n) != I:
        raise ValueError('multiple not verified')
    for d in divisors(M):
        if matpow(Q, d, n) == I:
            return d


def pisano_iter(n):
    if n == 1:
        return 1
    a, b, k = 0, 1, 0
    while True:
        a, b = b, (a + b) % n
        k += 1
        if (a, b) == (0, 1):
            return k


def legendre5(p):
    r = p % 5
    return 0 if r == 0 else (1 if r in (1, 4) else -1)


def main():
    # (1) Pisano
    N = 2000
    pis = {n: pisano_order(n) for n in range(1, N + 1)}
    mismatch = [n for n in range(1, 400) if pis[n] != pisano_iter(n)]
    six = [n for n in range(1, N + 1) if pis[n] == 6 * n]
    over = [n for n in range(1, N + 1) if pis[n] > 6 * n]
    print(f'(1) pi(n) by matrix order for n <= {N}; iteration cross-check mismatches (n < 400): {mismatch}')
    print(f'    n with pi(n) = 6n: {six}; n with pi(n) > 6n: {over}')
    print(f'    table: {dict((n, pis[n]) for n in (2, 3, 5, 7, 10, 11, 12, 13, 19, 20, 25, 49, 50, 55, 89, 121, 125))}')
    sp_ok = all((p - 1) % pis[p] == 0 for p in primerange(7, N + 1) if p % 5 in (1, 4))
    in_ok = all((2 * (p + 1)) % pis[p] == 0 and (p + 1) % pis[p] != 0 for p in primerange(7, N + 1) if p % 5 in (2, 3))
    # exact equality cases (for the owner's 'one below'): split primes with pi(p) = p - 1, inert with pi(p) = 2(p+1)
    sp_eq = [p for p in primerange(7, 200) if p % 5 in (1, 4) and pis[p] == p - 1]
    sp_lt = [p for p in primerange(7, 200) if p % 5 in (1, 4) and pis[p] < p - 1]
    print(f'    split p: pi(p) | p-1 for all p <= {N}: {sp_ok}; inert p: pi(p) | 2(p+1), not p+1: {in_ok}; pi(2), pi(3), pi(5) = {pis[2]}, {pis[3]}, {pis[5]}')
    print(f'    split p < 200 with pi(p) = p - 1: {sp_eq};  with pi(p) < p - 1: {sp_lt}')
    # prime powers: pi(p^2) = p pi(p) (equivalent to p not WSS) for p < 2000
    pp = [p for p in primerange(2, 45) if pisano_order(p * p) != p * pis[p]]
    print(f'    primes p < 45 with pi(p^2) != p pi(p): {pp}')
    # (2) Wall-Sun-Sun search to 2e6, two tests
    P = 2_000_000
    cnt = 0
    wss = []
    disagree = []
    for p in primerange(7, P + 1):
        cnt += 1
        m = p * p
        # Lucas L_p mod p^2 via Q^p: trace of Q^p = L_p
        M = matpow(Q, p, m)
        Lp = (M[0][0] + M[1][1]) % m
        e = p - legendre5(p)
        Fe = matpow(Q, e, m)[0][1]   # Q^e = [[F_(e+1), F_e], [F_e, F_(e-1)]]
        a = (Lp == 1)
        b = (Fe == 0)
        if a != b:
            disagree.append(p)
        if a or b:
            wss.append(p)
    print(f'(2) primes 7 <= p <= {P}: {cnt}; Wall-Sun-Sun (L_p == 1 or F_(p-(5/p)) == 0 mod p^2): {wss}; tests disagree at: {disagree}')
    # (3) THM-4657 orbit law at split p <= 199, golden ratio by brute force in Z/p^2
    rows = []
    bad = []
    for p in primerange(7, 200):
        if p % 5 not in (1, 4):
            continue
        m = p * p
        roots = [x for x in range(m) if (x * x - x - 1) % m == 0]
        assert len(roots) == 2
        for phi in roots:
            d = 1
            while pow(phi, d, p) != 1:
                d += 1
            dd = d
            while dd % 2 == 0:
                dd //= 2
            g = pow(phi, d // dd, m)
            eps = 0 if pow(phi, p - 1, m) == 1 else 1
            # orbits of <g> on Z/p^2 minus 0, and the negation pairing
            seen = [False] * m
            orbits = []
            for x in range(1, m):
                if not seen[x]:
                    o = []
                    y = x
                    while not seen[y]:
                        seen[y] = True
                        o.append(y)
                        y = y * g % m
                    orbits.append(frozenset(o))
            selfneg = sum(1 for o in orbits if (m - min(o)) % m in o)
            E = len(orbits) // 2
            pred = (p + 1) * (p - 1) // (2 * dd) if eps == 0 else (p - 1) // dd
            ordg = 1
            while pow(g, ordg, m) != 1:
                ordg += 1
            if selfneg or len(orbits) % 2 or E != pred or ordg != dd * p ** eps:
                bad.append((p, phi))
            rows.append((p, phi, d, dd, ordg, eps, E))
    print(f'(3) THM-4657 orbit law: {len(rows)} (p, phi) pairs with split p < 200; failures: {bad}; WSS among them: {[r for r in rows if r[5] == 0]}')
    print(f'    sample (p, phi mod p^2, ord_p phi, odd part d\', ord_(p^2) g, eps, E): {rows[:4]}')
    # (3b) explosion branch: units g of odd order d' with ord_(p^2)(g) = d' (Teichmueller-type), orbit count must be (p+1)(p-1)/(2d')
    bad2 = []
    tested = 0
    for p in primerange(3, 62):
        m = p * p
        # a primitive root mod p^2, then the subgroup of order p - 1 (elements with g^(p-1) = 1 mod p^2)
        for h in range(2, m):
            if gcd(h, p) == 1:
                o = 1
                y = h
                while y != 1:
                    y = y * h % m
                    o += 1
                if o == p * (p - 1):
                    break
        t = pow(h, p, m)   # generator of the order-(p-1) subgroup
        for dprime in divisors(p - 1):
            if dprime % 2 == 0:
                continue
            g = pow(t, (p - 1) // dprime, m)
            seen = [False] * m
            norb = 0
            for x in range(1, m):
                if not seen[x]:
                    norb += 1
                    y = x
                    while not seen[y]:
                        seen[y] = True
                        y = y * g % m
            tested += 1
            if norb != (p + 1) * (p - 1) // dprime:
                bad2.append((p, dprime, norb))
    print(f'    explosion branch: {tested} (p, d\') cases with ord_(p^2) g = d\' odd, p < 62: orbit count (p+1)(p-1)/d\' failures: {bad2}')
    # (4) Fibonacci Paley primes: p = 3 mod 4, split, square root of x^2-x-1 has order (p-1)/2
    fpp = []
    for p in primerange(7, 2000):
        if p % 4 != 3 or p % 5 not in (1, 4):
            continue
        roots = [x for x in range(p) if (x * x - x - 1) % p == 0]
        sq = [r for r in roots if pow(r, (p - 1) // 2, p) == 1]
        assert len(sq) == 1   # phi psi = -1 is a non-square when p = 3 mod 4
        r = sq[0]
        o = 1
        while pow(r, o, p) != 1:
            o += 1
        if o == (p - 1) // 2:
            fpp.append(p)
    print(f'(4) Fibonacci Paley primes below 2000: {len(fpp)}: {fpp[:12]} ...')
    # also: the other root psi = -1/phi is then a primitive root (so NO psi-invariant circulant tournament exists)
    print(f'    at p = 11: roots {[x for x in range(11) if (x * x - x - 1) % 11 == 0]}, <4> = {sorted({pow(4, k, 11) for k in range(5)})}, ord(8) = {[k for k in range(1, 11) if pow(8, k, 11) == 1][0]}')


if __name__ == '__main__':
    main()
