#!/usr/bin/env python3
"""Exact independent controls for the signed duplication SOS theorem.

The companion note gives an all-parameter proof. These finite controls check
the rational Gram identity, coefficient multiplicities, equality boundaries,
and the actual Laurent coefficient map, without importing repository proofs.
Only Python's standard library is used; every gate survives python -O.
"""

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import comb, factorial, gcd, prod
from random import Random


CHECKS = 0
TRACE = sha256()


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def record(*row):
    TRACE.update((repr(row) + "\n").encode())


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def multiply(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i+j] += x*y
    return out


def power(p, exponent):
    out = [1]
    for _ in range(exponent):
        out = multiply(out, p)
    return out


def roots_polynomial(roots):
    out = [1]
    for root in roots:
        out = multiply(out, [1, root])
    return out


def coefficient(p, k):
    return p[k] if 0 <= k < len(p) else 0


def square_coefficient(p, k):
    return sum(p[j]*coefficient(p, k-j) for j in range(len(p)))


def weight(k, s):
    return F(2**(2*k-2*s+1)*factorial(2*s-2)*factorial(k-s)*factorial(k),
             factorial(s-1)*factorial(2*k))


def integral_controls():
    for ell in range(41):
        literal_integral = 4**ell * sum(
            (F((-1)**j*comb(ell, j), ell+j+1) for j in range(ell+1)), F(0))
        increment = F(4**(ell+1), comb(2*ell+2, ell+1)) - F(4**ell, comb(2*ell, ell))
        require(increment == literal_integral, "beta increment from expanded polynomial")
    for k in range(1, 33):
        for s in range(1, k+1):
            literal = sum((F((-1)**j*comb(k-s, j), 2*s+2*j-1)
                           for j in range(k-s+1)), F(0))
            require(weight(k, s) == literal > 0, "SOS coefficient integral")
        require(weight(k, k) == F(1, 2*k-1), "terminal strictly positive weight")
    print("integrals: 41 distance increments; all 528 SOS weights k=1..32")


def gram_controls():
    entries = 0
    for n in range(1, 36):
        for k in range(1, n+1):
            D = comb(n, k)
            C = F(comb(2*n, 2*k), D*D)
            row = sum(comb(k, ell)*comb(n-k, ell)*F(4**ell, comb(2*ell, ell))
                      for ell in range(min(k, n-k)+1))
            require(row == C*D, "constant row sum via independent Johnson counts")
            for ell in range(min(k, n-k)+1):
                total = F(0)
                for s in range(1, k+1):
                    beta = F(comb(n-s, k-s), D)
                    entry = choose(k-ell, s) - 2*beta*comb(k, s) + beta*beta*comb(n, s)
                    total += weight(k, s)*entry
                require(total == C-F(4**ell, comb(2*ell, ell)), "every SOS Gram entry")
                entries += 1
                record("gram", n, k, ell, total)
    print(f"Gram identity: {entries} exact distance entries; n=1..35, every k=1..n")


def signed_controls():
    rows = 0
    equalities = 0
    for n in range(2, 8):
        for roots in combinations_with_replacement((-2, -1, 0, 1, 2), n):
            p = roots_polynomial(roots)
            for k in range(1, n):
                gap = F(comb(2*n, 2*k), comb(n, k)**2)*p[k]**2-square_coefficient(p, 2*k)
                expected_equal = sum(x != 0 for x in roots) < k or len(set(roots)) == 1
                require(gap >= 0, "signed sharp duplication inequality")
                require((gap == 0) == expected_equal, "complete zero-root equality classification")
                rows += 1
                equalities += int(gap == 0)
                record("signed", roots, k, gap)
    print(f"signed exhaustive head: {rows} rows, {equalities} equalities; n=2..7 multisets in -2..2")


def direct_sos_controls():
    rows = 0
    for n in range(2, 9):
        root_lists = [(0,)*n, (1,)*n, tuple((-1)**i for i in range(n)),
                      tuple(range(1, n+1)), tuple(i-2 for i in range(n))]
        for k in range(1, n):
            subsets = list(combinations(range(n), k))
            masks = [sum(1 << i for i in I) for I in subsets]
            for roots in root_lists:
                v = [prod(roots[i] for i in I) for I in subsets]
                p = roots_polynomial(roots)
                kernel = sum(v[i]*v[j]*F(4**((masks[i]^masks[j]).bit_count()//2),
                             comb((masks[i]^masks[j]).bit_count(),
                                  (masks[i]^masks[j]).bit_count()//2))
                             for i in range(len(v)) for j in range(len(v)))
                require(kernel == square_coefficient(p, 2*k), "literal coefficient versus subset kernel")
                require(sum(v) == p[k], "literal elementary coefficient")
                total = F(0)
                for s in range(1, k+1):
                    beta = F(comb(n-s, k-s), comb(n, k))
                    for S in combinations(range(n), s):
                        mask = sum(1 << i for i in S)
                        inner = sum(v[i] for i in range(len(v)) if mask & masks[i] == mask)
                        total += weight(k, s)*(inner-beta*p[k])**2
                gap = F(comb(2*n, 2*k), comb(n, k)**2)*p[k]**2-kernel
                require(gap == total, "direct subset-sum SOS versus expanded polynomial")
                rows += 1
    print(f"direct coefficient/kernel/SOS replay: {rows} signed vector rows; n=2..8")


def zero_coefficient_controls():
    rng = Random(4438)
    rows = 0
    for n in range(3, 17):
        for k in range(1, n):
            for _ in range(25):
                roots = [rng.choice((-7, -3, -2, -1, 1, 2, 3, 7)) for _ in range(n-1)]
                p = roots_polynomial(roots)
                den, num = p[k-1], -p[k]
                if den == 0 or num == 0:
                    continue
                H = multiply(p, [den, num])
                require(H[0] != 0 and H[-1] != 0 and H[k] == 0, "tuned exact real-rooted carrier")
                doubled = square_coefficient(H, 2*k)
                squared_roots = roots_polynomial([r*r for r in roots])
                terminal = den*den*coefficient(squared_roots, k)+num*num*coefficient(squared_roots, k-1)
                require(doubled < 0, "interior zero coefficient has strictly negative square coefficient")
                require(-(2*k-1)*doubled >= terminal > 0, "quantitative terminal-SOS margin")
                rows += 1
                record("tuned", n, k, roots, den, num, doubled)
    print(f"tuned rational-root carriers: {rows} exact zeros; n=3..16; strict quantitative margin")


def primitive(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    if not p:
        return []
    content = 0
    for value in p:
        content = gcd(content, value)
    if p[-1] < 0:
        content = -content
    return [value//content for value in p]


def polynomial_gcd(p, q):
    p, q = primitive(p), primitive(q)
    while q:
        remainder = p[:]
        while remainder and len(remainder) >= len(q):
            shift, lead = len(remainder)-len(q), remainder[-1]
            remainder = [q[-1]*value for value in remainder]
            for i, value in enumerate(q):
                remainder[i+shift] -= lead*value
            remainder = primitive(remainder)
        p, q = q, remainder
    return p


def channel(mass, target):
    return [factorial(mass)//(factorial(mass-target+j)*factorial(target-2*j)*factorial(j))
            for j in range(target//2+1)]


def literal_ap(a, spacing, mass):
    counts = Counter()
    fact = [factorial(j) for j in range(mass+1)]
    for x in range(mass+1):
        for y in range(mass-x+1):
            z = mass-x-y
            if -a*x+(spacing-a)*y+(2*spacing-a)*z == 0:
                counts[z] += fact[mass]//(fact[x]*fact[y]*fact[z])
    return [counts[j] for j in range(max(counts, default=-1)+1)]


def ap_controls():
    rows = 0
    for spacing in range(2, 17):
        for a in range(1, spacing):
            m0 = spacing//gcd(a, spacing)
            for rung in range(1, 4):
                mass = rung*m0
                target = a*mass//spacing
                p, q = channel(mass, target), channel(2*mass, 2*target)
                require(p == literal_ap(a, spacing, mass), "first literal Laurent charge fibre")
                require(q == literal_ap(a, spacing, 2*mass), "doubled literal Laurent charge fibre")
                require(len(polynomial_gcd(p, q)) == 1, "all-mass AP doubled-channel gcd")
                rows += 1
                record("AP", spacing, a, rung, p, q)
    for h in range(1, 21):
        a, spacing = 2*h, 2*h+1
        p = channel(spacing, a)
        require(len(p) == h+1 and all(v > 0 for v in p), "unbounded sharp-family complete first row")
        require(gcd(a, spacing) == 1 and 2*spacing == a+(2*spacing-a), "sharp width identity")
    require(channel(11, 10) == [11, 495, 4620, 11550, 6930, 462], "six-channel width-22 control")
    print(f"AP exact fibre/gcd controls: {rows} rows; spacing=2..16, all a, first three admissible masses")
    print("sharp unbounded-channel family: h=1..20; h=5 support=(-10,1,12), width=22")


def laurent_multiply(p, q):
    result = Counter()
    for i, x in p.items():
        for j, y in q.items():
            result[i+j] += x*y
    return {i: x for i, x in result.items() if x}


def core_controls():
    rows = 0
    cancellations = 0
    for N in range(1, 7):
        root_lists = [(1,)*N, tuple((-1)**j for j in range(N)), tuple(j+1 if j % 2 else -j-1 for j in range(N))]
        for roots in root_lists:
            R = roots_polynomial(roots)
            for spacing in range(1, 5):
                for a in range(1, spacing*N):
                    m0 = spacing//gcd(a, spacing)
                    f = {spacing*j-a: value for j, value in enumerate(R) if value}
                    running = {0: 1}
                    values = []
                    for mass in range(1, 2*m0+1):
                        running = laurent_multiply(running, f)
                        ct = running.get(0, 0)
                        expected = coefficient(power(R, mass), a*mass//spacing) if a*mass % spacing == 0 else 0
                        require(ct == expected, "literal Laurent versus ordinary carrier coefficient")
                        if mass % m0:
                            require(ct == 0, "all intermediate masses excluded by congruence")
                        values.append(ct)
                    first = next((j+1 for j, value in enumerate(values) if value), None)
                    require(first in (m0, 2*m0), "exact two-rung detection for real-rooted core")
                    if values[m0-1] == 0:
                        require(values[2*m0-1] < 0, "real-core normalized doubled sign")
                        cancellations += 1
                    if N == 1:
                        require(first == m0, "linear core has nonzero first congruent coefficient")
                    require(first <= spacing*N, "real-rooted-core linear width bound")
                    rows += 1
                    record("core", roots, spacing, a, values)
    print(f"literal arbitrary-support core controls: {rows} rows, {cancellations} first-rung cancellations; degree=1..6")


def hostile_controls():
    H = multiply(multiply([1, 1], [2, 1]), multiply([-10, 1], [-4, 1]))
    require(H == [80, 92, 0, -11, 1], "real-rooted individual-pair hostile")
    require(H[0]*H[4] == 80 > 0 and H[1]*H[3] == -1012, "positive and negative pairs coexist")
    require(square_coefficient(H, 4) == -1864, "grouped square is negative")
    require(square_coefficient([1, 0, 1], 2) == 2, "dropping real-rootedness gives strict positive hostile")
    require(square_coefficient([0, 0, 1], 2) == 0, "zero-root boundary defeats unjustified strictness")
    cubic = [1, 1, 0, -1]
    require(coefficient(power(cubic, 4), 3) == 0, "actual trinomial wall outside real-rooted-core scope")
    require(coefficient(power(cubic, 8), 6) == -224, "outside-scope wall still passes inherited two-channel theorem")
    require(-4*(-1)-27*(-1)**2 == -23, "compressed cubic has nonreal roots")
    require(polynomial_gcd([1, 2, 1], [1, 1]) == [1, 1], "gcd positive shared-factor control")
    require(polynomial_gcd([1, 0, 1], [1, 1]) == [1], "gcd coprime control")
    print("hostiles: mixed individual pair signs; nonreal quadratic +2; zero-root equality; cubic wall (0,-224), discriminant -23")


def hermite_controls():
    hermite = [[1], [0, 1]]
    for n in range(1, 64):
        p = [0] + hermite[n]
        for j, value in enumerate(hermite[n-1]):
            p[j] -= n*value
        hermite.append(p)
    for n in range(65):
        literal = [0]*(n+1)
        for j in range(n//2+1):
            literal[n-2*j] = (-1)**j*factorial(n)//(2**j*factorial(j)*factorial(n-2*j))
        require(hermite[n] == literal, "Hermite recurrence versus generating-function coefficient")
    for k in range(1, 33):
        scaled = [value*2**(j//2) if j % 2 == 0 else 0 for j, value in enumerate(hermite[2*k])]
        require(all(value == 0 for j, value in enumerate(hermite[2*k]) if j % 2), "even Hermite scaling has no irrational coefficients")
        rhs = [2**k*value for value in multiply(hermite[k], hermite[k])]
        for s in range(1, k+1):
            odd_factorial = prod(range(1, 2*s-2, 2))
            positive = comb(k, s)*2**(k-s)*odd_factorial
            require(positive > 0, "every Hermite negative-square coefficient is positive")
            for j, value in enumerate(multiply(hermite[k-s], hermite[k-s])):
                rhs[j] -= positive*value
        require(scaled == rhs, "full scaled Hermite duplication SOS coefficient identity")
        require(len(polynomial_gcd(hermite[k], scaled)) == 1, "scaled Hermite root gcd independent pseudoremainder")
        require(prod(range(1, 2*k-2, 2)) > 0, "strict nonzero constant in root-sign certificate")
        record("Hermite", k, scaled)
    print("Hermite scaled duplication: k=1..32 exact polynomial SOS/root-sign certificates and constant gcds; recurrence checked through64")


def main():
    print("SIGNED DUPLICATION SOS / REAL-ROOTED LAURENT CORE / ALL-MASS AP")
    print("status=FINITE-EXACT CHECKER of independently audited analytic proof; no repository proof imports")
    integral_controls()
    gram_controls()
    signed_controls()
    direct_sos_controls()
    zero_coefficient_controls()
    ap_controls()
    core_controls()
    hostile_controls()
    hermite_controls()
    print(f"trace_sha256={TRACE.hexdigest()}")
    print(f"PASS explicit_gates={CHECKS}")


if __name__ == "__main__":
    main()
