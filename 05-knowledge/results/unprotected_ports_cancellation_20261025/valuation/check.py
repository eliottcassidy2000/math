"""Independent arithmetic controls for the characteristic spectrum; stdlib only."""
from fractions import Fraction
from math import comb
import json


def need(value, message):
    if not value:
        raise RuntimeError(message)


def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    c[i + j] += x * y
    return c


def power(a, n):
    b = [1]
    while n:
        if n & 1:
            b = mul(b, a)
        n >>= 1
        if n:
            a = mul(a, a)
    return b


def vp(n, p):
    need(n != 0, "valuation of zero")
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def depth(a, ell=0):
    return next(i for i, x in enumerate(a) if x if not ell or x % ell)


def factors(n):
    ans = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            ans.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        ans.append(n)
    return ans


def check(p, r):
    q, h = p**r, p**(r - 1)
    P = []
    for k in range(q + 1):
        numerator = comb(q, k)
        if k % p == 0:
            numerator += (q - 1) * comb(h, k // p)
        need(numerator % q == 0, "nonintegral subset-orbit count")
        P.append(numerator // q)
    shifted = [-sum((-1)**k * comb(k, j) * P[k]
                    for k in range(j, q + 1)) for j in range(q + 1)]
    powered = power([0] + [comb(p, i) for i in range(1, p + 1)], h)
    numerator = [(q - 1) * x for x in powered]
    numerator[q] += 1
    need(all(x % q == 0 for x in numerator), "nonintegral positive formula")
    positive = [x // q for x in numerator]
    need(positive == shifted, "independent polynomial paths disagree")
    need(all(x >= 0 for x in shifted), "negative shifted coefficient")
    need(depth(shifted) == h, "wrong rational depth")
    need(shifted[h] == (q - 1) * p**(h - r), "wrong leading coefficient")
    C = [0] + [comb(p, i) // p for i in range(1, p)]
    for i in range(1, p):
        need(comb(p, i) % p == 0, "C not integral")
        need(C[i] % p == (-1)**(i - 1) * pow(i, -1, p) % p,
             "finite logarithm mismatch")
    residue = [0] * (q + 1)
    residue[q] = 1
    for i in range(1, p):
        residue[q - p + i] = -C[i] % p
    need([x % p for x in shifted] == residue, "wrong physical residue")
    for j in range(1, h + 1):
        need(vp(comb(h, j), p) == r - 1 - vp(j, p), "binomial valuation")
        a = Fraction(comb(h, j) * p**j, q)
        need(a.denominator == 1, "expansion scalar not integral")
        need(vp(a.numerator, p) == j - 1 - vp(j, p), "scalar valuation")
        need(j == 1 or a.numerator % p == 0, "higher term survives mod p")
    primes = sorted(set([2, 3, 5, 7, 11, 13, 17, 19, p] + factors(q - 1)))
    spectrum = {}
    for ell in primes:
        expected = q - p + 1 if ell == p else q if (q - 1) % ell == 0 else h
        actual = depth(shifted, ell)
        need(actual == expected, f"spectrum p={p} r={r} ell={ell}")
        if (q - 1) % ell == 0:
            need([x % ell for x in shifted] == [0] * q + [1], "pure monomial")
        spectrum[str(ell)] = actual
    return {"p": p, "r": r, "q": q, "rational_depth": h,
            "physical_depth": q - p + 1, "leading_coefficient": shifted[h],
            "checked_characteristics": spectrum}


def main():
    universe = [(p, r) for p, bound in [(3, 4), (5, 4), (7, 3), (11, 2),
                                       (13, 2), (19, 1), (23, 1)]
                for r in range(1, bound + 1)]
    rows = [check(p, r) for p, r in universe]
    cyclic = []
    for p, bound in [(3, 4), (5, 4), (7, 3), (11, 2), (13, 2), (19, 1), (23, 1)]:
        for r in range(1, bound + 1):
            q = p**r
            numerator = [0] * (q + 1)
            for a in range(r + 1):
                d = p**a
                phi = 1 if a == 0 else (p - 1) * p**(a - 1)
                term = power([0] + [comb(d, j) for j in range(1, d + 1)], q // d)
                for j, x in enumerate(term):
                    numerator[j] += phi * x
            need(all(x % q == 0 for x in numerator), "cyclic integrality")
            poly = [x // q for x in numerator]
            direct_numerators = [sum((1 if a == 0 else (p - 1) * p**(a - 1))
                            * comb(q // p**a, j // p**a)
                            for a in range(r + 1) if j % p**a == 0)
                        for j in range(q + 1)]
            need(all(x % q == 0 for x in direct_numerators), "cyclic orbit division")
            direct_P = [x // q for x in direct_numerators]
            literal_shift = [-sum((-1)**k * comb(k, j) * direct_P[k]
                                  for k in range(j, q + 1)) for j in range(q + 1)]
            need(poly == literal_shift, "cyclic independent shift")
            need(depth(poly) == 1, "cyclic rational depth")
            actual = depth(poly, p)
            need(actual == q // p, "cyclic modular depth")
            carry_residue = [0] * (q + 1)
            carry_residue[q] = 1
            for a in range(r):
                for j in range(1, p):
                    degree = q - p**(a + 1) + p**a * j
                    carry_residue[degree] -= comb(p, j) // p
            need([x % p for x in poly] == [x % p for x in carry_residue],
                 "cyclic all-scale modular formula")
            need(sum(x % p != 0 for x in poly) == 1 + r * (p - 1),
                 "cyclic scale-sparse monomial count")
            elementary_numerator = [(q - 1) * x for x in
                power([0] + [comb(p, i) for i in range(1, p + 1)], q // p)]
            elementary_numerator[q] += 1
            elementary = [x // q for x in elementary_numerator]
            need(all(j * (poly[j] - elementary[j]) % p == 0
                     for j in range(1, q + 1)), "common modular derivative")
            need(all((poly[j] - elementary[j]) % p == 0
                     for j in range(q + 1) if j % p or j in (0, q)),
                 "Frobenius difference support")
            cyclic.append({"p": p, "r": r, "q": q, "physical_depth": actual,
                           "nonzero_modular_monomials": 1 + r * (p - 1)})
    print(json.dumps({"status": "PASS", "cases": len(rows), "max_q": 625,
                      "universe": rows, "cyclic_finite_cases": cyclic}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
