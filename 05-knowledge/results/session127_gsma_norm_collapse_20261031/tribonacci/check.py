#!/usr/bin/env python3
"""Exact audit of the Tribonacci 305 and degree-3/4/5 CF/norm bridge.

Standard library only. Checks remain active under python -O.
Run: python3 check.py > receipt.json; python3 -O check.py > receipt.optimized.json
"""
from fractions import Fraction as Q
from itertools import permutations, product
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def form(k, p, q):
    return p**k - sum(p**j * q**(k-j) for j in range(k))


def poly(k, x):
    return x**k - sum(x**j for j in range(k))


def irreducible_mod_prime(k, prime):
    """Exhaust all monic proper possible factors of degree <= floor(k/2)."""
    source = [prime-1]*k+[1]
    for degree in range(1, k//2+1):
        for coeffs in product(range(prime), repeat=degree):
            divisor = list(coeffs)+[1]
            rem = source[:]
            for index in range(k, degree-1, -1):
                lead = rem[index]
                for j in range(degree+1):
                    rem[index-degree+j] = (rem[index-degree+j]-lead*divisor[j]) % prime
            if not any(rem):
                return False
    return True


def derivative(k, x):
    return k*x**(k-1) - sum(j*x**(j-1) for j in range(1, k))


def isolate(k, bits=512):
    lo, hi = Q(3, 2), Q(2)
    require(poly(k, lo) < 0 < poly(k, hi), "root bracket")
    for _ in range(bits):
        mid = (lo + hi)/2
        if poly(k, mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def cf_interval(lo, hi, count):
    result = []
    for _ in range(count):
        a, b = lo.numerator//lo.denominator, hi.numerator//hi.denominator
        require(a == b, "insufficient interval precision")
        result.append(a)
        lo, hi = 1/(hi-a), 1/(lo-a)
    return result


def cf_sign(k, count):
    """Independent CF extraction by integer signs of F_k at Mobius bounds.

    f_k has exactly one positive root (Descartes), so for all x>0 the
    sign of f_k(x) is the sign of x-alpha. No floating root enters.
    """
    digits = [1]
    p, q, r, s = 1, 1, 1, 0
    while len(digits) < count:
        orientation = p*s-r*q
        require(abs(orientation) == 1, "unimodular CF frame")
        def above(t):
            return orientation * form(k, p*t+r, q*t+s) > 0
        lo, hi = 1, 2
        require(not above(lo), "complete quotient >1")
        while not above(hi):
            hi *= 2
        while hi-lo > 1:
            mid = (lo+hi)//2
            if above(mid):
                hi = mid
            else:
                lo = mid
        digits.append(lo)
        p, q, r, s = p*lo+r, q*lo+s, p, q
    return digits


def convergents(digits):
    p, q, r, s = 1, 0, 0, 1
    out = []
    for a in digits:
        p, q, r, s = a*p+r, a*q+s, p, q
        out.append((p, q, r, s))
    return out


def determinant(m):
    n = len(m)
    result = 0
    for perm in permutations(range(n)):
        inv = sum(perm[i] > perm[j] for i in range(n) for j in range(i+1,n))
        prod = (-1)**inv
        for i in range(n):
            prod *= m[i][perm[i]]
        result += prod
    return result


def matrix_norm(k, p, q):
    # Multiplication by alpha in basis 1,alpha,...,alpha^(k-1).
    c = [[0]*k for _ in range(k)]
    for j in range(k-1):
        c[j+1][j] = 1
    for i in range(k):
        c[i][-1] = 1
    return determinant([[p*(i == j)-q*c[i][j] for j in range(k)] for i in range(k)])


def product_mod(k, a, b):
    v = [0]*(2*k-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            v[i+j] += x*y
    for i in range(len(v)-1,k-1,-1):
        coeff = v[i]
        v[i] = 0
        for j in range(k):
            v[i-k+j] += coeff
    return v[:k]


def power_mod(k, exponent):
    require(exponent >= 0, "positive power input")
    a, b = [1]+[0]*(k-1), [0,1]+[0]*(k-2)
    for _ in range(exponent):
        a = product_mod(k,a,b)
    return a


def inverse_power_mod(k, exponent):
    a = [1]+[0]*(k-1)
    inverse = [-1]*(k-1)+[1]
    for _ in range(exponent):
        a = product_mod(k,a,inverse)
    return a


def secant_interval(k, x, lo, hi):
    """Outward rational bounds for the divided difference D_k(x,alpha)."""
    lower = sum(x**(k-1-h)*lo**h for h in range(k))
    upper = sum(x**(k-1-h)*hi**h for h in range(k))
    for j in range(1,k):
        lower -= sum(x**(j-1-h)*hi**h for h in range(j))
        upper -= sum(x**(j-1-h)*lo**h for h in range(j))
    return lower, upper


def interval_json(lo, hi):
    # Export compact outward rational enclosures, not enormous intermediate
    # fractions. Exact fine bounds above still perform every logical check.
    scale = 10**12
    a, b = lo*scale, hi*scale
    lower = Q(a.numerator//a.denominator, scale)
    upper = Q(-(-b.numerator//b.denominator), scale)
    require(lower <= lo <= hi <= upper, "outward output rounding")
    return {"lower": str(lower), "upper": str(upper)}


def main():
    report = {"status":"PASS", "arithmetic":"integer/Fraction; floats display-only",
              "cf_terms":40, "degrees":{}}
    for k in (3,4,5):
        prime = {3:3,4:2,5:5}[k]
        require(irreducible_mod_prime(k,prime), "irreducibility certificate")
        lo, hi = isolate(k)
        digits = cf_interval(lo, hi, 40)
        require(digits == cf_sign(k,40), "independent CF path disagreement")
        convs = convergents(digits)
        records = []
        for index, (p,q,r,s) in enumerate(convs[:-1]):
            n = form(k,p,q)
            require(n == matrix_norm(k,p,q), "norm/determinant disagreement")
            dl, du = secant_interval(k,Q(p,q),lo,hi)
            require(dl > 0, "positive secant")
            tl = q**(k-2)*dl/abs(n)-Q(s,q)
            tu = q**(k-2)*du/abs(n)-Q(s,q)
            nxt = digits[index+1]
            require(nxt <= tl <= tu < nxt+1, "tail quotient bound")
            if index <= 5 or nxt >= 20:
                records.append({"index":index,"p":p,"q":q,"previous_p":r,
                                "previous_q":s,"norm":n,"next_digit":nxt,
                                "tail_bounds":interval_json(tl,tu)})
        require(inverse_power_mod(k,k) == [2,-1]+[0]*(k-2), "uniform sparse unit")
        # Universal unit-norm control; its next CF digit was checked separately.
        require(form(k,2,1) == 1, "universal norm-one control")
        report["degrees"][str(k)] = {"cf":digits, "selected_convergents":records,
            "uniform_sparse_unit_exponent":k,"irreducibility_prime":prime}

    require(inverse_power_mod(3,16) == [-103,56,0], "exceptional cubic unit")
    require(product_mod(3,power_mod(3,16),[-103,56,0]) == [1,0,0],
            "independent positive-power identity")
    require(form(3,103,56) == -1, "exceptional norm")
    require(form(3,31461,17105) == 83326, "305 upper boundary")
    require(form(3,31564,17161) == -10037, "306 lower boundary")
    lo, hi = isolate(3,48)
    tail_lo, tail_hi = (lo**16-25)/56, (hi**16-25)/56
    require(305 < tail_lo < tail_hi < 306, "305 tail proof")
    report["exceptional_unit"] = {"inverse_power":16,"coefficients":[-103,56,0],
        "norm_of_103_minus_56tau":-1,"tail": "(tau^16-25)/56",
        "tail_bounds": interval_json(tail_lo,tail_hi),
        "boundary_forms":{"31461/17105":83326,"31564/17161":-10037}}

    legacy = [1,1,5,4,2,305,1,8,2,1,4,6,17,5,1,5,4,6,3,7]
    exact = report["degrees"]["3"]["cf"]
    first = next(i for i,(a,b) in enumerate(zip(legacy,exact)) if a != b)
    require(first == 12 and exact[first] == 14, "legacy hostile mismatch")
    report["historical_correction"] = {"first_wrong_zero_based_index":first,
        "legacy_digit":legacy[first],"correct_digit":exact[first],
        "exact_20_digits":exact[:20]}

    plus_one = [(0,-1),(1,0),(2,1),(-103,-56)]
    require(all(form(3,p,q) == matrix_norm(3,p,q) == 1 for p,q in plus_one),
            "four cited global Thue solutions")
    # The discriminant of x^3+a*x^2+b*x+c.
    a=b=c=-1
    disc = a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c
    require(disc == -44, "cubic discriminant")
    report["cubic_thue"] = {"discriminant":disc,"norm_plus_one_solutions":plus_one,
        "completeness":"CITED Delone-Nagell via Knapp-Lemos-Neumann Lemma2.4, Example4.1; not established by finite enumeration",
        "positive_sparse_inverse_power_exponents":[3,16]}
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
