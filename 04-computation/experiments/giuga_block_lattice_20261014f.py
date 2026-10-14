"""Exact signed two-prime completion and retained star-lattice controls.

No cutoff on the two new primes: an index interval and divisor identity
make the universe finite. Trial division is exact, with no probabilistic
primality input. The mod-eight theorem is a cited theorem, not proved by
the finite controls here. All checks remain active with python -O.
"""
from fractions import Fraction
from itertools import combinations, product
from math import gcd, isqrt, prod
import json

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(label)


def exact_int(n, lower=1):
    if type(n) is not int or n < lower:
        raise ValueError("integer domain")
    return n


def sign(mu):
    if type(mu) is not int or mu not in (-1, 1):
        raise ValueError("sign domain")
    return mu


def prime(n):
    return type(n) is int and n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def factors(n):
    exact_int(n)
    out = []
    p = 2
    while p*p <= n:
        e = 0
        while n % p == 0:
            e += 1
            n //= p
        if e:
            out.append((p, e))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


def squarefree_support(n):
    exact_int(n, 2)
    fs = factors(n)
    if any(e != 1 for _, e in fs):
        raise ValueError("core must be squarefree")
    return tuple(p for p, _ in fs)


def divisors(n):
    result = [1]
    for p, e in factors(n):
        result = [a*p**b for a in result for b in range(e+1)]
    return sorted(result)


def first_allowed(n, count):
    result = []
    p = 2
    while len(result) < count:
        if n % p and prime(p):
            result.append(p)
        p += 1
    return result


def signed_index_from_support(ps, mu):
    return sum((Fraction(1, p) for p in ps), Fraction()) + Fraction(mu, prod(ps))


def compile_two(n, mu):
    ps = squarefree_support(n)
    sign(mu)
    S = sum((Fraction(1, p) for p in ps), Fraction())
    A0 = sum(n//p for p in ps)
    q0, r0 = first_allowed(n, 2)
    upper = S+Fraction(1, q0)+Fraction(1, r0)+Fraction(mu, n*q0*r0)
    answers = []
    carriers = []
    for k in range(S.numerator//S.denominator+1, upper.numerator//upper.denominator+1):
        a = k*n-A0
        D = n*n+mu*a
        if a <= 0 or D <= 0:
            raise ValueError("positive carrier theorem failed")
        carriers.append([k, a, D])
        for u in divisors(D):
            v = D//u
            if u >= v:
                break
            if (u+n) % a or (v+n) % a:
                continue
            q, r = (u+n)//a, (v+n)//a
            if q < r and prime(q) and prime(r) and n % q and n % r:
                answers.append({"index": k, "primes": [q, r], "number": n*q*r,
                                "factors": [u, v]})
    return {"core": n, "mu": mu, "core_prime_support": list(ps),
            "index_interval": [str(S), str(upper)],
            "carriers_k_a_D": carriers, "completions": answers}


def native_signed(ps, mu):
    n = prod(ps)
    return all((n//p+mu) % p == 0 for p in ps)


def native_korselt(ps):
    n = prod(ps)
    return len(ps) > 1 and all((n-1) % (p-1) == 0 for p in ps)


def independent_small_controls():
    bank = [p for p in range(2, 32) if prime(p)]
    total_cores = 0
    for n in range(2, 151):
        fs = factors(n)
        if any(e != 1 for _, e in fs):
            continue
        total_cores += 1
        ps = tuple(p for p, _ in fs)
        for mu in (-1, 1):
            result = compile_two(n, mu)
            seen = {tuple(row['primes']) for row in result['completions']}
            for q, r in combinations(bank, 2):
                if n % q == 0 or n % r == 0:
                    continue
                target = ps+(q, r)
                check(((q, r) in seen) == native_signed(target, mu),
                      "independent local-prime criterion versus compiler")
            for row in result['completions']:
                q, r = row['primes']
                target = ps+(q, r)
                check(native_signed(target, mu), "all-height completion native replay")
                check(signed_index_from_support(target, mu) == row['index'],
                      "independent reciprocal index")
    return {"squarefree_cores": total_cores, "core_interval": [2, 150],
            "literal_new_prime_bound": 31,
            "compiler_prime_bound": "none; exact finite divisor universe"}


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    det = Fraction(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return Fraction(0)
        if j != i:
            a[i], a[j] = a[j], a[i]
            det = -det
        pivot = a[i][i]
        det *= pivot
        for j in range(i+1, len(a)):
            scale = a[j][i]/pivot
            for k in range(i+1, len(a)):
                a[j][k] -= scale*a[i][k]
    return det


def lattice(ps, mu):
    if not ps or tuple(sorted(set(ps))) != tuple(ps) or not all(prime(p) for p in ps):
        raise ValueError("distinct sorted prime support required")
    sign(mu)
    n = prod(ps)
    k = signed_index_from_support(ps, mu)
    if k.denominator != 1:
        raise ValueError("integer central coefficient required")
    r = len(ps)
    B = [[int(k)]+[1]*r]
    B.extend([[1]+[p if i == j else 0 for j in range(r)] for i, p in enumerate(ps)])
    check(determinant(B) == mu, "independent star determinant")
    schur = Fraction(int(k))-sum((Fraction(1, p) for p in ps), Fraction())
    check(schur == Fraction(mu, n), "exact Schur complement and inertia")
    chars = [w for w in product((0, 1), repeat=r+1)
             if all((sum(B[i][j]*w[j] for j in range(r+1))-B[i][i]) % 2 == 0
                    for i in range(r+1))]
    check(len(chars) == 1, "one characteristic residue class")
    w = chars[0]
    norm = sum(w[i]*B[i][j]*w[j] for i in range(r+1) for j in range(r+1))
    check((norm-r-mu) % 8 == 0, "cited van der Blij consistency")
    if n % 2:
        check(w == (1,)+(0,)*r and (r+mu-int(k)) % 8 == 0,
              "odd star exact characteristic vector")
    return {"number": n, "support": list(ps), "mu": mu, "index": int(k),
            "determinant": mu, "signature": r+mu,
            "characteristic_vector": list(w), "characteristic_norm": norm}


def block_bound(n):
    ps = squarefree_support(n)
    k = signed_index_from_support(ps, -1)
    if k.denominator != 1 or len(ps) < 2:
        raise ValueError("Giuga root required")
    s = Fraction(1, n)
    previous = s
    t = 0
    p = 2
    while s <= 1:
        if n % p and prime(p):
            previous = s
            s += Fraction(1, p)
            t += 1
        p += 1
    check(previous <= 1 < s, "strict reciprocal block threshold")
    return {"root": n, "minimum_new_primes": t, "last_prime": p-1,
            "previous_upper_bound": str(previous), "first_upper_bound_above_one": str(s)}


def main():
    controls = independent_small_controls()
    roots = (2, 6, 30, 42, 561, 858, 1105, 1722, 1729, 1806, 47058)
    rows = [compile_two(n, mu) for n in roots for mu in (-1, 1)]
    row561 = next(r for r in rows if r['core'] == 561 and r['mu'] == -1)
    check([r['primes'] for r in row561['completions']] == [[2, 59]],
          "complete 561 Giuga two-prime catalog")
    check(native_korselt((3, 11, 17)) and not native_signed((3, 11, 17), -1),
          "Carmichael-only source")
    check(native_signed((2, 3, 11, 17, 59), -1)
          and not native_korselt((2, 3, 11, 17, 59)), "Giuga-only child")
    check(native_korselt((3, 11, 17, 401, 3041))
          and not native_signed((3, 11, 17, 401, 3041), -1),
          "inherited Carmichael two-prime repair remains separate")
    supports = [((2,), 1), ((2, 3), 1), ((2, 3, 7), 1), ((2, 3, 7, 43), 1),
                ((2, 3, 11, 23, 31), 1), ((2, 3, 5), -1),
                ((2, 3, 11, 13), -1), ((2, 3, 7, 41), -1),
                ((2, 3, 11, 17, 59), -1)]
    supports += [((p,), -1) for p in range(3, 102) if prime(p)]
    lattice_rows = [lattice(ps, mu) for ps, mu in supports]
    check(determinant([[1, 1, 1], [1, 3, 0], [1, 0, 5]]) == 7,
          "nonunimodular hostile determinant")
    check((1-3) % 8 != 0, "nonunimodular hostile signature mismatch")
    for ps, mu in [((3, 5), -1), ((2, 2), 1), ((2, 9), 1)]:
        try:
            lattice(ps, mu)
        except ValueError:
            check(True, "malformed lattice source rejected")
        else:
            check(False, "malformed lattice source accepted")
    for n, mu in [(True, 1), (4, -1), (6, True), (6, 0)]:
        try:
            compile_two(n, mu)
        except ValueError:
            check(True, "malformed completion source rejected")
        else:
            check(False, "malformed completion source accepted")
    bounds = [block_bound(n) for n in (30, 858, 1722, 66198, 2214408306)]
    check([r['minimum_new_primes'] for r in bounds] == [56, 55, 55, 54, 54],
          "source-specific block bounds")
    print(json.dumps({"status": "PROVED scoped algebra; CITED lattice theorem; FINITE-EXACT controls; Giuga/Collatz OPEN",
                      "controls": controls, "two_prime_catalog": rows,
                      "lattice_controls": lattice_rows, "proper_Giuga_multiple_bounds": bounds,
                      "checks": CHECKS}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
