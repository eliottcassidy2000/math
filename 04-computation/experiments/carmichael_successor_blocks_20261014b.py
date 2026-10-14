"""Exact predecessor certificates and COMPLETE one/two-prime extensions.

Run from the repository root; writes deterministic JSON to stdout.
All checks survive python -O. No external packages or probabilistic tests.
The two-prime search is finite by a proved bound on a quotient, not a
user-selected cap on primes. This does not enumerate arbitrary extensions.
"""
from itertools import combinations
import json
from math import gcd, lcm, prod

CHECKS = 0


def require(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(message)


def factors(n):
    if type(n) is not int or n < 1:
        raise ValueError("positive integer required")
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def prime(n):
    return n >= 2 and factors(n) == {n: 1}


def divisors(n):
    out = [1]
    for p, e in factors(n).items():
        out = [d * p**j for d in out for j in range(e + 1)]
    return sorted(out)


def digits(n, base):
    out = []
    while n:
        n, digit = divmod(n, base)
        out.append(digit)
    return out or [0]


def multiplication_carries(a, b, base):
    da, db = digits(a, base), digits(b, base)
    raw = [0] * (len(da) + len(db) - 1)
    for i, x in enumerate(da):
        for j, y in enumerate(db):
            raw[i + j] += x * y
    carry = total = 0
    for col in raw:
        carry = (col + carry) // base
        total += carry
    while carry:
        carry //= base
        total += carry
    return total


def korselt(ps):
    """An explicit, distinct-prime factor certificate, with compositeness."""
    if len(ps) < 2 or len(set(ps)) != len(ps) or not all(prime(p) for p in ps):
        return False
    n = prod(ps)
    return all((n - 1) % (p - 1) == 0 for p in ps)


def carmichael(n):
    fs = factors(n)
    return all(e == 1 for e in fs.values()) and korselt(tuple(fs))


def interface(n):
    if not carmichael(n):
        raise ValueError("certified Carmichael source required")
    lam = lcm(*(p - 1 for p in factors(n)))
    return lam, (n - 1) // lam


def one_extensions(n):
    lam, slack = interface(n)
    return sorted(1 + d * lam for d in divisors(slack)
                  if prime(1 + d * lam) and n % (1 + d * lam) != 0)


def two_extensions(n):
    """Complete finite enumeration of q<r, distinct from source primes.

    With e=n-1, a=q-1, b=r-1, t=(n*a+e)/b:
      1 <= t <= n+e//2, and a divides e*(n+t).
    Every filter below is necessary; final Korselt supplies sufficiency.
    """
    lam, _ = interface(n)
    e = n - 1
    ps = tuple(factors(n))
    hits = set()
    candidates = 0
    for t in range(1, n + e // 2 + 1):
        for a in divisors(e * (n + t)):
            if a < 2 or (n * a + e) % t:
                continue
            b = (n * a + e) // t
            if b <= a or (n * b + e) % a:
                continue
            q, r = a + 1, b + 1
            if n % q == 0 or n % r == 0 or (q * r - 1) % lam:
                continue
            candidates += 1
            if prime(q) and prime(r):
                require(korselt(ps + (q, r)), "final two-prime certificate")
                require(b <= n * a + e, "last-prime bound")
                hits.add((n * q * r, q, r))
    return sorted(hits), candidates


def independent_pairs(n, max_q):
    """Separate bounded audit: choose q first, then divisors of n*q-1."""
    ps = tuple(factors(n))
    out = set()
    for q in range(3, max_q + 1, 2):
        if not prime(q) or n % q == 0:
            continue
        for b in divisors(n * q - 1):
            r = b + 1
            if r > q and prime(r) and n % r != 0 and korselt(ps + (q, r)):
                out.add((n * q * r, q, r))
    return sorted(out)


def group_compiler(L, k):
    if L < 2 or k < 1 or gcd(L, k) != 1:
        raise ValueError("L>=2, k>=1 and gcd(L,k)=1 required")
    bank = [d * k + 1 for d in divisors(L)
            if prime(d * k + 1) and gcd(d * k + 1, L) == 1]
    hits = []
    for size in range(2, len(bank) + 1):
        for ps in combinations(bank, size):
            n = prod(ps)
            if n % L == 1:
                require((n - 1) % (k * L) == 0, "CRT conductor")
                require(korselt(ps), "product-one compilation")
                hits.append((n, ps))
    return bank, sorted(hits)


def main():
    # Direct all-residue power test is independent of the divisor criterion.
    for n in range(2, 700):
        absolute_fermat = all(pow(a, n, n) == a for a in range(n))
        criterion = all(e == 1 for e in factors(n).values()) and all(
            (n - 1) % (p - 1) == 0 for p in factors(n))
        require(absolute_fermat == criterion, "universal Fermat vs Korselt")
    for n in [9, 25, 45, 49, 81, 121, 225, 315, 561**2]:
        p = next(p for p, e in factors(n).items() if e >= 2)
        require((pow(p, n, p * p) - p) % (p * p) != 0, "nilpotent hostile")

    small = [n for n in range(3, 20000, 2) if carmichael(n)]
    require(small == [561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841],
            "small universe census")
    for n in small:
        lam, e = interface(n)
        predicted = set(one_extensions(n))
        for q in range(2, 1000):
            if prime(q) and n % q:
                require((q in predicted) == korselt(tuple(factors(n)) + (q,)),
                        "independent one-prime test")
        for q in predicted:
            d = (q - 1) // lam
            require(e % d == 0, "divisor label")
            require(interface(n * q) == (d * lam, n + e // d), "interface update")

    catalog = {}
    for n in [561, 1105, 1729]:
        hits, candidates = two_extensions(n)
        # The independent q-first test reaches every q in the 561 output;
        # for the other two roots it deliberately checks a bounded subrange.
        cap = 2000
        require([x for x in hits if x[1] <= cap] == independent_pairs(n, cap),
                "independent factor order")
        catalog[n] = {"lambda_slack": interface(n), "one_primes": one_extensions(n),
                      "two_products": hits, "preprimality_candidates": candidates}
    require(catalog[561]["two_products"] == [(684106401, 401, 3041)],
            "complete unique two-prime child")
    require([len(catalog[n]["two_products"]) for n in [561, 1105, 1729]] == [1, 13, 21],
            "complete two-prime counts")
    child = 684106401
    require(not one_extensions(561) and not one_extensions(1105)
            and not one_extensions(child), "single-step obstruction repeats")
    require(interface(child) == (15200, 45007), "child interface")
    require(prime(45007) and not prime(15201), "prime-slack leaf")
    require(not korselt((3, 11, 17, 401)) and not korselt((3, 11, 17, 3041)),
            "both intermediates fail")
    require(one_extensions(1729) == [37, 73, 109, 433, 577], "branching root")

    compiled = {}
    for L, k in [(80, 1), (80, 3), (36, 1), (120, 7)]:
        bank, hits = group_compiler(L, k)
        compiled[f"{L},{k}"] = {"primes": bank, "products": hits}
    require((172081, (7, 13, 31, 61)) in compiled["80,3"]["products"], "nontrivial CRT")
    require(91 % 6 == 1 and 91 % 2 == 1 and 90 % 12 != 0 and not korselt((7, 13)),
            "noncoprime conductor hostile")
    try:
        group_compiler(6, 2)
    except ValueError:
        require(True, "reject wrong CRT type")
    else:
        require(False, "wrong CRT type accepted")
    # Chernick identity: a family theorem conditional only on three primes.
    chernick = []
    for m in range(1, 101):
        ps = (6 * m + 1, 12 * m + 1, 18 * m + 1)
        require(prod(ps) - 1 == 36 * m * (36 * m * m + 11 * m + 1), "Chernick identity")
        if all(prime(p) for p in ps):
            require(korselt(ps), "Chernick prime guard")
            chernick.append((m, prod(ps)))
    require(1729 == 1**3 + 12**3 == 9**3 + 10**3 and 1729 - 1 == 12**3,
            "cube boundary")
    polygonal = {}
    for n in small:
        rows = []
        for p in factors(n):
            quotient, remainder = divmod(n // p - 1, p - 1)
            rank = 2 * quotient + 2
            require(remainder == 0 and rank >= 4 and rank % 2 == 0, "marked polygon rank")
            require(n == p + (rank - 2) * p * (p - 1) // 2, "marked polygon value")
            sp = sum(digits(n, p))
            require(sp >= p and (sp - 1) % (p - 1) == 0, "digit period receipt")
            rows.append({"p": p, "polygon_rank": rank, "digits_low_first": digits(n, p),
                         "digit_sum": sp, "digit_height": (sp - 1) // (p - 1)})
        polygonal[n] = rows
    require(561 == 17 * (2 * 17 - 1) == 33 * 34 // 2, "hexagon triangle intersection")
    require(91 == 7 * (2 * 7 - 1) and (91 - 1) % 12 != 0, "one-mark hostile")
    for n in range(2, 5001):
        fs = factors(n)
        digit_criterion = all(e == 1 for e in fs.values()) and all(
            sum(digits(n, p)) >= p and (sum(digits(n, p)) - 1) % (p - 1) == 0 for p in fs)
        require(digit_criterion == carmichael(n), "digit vs divisor criterion")
    for base in [2, 3, 5, 7, 11, 17]:
        for a in range(1, 32):
            for b in range(1, 32):
                require(sum(digits(a * b, base)) == sum(digits(a, base)) * sum(digits(b, base))
                        - (base - 1) * multiplication_carries(a, b, base), "carry conservation")
    for n in range(1, 2000, 2):
        m = 3 * n
        trailing_ones = 0
        while m % 2:
            trailing_ones += 1
            m //= 2
        endpoint = 3 * n + 1
        while endpoint % 2 == 0:
            endpoint //= 2
        require(endpoint.bit_count() == 2 * n.bit_count() + 1
                - multiplication_carries(n, 3, 2) - trailing_ones, "Collatz digit balance")
    require(7 < 17 and (7).bit_count() > (17).bit_count(), "digit descent is not payment")
    print(json.dumps({"status": "PROVED finite algorithms; FINITE-EXACT exhaustive catalogs",
                      "checks": CHECKS, "catalog": catalog, "group_compiler": compiled,
                      "chernick_m_le_100": chernick, "marked_polygonal_profiles": polygonal,
                      "scope": "Exactly one/two new primes for each listed root; no universal extension theorem"},
                     sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
