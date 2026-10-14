"""Full predecessor guards versus marked-prime supply, with literal Collatz receipts.

Production APIs never discover a ROOT suffix or a prime in a progression.
The bounded prime search in main() is an explicitly finite control only.
"""
from dataclasses import dataclass, replace
from itertools import product
from math import gcd, isqrt, lcm, prod
import hashlib
import json

import collatz_successor_holonomy_20261014b as native


def integer(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError("exact integer outside domain")
    return n


def prime(n):
    integer(n)
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n)+1, 2))


def v2(n):
    integer(n, 1)
    return (n & -n).bit_length()-1


@dataclass(frozen=True)
class KorseltSource:
    source: int
    primes: tuple


def audit_source(s):
    if type(s) is not KorseltSource:
        raise ValueError("exact source certificate required")
    native.odd(s.source)
    if type(s.primes) is not tuple or len(s.primes) < 2:
        raise ValueError("at least two distinct prime factors required")
    for p in s.primes:
        integer(p, 3)
        if not prime(p):
            raise ValueError("composite factor")
    if tuple(sorted(set(s.primes))) != s.primes:
        raise ValueError("canonical distinct factors required")
    if prod(s.primes) != s.source:
        raise ValueError("unchanged source factorization mismatch")
    if any((s.source-1) % (p-1) for p in s.primes):
        raise ValueError("full predecessor guard failed")
    return s


def predecessor_state(s):
    audit_source(s)
    conductor = lcm(*(p-1 for p in s.primes))
    slack = (s.source-1)//conductor
    return conductor, slack, v2(conductor)+v2(slack)


def descent(s):
    """Return an actual paid prefix, or None exactly for the first-rise branch."""
    _, _, height = predecessor_state(s)
    if height == 1:
        return None
    w = ((v2(3*s.source+1),) if height == 2
         else (2,)*((height-1)//2))
    r = native.receipt(s.source, w)
    if not 1 <= r.target < r.source:
        raise ValueError("declared descent did not pay")
    return r


def discharge(s, supplied_child_root_word):
    r = descent(s)
    if r is None:
        raise ValueError("the selected predecessor rule has no paid exit")
    return native.discharge(r, supplied_child_root_word)


@dataclass(frozen=True)
class PartialPlan:
    word: tuple
    prime_residue: int
    prime_modulus: int
    minimum_prime: int


def compile_partial(w):
    native.word(w)
    if not w:
        raise ValueError("nonempty native word required")
    r, modulus = native.native_cell(w)
    qdyadic = r*pow(33, -1, modulus) % modulus
    # CRT: q=257+330*j; divide the compatible even congruence by two.
    half = modulus//2
    j = ((qdyadic-257)//2)*pow(165, -1, half) % half
    period = 165*modulus
    residue = (257+330*j) % period
    lower = max(257, (half+32)//33)
    return PartialPlan(w, residue, period, lower)


def audit_plan(plan):
    if type(plan) is not PartialPlan:
        raise ValueError("exact partial plan required")
    native.word(plan.word)
    integer(plan.prime_residue)
    integer(plan.prime_modulus, 1)
    integer(plan.minimum_prime, 1)
    if compile_partial(plan.word) != plan:
        raise ValueError("noncanonical partial plan")
    return plan


def marked_guards(source, p):
    native.odd(source)
    integer(p, 3)
    if not prime(p) or source % p:
        raise ValueError("supplied marked prime must divide the source")
    return ((source//p-1) % p == 0, (source-1) % (p-1) == 0)


def realize_partial(plan, q):
    """Authenticate a supplied prime and the actual word; no full completion claim."""
    audit_plan(plan)
    integer(q, plan.minimum_prime)
    if not prime(q) or q % plan.prime_modulus != plan.prime_residue:
        raise ValueError("supplied prime outside the native phase")
    source = 33*q
    if any(marked_guards(source, p) != (True, True) for p in (3, 11)):
        raise ValueError("marked guard failure")
    if marked_guards(source, q) != (False, False):
        raise ValueError("new-prime obstruction was lost")
    return native.receipt(source, plan.word)


SOURCES = (
    (561, (3, 11, 17)), (1105, (5, 13, 17)),
    (1729, (7, 13, 19)), (2465, (5, 17, 29)),
    (2821, (7, 13, 31)), (6601, (7, 23, 41)),
    (8911, (7, 19, 67)), (10585, (5, 29, 73)),
    (15841, (7, 31, 73)), (29341, (13, 37, 61)),
    (41041, (7, 11, 13, 41)), (46657, (13, 37, 97)),
    (52633, (7, 73, 103)), (62745, (3, 5, 47, 89)),
    (63973, (7, 13, 19, 37)), (75361, (11, 13, 17, 31)),
)
# A discovered finite word is retained as data, never rediscovered by an API.
ROOT_561 = (2, 4, 1, 1, 1, 3, 4, 1, 3, 1, 2, 3, 4)
FIRST_DESCENT_8911 = (1, 1, 1, 3, 1, 1, 1, 5)
ROOT_3569 = (2, 5, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3,
             1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4)


def ray_descent(source):
    """Actual first descent on8911+32768*t, without a child ROOT assertion."""
    native.odd(source)
    if source < 8911 or (source-8911) % 32768:
        raise ValueError("source outside the retained first-descent ray")
    t = (source-8911)//32768
    result = native.receipt(source, FIRST_DESCENT_8911)
    if result.target != 3569+13122*t or result.source-result.target != 5342+19646*t:
        raise ValueError("native ray affine identity mismatch")
    return result


def main():
    checks = 0

    def check(condition):
        nonlocal checks
        checks += 1
        if not condition:
            raise RuntimeError("control failed at " + str(checks))

    def reject(fn):
        try:
            fn()
        except ValueError:
            check(True)
            return
        check(False)

    # Independent literal loop uses repeated division, not v2 or affine carriers.
    def literal(n, w):
        for a in w:
            if n == 1:
                raise RuntimeError("unexpected early ROOT")
            z, count = 3*n+1, 0
            while z % 2 == 0:
                z //= 2
                count += 1
            check(count == a)
            n = z
        return n

    # Exact bounded discovery only for the finite declared test universe.
    def proposal(plan):
        q = plan.prime_residue
        if q < plan.minimum_prime:
            q += ((plan.minimum_prime-q+plan.prime_modulus-1)//plan.prime_modulus)*plan.prime_modulus
        for _ in range(10000):
            if prime(q):
                return q
            q += plan.prime_modulus
        raise RuntimeError("finite test proposal cap exhausted")

    rows = []
    unresolved = []
    for n, ps in SOURCES:
        s = KorseltSource(n, ps)
        conductor, slack, height = predecessor_state(s)
        check(conductor*slack == n-1 and height == v2(n-1))
        first = v2(3*n+1)
        check(((3*n+1)//(1 << first) < n) == (n % 4 == 1))
        check((n % 4 == 3) == (len(ps) % 2 == 1 and all(p % 4 == 3 for p in ps)))
        if any(p % 8 == 1 for p in ps):
            check(first == 2)
        r = descent(s)
        if r is None:
            unresolved.append(n)
            check(first == 1)
        else:
            check(literal(n, r.word) == r.target < n)
            if height >= 3:
                k = len(r.word)
                check(r.target == 1+3**k*(n-1)//4**k)
                check((r.target % 4 == 3) if height % 2 else (r.target % 8 == 5))
        rows.append({"source": n, "prime_factors": ps, "conductor": conductor,
                     "slack": slack, "height": height,
                     "paid_word": None if r is None else r.word,
                     "endpoint": None if r is None else r.target})
    check(unresolved == [8911])
    # The all-source twos-run identity has no Korselt hypothesis.
    for n in range(3, 8192, 2):
        h = v2(n-1)
        if h >= 3:
            k = (h-1)//2
            y = native.replay(n, (2,)*k)
            check(1 < y < n and y == 1+3**k*(n-1)//4**k)
            check((v2(3*y+1) == 1) if h % 2 else (v2(3*y+1) >= 3))

    partial_rows = []
    words = [w for length in range(1, 4) for w in product(range(1, 5), repeat=length)]
    words += [(1,)*k for k in (4, 8, 12)]
    for w in words:
        plan = compile_partial(w)
        check(gcd(plan.prime_residue, plan.prime_modulus) == 1)
        check(plan.prime_residue % 330 == 257)
        r, mod = native.native_cell(w)
        check(33*plan.prime_residue % mod == r)
        q = proposal(plan)
        rec = realize_partial(plan, q)
        check(literal(rec.source, w) == rec.target)
        # Reverse affine authentication uses exact integer division at each step.
        z = rec.target
        for a in reversed(w):
            z, remainder = divmod((1 << a)*z-1, 3)
            check(remainder == 0 and z % 2 == 1)
        check(z == rec.source)
        check(q > 32 and 32 % q != 0 and 32 % (q-1) != 0)
        check(marked_guards(rec.source, 3) == marked_guards(rec.source, 11) == (True, True))
        reject(lambda q=q: audit_source(KorseltSource(33*q, tuple(sorted((3, 11, q))))))
        if len(w) > 3 or w in ((1,), (2,), (4,), (1, 2, 3)):
            partial_rows.append({"word": w, "prime_residue": plan.prime_residue,
                                 "prime_modulus": plan.prime_modulus,
                                 "supplied_prime": q, "source": rec.source,
                                 "endpoint": rec.target})

    full = discharge(KorseltSource(561, (3,11,17)), ROOT_561[1:])
    check(full.target == 1 and full.word == ROOT_561)
    check(literal(561, full.word) == 1)
    # Singleton supplied data: the generic R=1 dispatcher remains unchanged.
    first_8911 = native.receipt(8911, FIRST_DESCENT_8911)
    check(first_8911.target == 3569 < 8911)
    for length in range(1, len(FIRST_DESCENT_8911)):
        check(literal(8911, FIRST_DESCENT_8911[:length]) >= 8911)
    full_8911 = native.discharge(first_8911, ROOT_3569)
    check(literal(8911, full_8911.word) == 1)
    check((len(full_8911.word), sum(full_8911.word)) == (32, 64))
    first_15841 = descent(KorseltSource(15841, (7,31,73)))
    check(first_15841.target == full_8911.source)
    full_15841 = native.concatenate(first_15841, full_8911)
    check(literal(15841, full_15841.word) == 1)
    check((len(full_15841.word), sum(full_15841.word)) == (34, 68))
    check(native.carrier(FIRST_DESCENT_8911) == (6561, 16384, 9425))
    check(native.native_cell(FIRST_DESCENT_8911) == (8911, 32768))
    for length in range(1, len(FIRST_DESCENT_8911)):
        p, q, b = native.carrier(FIRST_DESCENT_8911[:length])
        check(p > q and b > 0)
    for t in range(128):
        rec = ray_descent(8911+32768*t)
        check(literal(rec.source, rec.word) == rec.target == 3569+13122*t)
        check(rec.source-rec.target == 5342+19646*t > 0)
    check(ray_descent(8911) == first_8911)
    check(8912 & (8912-1) != 0)
    check(pow(2, 15, 32768) == 0)
    for exponent in range(1, 41):
        check((pow(2, exponent, 32768)-1) % 32768 != 8911)
    for n in (8911.0, True, 8913, 8191, 2):
        reject(lambda n=n: ray_descent(n))
    root_rows = [{"source": rr.source, "odd_rank": len(rr.word), "cost": sum(rr.word),
                  "word": rr.word, "word_sha256": hashlib.sha256(bytes(rr.word)).hexdigest()}
                 for rr in (full, full_8911, full_15841)]
    reject(lambda: native.discharge(first_8911, ROOT_3569+(2,)))
    reject(lambda: native.discharge(first_8911, ROOT_561))
    reject(lambda: native.concatenate(first_15841, full))
    reject(lambda: discharge(KorseltSource(8911, (7,19,67)), ROOT_561))
    reject(lambda: discharge(KorseltSource(561, (3,11,17)), ROOT_561))
    reject(lambda: native.discharge(full, (2,)))
    plan = compile_partial((1,))
    for bad in (True, 3.0, -1):
        reject(lambda bad=bad: prime(bad))
    for s in (KorseltSource(561.0,(3,11,17)), KorseltSource(561,(3,11,17.0)),
              KorseltSource(561,(3,3,17)), KorseltSource(561,[3,11,17]),
              KorseltSource(561,(3,11,19)), KorseltSource(33*257,(3,11,257))):
        reject(lambda s=s: audit_source(s))
    for w in ([], (), (True,), (1.0,), (0,)):
        reject(lambda w=w: compile_partial(w))
    for p in (replace(plan, prime_residue=float(plan.prime_residue)),
              replace(plan, prime_modulus=plan.prime_modulus+1),
              replace(plan, minimum_prime=True), replace(plan, word=(2,))):
        reject(lambda p=p: audit_plan(p))
    reject(lambda: realize_partial(plan, True))
    reject(lambda: realize_partial(plan, 257))  # native (1) requires the other parity class mod4
    output = {"status": "FINITE-EXACT", "checks": checks,
              "universe": "16 inherited Carmichael sources; odd n=3..8191; 84 words on1..4 of length1..3 plus ones lengths4,8,12",
              "full_korselt_sources": rows, "unresolved_first_rise": unresolved,
              "partial_prime_proposals": len(words), "selected_partial_controls": partial_rows,
              "root_controls": root_rows,
              "singleton_first_descent": {"source": 8911, "endpoint": 3569,
                                           "word": FIRST_DESCENT_8911,
                                           "odd_rank": 8, "cost": 14},
              "native_first_descent_ray": {"source_base": 8911, "source_step": 32768,
                                           "child_base": 3569, "child_step": 13122,
                                           "literal_parameter_controls": [0, 127],
                                           "positive_mersenne_intersection": "empty for every exponent",
                                           "root_grounded_parameters_in_this_packet": [0]},
              "scope": "Marked-prime tests are not a full Giuga/Korselt certificate. Paid prefixes are not ROOT certificates without a supplied suffix."}
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
