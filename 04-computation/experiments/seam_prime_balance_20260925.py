"""Divisor-balance typing, safe quartic filter, and exact profile automaton.

Pure Python; explicit finite universes printed below. General completeness
rests on the elementary proof in the matching result note.
"""
from itertools import combinations_with_replacement
from math import isqrt, prod


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def factor(n):
    result = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def profile(n):
    return tuple(sorted(factor(n).values()))


def counts_profile(a):
    if not a:
        return (0, 0, 0)
    r = len(a)
    return (prod(x + 1 for x in a) - 2,
            2 ** r - 1 - int(max(a) == 1),
            r - int(a == (1,)))


def counts_direct(n):
    divisors = [d for d in range(2, n) if n % d == 0]
    squarefree = lambda d: all(d % (k * k) for k in range(2, isqrt(d) + 1))
    prime = lambda d: all(d % k for k in range(2, isqrt(d) + 1))
    return (len(divisors), sum(squarefree(d) for d in divisors),
            sum(prime(d) for d in divisors))


def balanced(a):
    f, s, u = counts_profile(a)
    return f == s + u


LIVE = {(), (1,), (2,), (3,), (1, 1), (1, 2), (1, 1, 1), (1, 1, 2)}
ACCEPT = {(), (1,), (3,), (1, 1, 2)}


def add_new(a):
    return tuple(sorted(a + (1,)))


def raise_exponent(a, e):
    b = list(a)
    b[b.index(e)] += 1
    return tuple(sorted(b))


def state(a):
    return a if a in LIVE else "REJECT"


def divides_balanced_profile(a):
    # A divisor of a balanced shape embeds into one of the three targets.
    # Sorted coordinatewise embedding into largest slots is sufficient.
    for target in ((1,), (3,), (1, 1, 2)):
        if len(a) <= len(target):
            slots = target[len(target) - len(a):] if a else ()
            if all(x <= y for x, y in zip(a, slots)):
                return True
    return False


def main():
    print("seam_prime_balance_20260925: divisor counts, not tournament sizes")
    for n in range(1, 1001):
        check(counts_direct(n) == counts_profile(profile(n)), ("direct counts", n))
    print("Independent direct divisor counts n=1..1000: PASS")
    count = 0
    for r in range(1, 9):
        for a in combinations_with_replacement(range(1, 9), r):
            count += 1
            check(balanced(a) == (a in ACCEPT), ("classification", a))
            check((a in LIVE) == divides_balanced_profile(a), ("live quotient", a))
            if max(a) >= 4:
                f, s, u = counts_profile(a)
                bound = 3 * 2 ** (r - 1) - r - 1
                check(f - s - u >= bound > 0, ("quartic bound", a))
            if a not in LIVE:
                check(add_new(a) not in LIVE, ("new prime reject", a))
                for e in set(a):
                    check(raise_exponent(a, e) not in LIVE, ("raise reject", a, e))
    print("Profiles:r=1..8,each exponent1..8,sorted profiles", count, "PASS")
    first_mixed = next(n for n in range(2, 1001) if profile(n) == (1, 1, 2))
    check(first_mixed == 60, "first mixed")
    print("First mixed balanced example:", first_mixed, factor(first_mixed), counts_direct(first_mixed))
    for n in (1, 2, 4, 6, 8, 12, 16, 24, 27, 32, 60, 81, 120):
        f, s, u = counts_direct(n)
        print("n", n, "profile", profile(n), "(F,S,U)", (f, s, u), "D", f-s-u)
    for a in sorted(LIVE, key=lambda x: (len(x), x)):
        transitions = [("new prime", state(add_new(a)))]
        transitions.extend(("raise exponent" + str(e), state(raise_exponent(a, e)))
                           for e in sorted(set(a)))
        print("automaton", a, "accept", a in ACCEPT, transitions)
    check(profile(6) == profile(15) and not balanced(profile(24))
          and balanced(profile(60)), "prime overlap hostile")
    print("Forgetting named support:6,15 same profile;*4 gives24 false,60 true")
    check(balanced(profile(60)) and sum(profile(60)) == 4, "total quartic hostile")
    check(not balanced(profile(16)) and balanced(profile(8)), "cap exponent hostile")
    check(not balanced(profile(32)) and balanced(profile(2)), "fourth-power quotient hostile")
    check((3 * 81 + 1) // 4 == 61 and not balanced(profile(81))
          and balanced(profile(61)), "Collatz static filter hostile")
    print("Static quartic filter is not Collatz closed:16->8;odd-accelerated81->61")
    check(counts_profile((4,)) == (3, 1, 1), "p4 defect")
    check(counts_profile((1, 1, 2)) == (10, 7, 3), "mixed degree4")
    # Odd function need not take odd values on odd integers.
    check(all((n ** 3 + n) % 2 == 0 for n in range(-101, 102)), "odd function parity")
    check(all((-n) ** 3 + (-n) == -(n ** 3 + n) for n in range(-101, 102)),
          "odd function symmetry")
    print("Odd-function hostile:f(x)=x^3+x is odd but integer values always even")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
